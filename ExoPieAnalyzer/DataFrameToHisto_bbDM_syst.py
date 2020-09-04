from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain,TGraphAsymmErrors, TMath, TH2D, TLorentzVector, AddressOf, gROOT, TNamed
import ROOT as ROOT
import os,traceback
import sys, optparse,argparse
from array import array
import math
import numpy as numpy
import pandas
from root_pandas import read_root
from pandas import  DataFrame, concat
from pandas import Series
import time
import glob


## ----- start of clock
start = time.clock()


## ----- command line argument
usage = "python DataframeToHist.py -F -inDir directoryName -D outputDir "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--inputfile",  dest="inputfile",default="myfiles.root")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="out.root")
parser.add_argument("-F", "--farmout", action="store_true",  dest="farmout")
parser.add_argument("-inDir", "--inputDir",  dest="inputDir",default=".")
parser.add_argument("-D", "--outputdir", dest="outputdir",default=".")

args = parser.parse_args()

if args.farmout==None:
    isfarmout = False
else:
    isfarmout = args.farmout

if args.inputDir and isfarmout:
    inDir=args.inputDir

outputdir = '.'
if args.outputdir:
    outputdir = str(args.outputdir)


infile  = args.inputfile


args = parser.parse_args()


filename = 'OutputFiles/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8.root'


def SetHist(HISTNAME,binning):
    h=TH1F()
    if len(binning) == 3:
        h = TH1F(HISTNAME, HISTNAME, binning[0], binning[1], binning[2])
    else:
        nBins = len(binning) -1
        #h = TH1F(HISTNAME, HISTNAME, binning[0], binning[1], binning[2])  ## make it variable binning histogram
        h = TH1F(HISTNAME, HISTNAME, nBins, array('d',binning))
    return h


def VarToHist(df_var,df_weight,df_weight_den,df_weight_num,HISTNAME,binning):

    #df_var    = df[varname]
    #df_weight = df["weight"]

    h_var  = SetHist(HISTNAME, binning)
    weight=1.0
    weightPU= 1.0
    btag = 1.0
    for ij in df_var.index:
        value = df_var[ij]
        weight= df_weight[ij]
        #prefire_wgt = df["weightPrefire"][ij]
        #print df
        numerator   = df_weight_num[ij]
        if 'weightJEC' in HISTNAME:
            denominator = 1
        else:
            denominator = df_weight_den[ij]
        scale = numerator/denominator
        if '_nPV' in HISTNAME: scale = 1/denominator
        if weight==0.0:scale=1.0
        if ApplyWeight: h_var.Fill(value, (weight*scale))
        if not ApplyWeight:h_var.Fill(value)

    return h_var

def getBinRange(nBins, xlow,xhigh):
    diff = float(xhigh - xlow)/float(nBins)
    binRange = [xlow+ij*diff for ij in range(nBins+1)]
    return binRange

#def HistWrtter(df, inFile,treeName, mode="UPDATE"):
def HistWrtter(df, outfilename, treeName,mode="UPDATE"):
    h_list = []
    if 'preselR' in treeName:
        reg = treeName.split('_')[1]
    else:
        reg = treeName.split('_')[1]+'_'+treeName.split('_')[2]
    if 'SR' in reg or 'preselR' in reg:
        #CENTRAL AND SYSTEMATICS FOR MET HISTOGRAM
        h_list.append(VarToHist(df["MET"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET",[200,250,350,500,1000]))
        #B-TAG SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightB"],df["weightB_up"],"h_reg_"+reg+"_MET_weightB_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightB"],df["weightB_down"],"h_reg_"+reg+"_MET_weightB_down",[200,250,350,500,1000]))
        #EWK SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightEWK"],df["weightEWK_up"],"h_reg_"+reg+"_MET_weightEWK_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightEWK"],df["weightEWK_down"],"h_reg_"+reg+"_MET_weightEWK_down",[200,250,350,500,1000]))
        #Top pT REWEIGHTING
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightTop"],df["weightTop_up"],"h_reg_"+reg+"_MET_weightTop_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightTop"],df["weightTop_down"],"h_reg_"+reg+"_MET_weightTop_down",[200,250,350,500,1000]))
        #MET Trigger SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightMET"],df["weightMET_up"],"h_reg_"+reg+"_MET_weightMET_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightMET"],df["weightMET_down"],"h_reg_"+reg+"_MET_weightMET_down",[200,250,350,500,1000]))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightEle"],df["weightEle_up"],"h_reg_"+reg+"_MET_weightEle_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightEle"],df["weightEle_down"],"h_reg_"+reg+"_MET_weightEle_down",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightMu"],df["weightMu_up"],"h_reg_"+reg+"_MET_weightMu_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightMu"],df["weightMu_down"],"h_reg_"+reg+"_MET_weightMu_down",[200,250,350,500,1000]))
        #pu WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightPU"],df["weightPU_up"],"h_reg_"+reg+"_MET_weightPU_down",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["weightPU"],df["weightPU_down"],"h_reg_"+reg+"_MET_weightPU_up",[200,250,350,500,1000]))
        #weightJEC SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["weight"],df["weightJEC_up"],"h_reg_"+reg+"_MET_weightJEC_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["weight"],df["weightJEC_down"],"h_reg_"+reg+"_MET_weightJEC_down",[200,250,350,500,1000]))
        #JER SYSTEMATICS
        h_list.append(VarToHist(df["MET_Res_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_Res_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET_Res_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_Res_down",[200,250,350,500,1000]))

        h_list.append(VarToHist(df["MET_En_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_En_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["MET_En_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_En_down",[200,250,350,500,1000]))

        h_list.append(VarToHist(df["Njets_PassID"],   df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df["Jet1Pt"],  df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Pt",[50,30,1000]))
        h_list.append(VarToHist(df["Jet1Eta"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet1Phi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet1deepCSV"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1deepCSV",[15,0,1.1]))
        h_list.append(VarToHist(df["Jet2Pt"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet2Pt",[15,30,800]))
        h_list.append(VarToHist(df["Jet2Eta"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet2Eta",[15,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet2Phi"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet2Phi",[15,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet2deepCSV"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet2deepCSV",[15,0,1.1]))
        h_list.append(VarToHist(df["nPV"], df["weight"], df["weightPU"], df["weight"], "h_reg_"+reg+"_nPV", [70, 0, 70]))
        h_list.append(VarToHist(df["nPV"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_PUnPV",[70,0,70]))
        h_list.append(VarToHist(df["dPhi_jetMET"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_min_dPhi",[15,0.5,3.2]))
        h_list.append(VarToHist(df["METPhi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_METPhi",[15,-3.14,3.14]))
    else:
        h_list.append(VarToHist(df["MET"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_MET",[30,0,1000]))
        #CENTRAL AND SYSTEMATICS FOR Recoil HISTOGRAM
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil",[200,250,350,500,1000]))
        #B-TAG SYSTEMATICS
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightB"],df["weightB_up"],"h_reg_"+reg+"_Recoil_weightB_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightB"],df["weightB_down"],"h_reg_"+reg+"_Recoil_weightB_down",[200,250,350,500,1000]))
        #EWK SYSTEMATICS
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightEWK"],df["weightEWK_up"],"h_reg_"+reg+"_Recoil_weightEWK_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightEWK"],df["weightEWK_down"],"h_reg_"+reg+"_Recoil_weightEWK_down",[200,250,350,500,1000]))
        #Top pT REWEIGHTING
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightTop"],df["weightTop_up"],"h_reg_"+reg+"_Recoil_weightTop_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightTop"],df["weightTop_down"],"h_reg_"+reg+"_Recoil_weightTop_down",[200,250,350,500,1000]))
        #MET Trigger SYSTEMATICS
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightRecoil"],df["weightRecoil_up"],"h_reg_"+reg+"_Recoil_weightRecoil_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightRecoil"],df["weightRecoil_down"],"h_reg_"+reg+"_Recoil_weightRecoil_down",[200,250,350,500,1000]))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightEle"],df["weightEle_up"],"h_reg_"+reg+"_Recoil_weightEle_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightEle"],df["weightEle_down"],"h_reg_"+reg+"_Recoil_weightEle_down",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightMu"],df["weightMu_up"],"h_reg_"+reg+"_Recoil_weightMu_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightMu"],df["weightMu_down"],"h_reg_"+reg+"_Recoil_weightMu_down",[200,250,350,500,1000]))
        #pu WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightPU"],df["weightPU_up"],"h_reg_"+reg+"_Recoil_weightPU_down",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weightPU"],df["weightPU_down"],"h_reg_"+reg+"_Recoil_weightPU_up",[200,250,350,500,1000]))
        #weightJEC SYSTEMATICS
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weight"],df["weightJEC_up"],"h_reg_"+reg+"_Recoil_weightJEC_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil"], df["weight"],df["weight"],df["weightJEC_down"],"h_reg_"+reg+"_Recoil_weightJEC_down",[200,250,350,500,1000]))
        #JER SYSTEMATICS
        h_list.append(VarToHist(df["Recoil_Res_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil_Res_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil_Res_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil_Res_down",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil_En_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil_En_up",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Recoil_En_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil_En_down",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Jet1Pt"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Pt",[50,30,1000]))
        h_list.append(VarToHist(df["Jet1Eta"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet1Phi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet1deepCSV"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1deepCSV",[15,0,1.1]))
        h_list.append(VarToHist(df["Njets_PassID"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df["NEle"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_NEle",[10,0,10]))
        h_list.append(VarToHist(df["NMu"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_NMu",[10,0,10]))
        h_list.append(VarToHist(df["Jet2Pt"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet2Pt",[15,30,800]))
        h_list.append(VarToHist(df["Jet2Eta"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet2Eta",[15,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet2Phi"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet2Phi",[15,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet2deepCSV"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet2deepCSV",[15,0,1.1]))
        h_list.append(VarToHist(df["nPV"],df["weight"],df["weightPU"],df["weight"], "h_reg_"+reg+"_nPV",[70,0,70]))
        h_list.append(VarToHist(df["nPV"],df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_PUnPV",[70,0,70]))
        h_list.append(VarToHist(df["METPhi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_METPhi",[15,-3.14,3.14]))
        h_list.append(VarToHist(df["RecoilPhi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_RecoilPhi",[15,-3.14,3.14]))
        h_list.append(VarToHist(df["dPhi_jetMET"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_min_dPhi",[15,0.5,3.2]))#min_dPhi)
        h_list.append(VarToHist(df["leadingLepPt"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep1_pT",[15,30,500]))
        h_list.append(VarToHist(df["leadingLepEta"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep1_eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["leadingLepPhi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep1_Phi",[30,-3.14,3.14]))
        if 'munu' in reg or 'enu' in reg:
            h_list.append(VarToHist(df["Wmass"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Wmass",[15,0,160]))
            h_list.append(VarToHist(df["WpT"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_WpT",[15,0,700]))
        if 'Zmumu' in reg or 'Zee' in reg:
            h_list.append(VarToHist(df["Zmass"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Zmass",[15,60,120]))
            h_list.append(VarToHist(df["ZpT"], df["weight"], df["weight"],df["weight"],"h_reg_"+reg+"_ZpT",[15,0,700]))
            h_list.append(VarToHist(df["subleadingLepPt"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep2_pT",[15,30,500]))
            h_list.append(VarToHist(df["subleadingLepEta"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep2_eta",[30,-2.5,2.5]))
            h_list.append(VarToHist(df["subleadingLepPhi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep2_Phi",[30,-3.14,3.14]))
    #outfilename = 'Output_'+inFile.split('/')[-1]
    fout = TFile(outfilename, mode)
    for ih in h_list: ih.Write()


def emptyHistWritter(treeName,outfilename,mode="UPDATE"):
    h_list = []
    if 'preselR' in treeName:
        reg = treeName.split('_')[1]
    else:
        reg = treeName.split('_')[1]+'_'+treeName.split('_')[2]
    if 'SR' in reg or 'preselR' in reg:
        h_list.append(SetHist("h_reg_"+reg+"_MET",[200,250,350,500,1000]))

        h_list.append(SetHist("h_reg_"+reg+"_MET_weightB_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightB_down",[200,250,350,500,1000]))
        #EWK SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightEWK_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightEWK_down",[200,250,350,500,1000]))
        #Top pT REWEIGHTING
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightTop_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightTop_down",[200,250,350,500,1000]))
        #MET Trigger SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightMET_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightMET_down",[200,250,350,500,1000]))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightEle_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightEle_down",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightMu_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightMu_down",[200,250,350,500,1000]))
        #pu WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightPU_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightPU_down",[200,250,350,500,1000]))
        #weightJEC SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightJEC_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_weightJEC_down",[200,250,350,500,1000]))
        #JER SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_Res_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_Res_down",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_En_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_En_down",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        #h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Pt",[50,30,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1deepCSV",[15,0,1.1]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Pt",[15,30,800]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Eta",[15,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Phi",[15,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2deepCSV",[15,0,1.1]))
        h_list.append(SetHist("h_reg_"+reg+"_nPV",[70,0,70]))
        h_list.append(SetHist("h_reg_"+reg+"_PUnPV",[70,0,70]))
        h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[15,0.5,3.2]))
        h_list.append(SetHist("h_reg_"+reg+"_METPhi",[15,-3.14,3.14]))

    else:
        h_list.append(SetHist("h_reg_"+reg+"_MET",   [30,0,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil",[200,250,350,500,1000]))
        #btag SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightB_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightB_down",[200,250,350,500,1000]))
        #EWK SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightEWK_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightEWK_down",[200,250,350,500,1000]))
        #Top pT REWEIGHTING
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightTop_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightTop_down",[200,250,350,500,1000]))
        #MET Trigger SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightRecoil_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightRecoil_down",[200,250,350,500,1000]))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightEle_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightEle_down",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightMu_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightMu_down",[200,250,350,500,1000]))
        #pu WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightPU_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightPU_down",[200,250,350,500,1000]))
        #weightJEC SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightJEC_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_weightJEC_down",[200,250,350,500,1000]))
        #JER SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_Res_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_Res_down",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_En_up",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_En_down",[200,250,350,500,1000]))

        h_list.append(SetHist("h_reg_"+reg+"_Jet1Pt",[50,30,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1deepCSV",[15,0,1.1]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Pt",[15,30,800]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Eta",[15,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Phi",[15,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2deepCSV",[15,0,1.1]))
        h_list.append(SetHist("h_reg_"+reg+"_nPV",[70,0,70]))
        h_list.append(SetHist("h_reg_"+reg+"_PUnPV",[70,0,70]))
        h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_NEle",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_NMu",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_METPhi",[15,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_RecoilPhi",[15,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[15,0.5,3.2]))#mini_dPhi)
        h_list.append(SetHist("h_reg_"+reg+"_lep1_pT",[15,30,500]))
        h_list.append(SetHist("h_reg_"+reg+"_lep1_eta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_lep1_Phi",[30,-3.14,3.14]))
        if 'Wmunu' in reg or 'Wenu' in reg:
            h_list.append(SetHist("h_reg_"+reg+"_Wmass",[15,0,160]))
            h_list.append(SetHist("h_reg_"+reg+"_WpT",[15,0,700]))
        if 'Zmumu' in reg or 'Zee' in reg:
            h_list.append(SetHist("h_reg_"+reg+"_Zmass",[15,60,120]))
            h_list.append(SetHist("h_reg_"+reg+"_ZpT",[15,0,700]))
            h_list.append(SetHist("h_reg_"+reg+"_lep2_pT",[15,30,500]))
            h_list.append(SetHist("h_reg_"+reg+"_lep2_eta",[30,-2.5,2.5]))
            h_list.append(SetHist("h_reg_"+reg+"_lep2_Phi",[30,-3.14,3.14]))
    #outfilename = 'Output_'+inFile.split('/')[-1]
    fout = TFile(outfilename, mode)
    for ih in h_list: ih.Write()


'''
---------------------------------------------------------------
START MAKING HISTOGRAMS
---------------------------------------------------------------
'''

trees =['bbDM_preselR','bbDM_SR_1b','bbDM_SR_2b','bbDM_ZeeCR_1b','bbDM_ZeeCR_2b','bbDM_ZmumuCR_1b','bbDM_ZmumuCR_2b','bbDM_WenuCR_1b','bbDM_WenuCR_2b','bbDM_WmunuCR_1b','bbDM_WmunuCR_2b','bbDM_TopenuCR_1b','bbDM_TopenuCR_2b','bbDM_TopmunuCR_1b','bbDM_TopmunuCR_2b']

#inputFilename=infile
filename=infile
ApplyWeight = True
def runFile(filename,trees):
    tf =  ROOT.TFile(filename)
    h_reg_preselR_cutFlow = tf.Get('h_reg_preselR_cutFlow')
    h_reg_SR_1b_cutFlow = tf.Get('h_reg_SR_1b_cutFlow')
    h_reg_SR_2b_cutFlow = tf.Get('h_reg_SR_2b_cutFlow')
    h_reg_ZeeCR_1b_cutFlow = tf.Get('h_reg_ZeeCR_1b_cutFlow')
    h_reg_ZeeCR_2b_cutFlow = tf.Get('h_reg_ZeeCR_2b_cutFlow')
    h_reg_ZmumuCR_1b_cutFlow = tf.Get('h_reg_ZmumuCR_1b_cutFlow')
    h_reg_ZmumuCR_2b_cutFlow = tf.Get('h_reg_ZmumuCR_2b_cutFlow')
    h_reg_WenuCR_1b_cutFlow = tf.Get('h_reg_WenuCR_1b_cutFlow')
    h_reg_WenuCR_2b_cutFlow = tf.Get('h_reg_WenuCR_2b_cutFlow')
    h_reg_WmunuCR_1b_cutFlow = tf.Get('h_reg_WmunuCR_1b_cutFlow')
    h_reg_WmunuCR_2b_cutFlow = tf.Get('h_reg_WmunuCR_2b_cutFlow')
    h_reg_TopenuCR_1b_cutFlow = tf.Get('h_reg_TopenuCR_1b_cutFlow')
    h_reg_TopenuCR_2b_cutFlow = tf.Get('h_reg_TopenuCR_2b_cutFlow')
    h_reg_TopmunuCR_1b_cutFlow = tf.Get('h_reg_TopmunuCR_1b_cutFlow')
    h_reg_TopmunuCR_2b_cutFlow = tf.Get('h_reg_TopmunuCR_2b_cutFlow')
    global ApplyWeight
    if ('SingleElectron' in filename) or ('MET' in filename) or ('EGamma' in filename): ApplyWeight = False
    else:ApplyWeight = True

    print 'ApplyWeight',ApplyWeight
    h_total = tf.Get('h_total')
    h_total_weight = tf.Get('h_total_mcweight')
    #print 'total',h_total_weight.Integral()
    outfilename = outputdir+'/'+'Output_'+filename.split('/')[-1]
    for index, tree in enumerate(trees):
        #print 'tree',tree
        tt = tf.Get(tree)
        nent = tt.GetEntries()

        if index==0: mode="RECREATE"
        if index>0: mode="UPDATE"

        if nent > 0:
            df = read_root(filename,tree)
            df = df[df.Jet1Pt > 50.0]
            #df = df[df.nJets <=2 ]
            HistWrtter(df, outfilename,tree,mode)
        else:
            emptyHistWritter(tree,outfilename,mode)
    f = TFile(outfilename, "UPDATE")
    h_reg_preselR_cutFlow.Write()
    h_reg_SR_1b_cutFlow.Write()
    h_reg_SR_2b_cutFlow.Write()
    h_reg_ZeeCR_1b_cutFlow.Write()
    h_reg_ZeeCR_2b_cutFlow.Write()
    h_reg_ZmumuCR_1b_cutFlow.Write()
    h_reg_ZmumuCR_2b_cutFlow.Write()
    h_reg_WenuCR_1b_cutFlow.Write()
    h_reg_WenuCR_2b_cutFlow.Write()
    h_reg_WmunuCR_1b_cutFlow.Write()
    h_reg_WmunuCR_2b_cutFlow.Write()
    h_reg_TopenuCR_1b_cutFlow.Write()
    h_reg_TopenuCR_2b_cutFlow.Write()
    h_reg_TopmunuCR_1b_cutFlow.Write()
    h_reg_TopmunuCR_2b_cutFlow.Write()
    h_total_weight.Write()
    h_total.Write()

if isfarmout:
    path=inDir
    files=glob.glob(path+'/*')
    for inputFile in files:
        print 'running code for file:  ',inputFile
	runFile(inputFile,trees)

if not isfarmout:
    filename=infile
    print 'running code for file:  ',filename
    runFile(filename,trees)


stop = time.clock()
print "%.4gs" % (stop-start)
