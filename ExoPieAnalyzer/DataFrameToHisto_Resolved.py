
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
    if ApplyWeight:print 'Filling Histogram',HISTNAME, 'for Background'
    if not ApplyWeight:print 'Filling Histogram',HISTNAME, 'for Data'

    h_var  = SetHist(HISTNAME, binning)
    weight=1.0
    puweight= 1.0
    btag = 1.0
    for ij in df_var.index:
        value = df_var[ij]
        weight= df_weight[ij]
        #print df
        numerator   = df_weight_num[ij]
        denominator = df_weight_den[ij]
        scale       = numerator/denominator

        if weight==0.0:scale=1.0
        if ApplyWeight: h_var.Fill(value, weight*scale)
        if not ApplyWeight:h_var.Fill(value)

    '''
    for index, value in enumerate(df_var):
	#print 'index',index,'value',value, 'weight',df_weight[index]
        weight = df_weight[index]
        h_var.Fill(value, weight)
        #print 'weight',weight
    '''
    return h_var

def getBinRange(nBins, xlow,xhigh):
    diff = float(xhigh - xlow)/float(nBins)
    binRange = [xlow+ij*diff for ij in range(nBins+1)]
    return binRange

#def HistWrtter(df, inFile,treeName, mode="UPDATE"):
def HistWrtter(df, outfilename, treeName,mode="UPDATE"):
    h_list = []
    reg=treeName.split('_')[1]
    '''
    FjetBins = getBinRange(15,200,1000)
    leppTbins = getBinRange(15,30,500)
    fjSDBins  = getBinRange(15,100,150) 
    '''
    if 'SBand' in reg or 'SR' in reg:


        #CENTRAL AND SYSTEMATICS FOR MET HISTOGRAM
        h_list.append(VarToHist(df["MET"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET",[200,270,345,480,1000]))
        #B-TAG SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["btagweight"],df["btagweight_up"],"h_reg_"+reg+"_MET_btagweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["btagweight"],df["btagweight_down"],"h_reg_"+reg+"_MET_btagweight_down",[200,270,345,480,1000]))
        #EWK SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["ewkweight"],df["ewkweight_up"],"h_reg_"+reg+"_MET_ewkweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["ewkweight"],df["ewkweight_down"],"h_reg_"+reg+"_MET_ewkweight_down",[200,270,345,480,1000]))
        #Top pT REWEIGHTING
        h_list.append(VarToHist(df["MET"], df["weight"],df["toppTweight"],df["toppTweight_up"],"h_reg_"+reg+"_MET_toppTweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["toppTweight"],df["toppTweight_down"],"h_reg_"+reg+"_MET_toppTweight_down",[200,270,345,480,1000]))
        #MET Trigger SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["METweight"],df["METweight_up"],"h_reg_"+reg+"_MET_metTrigweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["METweight"],df["METweight_down"],"h_reg_"+reg+"_MET_metTrigweight_down",[200,270,345,480,1000]))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["lepweight"],df["lepweight_up"],"h_reg_"+reg+"_MET_lepweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["lepweight"],df["lepweight_down"],"h_reg_"+reg+"_MET_lepweight_down",[200,270,345,480,1000]))
        #pu WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["puweight"],df["puweight_up"],"h_reg_"+reg+"_MET_puweight_down",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["puweight"],df["puweight_down"],"h_reg_"+reg+"_MET_puweight_up",[200,270,345,480,1000]))
        #JEC SYSTEMATICS
        h_list.append(VarToHist(df["MET"], df["weight"],df["jec"],df["jec_up"],"h_reg_"+reg+"_MET_jec_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["MET"], df["weight"],df["jec"],df["jec_down"],"h_reg_"+reg+"_MET_jec_down",[200,270,345,480,1000]))
        #JER SYSTEMATICS
        h_list.append(VarToHist(df["METRes_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_Res_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["METRes_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_Res_down",[200,270,345,480,1000]))

        h_list.append(VarToHist(df["METEn_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_En_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["METEn_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_En_down",[200,270,345,480,1000]))
 
	#h_list.append(VarToHist(df,df["MET"], df["weight"], "h_reg_"+reg+"_MET",[200,270,345,480,1000]))
        #h_list.append(VarToHist(df["FJetPt"], df["weight"], "h_reg_"+reg+"_FJetPt",[50,200,1000]))
        #h_list.append(VarToHist(df["FJetMass"], df["weight"],"h_reg_"+reg+"_FJetMass",[35,30,350]))#FJetMass
        #h_list.append(VarToHist(df["FJetEta"], df["weight"], "h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        #h_list.append(VarToHist(df["FJetPhi"], df["weight"], "h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        #h_list.append(VarToHist(df["FJetCSV"], df["weight"], "h_reg_"+reg+"_FJetCSV",[30,0,1]))
        h_list.append(VarToHist(df["nJets"],   df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_nJets",[10,0,10]))
        #h_list.append(VarToHist(df["min_dPhi"],   df["weight"], "h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
        h_list.append(VarToHist(df["Jet1Pt"],  df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Pt",[50,30,1000]))
        h_list.append(VarToHist(df["Jet1Eta"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet1Phi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14])) #DiJetMass
        h_list.append(VarToHist(df["DiJetMass"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_DiJetMass",[50,0,400]))#FJetMass
    else:
	h_list.append(VarToHist(df["MET"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_MET",[30,0,1000]))
        #CENTRAL AND SYSTEMATICS FOR RECOIL HISTOGRAM
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil",[200,270,345,480,1000]))
        #B-TAG SYSTEMATICS
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["btagweight"],df["btagweight_up"],"h_reg_"+reg+"_Recoil_btagweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["btagweight"],df["btagweight_down"],"h_reg_"+reg+"_Recoil_btagweight_down",[200,270,345,480,1000]))
        #EWK SYSTEMATICS
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["ewkweight"],df["ewkweight_up"],"h_reg_"+reg+"_Recoil_ewkweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["ewkweight"],df["ewkweight_down"],"h_reg_"+reg+"_Recoil_ewkweight_down",[200,270,345,480,1000]))
        #Top pT REWEIGHTING
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["toppTweight"],df["toppTweight_up"],"h_reg_"+reg+"_Recoil_toppTweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["toppTweight"],df["toppTweight_down"],"h_reg_"+reg+"_Recoil_toppTweight_down",[200,270,345,480,1000]))
        #MET Trigger SYSTEMATICS
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["recoilweight"],df["recoilweight_up"],"h_reg_"+reg+"_Recoil_metTrigweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["recoilweight"],df["recoilweight_down"],"h_reg_"+reg+"_Recoil_metTrigweight_down",[200,270,345,480,1000]))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["lepweight"],df["lepweight_up"],"h_reg_"+reg+"_Recoil_lepweight_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["lepweight"],df["lepweight_down"],"h_reg_"+reg+"_Recoil_lepweight_down",[200,270,345,480,1000]))
        #pu WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["puweight"],df["puweight_up"],"h_reg_"+reg+"_Recoil_puweight_down",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["puweight"],df["puweight_down"],"h_reg_"+reg+"_Recoil_puweight_up",[200,270,345,480,1000]))
        #JEC SYSTEMATICS
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["jec"],df["jec_up"],"h_reg_"+reg+"_Recoil_jec_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["RECOIL"], df["weight"],df["jec"],df["jec_down"],"h_reg_"+reg+"_Recoil_jec_down",[200,270,345,480,1000]))
        #JER SYSTEMATICS
        h_list.append(VarToHist(df["recoilRes_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil_Res_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["recoilRes_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil_Res_down",[200,270,345,480,1000]))

        h_list.append(VarToHist(df["recoilEn_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil_En_up",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["recoilEn_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Recoil_En_down",[200,270,345,480,1000]))

        #h_list.append(VarToHist(df["FJetPt"], df["weight"], "h_reg_"+reg+"_FJetPt",[15,200,1000]))
        #h_list.append(VarToHist(df["FJetMass"], df["weight"],"h_reg_"+reg+"_FJetMass",[15,100,150]))#FJetMass        
        #h_list.append(VarToHist(df["FJetEta"], df["weight"], "h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        #h_list.append(VarToHist(df["FJetPhi"], df["weight"], "h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        #h_list.append(VarToHist(df["FJetCSV"], df["weight"], "h_reg_"+reg+"_FJetCSV",[30,0,1]))
        h_list.append(VarToHist(df["Jet1Pt"],  df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Pt",[50,30,1000]))
	h_list.append(VarToHist(df["Jet1Eta"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet1Phi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(VarToHist(df["nJets"],   df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df["min_dPhi"],   df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)  "_met_Phi",[50,-4,4]))
        h_list.append(VarToHist(df["met_Phi"],   df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_met_Phi",[50,-4,4]))
        h_list.append(VarToHist(df["lep1_pT"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep1_pT",[15,30,500]))
        h_list.append(VarToHist(df["lep1_eta"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep1_eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["lep1_Phi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep1_Phi",[30,-3.14,3.14]))       
        if 'Zmumu' in reg or 'Zee' in reg:
	    h_list.append(VarToHist(df["Zmass"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_Zmass",[15,60,120]))
	    h_list.append(VarToHist(df["ZpT"], df["weight"], df["weight"],df["weight"],"h_reg_"+reg+"_ZpT",[15,0,700])) 
            h_list.append(VarToHist(df["lep2_pT"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep2_pT",[15,30,500]))
            h_list.append(VarToHist(df["lep2_eta"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep2_eta",[30,-2.5,2.5]))
            h_list.append(VarToHist(df["lep2_Phi"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_lep2_Phi",[30,-3.14,3.14]))
    #outfilename = 'Output_'+inFile.split('/')[-1]
    fout = TFile(outfilename, mode)
    for ih in h_list: ih.Write()


def emptyHistWritter(treeName,outfilename,mode="UPDATE"):
    h_list = []
    reg=treeName.split('_')[1]
    '''
    FjetBins = getBinRange(15,200,1000)
    leppTbins = getBinRange(15,30,500)
    fjSDBins  = getBinRange(15,100,150)
    '''
    if 'SBand' in reg or 'SR' in reg:
        h_list.append(SetHist("h_reg_"+reg+"_MET",[200,270,345,480,1000]))


        h_list.append(SetHist("h_reg_"+reg+"_MET_btagweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_btagweight_down",[200,270,345,480,1000]))
        #EWK SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_ewkweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_ewkweight_down",[200,270,345,480,1000]))
        #Top pT REWEIGHTING
        h_list.append(SetHist("h_reg_"+reg+"_MET_toppTweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_toppTweight_down",[200,270,345,480,1000]))
        #MET Trigger SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_metTrigweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_metTrigweight_down",[200,270,345,480,1000]))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_lepweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_lepweight_down",[200,270,345,480,1000]))
        #pu WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_puweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_puweight_down",[200,270,345,480,1000]))
        #JEC SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_jec_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_jec_down",[200,270,345,480,1000]))
        #JER SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_Res_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_Res_down",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_En_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_En_down",[200,270,345,480,1000]))

        #h_list.append(SetHist("h_reg_"+reg+"_FJetPt",[50,200,1000]))
        #h_list.append(SetHist("h_reg_"+reg+"_FJetMass",[35,30,350]))#FJetMass
        #h_list.append(SetHist("h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        #h_list.append(SetHist("h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        #h_list.append(SetHist("h_reg_"+reg+"_FJetCSV",[30,0,1]))
	h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        #h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Pt",[50,30,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_DiJetMass",[50,0,400]))
    else:
	h_list.append(SetHist("h_reg_"+reg+"_MET",   [30,0,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil",[200,270,345,480,1000]))

        h_list.append(SetHist("h_reg_"+reg+"_Recoil_btagweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_btagweight_down",[200,270,345,480,1000]))
        #EWK SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_ewkweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_ewkweight_down",[200,270,345,480,1000]))
        #Top pT REWEIGHTING
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_toppTweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_toppTweight_down",[200,270,345,480,1000]))
        #MET Trigger SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_metTrigweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_metTrigweight_down",[200,270,345,480,1000]))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_lepweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_lepweight_down",[200,270,345,480,1000]))
        #pu WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_puweight_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_puweight_down",[200,270,345,480,1000]))
        #JEC SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_jec_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_jec_down",[200,270,345,480,1000]))
        #JER SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_Res_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_Res_down",[200,270,345,480,1000]))

        h_list.append(SetHist("h_reg_"+reg+"_Recoil_En_up",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil_En_down",[200,270,345,480,1000]))

        #h_list.append(SetHist("h_reg_"+reg+"_FJetPt",[15,200,1000]))
        #h_list.append(SetHist("h_reg_"+reg+"_FJetMass",[15,100,150]))#FJetMass
        #h_list.append(SetHist("h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        #h_list.append(SetHist("h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        #h_list.append(SetHist("h_reg_"+reg+"_FJetCSV",[30,0,1]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Pt",[50,30,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
        h_list.append(SetHist("h_reg_"+reg+"_met_Phi",[50,-4,4]))
        h_list.append(SetHist("h_reg_"+reg+"_lep1_pT",[15,30,500]))
        h_list.append(SetHist("h_reg_"+reg+"_lep1_eta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_lep1_Phi",[30,-3.14,3.14]))
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

trees =['monoHbb_Tope_resolved','monoHbb_Topmu_resolved','monoHbb_We_resolved','monoHbb_Wmu_resolved','monoHbb_Zmumu_resolved','monoHbb_Zee_resolved','monoHbb_SR_resolved','monoHbb_SBand_resolved']
# 'monoHbb_Zmumu_boosted','monoHbb_Zee_boosted'

#inputFilename=infile
filename=infile
ApplyWeight = True
def runFile(filename,trees):
    tf =  ROOT.TFile(filename)
    global ApplyWeight
    if 'SingleElectron' in filename or 'MET' in filename: ApplyWeight = False
    else:ApplyWeight = True

    #print 'ApplyWeight',ApplyWeight
    h_total = tf.Get('h_total')
    h_total_weight = tf.Get('h_total_mcweight')
    
    '''
    h_reg_WenuCR_resolved_cutFlow = tf.Get('h_reg_WenuCR_resolved_cutFlow')
    h_reg_WmunuCR_resolved_cutFlow= tf.Get('h_reg_WmunuCR_resolved_cutFlow')
    h_reg_TopenuCR_resolved_cutFlow= tf.Get('h_reg_TopenuCR_resolved_cutFlow')
    h_reg_TopmunuCR_resolved_cutFlow= tf.Get('h_reg_TopmunuCR_resolved_cutFlow')
    h_reg_SBand_resolved_cutFlow= tf.Get('h_reg_SBand_resolved_cutFlow')
    h_reg_SBand_boosted_cutFlow= tf.Get('h_reg_SBand_boosted_cutFlow')
    '''
    #print 'total',h_total_weight.Integral()
    outfilename = outputdir+'/'+'Output_'+filename.split('/')[-1]
    for index, tree in enumerate(trees):
        print 'tree',tree
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
    #f.cd()
    h_total_weight.Write()
    h_total.Write()

    '''
    h_reg_WenuCR_resolved_cutFlow.Write()
    h_reg_WmunuCR_resolved_cutFlow.Write()
    h_reg_TopenuCR_resolved_cutFlow.Write()
    h_reg_TopmunuCR_resolved_cutFlow.Write()
    h_reg_SBand_resolved_cutFlow.Write()
    h_reg_SBand_boosted_cutFlow.Write()
    '''

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
