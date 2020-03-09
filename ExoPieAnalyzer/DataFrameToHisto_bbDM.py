
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


def VarToHist(df_var,df,df_weight, HISTNAME,binning):
    #df_var    = df[varname]
    #df_weight = df,df["weight"]
    df_weightEWK = df["weightEWK"]
    h_var  = SetHist(HISTNAME, binning)
    weight=1.0
    for index, value in enumerate(df_var):
        #print 'index',index,'value',value, 'weight',df_weight[index]
        weight = df_weight[index]
        weightEWK = df_weightEWK[index]
        if '_nPV' in HISTNAME:
            h_var.Fill(value, 1)
        else:
            h_var.Fill(value, weight)
            #h_var.Fill(value, weight/weightEWK)
    return h_var

def getBinRange(nBins, xlow,xhigh):
    diff = float(xhigh - xlow)/float(nBins)
    binRange = [xlow+ij*diff for ij in range(nBins+1)]
    return binRange

#def HistWrtter(df, inFile,treeName, mode="UPDATE"):
def HistWrtter(df, outfilename, treeName,mode="UPDATE"):
    h_list = []
    reg=treeName.split('_')[1]+'_'+treeName.split('_')[2]

    if 'SR' in reg:
        h_list.append(VarToHist(df["MET"], df,df["weight"], "h_reg_"+reg+"_MET",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Jet1Pt"], df,df["weight"], "h_reg_"+reg+"_Jet1Pt",[15,30,800]))
        h_list.append(VarToHist(df["Jet1Eta"], df,df["weight"], "h_reg_"+reg+"_Jet1Eta",[15,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet1Phi"], df,df["weight"], "h_reg_"+reg+"_Jet1Phi",[15,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet1deepCSV"], df,df["weight"], "h_reg_"+reg+"_Jet1deepCSV",[15,0,1.1]))
        h_list.append(VarToHist(df["Jet2Pt"], df,df["weight"], "h_reg_"+reg+"_Jet2Pt",[15,30,800]))
        h_list.append(VarToHist(df["Jet2Eta"], df,df["weight"], "h_reg_"+reg+"_Jet2Eta",[15,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet2Phi"], df,df["weight"], "h_reg_"+reg+"_Jet2Phi",[15,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet2deepCSV"], df,df["weight"], "h_reg_"+reg+"_Jet2deepCSV",[15,0,1.1]))
        h_list.append(VarToHist(df["Njets_PassID"],   df,df["weight"], "h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df["dPhi_jetMET"],   df,df["weight"], "h_reg_"+reg+"_min_dPhi",[15,0.5,3.2]))#mini_dPhi)
        h_list.append(VarToHist(df["nPV"],   df,df["weight"], "h_reg_"+reg+"_nPV",[70,0,70]))
        h_list.append(VarToHist(df["nPV"],   df,df["weightPU"], "h_reg_"+reg+"_PUnPV",[70,0,70]))
    else:
        h_list.append(VarToHist(df["MET"], df,df["weight"], "h_reg_"+reg+"_MET",[15,0,700]))
        h_list.append(VarToHist(df["Recoil"], df,df["weight"], "h_reg_"+reg+"_Recoil",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Jet1Pt"], df,df["weight"], "h_reg_"+reg+"_Jet1Pt",[15,30,800]))
        h_list.append(VarToHist(df["Jet1Eta"], df,df["weight"], "h_reg_"+reg+"_Jet1Eta",[15,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet1Phi"], df,df["weight"], "h_reg_"+reg+"_Jet1Phi",[15,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet1deepCSV"], df,df["weight"], "h_reg_"+reg+"_Jet1deepCSV",[15,0,1.1]))
        h_list.append(VarToHist(df["Jet2Pt"],  df,df["weight"], "h_reg_"+reg+"_Jet2Pt",[15,30,800]))
        h_list.append(VarToHist(df["Jet2Eta"], df,df["weight"], "h_reg_"+reg+"_Jet2Eta",[15,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet2Phi"], df,df["weight"], "h_reg_"+reg+"_Jet2Phi",[15,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet2deepCSV"], df,df["weight"], "h_reg_"+reg+"_Jet2deepCSV",[15,0,1.1]))
        h_list.append(VarToHist(df["Njets_PassID"],   df,df["weight"], "h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df["dPhi_jetMET"],   df,df["weight"], "h_reg_"+reg+"_min_dPhi",[15,0.5,3.2]))#mini_dPhi)
        h_list.append(VarToHist(df["leadingLepPt"], df,df["weight"], "h_reg_"+reg+"_lep1_pT",[15,30,500]))
        h_list.append(VarToHist(df["leadingLepEta"], df,df["weight"], "h_reg_"+reg+"_lep1_eta",[15,-2.5,2.5]))
        h_list.append(VarToHist(df["leadingLepPhi"], df,df["weight"], "h_reg_"+reg+"_lep1_Phi",[15,-3.14,3.14]))
        h_list.append(VarToHist(df["nPV"],   df,df["weight"], "h_reg_"+reg+"_nPV",[70,0,70]))
        h_list.append(VarToHist(df["nPV"],   df,df["weightPU"], "h_reg_"+reg+"_PUnPV",[70,0,70]))
        if 'munu' in reg or 'enu' in reg:
            h_list.append(VarToHist(df["Wmass"], df,df["weight"],"h_reg_"+reg+"_Wmass",[15,0,160]))
            h_list.append(VarToHist(df["WpT"], df,df["weight"], "h_reg_"+reg+"_WpT",[15,0,700]))
        if 'Zmumu' in reg or 'Zee' in reg:
            h_list.append(VarToHist(df["Zmass"], df,df["weight"],"h_reg_"+reg+"_Zmass",[15,60,120]))
            h_list.append(VarToHist(df["ZpT"], df,df["weight"], "h_reg_"+reg+"_ZpT",[15,0,700]))
            h_list.append(VarToHist(df["subleadingLepPt"], df,df["weight"], "h_reg_"+reg+"_lep2_pT",[15,10,500]))
            h_list.append(VarToHist(df["subleadingLepEta"], df,df["weight"], "h_reg_"+reg+"_lep2_eta",[15,-2.5,2.5]))
            h_list.append(VarToHist(df["subleadingLepPhi"], df,df["weight"], "h_reg_"+reg+"_lep2_Phi",[15,-3.14,3.14]))
    fout = TFile(outfilename, mode)
    for ih in h_list: ih.Write()


def emptyHistWritter(treeName,outfilename,mode="UPDATE"):
    h_list = []
    reg=treeName.split('_')[1]+'_'+treeName.split('_')[2]
    if 'SR' in reg:
        h_list.append(SetHist("h_reg_"+reg+"_MET",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Pt",[15,30,800]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Eta",[15,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Phi",[15,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1deepCSV",[15,0,1.1]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Pt",[15,30,800]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Eta",[15,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Phi",[15,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2deepCSV",[15,0,1.1]))
        h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[15,0.5,3.2]))#mini_dPhi)
        h_list.append(SetHist("h_reg_"+reg+"_nPV",[70,0,70]))
        h_list.append(SetHist("h_reg_"+reg+"_PUnPV",[70,0,70]))
    else:
        h_list.append(SetHist("h_reg_"+reg+"_MET",[15,0,700]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Pt",[15,30,800]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Eta",[15,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Phi",[15,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1deepCSV",[15,0,1.1]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Pt",[15,30,800]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Eta",[15,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Phi",[15,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2deepCSV",[15,0,1.1]))
        h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[15,0.5,3.2]))#mini_dPhi)
        h_list.append(SetHist("h_reg_"+reg+"_lep1_pT",[15,30,500]))
        h_list.append(SetHist("h_reg_"+reg+"_lep1_eta",[15,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_lep1_Phi",[15,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_nPV",[70,0,70]))
        h_list.append(SetHist("h_reg_"+reg+"_PUnPV",[70,0,70]))
        if 'munu' in reg or 'enu' in reg:
            h_list.append(SetHist("h_reg_"+reg+"_Wmass",[15,0,160]))
            h_list.append(SetHist("h_reg_"+reg+"_WpT",[15,0,700]))
        if 'Zmumu' in reg or 'Zee' in reg:
            h_list.append(SetHist("h_reg_"+reg+"_Zmass",[15,60,120]))
            h_list.append(SetHist("h_reg_"+reg+"_ZpT",[15,0,700]))
            h_list.append(SetHist("h_reg_"+reg+"_lep2_pT",[15,10,500]))
            h_list.append(SetHist("h_reg_"+reg+"_lep2_eta",[15,-2.5,2.5]))
            h_list.append(SetHist("h_reg_"+reg+"_lep2_Phi",[15,-3.14,3.14]))
    #outfilename = 'Output_'+inFile.split('/')[-1]
    fout = TFile(outfilename, mode)
    for ih in h_list: ih.Write()

'''
---------------------------------------------------------------
START MAKING HISTOGRAMS
---------------------------------------------------------------
'''

trees =['bbDM_SR_1b','bbDM_SR_2b','bbDM_ZeeCR_1b','bbDM_ZeeCR_2b','bbDM_ZmumuCR_1b','bbDM_ZmumuCR_2b','bbDM_WenuCR_1b','bbDM_WenuCR_2b','bbDM_WmunuCR_1b','bbDM_WmunuCR_2b','bbDM_TopenuCR_1b','bbDM_TopenuCR_2b','bbDM_TopmunuCR_1b','bbDM_TopmunuCR_2b']

#inputFilename=infile
filename=infile

def runFile(filename,trees):
    tf =  ROOT.TFile(filename)
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
    h_total = tf.Get('h_total')
    h_total_weight = tf.Get('h_total_mcweight')
    outfilename = outputdir+'/'+'Output_'+filename.split('/')[-1]
    for index, tree in enumerate(trees):
        tt = tf.Get(tree)
        nent = tt.GetEntries()

        if index==0: mode="RECREATE"
        if index>0: mode="UPDATE"

        if nent > 0:
            df = read_root(filename,tree)
            HistWrtter(df, outfilename,tree,mode)
        else:
            emptyHistWritter(tree,outfilename,mode)

    f = TFile(outfilename, "UPDATE")
    #f.cd()
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
