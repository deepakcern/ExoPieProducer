
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


def VarToHist(df_var,df_weight, HISTNAME,binning):

    #df_var    = df[varname]
    #df_weight = df["weight"]

    h_var  = SetHist(HISTNAME, binning)
    weight=1.0
    for index, value in enumerate(df_var):
        #print 'index',index,'value',value, 'weight',df_weight[index]
        weight = df_weight[index]
        h_var.Fill(value, weight)

    return h_var

def getBinRange(nBins, xlow,xhigh):
    diff = float(xhigh - xlow)/float(nBins)
    binRange = [xlow+ij*diff for ij in range(nBins+1)]
    return binRange

#def HistWrtter(df, inFile,treeName, mode="UPDATE"):
def HistWrtter(df, outfilename, treeName,mode="UPDATE"):
    h_list = []
    reg=treeName.split('_')[1]+'_'+treeName.split('_')[2]
    '''
    FjetBins = getBinRange(15,200,1000)
    leppTbins = getBinRange(15,30,500)
    fjSDBins  = getBinRange(15,100,150)
    '''
    if 'SR' in reg:
        h_list.append(VarToHist(df["MET"], df["weight"], "h_reg_"+reg+"_MET",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Jet1Pt"], df["weight"], "h_reg_"+reg+"_Jet1Pt",[10,30,800]))
        h_list.append(VarToHist(df["Jet1Eta"], df["weight"], "h_reg_"+reg+"_Jet1Eta",[10,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet1Phi"], df["weight"], "h_reg_"+reg+"_Jet1Phi",[10,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet1deepCSV"], df["weight"], "h_reg_"+reg+"_Jet1deepCSV",[10,0,1.2]))
        h_list.append(VarToHist(df["Njets_PassID"],   df["weight"], "h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df["dPhi_jetMET"],   df["weight"], "h_reg_"+reg+"_min_dPhi",[20,0,4]))#mini_dPhi)
    else:
        h_list.append(VarToHist(df["MET"], df["weight"], "h_reg_"+reg+"_MET",[10,0,700]))
        h_list.append(VarToHist(df["Recoil"], df["weight"], "h_reg_"+reg+"_Recoil",[200,250,350,500,1000]))
        h_list.append(VarToHist(df["Jet1Pt"], df["weight"], "h_reg_"+reg+"_Jet1Pt",[10,30,800]))
        h_list.append(VarToHist(df["Jet1Eta"], df["weight"], "h_reg_"+reg+"_Jet1Eta",[10,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet1Phi"], df["weight"], "h_reg_"+reg+"_Jet1Phi",[10,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet1deepCSV"], df["weight"], "h_reg_"+reg+"_Jet1deepCSV",[10,0,1.2]))
        h_list.append(VarToHist(df["Jet2Pt"],  df["weight"], "h_reg_"+reg+"_Jet2Pt",[10,30,800]))
        h_list.append(VarToHist(df["Jet2Eta"], df["weight"], "h_reg_"+reg+"_Jet2Eta",[10,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet2Phi"], df["weight"], "h_reg_"+reg+"_Jet2Phi",[10,-3.14,3.14]))
        h_list.append(VarToHist(df["Jet2deepCSV"], df["weight"], "h_reg_"+reg+"_Jet2deepCSV",[10,0,1.2]))
        h_list.append(VarToHist(df["Njets_PassID"],   df["weight"], "h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df["dPhi_jetMET"],   df["weight"], "h_reg_"+reg+"_min_dPhi",[20,0,4]))#mini_dPhi)
        h_list.append(VarToHist(df["leadingLepPt"], df["weight"], "h_reg_"+reg+"_lep1_pT",[10,30,500]))
        h_list.append(VarToHist(df["leadingLepEta"], df["weight"], "h_reg_"+reg+"_lep1_eta",[10,-2.5,2.5]))
        h_list.append(VarToHist(df["leadingLepPhi"], df["weight"], "h_reg_"+reg+"_lep1_Phi",[10,-3.14,3.14]))
        if 'munu' in reg or 'enu' in reg:
            h_list.append(VarToHist(df["Wmass"], df["weight"],"h_reg_"+reg+"_Wmass",[15,0,160]))
            h_list.append(VarToHist(df["WpT"], df["weight"], "h_reg_"+reg+"_WpT",[15,0,700]))
        if 'Zmumu' in reg or 'Zee' in reg:
            h_list.append(VarToHist(df["Zmass"], df["weight"],"h_reg_"+reg+"_Zmass",[15,60,120]))
            h_list.append(VarToHist(df["ZpT"], df["weight"], "h_reg_"+reg+"_ZpT",[15,0,700]))
            h_list.append(VarToHist(df["subleadingLepPt"], df["weight"], "h_reg_"+reg+"_lep2_pT",[15,30,500]))
            h_list.append(VarToHist(df["subleadingLepEta"], df["weight"], "h_reg_"+reg+"_lep2_eta",[10,-2.5,2.5]))
            h_list.append(VarToHist(df["subleadingLepPhi"], df["weight"], "h_reg_"+reg+"_lep2_Phi",[10,-3.14,3.14]))
    #outfilename = 'Output_'+inFile.split('/')[-1]
    fout = TFile(outfilename, mode)
    for ih in h_list: ih.Write()


def emptyHistWritter(treeName,outfilename,mode="UPDATE"):
    h_list = []
    reg=treeName.split('_')[1]+'_'+treeName.split('_')[2]
    '''
    FjetBins = getBinRange(15,200,1000)
    leppTbins = getBinRange(15,30,500)
    fjSDBins  = getBinRange(15,100,150)
    '''
    if 'SR' in reg:
        h_list.append(SetHist("h_reg_"+reg+"_MET",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Pt",[10,30,800]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Eta",[10,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Phi",[10,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1deepCSV",[10,0,1.2]))
        h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[20,0,4]))#mini_dPhi)
    else:
        h_list.append(SetHist("h_reg_"+reg+"_MET",[20,0,700]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil",[200,250,350,500,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Pt",[10,30,800]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Eta",[10,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Phi",[10,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1deepCSV",[10,0,1.2]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Pt",[10,30,800]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Eta",[10,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2Phi",[10,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet2deepCSV",[10,0,1.2]))
        h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[20,0,4]))#mini_dPhi)
        h_list.append(SetHist("h_reg_"+reg+"_lep1_pT",[10,30,500]))
        h_list.append(SetHist("h_reg_"+reg+"_lep1_eta",[10,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_lep1_Phi",[10,-3.14,3.14]))
        if 'Wmunu' in reg or 'Wenu' in reg:
            h_list.append(SetHist("h_reg_"+reg+"_Wmass",[30,0,160]))
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

trees =['bbDM_SR_1b','bbDM_SR_2b','bbDM_ZeeCR_1b','bbDM_ZeeCR_2b','bbDM_ZmumuCR_1b','bbDM_ZmumuCR_2b','bbDM_WenuCR_1b','bbDM_WenuCR_2b','bbDM_WmunuCR_1b','bbDM_WmunuCR_2b','bbDM_TopenuCR_1b','bbDM_TopenuCR_2b','bbDM_TopmunuCR_1b','bbDM_TopmunuCR_2b']
# 'monoHbb_Zmumu_boosted','monoHbb_Zee_boosted'

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
