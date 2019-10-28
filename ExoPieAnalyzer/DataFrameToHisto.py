
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
        #print 'weight',weight
    
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
	h_list.append(VarToHist(df["MET"], df["weight"], "h_reg_"+reg+"_MET",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["FJetPt"], df["weight"], "h_reg_"+reg+"_FJetPt",[50,200,1000]))
        h_list.append(VarToHist(df["FJetMass"], df["weight"],"h_reg_"+reg+"_FJetMass",[35,30,350]))#FJetMass
        h_list.append(VarToHist(df["FJetEta"], df["weight"], "h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["FJetPhi"], df["weight"], "h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        h_list.append(VarToHist(df["FJetCSV"], df["weight"], "h_reg_"+reg+"_FJetCSV",[30,0,1]))
        h_list.append(VarToHist(df["nJets"],   df["weight"], "h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df["min_dPhi"],   df["weight"], "h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
    else:
	h_list.append(VarToHist(df["MET"], df["weight"], "h_reg_"+reg+"_MET",[200,270,345,480,1000]))
	h_list.append(VarToHist(df["RECOIL"], df["weight"], "h_reg_"+reg+"_Recoil",[200,270,345,480,1000]))
        h_list.append(VarToHist(df["FJetPt"], df["weight"], "h_reg_"+reg+"_FJetPt",[15,200,1000]))
        h_list.append(VarToHist(df["FJetMass"], df["weight"],"h_reg_"+reg+"_FJetMass",[15,100,150]))#FJetMass        
        h_list.append(VarToHist(df["FJetEta"], df["weight"], "h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["FJetPhi"], df["weight"], "h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        h_list.append(VarToHist(df["FJetCSV"], df["weight"], "h_reg_"+reg+"_FJetCSV",[30,0,1]))
        h_list.append(VarToHist(df["Jet1Pt"],  df["weight"], "h_reg_"+reg+"_Jet1Pt",[50,200,1000]))
	h_list.append(VarToHist(df["Jet1Eta"], df["weight"], "h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["Jet1Phi"], df["weight"], "h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(VarToHist(df["nJets"],   df["weight"], "h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df["min_dPhi"],   df["weight"], "h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
        h_list.append(VarToHist(df["lep1_pT"], df["weight"], "h_reg_"+reg+"_lep1_pT",[15,30,500]))
        h_list.append(VarToHist(df["lep1_eta"], df["weight"], "h_reg_"+reg+"_lep1_eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["lep1_Phi"], df["weight"], "h_reg_"+reg+"_lep1_Phi",[30,-3.14,3.14]))       
        if 'Zmumu' in reg or 'Zee' in reg:
	    h_list.append(VarToHist(df["Zmass"], df["weight"],"h_reg_"+reg+"_Zmass",[15,60,120]))
	    h_list.append(VarToHist(df["ZpT"], df["weight"], "h_reg_"+reg+"_ZpT",[15,0,700])) 
            h_list.append(VarToHist(df["lep2_pT"], df["weight"], "h_reg_"+reg+"_lep2_pT",[15,30,500]))
            h_list.append(VarToHist(df["lep2_eta"], df["weight"], "h_reg_"+reg+"_lep2_eta",[30,-2.5,2.5]))
            h_list.append(VarToHist(df["lep2_Phi"], df["weight"], "h_reg_"+reg+"_lep2_Phi",[30,-3.14,3.14]))
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
        h_list.append(SetHist("h_reg_"+reg+"_FJetPt",[50,200,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetMass",[35,30,350]))#FJetMass
        h_list.append(SetHist("h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetCSV",[30,0,1]))
	h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
    else:
	h_list.append(SetHist("h_reg_"+reg+"_MET",   [200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Recoil",[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetPt",[15,200,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetMass",[15,100,150]))#FJetMass
        h_list.append(SetHist("h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetCSV",[30,0,1]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Pt",[50,200,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
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

trees =['monoHbb_SR_boosted','monoHbb_Tope_boosted','monoHbb_Topmu_boosted','monoHbb_We_boosted','monoHbb_Wmu_boosted','monoHbb_TopWmu_boosted','monoHbb_TopWe_boosted','monoHbb_Zmumu_boosted','monoHbb_Zee_boosted','monoHbb_SBand_boosted']
# 'monoHbb_Zmumu_boosted','monoHbb_Zee_boosted'

#inputFilename=infile
filename=infile

def runFile(filename,trees):
    tf =  ROOT.TFile(filename)
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
