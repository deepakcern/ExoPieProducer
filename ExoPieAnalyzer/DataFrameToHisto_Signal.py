
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


def VarToHist(df,df_var,df_weight,df_weight_den,df_weight_num,HISTNAME,binning,boost=False):#(df,df_var,df_weight, HISTNAME,binning,boost=False):

    #df_var    = df[varname]
    #df_weight = df["weight"]

    h_var  = SetHist(HISTNAME, binning)
    weight=1.0
    if not boost:
	df["FJetPt"]=1.0#
	print "filling hist for boosted"
    else:print "filling hist for resolved"

    scale=1.0
    for value,weight,numerator,denominator,fjetpT in zip(df_var,df_weight,df_weight_num,df_weight_den,df["FJetPt"]):
        #value = df_var[ij]
        #weight= df_weight[ij]

        #numerator   = df_weight_num[ij]
        #denominator = df_weight_den[ij]
        if denominator!=0:scale       = numerator/denominator


        if boost:
            #fjetpT = df["FJetPt"][ij]
	    #print "filling hist for boosted"
            if fjetpT < 350 and boost:
                h_var.Fill(value, weight*0.82*scale)
            if fjetpT > 350 and boost:
	        h_var.Fill(value, weight*0.77*scale)
        else:
            #print "filling hist for resolved"
            h_var.Fill(value, weight*scale)

    return h_var

def getBinRange(nBins, xlow,xhigh):
    diff = float(xhigh - xlow)/float(nBins)
    binRange = [xlow+ij*diff for ij in range(nBins+1)]
    return binRange

#def HistWrtter(df, inFile,treeName, mode="UPDATE"):
def HistWrtter(df, outfilename, treeName,mode="UPDATE"):
    h_list = []
    reg=treeName.split('_')[1]
    cat = treeName.split('_')[-1]
    boost=False
    if 'boosted' in cat:boost=True
    '''
    FjetBins = getBinRange(15,200,1000)
    leppTbins = getBinRange(15,30,500)
    fjSDBins  = getBinRange(15,100,150)
    '''
    if 'SBand' in reg or 'SR' in reg:
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_MET_"+cat,[200,270,345,480,1000],boost))
        if boost:h_list.append(VarToHist(df,df["FJetN2b1"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_FJetN2b1_"+cat,[40,-0.5,0.5]))
        if boost:h_list.append(VarToHist(df,df["N2DDT"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_N2DDT_"+cat,[40,-0.5,0.5]))
        #B-TAG SYSTEMATICS
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["btagweight"],df["btagweight_up"],"h_reg_"+reg+"_MET_btagweight_up_"+cat,[200,270,345,480,1000],boost))
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["btagweight"],df["btagweight_down"],"h_reg_"+reg+"_MET_btagweight_down_"+cat,[200,270,345,480,1000],boost))
        #EWK SYSTEMATICS
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["ewkweight"],df["ewkweight_up"],"h_reg_"+reg+"_MET_ewkweight_up_"+cat,[200,270,345,480,1000],boost))
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["ewkweight"],df["ewkweight_down"],"h_reg_"+reg+"_MET_ewkweight_down_"+cat,[200,270,345,480,1000],boost))
        #Top pT REWEIGHTING
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["toppTweight"],df["toppTweight_up"],"h_reg_"+reg+"_MET_toppTweight_up_"+cat,[200,270,345,480,1000],boost))
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["toppTweight"],df["toppTweight_down"],"h_reg_"+reg+"_MET_toppTweight_down_"+cat,[200,270,345,480,1000],boost))
        #MET Trigger SYSTEMATICS
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["METweight"],df["METweight_up"],"h_reg_"+reg+"_MET_metTrigweight_up_"+cat,[200,270,345,480,1000],boost))
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["METweight"],df["METweight_down"],"h_reg_"+reg+"_MET_metTrigweight_down_"+cat,[200,270,345,480,1000],boost))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["lepweight"],df["lepweight_up"],"h_reg_"+reg+"_MET_lepweight_up_"+cat,[200,270,345,480,1000],boost))
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["lepweight"],df["lepweight_down"],"h_reg_"+reg+"_MET_lepweight_down_"+cat,[200,270,345,480,1000],boost))
        #pu WEIGHT SYSTEMATICS
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["puweight"],df["puweight_up"],"h_reg_"+reg+"_MET_puweight_down_"+cat,[200,270,345,480,1000],boost))
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["puweight"],df["puweight_down"],"h_reg_"+reg+"_MET_puweight_up_"+cat,[200,270,345,480,1000],boost))
        #JEC SYSTEMATICS
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["jec"],df["jec_up"],"h_reg_"+reg+"_MET_jec_up_"+cat,[200,270,345,480,1000],boost))
        h_list.append(VarToHist(df,df["MET"], df["weight"],df["jec"],df["jec_down"],"h_reg_"+reg+"_MET_jec_down_"+cat,[200,270,345,480,1000],boost))
        #JER SYSTEMATICS
        h_list.append(VarToHist(df,df["METRes_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_Res_up_"+cat,[200,270,345,480,1000],boost))
        h_list.append(VarToHist(df,df["METRes_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_Res_down_"+cat,[200,270,345,480,1000],boost))

        h_list.append(VarToHist(df,df["METEn_up"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_En_up_"+cat,[200,270,345,480,1000],boost))
        h_list.append(VarToHist(df,df["METEn_down"], df["weight"],df["weight"],df["weight"],"h_reg_"+reg+"_MET_En_down_"+cat,[200,270,345,480,1000],boost))

        '''
        h_list.append(VarToHist(df["FJetPt"], df["weight"], "h_reg_"+reg+"_FJetPt",[50,200,1000]))
        h_list.append(VarToHist(df["FJetMass"], df["weight"],"h_reg_"+reg+"_FJetMass",[35,30,350]))#FJetMass
        h_list.append(VarToHist(df["FJetEta"], df["weight"], "h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df["FJetPhi"], df["weight"], "h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        h_list.append(VarToHist(df["FJetCSV"], df["weight"], "h_reg_"+reg+"_FJetCSV",[30,0,1]))
        '''
        h_list.append(VarToHist(df,df["nJets"],   df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_nJets_"+cat,[10,0,10],boost))
        #h_list.append(VarToHist(df["min_dPhi"],   df["weight"], "h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)  min_dphi_jets
        if 'boosted' in cat and 'SR' in reg:
		h_list.append(VarToHist(df,df["min_dphi_jets"],   df["weight"],df["weight"],df["weight"], "h_reg_"+reg+"_min_dphi_jets_"+cat,[50,0,4],boost))#mini_dPhi)
		#h_list.append(VarToHist(df,df["max_dphi_jets"],   df["weight"], "h_reg_"+reg+"_max_dphi_jets_"+cat,[50,0,4],boost))
    else:
        print ''
        '''
	h_list.append(VarToHist(df,df["MET"], df["weight"], "h_reg_"+reg+"_MET",[200,270,345,480,1000]))
	h_list.append(VarToHist(df,df["RECOIL"], df["weight"], "h_reg_"+reg+"_Recoil",[200,270,345,480,1000]))
        h_list.append(VarToHist(df,df["FJetPt"], df["weight"], "h_reg_"+reg+"_FJetPt",[15,200,1000]))
        h_list.append(VarToHist(df,df["FJetMass"], df["weight"],"h_reg_"+reg+"_FJetMass",[15,100,150]))#FJetMass
        h_list.append(VarToHist(df,df["FJetEta"], df["weight"], "h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df,df["FJetPhi"], df["weight"], "h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        h_list.append(VarToHist(df,df["FJetCSV"], df["weight"], "h_reg_"+reg+"_FJetCSV",[30,0,1]))
        h_list.append(VarToHist(df,df["Jet1Pt"],  df["weight"], "h_reg_"+reg+"_Jet1Pt",[50,200,1000]))
	h_list.append(VarToHist(df,df["Jet1Eta"], df["weight"], "h_reg_"+reg+"_Jet1Eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df,df["Jet1Phi"], df["weight"], "h_reg_"+reg+"_Jet1Phi",[30,-3.14,3.14]))
        h_list.append(VarToHist(df,df["nJets"],   df["weight"], "h_reg_"+reg+"_nJets",[10,0,10]))
        h_list.append(VarToHist(df,df["min_dPhi"],   df["weight"], "h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
        h_list.append(VarToHist(df,df["lep1_pT"], df["weight"], "h_reg_"+reg+"_lep1_pT",[15,30,500]))
        h_list.append(VarToHist(df,df["lep1_eta"], df["weight"], "h_reg_"+reg+"_lep1_eta",[30,-2.5,2.5]))
        h_list.append(VarToHist(df,df["lep1_Phi"], df["weight"], "h_reg_"+reg+"_lep1_Phi",[30,-3.14,3.14]))
        if 'Zmumu' in reg or 'Zee' in reg:
	    h_list.append(VarToHist(df,df["Zmass"], df["weight"],"h_reg_"+reg+"_Zmass",[15,60,120]))
	    h_list.append(VarToHist(df,df["ZpT"], df["weight"], "h_reg_"+reg+"_ZpT",[15,0,700]))
            h_list.append(VarToHist(df,df["lep2_pT"], df["weight"], "h_reg_"+reg+"_lep2_pT",[15,30,500]))
            h_list.append(VarToHist(df,df["lep2_eta"], df["weight"], "h_reg_"+reg+"_lep2_eta",[30,-2.5,2.5]))
            h_list.append(VarToHist(df,df["lep2_Phi"], df["weight"], "h_reg_"+reg+"_lep2_Phi",[30,-3.14,3.14]))
        '''
    #outfilename = 'Output_'+inFile.split('/')[-1]
    fout = TFile(outfilename, mode)
    for ih in h_list: ih.Write()


def emptyHistWritter(treeName,outfilename,mode="UPDATE"):
    h_list = []
    reg=treeName.split('_')[1]
    cat = treeName.split('_')[-1]
    '''
    FjetBins = getBinRange(15,200,1000)
    leppTbins = getBinRange(15,30,500)
    fjSDBins  = getBinRange(15,100,150)
    '''
    if 'SBand' in reg or 'SR' in reg:
        h_list.append(SetHist("h_reg_"+reg+"_MET_"+cat,[200,270,345,480,1000]))#200,270,345,480,1000 # [50,0,700]  FJetN2b1
        if boost:h_list.append(SetHist("h_reg_"+reg+"_FJetN2b1_"+cat,[40,-0.5,0.5]))
        if boost:h_list.append(SetHist("h_reg_"+reg+"_N2DDT_"+cat,[40,-0.5,0.5]))
        '''
        h_list.append(SetHist("h_reg_"+reg+"_FJetPt",[50,200,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetMass",[35,30,350]))#FJetMass
        h_list.append(SetHist("h_reg_"+reg+"_FJetEta",[30,-2.5,2.5]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetPhi",[30,-3.14,3.14]))
        h_list.append(SetHist("h_reg_"+reg+"_FJetCSV",[30,0,1]))
        '''
	h_list.append(SetHist("h_reg_"+reg+"_nJets_"+cat,[10,0,10]))

        h_list.append(SetHist("h_reg_"+reg+"_MET_btagweight_up_"+cat,[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_btagweight_down_"+cat,[200,270,345,480,1000]))
        #EWK SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_ewkweight_up_"+cat,[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_ewkweight_down_"+cat,[200,270,345,480,1000]))
        #Top pT REWEIGHTING
        h_list.append(SetHist("h_reg_"+reg+"_MET_toppTweight_up_"+cat,[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_toppTweight_down_"+cat,[200,270,345,480,1000]))
        #MET Trigger SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_metTrigweight_up_"+cat,[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_metTrigweight_down_"+cat,[200,270,345,480,1000]))
        #LEPTON WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_lepweight_up_"+cat,[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_lepweight_down_"+cat,[200,270,345,480,1000]))
        #pu WEIGHT SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_puweight_up_"+cat,[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_puweight_down_"+cat,[200,270,345,480,1000]))
        #JEC SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_jec_up_"+cat,[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_jec_down_"+cat,[200,270,345,480,1000]))
        #JER SYSTEMATICS
        h_list.append(SetHist("h_reg_"+reg+"_MET_Res_up_"+cat,[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_Res_down_"+cat,[200,270,345,480,1000]))

        h_list.append(SetHist("h_reg_"+reg+"_MET_En_up_"+cat,[200,270,345,480,1000]))
        h_list.append(SetHist("h_reg_"+reg+"_MET_En_down_"+cat,[200,270,345,480,1000]))

        #h_list.append(SetHist("h_reg_"+reg+"_min_dPhi",[50,0,4]))#mini_dPhi)
        if 'boosted' in cat and 'SR' in reg:
		h_list.append(SetHist("h_reg_"+reg+"_min_dphi_jets_"+cat,[50,0,4]))
		h_list.append(SetHist("h_reg_"+reg+"_max_dphi_jets_"+cat,[50,0,4]))
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

#trees =['monoHbb_SR_boosted','monoHbb_Tope_boosted','monoHbb_Topmu_boosted','monoHbb_We_boosted','monoHbb_Wmu_boosted','monoHbb_TopWmu_boosted','monoHbb_TopWe_boosted','monoHbb_Zmumu_boosted','monoHbb_Zee_boosted','monoHbb_SBand_boosted']
trees =[ 'monoHbb_SR_boosted','monoHbb_SR_resolved']

#inputFilename=infile
filename=infile

ApplyN2DDT = False

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
	    #df = read_root(filename,tree)
            #df = df[df.isak4JetBasedHemEvent==0]
            #df = df[df.isak8JetBasedHemEvent==0]
            #df = df[df.ismetphiBasedHemEvent1==0]
            #df = df[df.nJets <=2 ]
            if 'resolved' in tree:
                df = df[df.Jet1Pt > 50.0]
                print 'applying jet1 pT for tree', tree
	    #if 'boosted' in tree and ApplyN2DDT:
	#	df = df[df.N2DDT < 0]
	#	df = df[df.N2DDT > -1]
	#	print 'applying N2DDT cut for tree', tree
            HistWrtter(df, outfilename,tree,mode)
        else:
            #print 'writing empty histograms for tree : ',tree
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
