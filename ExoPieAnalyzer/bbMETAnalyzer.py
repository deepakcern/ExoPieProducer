#!/usr/bin/env python
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

## for parallel threads in interactive running
from multiprocessing import Process
import multiprocessing as mp


isCondor = True
runInteractive = False
testing=True
## from commonutils
if isCondor:sys.path.append('ExoPieUtils/commonutils/')
else:sys.path.append('../../ExoPieUtils/commonutils/')
#import MathUtils as mathutil
from MathUtils import *
import BooleanUtils as boolutil


## from analysisutils
if isCondor:sys.path.append('ExoPieUtils/analysisutils/')
else:sys.path.append('../../ExoPieUtils/analysisutils/')
import analysis_utils as anautil


sys.path.append('configs')
import variables as var
import outvars_bbDM as out
import getRecoil as getRecoil

## from analysisutils
if isCondor:sys.path.append('ExoPieUtils/scalefactortools/')
else:sys.path.append('../../ExoPieUtils/scalefactortools/')


year_file= open("Year.py","w")


######################################################################################################
## All import are done before this
######################################################################################################

## ----- start of clock
start = time.clock()



## ----- command line argument
usage = "analyzer for bb+DM (istestging) "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--inputfile",  dest="inputfile",default="myfiles.txt")
parser.add_argument("-inDir", "--inputDir",  dest="inputDir",default=".")
parser.add_argument("-runOnTXT", "--runOnTXT",action="store_true", dest="runOnTXT")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="out.root")
parser.add_argument("-D", "--outputdir", dest="outputdir")
parser.add_argument("-F", "--farmout", action="store_true",  dest="farmout")
parser.add_argument("-T", "--testing", action="store_true",  dest="testing")
parser.add_argument("-y", "--year", dest="year", default="Year")

args = parser.parse_args()

if args.farmout==None:
    isfarmout = False
else:
    isfarmout = args.farmout

if args.testing==None:
    istest = False
else:
    istest = args.testing

if args.inputDir and isfarmout:
    dirName=args.inputDir

runOnTxt=False
if args.runOnTXT:
    runOnTxt = True

if args.year=='2016':
    print 'code is running for 2016'
    year_file.write('era="2016"')
elif args.year=='2017':
    print 'code is running for 2017'
    year_file.write('era="2017"')
elif args.year=='2018':
    print 'code is running for 2018'
    year_file.write('era="2018"')
else:
    print 'please provide year'
    sys.exit()
year_file.close()

if isfarmout:
    infile  = args.inputfile

else: print "No file is provided for farmout"


import ana_weight as wgt
from Year import era

outputdir = '.'
if args.outputdir:
    outputdir = str(args.outputdir)

infilename = "ExoPieTuples.root"

outDir=outputdir

def TextToList(textfile):
    return([iline.rstrip() for iline in open(textfile)])

## the input file list and key is caught in one variable as a python list,
#### first element is the list of rootfiles
#### second element is the key, user to name output.root
def getJECWeight(ep_THINjetCorrUnc):
    JECWeight_up   = 1.0
    JECWeight_down = 1.0
    for corr in ep_THINjetCorrUnc:
        JECWeight_up   *= (1+corr)
        JECWeight_down *= (1-corr)
    return JECWeight_up, JECWeight_down

def weight_(common_weight,ep_pfMetCorrPt,ep_ZmumuRecoil,ep_WmunuRecoil,nEle,ep_elePt,ep_eleEta,nMu,ep_muPt,ep_muEta):
    tot_weight = 1.0; weightMET = 1.0; weightEle = 1.0; weightMu = 1.0; weightRecoil = 1.0; weightEleTrig=1.0
    weightMET_up = 1.0; weightEle_up = 1.0; weightMu_up = 1.0; weightRecoil_up = 1.0; weightEleTrig_up = 1.0
    weightMET_down = 1.0; weightEle_down = 1.0; weightMu_down = 1.0; weightRecoil_down = 1.0; weightEleTrig_down = 1.0
    if (nEle==0 and nMu==0):
        if ep_pfMetCorrPt > 200:
            weightMET,weightMET_up,weightMET_down=wgt.getMETtrig_First(ep_pfMetCorrPt,'R')
        tot_weight = weightMET*common_weight

    if (nEle>0 and nMu==0):
        weightEleTrig = wgt.eletrig_weight(ep_elePt[0],ep_eleEta[0])[0]
        if (nEle==1):
            weightEle,weightEle_up,weightEle_down=wgt.ele_weight(ep_elePt[0],ep_eleEta[0],'T')
            weightEleTrig,weightEleTrig_up,weightEleTrig_down = wgt.eletrig_weight(ep_elePt[0],ep_eleEta[0])
            tot_weight = weightEle*common_weight*weightEleTrig
        if (nEle==2):
            weightEle = wgt.ele_weight(ep_elePt[0],ep_eleEta[0],'T')[0] * wgt.ele_weight(ep_elePt[1],ep_eleEta[1],'L')[0]
            weightEleTrig = wgt.eletrig_weight(ep_elePt[0],ep_eleEta[0])[0]

            weightEle_up = wgt.ele_weight(ep_elePt[0],ep_eleEta[0],'T')[1] * wgt.ele_weight(ep_elePt[1],ep_eleEta[1],'L')[1]
            weightEleTrig_up = wgt.eletrig_weight(ep_elePt[0],ep_eleEta[0])[1]

            weightEle_down = wgt.ele_weight(ep_elePt[0],ep_eleEta[0],'T')[2] * wgt.ele_weight(ep_elePt[1],ep_eleEta[1],'L')[2]
            weightEleTrig_down = wgt.eletrig_weight(ep_elePt[0],ep_eleEta[0])[2]

            tot_weight = weightEle*common_weight*weightEleTrig

    if (nEle==0 and nMu==1):
        weightMu,weightMu_up,weightMu_down=wgt.mu_weight(ep_muPt[0],ep_muEta[0],'T')
        if ep_WmunuRecoil>200:
            weightRecoil,weightRecoil_up,weightRecoil_down=wgt.getMETtrig_First(ep_WmunuRecoil,'R')
        tot_weight = weightMu*common_weight*weightRecoil
    if (nEle==0 and nMu==2):
        weightMu=wgt.mu_weight(ep_muPt[0],ep_muEta[0],'T')[0]*wgt.mu_weight(ep_muPt[1],ep_muEta[1],'L')[0]
        weightMu_up=wgt.mu_weight(ep_muPt[0],ep_muEta[0],'T')[1]*wgt.mu_weight(ep_muPt[1],ep_muEta[1],'L')[1]
        weightMu_down=wgt.mu_weight(ep_muPt[0],ep_muEta[0],'T')[2]*wgt.mu_weight(ep_muPt[1],ep_muEta[1],'L')[2]
        if ep_ZmumuRecoil>200:
            weightRecoil,weightRecoil_up,weightRecoil_down=wgt.getMETtrig_First(ep_ZmumuRecoil,'R')
        tot_weight = weightMu*common_weight*weightRecoil

    ele_wgt = [weightEle*weightEleTrig,weightEle_up*weightEleTrig_up,weightEle_down*weightEleTrig_down]
    mu_wgt = [weightMu,weightMu_up,weightMu_down]
    met_wgt = [weightMET,weightMET_up,weightMET_down]
    recoil_wgt = [weightRecoil,weightRecoil_up,weightRecoil_down]
    return tot_weight,weightEleTrig,ele_wgt,mu_wgt,recoil_wgt,met_wgt


dummy = -9999.0
def runbbdm(txtfile):

    print "in main function"

    infile_=[]
    outfilename=""
    prefix="Skimmed_"
    ikey_ = ""

    if  runInteractive:
        print "running for ", txtfile[0]
        infile_  = TextToList(txtfile[0])
        key_=txtfile[1]
        outfilename= txtfile[0].split('/')[-1].replace('.root.txt','.root')#prefix+key_+".root"
        # print "running for ", txtfile
        # infile_  = TextToList(txtfile)
        # outfilename= outDir+'/'+txtfile.split('/')[-1].replace('.txt','.root')#prefix+key_+".root"

    if not runInteractive:
        infile_=TextToList(txtfile)
        prefix_ = '' #'/eos/cms/store/group/phys_exotica/bbMET/2017_skimmedFiles/locallygenerated/'
        if outputdir!='.': prefix_ = outputdir+'/'
        print "prefix_", prefix_
        outfilename = prefix_+txtfile.split('/')[-1].replace('.txt','.root')#"SkimmedTree.root"
        print 'outfilename',  outfilename


    ## define global dataframes
    df_out_SR_1b = out.df_out_SR_1b
    df_out_SR_2b = out.df_out_SR_2b

    df_out_ZeeCR_1b = out.df_out_ZeeCR_1b
    df_out_ZeeCR_2b = out.df_out_ZeeCR_2b
    df_out_ZmumuCR_1b = out.df_out_ZmumuCR_1b
    df_out_ZmumuCR_2b = out.df_out_ZmumuCR_2b

    df_out_WenuCR_1b = out.df_out_WenuCR_1b
    df_out_WenuCR_2b = out.df_out_WenuCR_2b
    df_out_WmunuCR_1b = out.df_out_WmunuCR_1b
    df_out_WmunuCR_2b = out.df_out_WmunuCR_2b

    df_out_TopenuCR_1b = out.df_out_TopenuCR_1b
    df_out_TopenuCR_2b = out.df_out_TopenuCR_2b
    df_out_TopmunuCR_1b = out.df_out_TopmunuCR_1b
    df_out_TopmunuCR_2b = out.df_out_TopmunuCR_2b

    h_total = TH1F('h_total','h_total',2,0,2)
    h_total_mcweight = TH1F('h_total_mcweight','h_total_mcweight',2,0,2)
    h_eventCounter = TH1F('h_eventCounter','h_eventCounter',2,0.5,2.5)

    h_reg_SR_1b_cutFlow  = TH1F("h_reg_SR_1b_cutFlow", "h_reg_SR_1b_cutFlow", 7,0,7)
    h_reg_SR_2b_cutFlow  = TH1F("h_reg_SR_2b_cutFlow", "h_reg_SR_2b_cutFlow", 7,0,7)

    h_reg_ZeeCR_1b_cutFlow = TH1F("h_reg_ZeeCR_1b_cutFlow", "h_reg_ZeeCR_1b_cutFlow", 9,0,9)
    h_reg_ZeeCR_2b_cutFlow = TH1F("h_reg_ZeeCR_2b_cutFlow", "h_reg_ZeeCR_2b_cutFlow", 9,0,9)
    h_reg_ZmumuCR_1b_cutFlow = TH1F("h_reg_ZmumuCR_1b_cutFlow", "h_reg_ZmumuCR_1b_cutFlow", 9,0,9)
    h_reg_ZmumuCR_2b_cutFlow = TH1F("h_reg_ZmumuCR_2b_cutFlow", "h_reg_ZmumuCR_2b_cutFlow", 9,0,9)

    h_reg_WenuCR_1b_cutFlow = TH1F("h_reg_WenuCR_1b_cutFlow", "h_reg_WenuCR_1b_cutFlow", 9,0,9)
    h_reg_WenuCR_2b_cutFlow = TH1F("h_reg_WenuCR_2b_cutFlow", "h_reg_WenuCR_2b_cutFlow", 9,0,9)
    h_reg_WmunuCR_1b_cutFlow = TH1F("h_reg_WmunuCR_1b_cutFlow", "h_reg_WmunuCR_1b_cutFlow", 9,0,9)
    h_reg_WmunuCR_2b_cutFlow = TH1F("h_reg_WmunuCR_2b_cutFlow", "h_reg_WmunuCR_2b_cutFlow", 9,0,9)

    h_reg_TopenuCR_1b_cutFlow = TH1F("h_reg_TopenuCR_1b_cutFlow", "h_reg_TopenuCR_1b_cutFlow", 9,0,9)
    h_reg_TopenuCR_2b_cutFlow = TH1F("h_reg_TopenuCR_2b_cutFlow", "h_reg_TopenuCR_2b_cutFlow", 9,0,9)
    h_reg_TopmunuCR_1b_cutFlow = TH1F("h_reg_TopmunuCR_1b_cutFlow", "h_reg_TopmunuCR_1b_cutFlow", 9,0,9)
    h_reg_TopmunuCR_2b_cutFlow = TH1F("h_reg_TopmunuCR_2b_cutFlow", "h_reg_TopmunuCR_2b_cutFlow", 9,0,9)

    for infl in infile_:
        f_tmp = TFile.Open(infl,'READ')
        h_tmp = f_tmp.Get('h_total')
        h_tmp_weight = f_tmp.Get('h_total_mcweight')
        h_tmp_count = f_tmp.Get('h_eventCounter')
        h_total.Add(h_tmp)
        h_total_mcweight.Add(h_tmp_weight)
        h_eventCounter.Add(h_tmp_count)

    filename = infile_
    ieve = 0;icount = 0

    SR1bcount = 0.0; SR2bcount = 0.0
    ZeeCR1bcount=0.0; ZeeCR2bcount=0.0; ZmumuCR1count=0.0; ZmumuCR2count=0.0
    WenuCR1bcount=0.0; WenuCR2bcount=0.0; WmunuCR1bcount=0.0; WmunuCR2bcount=0.0
    TopenuCR1bcount=0.0; TopenuCR2bcount=0.0; TopmunuCR1bcount=0.0; TopmunuCR2bcount=0.0


    for df in read_root(filename, 'outTree', columns=var.allvars_bbDM, chunksize=125000):
        for ep_runId, ep_lumiSection, ep_eventId, \
        ep_prefiringweight,ep_prefiringweightup,ep_prefiringweightdown,\
        ep_pfMetCorrPt, ep_pfMetCorrPhi, ep_pfMetUncJetResUp, ep_pfMetUncJetResDown, ep_pfMetUncJetEnUp, ep_pfMetUncJetEnDown, \
        ep_WenuPhi, ep_WmunuPhi, ep_ZeePhi, ep_ZmumuPhi, \
        ep_ZeeRecoil, ep_ZmumuRecoil, ep_WenuRecoil, ep_WmunuRecoil, \
        ep_Zeemass, ep_Zmumumass, ep_Wenumass, ep_Wmunumass, \
        ep_isData, \
        ep_THINnJet, ep_THINjetPx, ep_THINjetPy, ep_THINjetPz, ep_THINjetEnergy, \
        ep_THINjetDeepCSV, ep_THINjetHadronFlavor, ep_THINjetNPV, \
        ep_THINjetCorrUnc, \
        ep_nEle, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy, \
        ep_eleIsPassTight, ep_eleIsPassLoose,ep_eleCharge, \
        ep_nPho, ep_phoIsPassTight, ep_phoPx, ep_phoPy, ep_phoPz, ep_phoEnergy, \
        ep_nMu, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy, ep_isTightMuon,ep_muCharge, \
        ep_nTau_DRBased_EleMuVeto,ep_nTau_discBased_looseElelooseMuVeto,ep_nTau_discBased_looseEleTightMuVeto,ep_nTau_discBased_mediumElelooseMuVeto,ep_nTau_discBased_TightEleTightMuVeto,\
        ep_pu_nTrueInt, ep_pu_nPUVert, \
        ep_THINjetNPV, \
        ep_mcweight, ep_genParPt, ep_genParSample, eletrigdecision, mutrigdecision, mettrigdecision \
        in zip(df.st_runId, df.st_lumiSection, df.st_eventId, \
               df.st_prefiringweight,df.st_prefiringweightup,df.st_prefiringweightdown,\
               df.st_pfMetCorrPt, df.st_pfMetCorrPhi, df.st_pfMetUncJetResUp, df.st_pfMetUncJetResDown, df.st_pfMetUncJetEnUp, df.st_pfMetUncJetEnDown, \
               df.WenuPhi, df.WmunuPhi, df.ZeePhi, df.ZmumuPhi, \
               df.ZeeRecoil, df.ZmumuRecoil, df.WenuRecoil, df.WmunuRecoil, \
               df.ZeeMass, df.ZmumuMass, df.Wenumass, df.Wmunumass, \
               df.st_isData, \
               df.st_THINnJet, df.st_THINjetPx, df.st_THINjetPy, df.st_THINjetPz, df.st_THINjetEnergy, \
               df.st_THINjetDeepCSV, df.st_THINjetHadronFlavor, df.st_THINjetNPV, \
               df.st_THINjetCorrUnc, \
               df.st_nEle, df.st_elePx, df.st_elePy, df.st_elePz, df.st_eleEnergy, \
               df.st_eleIsPassTight, df.st_eleIsPassLoose,df.st_eleCharge, \
               df.st_nPho, df.st_phoIsPassTight, df.st_phoPx, df.st_phoPy, df.st_phoPz, df.st_phoEnergy, \
               df.st_nMu, df.st_muPx, df.st_muPy, df.st_muPz, df.st_muEnergy, df.st_isTightMuon,df.st_muCharge, \
               df.st_nTau_DRBased_EleMuVeto,df.st_nTau_discBased_looseElelooseMuVeto,df.st_nTau_discBased_looseEleTightMuVeto,df.st_nTau_discBased_mediumElelooseMuVeto,df.st_nTau_discBased_TightEleTightMuVeto,\
               df.st_pu_nTrueInt, df.st_pu_nPUVert, \
               df.st_THINjetNPV, \
               df.mcweight, df.st_genParPt, df.st_genParSample, df.st_eletrigdecision, df.st_mutrigdecision, df.st_mettrigdecision):
            ieve = ieve + 1
            if ieve%10000==0: print "Processed",ieve,"Events"
            #if (ep_pfMetCorrPt <= 200.0) and (ep_ZeeRecoil <= 200.0) and (ep_ZmumuRecoil <= 200.0) and (ep_WenuRecoil <= 200.0) and (ep_WmunuRecoil <= 200.0) : continue
            isSR1b=False
            is1bCRWenu=False
            is1bCRWmunu=False
            is1bCRZee=False
            is1bCRZmumu=False
            is1bCRTopenu=False
            is1bCRTopmunu=False

            isSR2b=False
            is2bCRWenu=False
            is2bCRWmunu=False
            is2bCRZee=False
            is2bCRZmumu=False
            is2bCRTopenu=False
            is2bCRTopmunu=False

            #deepCSV_Med = 0.8484  # for old DMsimp sample, this deepcsv means CSVv2
            if era=='2016':
                deepCSV_Med = 0.6321
            elif era=='2017':
                deepCSV_Med = 0.4941
            elif era=='2018':
                deepCSV_Med = 0.4184

            '''
            -------------------------------------------------------------------------------
            electron VARS
            -------------------------------------------------------------------------------
            '''
            ep_nEle_       = [ij for ij in range(ep_nEle) if (ep_eleIsPassLoose[ij])]
            ep_nEle_index  = len(ep_nEle_)
            ep_elePt       = [getPt(ep_elePx[ij], ep_elePy[ij]) for ij in ep_nEle_]
            ep_eleEta      = [getEta(ep_elePx[ij], ep_elePy[ij], ep_elePz[ij]) for ij in ep_nEle_]
            ep_elePhi      = [getPhi(ep_elePx[ij], ep_elePy[ij]) for ij in ep_nEle_]
            ep_eleIsPTight = [ep_eleIsPassTight[ij] for ij in ep_nEle_]
            if era=="2016":
                minElePt   = 30.0
            else:
                minElePt   = 35.0

            '''
            -------------------------------------------------------------------------------
            muon VARS
            -------------------------------------------------------------------------------
            '''
            ep_muPt  = [getPt(ep_muPx[ij], ep_muPy[ij]) for ij in range(ep_nMu)]
            ep_muEta = [getEta(ep_muPx[ij], ep_muPy[ij], ep_muPz[ij]) for ij in range(ep_nMu)]
            ep_muPhi = [getPhi(ep_muPx[ij], ep_muPy[ij]) for ij in range(ep_nMu)]
            minMuPt  = 30.0
            '''

            -------------------------------------------------------------------------------
            photon VARS
            -------------------------------------------------------------------------------
            '''
            ep_phoPt = [getPt(ep_phoPx[ij], ep_phoPy[ij]) for ij in range(ep_nPho)]
            ep_phoEta = [getEta(ep_phoPx[ij], ep_phoPy[ij], ep_phoPz[ij]) for ij in range(ep_nPho)]
            ep_phoPhi = [getPhi(ep_phoPx[ij], ep_phoPy[ij]) for ij in range(ep_nPho)]

            myphotons = [True for ij in range(ep_nPho)]
            myeleBooleans = [True for ij in range(ep_nEle_index)]
            mymuBooleans = [True for ij in range(ep_nMu)]
            cleanedPho_ag_ele = []; cleanedPho_ag_mu = [];pass_pho_index_cleaned=[]
            if ep_nPho > 0: #and ep_nEle > 0:
                cleanedPho_ag_ele = anautil.jetcleaning(myphotons, myeleBooleans, ep_phoEta, ep_eleEta, ep_phoPhi, ep_elePhi, 0.4)
                cleanedPho_ag_mu  = anautil.jetcleaning(myphotons, mymuBooleans, ep_phoEta, ep_muEta, ep_phoPhi, ep_muPhi, 0.4)
                cleanedPhoton     = boolutil.logical_AND_List2(cleanedPho_ag_ele,cleanedPho_ag_mu)
                pass_pho_index_cleaned = boolutil.WhereIsTrue(cleanedPhoton)
                #print 'cleanedPho_ag_ele',cleanedPho_ag_ele, 'cleanedPho_ag_mu', cleanedPho_ag_mu
            nPho = len(pass_pho_index_cleaned)


            '''
            -------------------------------------------------------------------------------
            THIN JET VARS
            -------------------------------------------------------------------------------
            '''
            ep_THINjetPt = [getPt(ep_THINjetPx[ij], ep_THINjetPy[ij]) for ij in range(ep_THINnJet)]
            ep_THINjetEta = [getEta(ep_THINjetPx[ij], ep_THINjetPy[ij], ep_THINjetPz[ij]) for ij in range(ep_THINnJet)]
            ep_THINjetPhi = [getPhi(ep_THINjetPx[ij], ep_THINjetPy[ij]) for ij in range(ep_THINnJet)]

            if era=='2016':
                ep_THINbjets_index = [ij for ij in range(ep_THINnJet) if (ep_THINjetDeepCSV[ij] > deepCSV_Med and abs(ep_THINjetEta[ij]) < 2.4)]
            else:
                ep_THINbjets_index = [ij for ij in range(ep_THINnJet) if (ep_THINjetDeepCSV[ij] > deepCSV_Med and abs(ep_THINjetEta[ij]) < 2.5)]
            nBjets = len(ep_THINbjets_index)

            if len(ep_THINjetPt)==0 : continue

            min_dPhi_jet_MET = min([DeltaPhi(jet_phi,ep_pfMetCorrPhi) for jet_phi in ep_THINjetPhi])

            Jet2Pt  = dummy;Jet2Eta     = dummy
            Jet2Phi = dummy;Jet2deepCSV = dummy
            Jet3Pt  = dummy;Jet3Eta     = dummy
            Jet3Phi = dummy;Jet3deepCSV = dummy

            '''
            -------------------------------------------------------------------------------
            HADRONIC RECOIL
            -------------------------------------------------------------------------------
            '''
            #======   usage: ZRecoil_Phi_Zmass(nEle, eleCharge_, elepx_, elepy_, elepz_, elee_,met_,metphi_)=====
            ep_ZeeRecoil,ep_ZeeRecoil_dPhi,ep_Zeemass = getRecoil.ZRecoil_Phi_Zmass(ep_nEle_index, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetCorrPt, ep_pfMetCorrPhi)
            ep_ZeeRecoilResUp,ep_ZeeRecoil_dPhiResUp,ep_ZeemassResUp = getRecoil.ZRecoil_Phi_Zmass(ep_nEle_index, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetUncJetResUp, ep_pfMetCorrPhi)
            ep_ZeeRecoilResDown,ep_ZeeRecoil_dPhiResDown,ep_ZeemassResDown = getRecoil.ZRecoil_Phi_Zmass(ep_nEle_index, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetUncJetResDown, ep_pfMetCorrPhi)
            ep_ZeeRecoilEnUp,ep_ZeeRecoil_dPhiEnUp,ep_ZeemassEnUp = getRecoil.ZRecoil_Phi_Zmass(ep_nEle_index, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetUncJetEnUp, ep_pfMetCorrPhi)
            ep_ZeeRecoilEnDown,ep_ZeeRecoil_dPhiEnDown,ep_ZeemassEnDown = getRecoil.ZRecoil_Phi_Zmass(ep_nEle_index, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetUncJetEnDown, ep_pfMetCorrPhi)

            ep_ZmumuRecoil,ep_ZmumuRecoil_dPhi,ep_Zmumumass = getRecoil.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetCorrPt, ep_pfMetCorrPhi)
            ep_ZmumuRecoilResUp,ep_ZmumuRecoil_dPhiResUp,ep_ZmumumassResUp = getRecoil.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetUncJetResUp,ep_pfMetCorrPhi)
            ep_ZmumuRecoilResDown,ep_ZmumuRecoil_dPhiResDown,ep_ZmumumassResDown = getRecoil.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetUncJetResDown,ep_pfMetCorrPhi)
            ep_ZmumuRecoilEnUp,ep_ZmumuRecoil_dPhiEnUp,ep_ZmumumassEnUp = getRecoil.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetUncJetEnUp,ep_pfMetCorrPhi)
            ep_ZmumuRecoilEnDown,ep_ZmumuRecoil_dPhiEnDown,ep_ZmumumassEnDown = getRecoil.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetUncJetEnDown,ep_pfMetCorrPhi)

            ep_WenuRecoil, ep_WenuRecoil_dPhi, ep_Wenumass = getRecoil.WRecoil_Phi_Wmass(ep_nEle_index,ep_elePt,ep_elePhi,ep_elePx,ep_elePy,ep_pfMetCorrPt,ep_pfMetCorrPhi)
            ep_WenuRecoilResUp, ep_WenuRecoil_dPhiResUp, ep_WenumassResUp = getRecoil.WRecoil_Phi_Wmass(ep_nEle_index,ep_elePt,ep_elePhi,ep_elePx,ep_elePy,ep_pfMetUncJetResUp,ep_pfMetCorrPhi)
            ep_WenuRecoilResDown, ep_WenuRecoil_dPhiResDown, ep_WenumassResDown = getRecoil.WRecoil_Phi_Wmass(ep_nEle_index,ep_elePt,ep_elePhi,ep_elePx,ep_elePy,ep_pfMetUncJetResDown,ep_pfMetCorrPhi)
            ep_WenuRecoilEnUp, ep_WenuRecoil_dPhiEnUp, ep_WenumassEnUp = getRecoil.WRecoil_Phi_Wmass(ep_nEle_index,ep_elePt,ep_elePhi,ep_elePx,ep_elePy,ep_pfMetUncJetEnUp,ep_pfMetCorrPhi)
            ep_WenuRecoilEnDown, ep_WenuRecoil_dPhiEnDown, ep_WenumassEnDown = getRecoil.WRecoil_Phi_Wmass(ep_nEle_index,ep_elePt,ep_elePhi,ep_elePx,ep_elePy,ep_pfMetUncJetEnDown,ep_pfMetCorrPhi)

            ep_WmunuRecoil, ep_WmunuRecoil_dPhi, ep_Wmunumass = getRecoil.WRecoil_Phi_Wmass(ep_nMu,ep_muPt,ep_muPhi,ep_muPx,ep_muPy,ep_pfMetCorrPt,ep_pfMetCorrPhi)
            ep_WmunuRecoilResUp, ep_WmunuRecoil_dPhiResUp, ep_WmunumassResUp = getRecoil.WRecoil_Phi_Wmass(ep_nMu,ep_muPt,ep_muPhi,ep_muPx,ep_muPy,ep_pfMetUncJetResUp, ep_pfMetCorrPhi)
            ep_WmunuRecoilResDown, ep_WmunuRecoil_dPhiResDown, ep_WmunumassResDown = getRecoil.WRecoil_Phi_Wmass(ep_nMu,ep_muPt,ep_muPhi,ep_muPx, ep_muPy,ep_pfMetUncJetResDown,ep_pfMetCorrPhi)
            ep_WmunuRecoilEnUp, ep_WmunuRecoil_dPhiEnUp, ep_WmunumassEnUp = getRecoil.WRecoil_Phi_Wmass(ep_nMu,ep_muPt,ep_muPhi, ep_muPx, ep_muPy,ep_pfMetUncJetEnUp, ep_pfMetCorrPhi)
            ep_WmunuRecoilEnDown, ep_WmunuRecoil_dPhiEnDown, ep_WmunumassEnDown = getRecoil.WRecoil_Phi_Wmass(ep_nMu,ep_muPt,ep_muPhi, ep_muPx, ep_muPy,ep_pfMetUncJetEnDown, ep_pfMetCorrPhi)

            if (ep_pfMetCorrPt <= 200.0) and (ep_ZeeRecoil <= 200.0) and (ep_ZmumuRecoil <= 200.0) and (ep_WenuRecoil <= 200.0) and (ep_WmunuRecoil <= 200.0) : continue

            '''
            -------------------------------------------------------------------------------
            CR VARS
            -------------------------------------------------------------------------------
            '''
            if ep_nEle_index ==2: ZpT_ee = math.sqrt( (ep_elePx[0]+ep_elePx[1])*(ep_elePx[0]+ep_elePx[1]) + (ep_elePy[0]+ep_elePy[1])*(ep_elePy[0]+ep_elePy[1]))
            if ep_nMu==2: ZpT_mumu = math.sqrt( (ep_muPx[0]+ep_muPx[1])*(ep_muPx[0]+ep_muPx[1]) + (ep_muPy[0]+ep_muPy[1])*(ep_muPy[0]+ep_muPy[1]))

            if ep_nEle_index ==1: WpT_enu = math.sqrt((ep_pfMetCorrPt*math.cos(ep_pfMetCorrPhi) + ep_elePx[0])**2 + ( ep_pfMetCorrPt*math.sin(ep_pfMetCorrPhi) + ep_elePy[0])**2)
            if ep_nMu==1: WpT_munu = math.sqrt((ep_pfMetCorrPt*math.cos(ep_pfMetCorrPhi) + ep_muPx[0])**2 + ( ep_pfMetCorrPt*math.sin(ep_pfMetCorrPhi) + ep_muPy[0])**2)

            '''
            --------------------------------------------------------------------------------
            COMMON WEIGHT CALCULATION FOR ALL REGIONS
            --------------------------------------------------------------------------------
            '''
            weight = presel_weight = weightPU = weightB = weightEWK = weightQCD = weightTop = weightEleTrig = weightEle = weightMu = weightMET = weightRecoil = weightPrefire = -999.0
            weightB_up =  weightEWK_up =  weightQCD_up =  weightTop_up = weightJEC_up = weightEleTrig_up =  weightEle_up =  weightMu_up =  weightMET_up =  weightRecoil_up = weightPU_up = weightJEC_up = weightPrefire_up =  1.0
            weightB_down =  weightEWK_down =  weightQCD_down =  weightTop_down = weightJEC_down = weightEleTrig_down =  weightEle_down =  weightMu_down =  weightMET_down =  weightRecoil_down = weightPU_down = weightJEC_down = weightPrefire_down =  1.0
            if ep_isData:
                weight = presel_weight = weightPU = weightB = weightEWK = weightQCD = weightTop = weightEleTrig = weightEle = weightMu = weightMET = weightRecoil = weightPrefire = 1.0
                weightB_up =  weightEWK_up =  weightQCD_up =  weightTop_up = weightJEC_up = weightEleTrig_up =  weightEle_up =  weightMu_up =  weightMET_up =  weightRecoil_up = weightPU_up = weightJEC_up = weightPrefire_up =  1.0
                weightB_down =  weightEWK_down =  weightQCD_down =  weightTop_down = weightJEC_down = weightEleTrig_down =  weightEle_down =  weightMu_down =  weightMET_down =  weightRecoil_down = weightPU_down = weightJEC_down = weightPrefire_down = 1.0
            else:
                weightB,weightB_up,weightB_down = wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV,'MWP')
                weightPU = wgt.puweight(ep_pu_nTrueInt)[0]
                weightEWK = 1.0; weightQCD = 1.0; weightTop = 1.0
                if ep_genParSample == 23 and len(ep_genParPt) > 0 :
                    weightEWK = wgt.getEWKZ(ep_genParPt[0])
                    weightEWK_up = wgt.getEWKZ(ep_genParPt[0])*1.5; weightEWK_down = wgt.getEWKZ(ep_genParPt[0])*0.5
                    weightQCD = wgt.getQCDZ(ep_genParPt[0])
                    weightQCD_up = wgt.getQCDZ(ep_genParPt[0]); weightQCD_down = wgt.getQCDZ(ep_genParPt[0])
                if ep_genParSample == 24 and len(ep_genParPt) > 0 :
                    weightEWK = wgt.getEWKW(ep_genParPt[0])
                    weightEWK_up = wgt.getEWKW(ep_genParPt[0])*1.5; weightEWK_down = wgt.getEWKW(ep_genParPt[0])*0.5
                    weightQCD = wgt.getQCDW(ep_genParPt[0])
                    weightQCD_up = wgt.getQCDW(ep_genParPt[0]); weightQCD_down = wgt.getQCDW(ep_genParPt[0])
                if ep_genParSample == 6 and len(ep_genParPt) > 0:
                    weightTop,weightTop_up,weightTop_down = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                weightPrefire = ep_prefiringweight
                weightPrefire_up = ep_prefiringweightup
                weightPrefire_down = ep_prefiringweightdown
                common_weight = weightB * weightEWK * weightQCD * weightTop * weightPU * weightPrefire
                presel_weight = weightEWK * weightQCD * weightTop * weightPU * weightPrefire
                weight,weightEleTrig,ele_wgt,mu_wgt,recoil_wgt,met_wgt = weight_(common_weight,ep_pfMetCorrPt,ep_ZmumuRecoil,ep_WmunuRecoil,ep_nEle_index,ep_elePt,ep_eleEta,ep_nMu,ep_muPt,ep_muEta)
                weightEle=ele_wgt[0];weightMu=mu_wgt[0];weightRecoil=recoil_wgt[0];weightMET=met_wgt[0]
                weightEle_up=ele_wgt[1];weightMu_up=mu_wgt[1];weightRecoil_up=recoil_wgt[1];weightMET_up=met_wgt[1]
                weightEle_down=ele_wgt[2];weightMu_down=mu_wgt[2];weightRecoil_down=recoil_wgt[2];weightMET_down=met_wgt[2]
                weightJEC_up = getJECWeight(ep_THINjetCorrUnc)[0]; weightJEC_down = getJECWeight(ep_THINjetCorrUnc)[1]

            '''
            --------------------------------------------------------------------------------
            SIGNAL REGION
            --------------------------------------------------------------------------------
            '''
            h_reg_SR_1b_cutFlow.AddBinContent(1, presel_weight)
            h_reg_SR_2b_cutFlow.AddBinContent(1, presel_weight)
            if mettrigdecision:
                h_reg_SR_1b_cutFlow.AddBinContent(2, presel_weight*weightMET)
                h_reg_SR_2b_cutFlow.AddBinContent(2, presel_weight*weightMET)
                if (ep_pfMetCorrPt > 200.):
                   h_reg_SR_1b_cutFlow.AddBinContent(3, presel_weight*weightMET)
                   h_reg_SR_2b_cutFlow.AddBinContent(3, presel_weight*weightMET)
                   if (ep_nEle_index == 0) and (ep_nMu == 0) and (nPho ==0) and (ep_nTau_discBased_TightEleTightMuVeto==0):
                       h_reg_SR_1b_cutFlow.AddBinContent(4, presel_weight*weightMET)
                       h_reg_SR_2b_cutFlow.AddBinContent(4, presel_weight*weightMET)
                       if (min_dPhi_jet_MET > 0.5):
                           h_reg_SR_1b_cutFlow.AddBinContent(5, presel_weight*weightMET)
                           h_reg_SR_2b_cutFlow.AddBinContent(5, presel_weight*weightMET)
                           if (ep_THINnJet <= 2) and (ep_THINjetPt[0] > 50.) :
                               h_reg_SR_1b_cutFlow.AddBinContent(6, presel_weight*weightMET)
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (nBjets==1):
                                   h_reg_SR_1b_cutFlow.AddBinContent(7, weight)
                                   isSR1b=True
                                   SR1bcount+=1
                                   if ep_THINnJet==2:
                                       Jet2Pt  = ep_THINjetPt[1]; Jet2Eta     = ep_THINjetEta[1]
                                       Jet2Phi = ep_THINjetPhi[1];Jet2deepCSV = ep_THINjetDeepCSV[1]
                           if (ep_THINnJet <=3 and ep_THINnJet > 1) and (ep_THINjetPt[0] > 50.) :
                               h_reg_SR_2b_cutFlow.AddBinContent(6, presel_weight*weightMET)
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and (nBjets==2):
                                   h_reg_SR_2b_cutFlow.AddBinContent(7, weight)
                                   isSR2b=True
                                   SR2bcount+=1
                                   if ep_THINnJet==3:
                                       Jet3Pt  = ep_THINjetPt[2]; Jet3Eta     = ep_THINjetEta[2]
                                       Jet3Phi = ep_THINjetPhi[2];Jet3deepCSV = ep_THINjetDeepCSV[2]
            '''
            --------------------------------------------------------------------------------
            ZEE CONTROL REGION
            --------------------------------------------------------------------------------
            '''
            h_reg_ZeeCR_1b_cutFlow.AddBinContent(1, presel_weight)
            h_reg_ZeeCR_2b_cutFlow.AddBinContent(1, presel_weight)
            if eletrigdecision:
                h_reg_ZeeCR_1b_cutFlow.AddBinContent(2, presel_weight*weightEleTrig)
                h_reg_ZeeCR_2b_cutFlow.AddBinContent(2, presel_weight*weightEleTrig)
                if (ep_pfMetCorrPt > 0.):
                    h_reg_ZeeCR_1b_cutFlow.AddBinContent(3, presel_weight*weightEleTrig)
                    h_reg_ZeeCR_2b_cutFlow.AddBinContent(3, presel_weight*weightEleTrig)
                    if (ep_nEle_index == 2) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_elePt[0] > minElePt) and (ep_eleIsPTight[0]) and nPho ==0:
                       h_reg_ZeeCR_1b_cutFlow.AddBinContent(4, presel_weight*weightEleTrig*weightEle)
                       h_reg_ZeeCR_2b_cutFlow.AddBinContent(4, presel_weight*weightEleTrig*weightEle)
                       if (ep_ZeeRecoil > 200.):
                           h_reg_ZeeCR_1b_cutFlow.AddBinContent(5, presel_weight*weightEleTrig*weightEle)
                           h_reg_ZeeCR_2b_cutFlow.AddBinContent(5, presel_weight*weightEleTrig*weightEle)
                           if (min_dPhi_jet_MET > 0.5):
                               h_reg_ZeeCR_1b_cutFlow.AddBinContent(6, presel_weight*weightEleTrig*weightEle)
                               h_reg_ZeeCR_2b_cutFlow.AddBinContent(6, presel_weight*weightEleTrig*weightEle)
                               if (ep_Zeemass >= 70 and ep_Zeemass <= 110):
                                   h_reg_ZeeCR_1b_cutFlow.AddBinContent(7, presel_weight*weightEleTrig*weightEle)
                                   h_reg_ZeeCR_2b_cutFlow.AddBinContent(7, presel_weight*weightEleTrig*weightEle)
                                   if (ep_THINnJet <= 2) and (ep_THINjetPt[0] > 50.):
                                       h_reg_ZeeCR_1b_cutFlow.AddBinContent(8, presel_weight*weightEleTrig*weightEle)
                                       if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (nBjets==1):
                                           h_reg_ZeeCR_1b_cutFlow.AddBinContent(9, weight)
                                           ZeeCR1bcount+=1
                                           is1bCRZee=True
                                   if (ep_THINnJet <= 3 and ep_THINnJet > 1) and (ep_THINjetPt[0] > 50.):
                                       h_reg_ZeeCR_2b_cutFlow.AddBinContent(8, presel_weight*weightEleTrig*weightEle)
                                       if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and (nBjets==2):
                                           h_reg_ZeeCR_2b_cutFlow.AddBinContent(9, weight)
                                           ZeeCR2bcount+=1
                                           is2bCRZee=True

            '''
            --------------------------------------------------------------------------------
            ZMUMU CONTROL REGION
            --------------------------------------------------------------------------------
            '''
            h_reg_ZmumuCR_1b_cutFlow.AddBinContent(1, presel_weight)
            h_reg_ZmumuCR_2b_cutFlow.AddBinContent(1, presel_weight)
            if mettrigdecision:
                h_reg_ZmumuCR_1b_cutFlow.AddBinContent(2, presel_weight*weightRecoil)
                h_reg_ZmumuCR_2b_cutFlow.AddBinContent(2, presel_weight*weightRecoil)
                if (ep_pfMetCorrPt > 0.):
                    h_reg_ZmumuCR_1b_cutFlow.AddBinContent(3, presel_weight*weightRecoil)
                    h_reg_ZmumuCR_2b_cutFlow.AddBinContent(3, presel_weight*weightRecoil)
                    if (ep_nEle_index == 0) and (ep_nMu == 2) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_muPt[0] > minMuPt) and (ep_isTightMuon[0]) and nPho ==0:
                        h_reg_ZmumuCR_1b_cutFlow.AddBinContent(4, presel_weight*weightRecoil*weightMu)
                        h_reg_ZmumuCR_2b_cutFlow.AddBinContent(4, presel_weight*weightRecoil*weightMu)
                        if (ep_ZmumuRecoil > 200. ) :
                            h_reg_ZmumuCR_1b_cutFlow.AddBinContent(5, presel_weight*weightRecoil*weightMu)
                            h_reg_ZmumuCR_2b_cutFlow.AddBinContent(5, presel_weight*weightRecoil*weightMu)
                            if (min_dPhi_jet_MET > 0.5):
                                h_reg_ZmumuCR_1b_cutFlow.AddBinContent(6, presel_weight*weightRecoil*weightMu)
                                h_reg_ZmumuCR_2b_cutFlow.AddBinContent(6, presel_weight*weightRecoil*weightMu)
                                if (ep_Zmumumass >= 70 and ep_Zmumumass <= 110):
                                    h_reg_ZmumuCR_1b_cutFlow.AddBinContent(7, presel_weight*weightRecoil*weightMu)
                                    h_reg_ZmumuCR_2b_cutFlow.AddBinContent(7, presel_weight*weightRecoil*weightMu)
                                    if (ep_THINnJet <=2) and (ep_THINjetPt[0] > 50.):
                                        h_reg_ZmumuCR_1b_cutFlow.AddBinContent(8, presel_weight*weightRecoil*weightMu)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (nBjets==1):
                                            h_reg_ZmumuCR_1b_cutFlow.AddBinContent(9, weight)
                                            ZmumuCR1count+=1
                                            is1bCRZmumu=True
                                    if (ep_THINnJet <= 3 and ep_THINnJet > 1) and (ep_THINjetPt[0] > 50.):
                                        h_reg_ZmumuCR_2b_cutFlow.AddBinContent(8, presel_weight*weightRecoil*weightMu)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and (nBjets==2):
                                            h_reg_ZmumuCR_2b_cutFlow.AddBinContent(9, weight)
                                            ZmumuCR2count+=1
                                            is2bCRZmumu=True
            '''
            --------------------------------------------------------------------------------
            WENU CONTROL REGION
            --------------------------------------------------------------------------------
            '''
            h_reg_WenuCR_1b_cutFlow.AddBinContent(1, presel_weight)
            h_reg_WenuCR_2b_cutFlow.AddBinContent(1, presel_weight)
            if eletrigdecision:
                h_reg_WenuCR_1b_cutFlow.AddBinContent(2, presel_weight*weightEleTrig)
                h_reg_WenuCR_2b_cutFlow.AddBinContent(2, presel_weight*weightEleTrig)
                if (ep_pfMetCorrPt > 0.):
                    h_reg_WenuCR_1b_cutFlow.AddBinContent(3, presel_weight*weightEleTrig*weightEle)
                    h_reg_WenuCR_2b_cutFlow.AddBinContent(3, presel_weight*weightEleTrig*weightEle)
                    if (ep_nEle_index == 1) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_elePt[0] > minElePt) and (ep_eleIsPTight[0]) and nPho ==0:
                        h_reg_WenuCR_1b_cutFlow.AddBinContent(4, presel_weight*weightEleTrig*weightEle)
                        h_reg_WenuCR_2b_cutFlow.AddBinContent(4, presel_weight*weightEleTrig*weightEle)
                        if (ep_WenuRecoil > 200.) :
                            h_reg_WenuCR_1b_cutFlow.AddBinContent(5, presel_weight*weightEleTrig*weightEle)
                            h_reg_WenuCR_2b_cutFlow.AddBinContent(5, presel_weight*weightEleTrig*weightEle)
                            if (min_dPhi_jet_MET > 0.5):
                                h_reg_WenuCR_1b_cutFlow.AddBinContent(6, presel_weight*weightEleTrig*weightEle)
                                h_reg_WenuCR_2b_cutFlow.AddBinContent(6, presel_weight*weightEleTrig*weightEle)
                                if (ep_Wenumass >= 0 and ep_Wenumass <= 160):
                                    h_reg_WenuCR_1b_cutFlow.AddBinContent(7, presel_weight*weightEleTrig*weightEle)
                                    h_reg_WenuCR_2b_cutFlow.AddBinContent(7, presel_weight*weightEleTrig*weightEle)
                                    if (ep_THINnJet ==1) and (ep_THINjetPt[0] > 50.):
                                        h_reg_WenuCR_1b_cutFlow.AddBinContent(8, presel_weight*weightEleTrig*weightEle)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (nBjets==1):
                                            h_reg_WenuCR_1b_cutFlow.AddBinContent(9, weight)
                                            WenuCR1bcount+=1
                                            is1bCRWenu=True
                                    if (ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.):
                                        h_reg_WenuCR_2b_cutFlow.AddBinContent(8, weight)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and (nBjets==2):
                                            h_reg_WenuCR_2b_cutFlow.AddBinContent(9, weight)
                                            WenuCR2bcount+=1
                                            is2bCRWenu=True
            '''
            --------------------------------------------------------------------------------
            WMUNU CONTROL REGION
            --------------------------------------------------------------------------------
            '''
            h_reg_WmunuCR_1b_cutFlow.AddBinContent(1, presel_weight)
            h_reg_WmunuCR_2b_cutFlow.AddBinContent(1, presel_weight*weightRecoil)
            if mettrigdecision:
                h_reg_WmunuCR_1b_cutFlow.AddBinContent(2, presel_weight*weightRecoil)
                h_reg_WmunuCR_2b_cutFlow.AddBinContent(2, presel_weight*weightRecoil)
                if (ep_pfMetCorrPt > 0.):
                    h_reg_WmunuCR_1b_cutFlow.AddBinContent(3, presel_weight*weightRecoil)
                    h_reg_WmunuCR_2b_cutFlow.AddBinContent(3, presel_weight*weightRecoil)
                    if (ep_nEle_index == 0) and (ep_nMu == 1) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_muPt[0] > minMuPt) and (ep_isTightMuon[0]) and nPho ==0:
                        h_reg_WmunuCR_1b_cutFlow.AddBinContent(4, presel_weight*weightRecoil*weightMu)
                        h_reg_WmunuCR_2b_cutFlow.AddBinContent(4, presel_weight*weightRecoil*weightMu)
                        if (ep_WmunuRecoil > 200.) :
                            h_reg_WmunuCR_1b_cutFlow.AddBinContent(5, presel_weight*weightRecoil*weightMu)
                            h_reg_WmunuCR_2b_cutFlow.AddBinContent(5, presel_weight*weightRecoil*weightMu)
                            if (min_dPhi_jet_MET > 0.5):
                                h_reg_WmunuCR_1b_cutFlow.AddBinContent(6, presel_weight*weightRecoil*weightMu)
                                h_reg_WmunuCR_2b_cutFlow.AddBinContent(6, presel_weight*weightRecoil*weightMu)
                                if (ep_Wmunumass >= 0 and ep_Wmunumass <= 160):
                                    h_reg_WmunuCR_1b_cutFlow.AddBinContent(7, presel_weight*weightRecoil*weightMu)
                                    h_reg_WmunuCR_2b_cutFlow.AddBinContent(7, presel_weight*weightRecoil*weightMu)
                                    if (ep_THINnJet ==1) and (ep_THINjetPt[0] > 50.):
                                        h_reg_WmunuCR_1b_cutFlow.AddBinContent(8, presel_weight*weightRecoil*weightMu)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (nBjets==1):
                                            h_reg_WmunuCR_1b_cutFlow.AddBinContent(9, weight)
                                            WmunuCR1bcount+=1
                                            is1bCRWmunu=True
                                    if (ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.):
                                        h_reg_WmunuCR_2b_cutFlow.AddBinContent(8, weight)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and (nBjets==2):
                                            h_reg_WmunuCR_2b_cutFlow.AddBinContent(9, weight)
                                            WmunuCR2bcount+=1
                                            is2bCRWmunu=True
            '''
            --------------------------------------------------------------------------------
            TOPENU CONTROL REGION
            --------------------------------------------------------------------------------
            '''
            h_reg_TopenuCR_1b_cutFlow.AddBinContent(1, presel_weight)
            h_reg_TopenuCR_2b_cutFlow.AddBinContent(1, presel_weight)
            if eletrigdecision:
                h_reg_TopenuCR_1b_cutFlow.AddBinContent(2, presel_weight*weightEleTrig)
                h_reg_TopenuCR_2b_cutFlow.AddBinContent(2, presel_weight*weightEleTrig)
                if (ep_pfMetCorrPt > 0.):
                    h_reg_TopenuCR_1b_cutFlow.AddBinContent(3, presel_weight*weightEleTrig*weightEle)
                    h_reg_TopenuCR_2b_cutFlow.AddBinContent(3, presel_weight*weightEleTrig*weightEle)
                    if (ep_nEle_index == 1) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_elePt[0] > minElePt) and (ep_eleIsPTight[0]) and nPho ==0:
                        h_reg_TopenuCR_1b_cutFlow.AddBinContent(4, presel_weight*weightEleTrig*weightEle)
                        h_reg_TopenuCR_2b_cutFlow.AddBinContent(4, presel_weight*weightEleTrig*weightEle)
                        if (ep_WenuRecoil > 200. ) :
                            h_reg_TopenuCR_1b_cutFlow.AddBinContent(5, presel_weight*weightEleTrig*weightEle)
                            h_reg_TopenuCR_2b_cutFlow.AddBinContent(5, presel_weight*weightEleTrig*weightEle)
                            if (min_dPhi_jet_MET > 0.5):
                                h_reg_TopenuCR_1b_cutFlow.AddBinContent(6, presel_weight*weightEleTrig*weightEle)
                                h_reg_TopenuCR_2b_cutFlow.AddBinContent(6, presel_weight*weightEleTrig*weightEle)
                                if (ep_Wenumass >= 0 and ep_Wenumass <= 160):
                                    h_reg_TopenuCR_1b_cutFlow.AddBinContent(7, presel_weight*weightEleTrig*weightEle)
                                    h_reg_TopenuCR_2b_cutFlow.AddBinContent(7, presel_weight*weightEleTrig*weightEle)
                                    if (ep_THINnJet >1) and (ep_THINjetPt[0] > 50.):
                                        h_reg_TopenuCR_1b_cutFlow.AddBinContent(8, presel_weight*weightEleTrig*weightEle)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (nBjets==1):
                                            h_reg_TopenuCR_1b_cutFlow.AddBinContent(9, weight)
                                            TopenuCR1bcount+=1
                                            is1bCRTopenu=True
                                    if (ep_THINnJet > 2) and (ep_THINjetPt[0] > 50.):
                                        h_reg_TopenuCR_2b_cutFlow.AddBinContent(8, weight)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and (nBjets==2):
                                            h_reg_TopenuCR_2b_cutFlow.AddBinContent(9, weight)
                                            TopenuCR2bcount+=1
                                            is2bCRTopenu=True
            '''
            --------------------------------------------------------------------------------
            TOPMUNU CONTROL REGION
            --------------------------------------------------------------------------------
            '''
            h_reg_TopmunuCR_1b_cutFlow.AddBinContent(1, presel_weight)
            h_reg_TopmunuCR_2b_cutFlow.AddBinContent(1, presel_weight)
            if mettrigdecision:
                h_reg_TopmunuCR_1b_cutFlow.AddBinContent(2, presel_weight*weightRecoil)
                h_reg_TopmunuCR_2b_cutFlow.AddBinContent(2, presel_weight*weightRecoil)
                if (ep_pfMetCorrPt > 0.):
                    h_reg_TopmunuCR_1b_cutFlow.AddBinContent(3, presel_weight*weightRecoil*weightMu)
                    h_reg_TopmunuCR_2b_cutFlow.AddBinContent(3, presel_weight*weightRecoil*weightMu)
                    if (ep_nEle_index == 0) and (ep_nMu == 1) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_muPt[0] > minMuPt) and (ep_isTightMuon[0]) and nPho ==0:
                        h_reg_TopmunuCR_1b_cutFlow.AddBinContent(4, presel_weight*weightRecoil*weightMu)
                        h_reg_TopmunuCR_2b_cutFlow.AddBinContent(4, presel_weight*weightRecoil*weightMu)
                        if (ep_WmunuRecoil > 200. ) :
                            h_reg_TopmunuCR_1b_cutFlow.AddBinContent(5, presel_weight*weightRecoil*weightMu)
                            h_reg_TopmunuCR_2b_cutFlow.AddBinContent(5, presel_weight*weightRecoil*weightMu)
                            if (min_dPhi_jet_MET > 0.5):
                                h_reg_TopmunuCR_1b_cutFlow.AddBinContent(6, presel_weight*weightRecoil*weightMu)
                                h_reg_TopmunuCR_2b_cutFlow.AddBinContent(6, presel_weight*weightRecoil*weightMu)
                                if (ep_Wmunumass >= 0 and ep_Wmunumass <= 160):
                                    h_reg_TopmunuCR_1b_cutFlow.AddBinContent(7, presel_weight*weightRecoil*weightMu)
                                    h_reg_TopmunuCR_2b_cutFlow.AddBinContent(7, presel_weight*weightRecoil*weightMu)
                                    if (ep_THINnJet >1) and (ep_THINjetPt[0] > 50.) :
                                        h_reg_TopmunuCR_1b_cutFlow.AddBinContent(8, presel_weight*weightRecoil*weightMu)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (nBjets==1):
                                            h_reg_TopmunuCR_1b_cutFlow.AddBinContent(9, weight)
                                            TopmunuCR1bcount+=1
                                            is1bCRTopmunu=True
                                    if (ep_THINnJet >2) and (ep_THINjetPt[0] > 50.) :
                                        h_reg_TopmunuCR_2b_cutFlow.AddBinContent(8, presel_weight*weightRecoil*weightMu)
                                        if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and (nBjets==2):
                                            h_reg_TopmunuCR_2b_cutFlow.AddBinContent(9, weight)
                                            TopmunuCR2bcount+=1
                                            is2bCRTopmunu=True

            if isSR1b:
                df_out_SR_1b = df_out_SR_1b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]),'Jet1Eta':float(ep_THINjetEta[0]),'Jet1Phi':float(ep_THINjetPhi[0]),'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(Jet2Pt),'Jet2Eta':float(Jet2Eta),'Jet2Phi':float(Jet2Phi),'Jet2deepCSV':float(Jet2deepCSV),
                                                    'Jet3Pt':float(dummy),'Jet3Eta':float(dummy),'Jet3Phi':float(dummy),'Jet3deepCSV':float(dummy),
                                                    'weight':float(weight),'weightMET':float(weightMET),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightMET_up':float(weightMET_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'MET_Res_up':float(ep_pfMetUncJetResUp),'MET_En_up':float(ep_pfMetUncJetEnUp),'MET_En_down':float(ep_pfMetUncJetEnDown),'MET_Res_down':float(ep_pfMetUncJetResDown),'weightJEC_down':float(weightJEC_down),'weightMET_down':float(weightMET_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('isSR1b')
            if isSR2b:
                df_out_SR_2b = df_out_SR_2b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]), 'Jet1Eta':float(ep_THINjetEta[0]), 'Jet1Phi':float(ep_THINjetPhi[0]), 'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(ep_THINjetPt[1]), 'Jet2Eta':float(ep_THINjetEta[1]), 'Jet2Phi':float(ep_THINjetPhi[1]), 'Jet2deepCSV':float(ep_THINjetDeepCSV[1]),
                                                    'Jet3Pt':float(Jet3Pt), 'Jet3Eta':float(Jet3Eta), 'Jet3Phi':float(Jet3Phi), 'Jet3deepCSV':float(Jet3deepCSV),
                                                    'weight':float(weight),'weightMET':float(weightMET),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightMET_up':float(weightMET_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'MET_Res_up':float(ep_pfMetUncJetResUp),'MET_En_up':float(ep_pfMetUncJetEnUp),'MET_En_down':float(ep_pfMetUncJetEnDown),'MET_Res_down':float(ep_pfMetUncJetResDown),'weightJEC_down':float(weightJEC_down),'weightMET_down':float(weightMET_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('isSR2b')

            if is1bCRZee:
                df_out_ZeeCR_1b = df_out_ZeeCR_1b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_ZeeRecoil ),'RecoilPhi':float(ep_ZeeRecoil_dPhi ),'Zmass':float(ep_Zeemass),'ZpT':float(ZpT_ee),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]),'Jet1Eta':float(ep_THINjetEta[0]),'Jet1Phi':float(ep_THINjetPhi[0]),'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(Jet2Pt),'Jet2Eta':float(Jet2Eta),'Jet2Phi':float(Jet2Phi),'Jet2deepCSV':float(Jet2deepCSV),
                                                    'Jet3Pt':float(dummy),'Jet3Eta':float(dummy),'Jet3Phi':float(dummy),'Jet3deepCSV':float(dummy),
                                                    'leadingLepPt':float(ep_elePt[0]),'leadingLepEta':float(ep_eleEta[0]),'leadingLepPhi':float(ep_elePhi[0]),
                                                    'subleadingLepPt':float(ep_elePt[1]),'subleadingLepEta':float(ep_eleEta[1]),'subleadingLepPhi':float(ep_elePhi[1]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_ZeeRecoilResUp),'Recoil_En_up':float(ep_ZeeRecoilEnUp),'Recoil_En_down':float(ep_ZeeRecoilEnDown),'Recoil_Res_down':float(ep_ZeeRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRZee')
            if is2bCRZee:
                df_out_ZeeCR_2b = df_out_ZeeCR_2b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_ZeeRecoil ),'RecoilPhi':float(ep_ZeeRecoil_dPhi ),'Zmass':float(ep_Zeemass),'ZpT':float(ZpT_ee),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]), 'Jet1Eta':float(ep_THINjetEta[0]), 'Jet1Phi':float(ep_THINjetPhi[0]), 'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(ep_THINjetPt[1]), 'Jet2Eta':float(ep_THINjetEta[1]), 'Jet2Phi':float(ep_THINjetPhi[1]), 'Jet2deepCSV':float(ep_THINjetDeepCSV[1]),
                                                    'Jet3Pt':float(Jet3Pt), 'Jet3Eta':float(Jet3Eta), 'Jet3Phi':float(Jet3Phi), 'Jet3deepCSV':float(Jet3deepCSV),
                                                    'leadingLepPt':float(ep_elePt[0]),'leadingLepEta':float(ep_eleEta[0]),'leadingLepPhi':float(ep_elePhi[0]),
                                                    'subleadingLepPt':float(ep_elePt[1]),'subleadingLepEta':float(ep_eleEta[1]),'subleadingLepPhi':float(ep_elePhi[1]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_ZeeRecoilResUp),'Recoil_En_up':float(ep_ZeeRecoilEnUp),'Recoil_En_down':float(ep_ZeeRecoilEnDown),'Recoil_Res_down':float(ep_ZeeRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRZee')

            if is1bCRZmumu:
                df_out_ZmumuCR_1b = df_out_ZmumuCR_1b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_ZmumuRecoil ),'RecoilPhi':float(ep_ZmumuRecoil_dPhi ),'Zmass':float(ep_Zmumumass),'ZpT':float(ZpT_mumu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]),'Jet1Eta':float(ep_THINjetEta[0]),'Jet1Phi':float(ep_THINjetPhi[0]),'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(Jet2Pt),'Jet2Eta':float(Jet2Eta),'Jet2Phi':float(Jet2Phi),'Jet2deepCSV':float(Jet2deepCSV),
                                                    'Jet3Pt':float(dummy),'Jet3Eta':float(dummy),'Jet3Phi':float(dummy),'Jet3deepCSV':float(dummy),
                                                    'leadingLepPt':float(ep_muPt[0]),'leadingLepEta':float(ep_muEta[0]),'leadingLepPhi':float(ep_muPhi[0]),
                                                    'subleadingLepPt':float(ep_muPt[1]),'subleadingLepEta':float(ep_muEta[1]),'subleadingLepPhi':float(ep_muPhi[1]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_ZmumuRecoilResUp),'Recoil_En_up':float(ep_ZmumuRecoilEnUp),'Recoil_En_down':float(ep_ZmumuRecoilEnDown),'Recoil_Res_down':float(ep_ZmumuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRZmumu')
            if is2bCRZmumu:
                df_out_ZmumuCR_2b = df_out_ZmumuCR_2b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_ZmumuRecoil ),'RecoilPhi':float(ep_ZmumuRecoil_dPhi ),'Zmass':float(ep_Zmumumass),'ZpT':float(ZpT_mumu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]), 'Jet1Eta':float(ep_THINjetEta[0]), 'Jet1Phi':float(ep_THINjetPhi[0]), 'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(ep_THINjetPt[1]), 'Jet2Eta':float(ep_THINjetEta[1]), 'Jet2Phi':float(ep_THINjetPhi[1]), 'Jet2deepCSV':float(ep_THINjetDeepCSV[1]),
                                                    'Jet3Pt':float(Jet3Pt), 'Jet3Eta':float(Jet3Eta), 'Jet3Phi':float(Jet3Phi), 'Jet3deepCSV':float(Jet3deepCSV),
                                                    'leadingLepPt':float(ep_muPt[0]),'leadingLepEta':float(ep_muEta[0]),'leadingLepPhi':float(ep_muPhi[0]),
                                                    'subleadingLepPt':float(ep_muPt[1]),'subleadingLepEta':float(ep_muEta[1]),'subleadingLepPhi':float(ep_muPhi[1]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_ZmumuRecoilResUp),'Recoil_En_up':float(ep_ZmumuRecoilEnUp),'Recoil_En_down':float(ep_ZmumuRecoilEnDown),'Recoil_Res_down':float(ep_ZmumuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRZmumu')
            if is1bCRWenu:
                df_out_WenuCR_1b = df_out_WenuCR_1b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_WenuRecoil ),'RecoilPhi':float(ep_WenuRecoil_dPhi ),'Wmass':float(ep_Wenumass),'WpT':float(WpT_enu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]),'Jet1Eta':float(ep_THINjetEta[0]),'Jet1Phi':float(ep_THINjetPhi[0]),'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(Jet2Pt),'Jet2Eta':float(Jet2Eta),'Jet2Phi':float(Jet2Phi),'Jet2deepCSV':float(Jet2deepCSV),
                                                    'Jet3Pt':float(dummy),'Jet3Eta':float(dummy),'Jet3Phi':float(dummy),'Jet3deepCSV':float(dummy),
                                                    'leadingLepPt':float(ep_elePt[0]),'leadingLepEta':float(ep_eleEta[0]),'leadingLepPhi':float(ep_elePhi[0]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_WenuRecoilResUp),'Recoil_En_up':float(ep_WenuRecoilEnUp),'Recoil_En_down':float(ep_WenuRecoilEnDown),'Recoil_Res_down':float(ep_WenuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRWenu')
            if is2bCRWenu:
                df_out_WenuCR_2b = df_out_WenuCR_2b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_WenuRecoil ),'RecoilPhi':float(ep_WenuRecoil_dPhi ),'Wmass':float(ep_Wenumass),'WpT':float(WpT_enu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]), 'Jet1Eta':float(ep_THINjetEta[0]), 'Jet1Phi':float(ep_THINjetPhi[0]), 'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(ep_THINjetPt[1]), 'Jet2Eta':float(ep_THINjetEta[1]), 'Jet2Phi':float(ep_THINjetPhi[1]), 'Jet2deepCSV':float(ep_THINjetDeepCSV[1]),
                                                    'Jet3Pt':float(Jet3Pt), 'Jet3Eta':float(Jet3Eta), 'Jet3Phi':float(Jet3Phi), 'Jet3deepCSV':float(Jet3deepCSV),
                                                    'leadingLepPt':float(ep_elePt[0]),'leadingLepEta':float(ep_eleEta[0]),'leadingLepPhi':float(ep_elePhi[0]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_WenuRecoilResUp),'Recoil_En_up':float(ep_WenuRecoilEnUp),'Recoil_En_down':float(ep_WenuRecoilEnDown),'Recoil_Res_down':float(ep_WenuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRWenu')

            if is1bCRWmunu:
                df_out_WmunuCR_1b = df_out_WmunuCR_1b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_WmunuRecoil ),'RecoilPhi':float(ep_WmunuRecoil_dPhi ),'Wmass':float(ep_Wmunumass),'WpT':float(WpT_munu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]),'Jet1Eta':float(ep_THINjetEta[0]),'Jet1Phi':float(ep_THINjetPhi[0]),'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(Jet2Pt),'Jet2Eta':float(Jet2Eta),'Jet2Phi':float(Jet2Phi),'Jet2deepCSV':float(Jet2deepCSV),
                                                    'Jet3Pt':float(dummy),'Jet3Eta':float(dummy),'Jet3Phi':float(dummy),'Jet3deepCSV':float(dummy),
                                                    'leadingLepPt':float(ep_muPt[0]),'leadingLepEta':float(ep_muEta[0]),'leadingLepPhi':float(ep_muPhi[0]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_WmunuRecoilResUp),'Recoil_En_up':float(ep_WmunuRecoilEnUp),'Recoil_En_down':float(ep_WmunuRecoilEnDown),'Recoil_Res_down':float(ep_WmunuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRWmunu')
            if is2bCRWmunu:
                df_out_WmunuCR_2b = df_out_WmunuCR_2b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_WmunuRecoil ),'RecoilPhi':float(ep_WmunuRecoil_dPhi ),'Wmass':float(ep_Wmunumass),'WpT':float(WpT_munu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]), 'Jet1Eta':float(ep_THINjetEta[0]), 'Jet1Phi':float(ep_THINjetPhi[0]), 'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(ep_THINjetPt[1]), 'Jet2Eta':float(ep_THINjetEta[1]), 'Jet2Phi':float(ep_THINjetPhi[1]), 'Jet2deepCSV':float(ep_THINjetDeepCSV[1]),
                                                    'Jet3Pt':float(Jet3Pt), 'Jet3Eta':float(Jet3Eta), 'Jet3Phi':float(Jet3Phi), 'Jet3deepCSV':float(Jet3deepCSV),
                                                    'leadingLepPt':float(ep_muPt[0]),'leadingLepEta':float(ep_muEta[0]),'leadingLepPhi':float(ep_muPhi[0]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_WmunuRecoilResUp),'Recoil_En_up':float(ep_WmunuRecoilEnUp),'Recoil_En_down':float(ep_WmunuRecoilEnDown),'Recoil_Res_down':float(ep_WmunuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRWmunu')
            if is1bCRTopenu:
                df_out_TopenuCR_1b = df_out_TopenuCR_1b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_WenuRecoil ),'RecoilPhi':float(ep_WenuRecoil_dPhi ),'Wmass':float(ep_Wenumass),'WpT':float(WpT_enu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]),'Jet1Eta':float(ep_THINjetEta[0]),'Jet1Phi':float(ep_THINjetPhi[0]),'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(Jet2Pt),'Jet2Eta':float(Jet2Eta),'Jet2Phi':float(Jet2Phi),'Jet2deepCSV':float(Jet2deepCSV),
                                                    'Jet3Pt':float(dummy),'Jet3Eta':float(dummy),'Jet3Phi':float(dummy),'Jet3deepCSV':float(dummy),
                                                    'leadingLepPt':float(ep_elePt[0]),'leadingLepEta':float(ep_eleEta[0]),'leadingLepPhi':float(ep_elePhi[0]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_WenuRecoilResUp),'Recoil_En_up':float(ep_WenuRecoilEnUp),'Recoil_En_down':float(ep_WenuRecoilEnDown),'Recoil_Res_down':float(ep_WenuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRTopenu')
            if is2bCRTopenu:
                df_out_TopenuCR_2b = df_out_TopenuCR_2b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_WenuRecoil ),'RecoilPhi':float(ep_WenuRecoil_dPhi ),'Wmass':float(ep_Wenumass),'WpT':float(WpT_enu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]), 'Jet1Eta':float(ep_THINjetEta[0]), 'Jet1Phi':float(ep_THINjetPhi[0]), 'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(ep_THINjetPt[1]), 'Jet2Eta':float(ep_THINjetEta[1]), 'Jet2Phi':float(ep_THINjetPhi[1]), 'Jet2deepCSV':float(ep_THINjetDeepCSV[1]),
                                                    'Jet3Pt':float(Jet3Pt), 'Jet3Eta':float(Jet3Eta), 'Jet3Phi':float(Jet3Phi), 'Jet3deepCSV':float(Jet3deepCSV),
                                                    'leadingLepPt':float(ep_elePt[0]),'leadingLepEta':float(ep_eleEta[0]),'leadingLepPhi':float(ep_elePhi[0]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_WenuRecoilResUp),'Recoil_En_up':float(ep_WenuRecoilEnUp),'Recoil_En_down':float(ep_WenuRecoilEnDown),'Recoil_Res_down':float(ep_WenuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRTopenu')
            if is1bCRTopmunu:
                df_out_TopmunuCR_1b = df_out_TopmunuCR_1b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_WmunuRecoil ),'RecoilPhi':float(ep_WmunuRecoil_dPhi ),'Wmass':float(ep_Wmunumass),'WpT':float(WpT_munu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]),'Jet1Eta':float(ep_THINjetEta[0]),'Jet1Phi':float(ep_THINjetPhi[0]),'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(Jet2Pt),'Jet2Eta':float(Jet2Eta),'Jet2Phi':float(Jet2Phi),'Jet2deepCSV':float(Jet2deepCSV),
                                                    'Jet3Pt':float(dummy),'Jet3Eta':float(dummy),'Jet3Phi':float(dummy),'Jet3deepCSV':float(dummy),
                                                    'leadingLepPt':float(ep_muPt[0]),'leadingLepEta':float(ep_muEta[0]),'leadingLepPhi':float(ep_muPhi[0]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_WmunuRecoilResUp),'Recoil_En_up':float(ep_WmunuRecoilEnUp),'Recoil_En_down':float(ep_WmunuRecoilEnDown),'Recoil_Res_down':float(ep_WmunuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRTopmunu')
            if is2bCRTopmunu:
                df_out_TopmunuCR_2b = df_out_TopmunuCR_2b.append({'run':float(ep_runId), 'lumi':float(ep_lumiSection), 'event':float(ep_eventId),'nPV':float(ep_THINjetNPV),
                                                    'MET':float(ep_pfMetCorrPt),'METPhi':float(ep_pfMetCorrPhi),'Recoil':float(ep_WmunuRecoil ),'RecoilPhi':float(ep_WmunuRecoil_dPhi ),'Wmass':float(ep_Wmunumass),'WpT':float(WpT_munu),
                                                    'dPhi_jetMET':float(min_dPhi_jet_MET),
                                                    'NTau':float(ep_nTau_discBased_TightEleTightMuVeto),'NEle':float(ep_nEle_index),'NMu':float(ep_nMu), 'nPho':float(nPho),
                                                    'Njets_PassID':float(ep_THINnJet),'Nbjets_PassID':float(nBjets),
                                                    'Jet1Pt':float(ep_THINjetPt[0]), 'Jet1Eta':float(ep_THINjetEta[0]), 'Jet1Phi':float(ep_THINjetPhi[0]), 'Jet1deepCSV':float(ep_THINjetDeepCSV[0]),
                                                    'Jet2Pt':float(ep_THINjetPt[1]), 'Jet2Eta':float(ep_THINjetEta[1]), 'Jet2Phi':float(ep_THINjetPhi[1]), 'Jet2deepCSV':float(ep_THINjetDeepCSV[1]),
                                                    'Jet3Pt':float(Jet3Pt), 'Jet3Eta':float(Jet3Eta), 'Jet3Phi':float(Jet3Phi), 'Jet3deepCSV':float(Jet3deepCSV),
                                                    'leadingLepPt':float(ep_muPt[0]),'leadingLepEta':float(ep_muEta[0]),'leadingLepPhi':float(ep_muPhi[0]),
                                                    'weight':float(weight),'weightRecoil':float(weightRecoil),'weightEle':float(weightEle),'weightMu':float(weightMu),'weightB':float(weightB),'weightEWK':float(weightEWK),'weightQCD':float(weightQCD),'weightTop':float(weightTop),'weightPU':float(weightPU),'weightPrefire':float(weightPrefire),'weightRecoil_up':float(weightRecoil_up),'weightEle_up':float(weightEle_up),'weightMu_up':float(weightMu_up),'weightB_up':float(weightB_up),'weightEWK_up':float(weightEWK_up),'weightQCD_up':float(weightQCD_up),'weightTop_up':float(weightTop_up),'weightPU_up':float(weightPU_up),'weightPrefire_up':float(weightPrefire_up),'weightJEC_up':float(weightJEC_up),'Recoil_Res_up':float(ep_WmunuRecoilResUp),'Recoil_En_up':float(ep_WmunuRecoilEnUp),'Recoil_En_down':float(ep_WmunuRecoilEnDown),'Recoil_Res_down':float(ep_WmunuRecoilResDown),'weightJEC_down':float(weightJEC_down),'weightRecoil_down':float(weightRecoil_down),'weightEle_down':float(weightEle_down),'weightMu_down':float(weightMu_down),'weightB_down':float(weightB_down),'weightEWK_down':float(weightEWK_down),'weightQCD_down':float(weightQCD_down),'weightTop_down':float(weightTop_down),'weightPU_down':float(weightPU_down),'weightPrefire_down':float(weightPrefire_down)
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRTopmunu')

    outfilenameis=outfilename
    df_out_SR_1b.to_root(outfilenameis, key='bbDM_SR_1b',mode='w')
    df_out_SR_2b.to_root(outfilenameis, key='bbDM_SR_2b',mode='a')

    df_out_ZeeCR_1b.to_root(outfilenameis, key='bbDM_ZeeCR_1b',mode='a')
    df_out_ZeeCR_2b.to_root(outfilenameis, key='bbDM_ZeeCR_2b',mode='a')
    df_out_ZmumuCR_1b.to_root(outfilenameis, key='bbDM_ZmumuCR_1b',mode='a')
    df_out_ZmumuCR_2b.to_root(outfilenameis, key='bbDM_ZmumuCR_2b',mode='a')

    df_out_WenuCR_1b.to_root(outfilenameis, key='bbDM_WenuCR_1b',mode='a')
    df_out_WenuCR_2b.to_root(outfilenameis, key='bbDM_WenuCR_2b',mode='a')
    df_out_WmunuCR_1b.to_root(outfilenameis, key='bbDM_WmunuCR_1b',mode='a')
    df_out_WmunuCR_2b.to_root(outfilenameis, key='bbDM_WmunuCR_2b',mode='a')

    df_out_TopenuCR_1b.to_root(outfilenameis, key='bbDM_TopenuCR_1b',mode='a')
    df_out_TopenuCR_2b.to_root(outfilenameis, key='bbDM_TopenuCR_2b',mode='a')
    df_out_TopmunuCR_1b.to_root(outfilenameis, key='bbDM_TopmunuCR_1b',mode='a')
    df_out_TopmunuCR_2b.to_root(outfilenameis, key='bbDM_TopmunuCR_2b',mode='a')

    print ('\n============SR cutflow============')
    print ('SR1bcount',SR1bcount,'SR2bcount',SR2bcount)
    print ('============SR cutflow============')

    print ('\n============Z cutflow============')
    print('ZeeCR1bcount',ZeeCR1bcount,'ZeeCR2bcount',ZeeCR2bcount)
    print ('ZmumuCR1count', ZmumuCR1count,'ZmumuCR2count', ZmumuCR2count)

    print ('\n============W cutflow============')
    print ('WenuCR1bcount', WenuCR1bcount,'WenuCR2bcount', WenuCR2bcount)
    print ('WmunuCR1bcount', WmunuCR1bcount,'WmunuCR2bcount', WmunuCR2bcount)

    print ('\n============Top cutflow============')
    print ('TopenuCR1bcount', TopenuCR1bcount,'TopenuCR2bcount', TopenuCR2bcount)
    print ('TopmunuCR1bcount', TopmunuCR1bcount,'TopmunuCR2bcount', TopmunuCR2bcount)

    cfsr_list = {1:'presel',2:'trigger',3:'MET',4:'nLep',5:'min_dPhi',6:'nJet',7:'nBjets'}
    cfcr_list = {1:'presel',2:'trigger',3:'MET',4:'nLep',5:'Recoil',6:'min_dPhi',7:'Z/W mass',8:'nJet',9:'nBjets'}
    for i in [1,2,3,4,5,6,7]:
        h_reg_SR_1b_cutFlow.GetXaxis().SetBinLabel(i,cfsr_list[i])
        h_reg_SR_2b_cutFlow.GetXaxis().SetBinLabel(i,cfsr_list[i])
    for i in [1,2,3,4,5,6,7,8,9]:
        h_reg_ZeeCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_ZeeCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_ZmumuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_ZmumuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_WenuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_WenuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_WmunuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_WmunuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_TopenuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_TopenuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_TopmunuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
        h_reg_TopmunuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cfcr_list[i])
    h_reg_SR_1b_cutFlow.SetEntries(1)
    h_reg_SR_2b_cutFlow.SetEntries(1)
    h_reg_ZeeCR_1b_cutFlow.SetEntries(1)
    h_reg_ZeeCR_2b_cutFlow.SetEntries(1)
    h_reg_ZmumuCR_1b_cutFlow.SetEntries(1)
    h_reg_ZmumuCR_2b_cutFlow.SetEntries(1)
    h_reg_WenuCR_1b_cutFlow.SetEntries(1)
    h_reg_WenuCR_2b_cutFlow.SetEntries(1)
    h_reg_WmunuCR_1b_cutFlow.SetEntries(1)
    h_reg_WmunuCR_2b_cutFlow.SetEntries(1)
    h_reg_TopenuCR_1b_cutFlow.SetEntries(1)
    h_reg_TopenuCR_2b_cutFlow.SetEntries(1)
    h_reg_TopmunuCR_1b_cutFlow.SetEntries(1)
    h_reg_TopmunuCR_2b_cutFlow.SetEntries(1)
    print ('===============================\n')
    print ("output written to ", outfilename)
    outfile = TFile(outfilenameis,'UPDATE')
    outfile.cd()
    h_total_mcweight.Write()
    h_total.Write()
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
    outfile.Write()
    outfile.Close()

    end = time.clock()
    print "%.4gs" % (end-start)



if __name__ == '__main__':
    if not runInteractive:
        txtFile=infile
        runbbdm(txtFile)

    if runInteractive and runOnTxt:
        filesPath = dirName+'/*txt'
        files     = glob.glob(filesPath)
        n = 8 #submit n txt files at a time, make equal to cores
        final = [files[i * n:(i + 1) * n] for i in range((len(files) + n - 1) // n )]
        if istest:
            runbbdm,final[0]
        else:
            for i in range(len(final)):
                try:
                    pool = mp.Pool(8)
                    pool.map(runbbdm,final[i])
                    pool.close()
                    pool.join()
                except Exception as e:
                    print e
                    print "Corrupt file inside input txt file is detected! Skipping this txt file:  ", final[i]
                    continue

    if runInteractive and not runOnTxt:
        ''' following part is for interactive running. This is still under testing because output file name can't be changed at this moment '''
        inputpath="/home/ptiwari/t3store3/CMSSW_10_3_0/src/test_syst/ExoPieProducer/ExoPieAnalyzer/testfile"

        os.system('rm dirlist.txt')
        os.system("ls -1 "+inputpath+" > dirlist.txt")

        allkeys=[idir.rstrip() for idir in open('dirlist.txt')]
        alldirs=[inputpath+"/"+idir.rstrip() for idir in open('dirlist.txt')]

        pool = mp.Pool(2)
        allsample=[]
        for ikey in allkeys:
            dirpath=inputpath+"/"+ikey
            txtfile=ikey+".txt"
            os.system ("find "+dirpath+"  -name \"*.root\" | grep -v \"failed\"  > "+txtfile)
            fileList=TextToList(txtfile)
            ## this is the list, first element is txt file with all the files and second element is the ikey (kind of sample name identifier)
            sample_  = [txtfile, ikey]
            ## push information about one sample into global list.
            allsample.append(sample_)
        print allsample
        if istest:
            runbbdm(allsample[0])
        else:
            pool.map(runbbdm, allsample)
        ## this works fine but the output file name get same value becuase it is done via a text file at the moment, need to find a better way,
