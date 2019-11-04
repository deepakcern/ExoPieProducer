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


isCondor = False
runInteractive = True
testing=True
## from commonutils
if isCondor:sys.path.append('ExoPieUtils/commonutils/')
else:sys.path.append('../../ExoPieUtils/commonutils/')
import MathUtils as mathutil
from MathUtils import *
import BooleanUtils as boolutil


## from analysisutils
if isCondor:sys.path.append('ExoPieUtils/analysisutils/')
else:sys.path.append('../../ExoPieUtils/analysisutils/')
import analysis_utils as anautil


sys.path.append('configs')
import variables as var
import outvars_bbDM as out

## from analysisutils
if isCondor:sys.path.append('ExoPieUtils/scalefactortools/')
else:sys.path.append('../../ExoPieUtils/scalefactortools/')

##please change the era accordingly
year_file= open("Year.py","w")
year_file.write('era="2016"')
year_file.close()
import ana_weight as wgt



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


if isfarmout:
    infile  = args.inputfile

else: print "No file is provided for farmout"


outputdir = '.'
if args.outputdir:
    outputdir = str(args.outputdir)

infilename = "NCUGlobalTuples.root"

outDir=outputdir

def TextToList(textfile):
    return([iline.rstrip() for iline in open(textfile)])

## the input file list and key is caught in one variable as a python list,
#### first element is the list of rootfiles
#### second element is the key, user to name output.root

dummy = -9999.0
def runbbdm(txtfile):

    print "in main function"

    infile_=[]
    outfilename=""
    prefix="Skimmed_"
    ikey_ = ""

    if  runInteractive:
        # print "running for ", txtfile[0]
        # infile_  = TextToList(txtfile[0])
        # key_=txtfile[1]
        # outfilename= txtfile[0].split('/')[-1].replace('.root.txt','.root')#prefix+key_+".root"
        print "running for ", txtfile
        infile_  = TextToList(txtfile)
        outfilename= outDir+'/'+txtfile.split('/')[-1].replace('.txt','.root')#prefix+key_+".root"

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

    #outputfilename = args.outputfile

    h_total = TH1F('h_total','h_total',2,0,2)
    h_total_mcweight = TH1F('h_total_mcweight','h_total_mcweight',2,0,2)

    h_reg_SR_1b_cutFlow  = TH1F("h_reg_SR_1b_cutFlow", "h_reg_SR_1b_cutFlow", 5,0,5)
    h_reg_SR_2b_cutFlow  = TH1F("h_reg_SR_2b_cutFlow", "h_reg_SR_2b_cutFlow", 5,0,5)

    h_reg_ZeeCR_1b_cutFlow = TH1F("h_reg_ZeeCR_1b_cutFlow", "h_reg_ZeeCR_1b_cutFlow", 6,0,6)
    h_reg_ZeeCR_2b_cutFlow = TH1F("h_reg_ZeeCR_2b_cutFlow", "h_reg_ZeeCR_2b_cutFlow", 6,0,6)
    h_reg_ZmumuCR_1b_cutFlow = TH1F("h_reg_ZmumuCR_1b_cutFlow", "h_reg_ZmumuCR_1b_cutFlow", 6,0,6)
    h_reg_ZmumuCR_2b_cutFlow = TH1F("h_reg_ZmumuCR_2b_cutFlow", "h_reg_ZmumuCR_2b_cutFlow", 6,0,6)

    h_reg_WenuCR_1b_cutFlow = TH1F("h_reg_WenuCR_1b_cutFlow", "h_reg_WenuCR_1b_cutFlow", 6,0,6)
    h_reg_WenuCR_2b_cutFlow = TH1F("h_reg_WenuCR_2b_cutFlow", "h_reg_WenuCR_2b_cutFlow", 6,0,6)
    h_reg_WmunuCR_1b_cutFlow = TH1F("h_reg_WmunuCR_1b_cutFlow", "h_reg_WmunuCR_1b_cutFlow", 6,0,6)
    h_reg_WmunuCR_2b_cutFlow = TH1F("h_reg_WmunuCR_2b_cutFlow", "h_reg_WmunuCR_2b_cutFlow", 6,0,6)

    h_reg_TopenuCR_1b_cutFlow = TH1F("h_reg_TopenuCR_1b_cutFlow", "h_reg_TopenuCR_1b_cutFlow", 6,0,6)
    h_reg_TopenuCR_2b_cutFlow = TH1F("h_reg_TopenuCR_2b_cutFlow", "h_reg_TopenuCR_2b_cutFlow", 6,0,6)
    h_reg_TopmunuCR_1b_cutFlow = TH1F("h_reg_TopmunuCR_1b_cutFlow", "h_reg_TopmunuCR_1b_cutFlow", 6,0,6)
    h_reg_TopmunuCR_2b_cutFlow = TH1F("h_reg_TopmunuCR_2b_cutFlow", "h_reg_TopmunuCR_2b_cutFlow", 6,0,6)

    for infl in infile_:
        f_tmp = TFile.Open(infl,'READ')
        h_tmp = f_tmp.Get('h_total')
        h_tmp_weight = f_tmp.Get('h_total_mcweight')
        h_total.Add(h_tmp)
        h_total_mcweight.Add(h_tmp_weight)

    filename = infile_
    ieve = 0;icount = 0
    cut_ep_THINnJet_1b = 0.; cut_ep_THINjetDeepCSV_1b = 0.
    cut_ep_THINnJet_2b = 0.; cut_ep_THINjetDeepCSV_2b = 0.
    cut_ep_nLep = 0.0; cut_ep_pfMetCorrPt = 0.; cut_min_dPhi = 0.0

    Zee_cut_ep_THINnJet_1b = 0.; Zee_cut_ep_THINjetDeepCSV_1b = 0.
    Zee_cut_ep_THINnJet_2b = 0.; Zee_cut_ep_THINjetDeepCSV_2b = 0.
    Zee_cut_ep_nLep = 0.0; Zee_cut_min_dPhi = 0.0; Zee_cut_ep_Recoil=0
    Zee_cut_ep_Zeemass = 0.0; ZeeCR1bcount=0.0; ZeeCR2bcount=0.0

    Zmumu_cut_ep_THINnJet_1b = 0.; Zmumu_cut_ep_THINjetDeepCSV_1b = 0.
    Zmumu_cut_ep_THINnJet_2b = 0.; Zmumu_cut_ep_THINjetDeepCSV_2b = 0.
    Zmumu_cut_ep_nLep = 0.0; Zmumu_cut_min_dPhi = 0.0; Zmumu_cut_ep_Recoil=0
    Zmumu_cut_ep_Zmumumass = 0.0; ZmumuCR1count=0.0; ZmumuCR2count=0.0

    Wenu_cut_ep_THINnJet_1b = 0.; Wenu_cut_ep_THINjetDeepCSV_1b = 0.
    Wenu_cut_ep_THINnJet_2b = 0.; Wenu_cut_ep_THINjetDeepCSV_2b = 0.
    Wenu_cut_ep_nLep = 0.0; Wenu_cut_min_dPhi = 0.0; Wenu_cut_ep_Recoil=0
    Wenu_cut_ep_Wenumass = 0.0; WenuCR1bcount=0.0; WenuCR2bcount=0.0

    Wmunu_cut_ep_THINnJet_1b = 0.; Wmunu_cut_ep_THINjetDeepCSV_1b = 0.
    Wmunu_cut_ep_THINnJet_2b = 0.; Wmunu_cut_ep_THINjetDeepCSV_2b = 0.
    Wmunu_cut_ep_nLep = 0.0; Wmunu_cut_min_dPhi = 0.0; Wmunu_cut_ep_Recoil=0
    Wmunu_cut_ep_Wmunumass = 0.0; WmunuCR1bcount=0.0; WmunuCR2bcount=0.0

    Topenu_cut_ep_THINnJet_1b = 0.; Topenu_cut_ep_THINjetDeepCSV_1b = 0.
    Topenu_cut_ep_THINnJet_2b = 0.; Topenu_cut_ep_THINjetDeepCSV_2b = 0.
    Topenu_cut_ep_nLep = 0.0; Topenu_cut_min_dPhi = 0.0; Topenu_cut_ep_Recoil=0
    Topenu_cut_ep_Topenumass = 0.0; TopenuCR1bcount=0.0; TopenuCR2bcount=0.0

    Topmunu_cut_ep_THINnJet_1b = 0.; Topmunu_cut_ep_THINjetDeepCSV_1b = 0.
    Topmunu_cut_ep_THINnJet_2b = 0.; Topmunu_cut_ep_THINjetDeepCSV_2b = 0.
    Topmunu_cut_ep_nLep = 0.0; Topmunu_cut_min_dPhi = 0.0; Topmunu_cut_ep_Recoil=0
    Topmunu_cut_ep_Topmunumass = 0.0; TopmunuCR1bcount=0.0; TopmunuCR2bcount=0.0

    SR1bcount = 0.0; SR2bcount = 0.0

    for df in read_root(filename, 'outTree', columns=var.allvars_bbDM, chunksize=125000):
        for ep_runId, ep_lumiSection, ep_eventId, \
            ep_pfMetCorrPt, ep_pfMetCorrPhi, ep_pfMetUncJetResUp, ep_pfMetUncJetResDown, ep_pfMetUncJetEnUp, ep_pfMetUncJetEnDown, \
            ep_WenuPhi, ep_WmunuPhi, ep_ZeePhi, ep_ZmumuPhi, \
            ep_ZeeRecoil, ep_ZmumuRecoil, ep_WenuRecoil, ep_WmunuRecoil, \
            ep_Zeemass, ep_Zmumumass, ep_Wenumass, ep_Wmunumass, \
            ep_isData, \
            ep_THINnJet, ep_THINjetPx, ep_THINjetPy, ep_THINjetPz, ep_THINjetEnergy, \
            ep_THINjetDeepCSV, ep_THINjetHadronFlavor, \
            ep_THINjetNHadEF, ep_THINjetCHadEF, ep_THINjetCEmEF, ep_THINjetPhoEF, ep_THINjetEleEF, ep_THINjetMuoEF, \
            ep_THINjetCorrUnc, \
            ep_nEle, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy, \
            ep_eleIsPassTight, ep_eleIsPassLoose, \
            ep_nPho, ep_phoIsPassTight, ep_phoPx, ep_phoPy, ep_phoPz, ep_phoEnergy, \
            ep_nMu, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy, ep_isTightMuon, \
            ep_nTau_discBased_looseElelooseMuVeto,ep_nTau_discBased_looseEleTightMuVeto,ep_nTau_discBased_mediumElelooseMuVeto,ep_nTau_discBased_TightEleTightMuVeto,\
            ep_pu_nTrueInt, ep_pu_nPUVert, \
            ep_THINjetNPV, \
            ep_mcweight, ep_genParPt, ep_genParSample, eletrigdecision, mutrigdecision, mettrigdecision \
            in zip(df.st_runId, df.st_lumiSection, df.st_eventId, \
                   df.st_pfMetCorrPt, df.st_pfMetCorrPhi, df.st_pfMetUncJetResUp, df.st_pfMetUncJetResDown, df.st_pfMetUncJetEnUp, df.st_pfMetUncJetEnDown, \
                   df.WenuPhi, df.WmunuPhi, df.ZeePhi, df.ZmumuPhi, \
                   df.ZeeRecoil, df.ZmumuRecoil, df.WenuRecoil, df.WmunuRecoil, \
                   df.ZeeMass, df.ZmumuMass, df.Wenumass, df.Wmunumass, \
                   df.st_isData, \
                   df.st_THINnJet, df.st_THINjetPx, df.st_THINjetPy, df.st_THINjetPz, df.st_THINjetEnergy, \
                   df.st_THINjetDeepCSV, df.st_THINjetHadronFlavor, \
                   df.st_THINjetNHadEF, df.st_THINjetCHadEF, df.st_THINjetCEmEF, df.st_THINjetPhoEF, df.st_THINjetEleEF, df.st_THINjetMuoEF, \
                   df.st_THINjetCorrUnc, \
                   df.st_nEle, df.st_elePx, df.st_elePy, df.st_elePz, df.st_eleEnergy, \
                   df.st_eleIsPassTight, df.st_eleIsPassLoose, \
                   df.st_nPho, df.st_phoIsPassTight, df.st_phoPx, df.st_phoPy, df.st_phoPz, df.st_phoEnergy, \
                   df.st_nMu, df.st_muPx, df.st_muPy, df.st_muPz, df.st_muEnergy, df.st_isTightMuon, \
                   df.st_nTau_discBased_looseElelooseMuVeto,df.st_nTau_discBased_looseEletightMuVeto,df.st_nTau_discBased_mediumElelooseMuVeto,df.st_nTau_discBased_tightEletightMuVeto,\
                   df.st_pu_nTrueInt, df.st_pu_nPUVert, \
                   df.st_THINjetNPV, \
                   df.mcweight, df.st_genParPt, df.st_genParSample, df.st_eletrigdecision, df.st_mutrigdecision, df.st_mettrigdecision \
                   ):

            ieve = ieve + 1
            if ieve%10000==0: print "Processed",ieve,"Events"

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
            deepCSV_Med = 0.6321

            '''
            -------------------------------------------------------------------------------
            electron VARS
            -------------------------------------------------------------------------------
            '''
            ep_nEle_ = [ij for ij in range(ep_nEle) if (ep_eleIsPassLoose[ij])]
            ep_nEle_index = len(ep_nEle_)
            ep_elePt  = [getPt(ep_elePx[ij], ep_elePy[ij]) for ij in ep_nEle_]
            ep_eleEta = [getEta(ep_elePx[ij], ep_elePy[ij], ep_elePz[ij]) for ij in ep_nEle_]
            ep_elePhi = [getPhi(ep_elePx[ij], ep_elePy[ij]) for ij in ep_nEle_]
            ep_eleIsTight  = [ep_eleIsPassTight[ij] for ij in ep_nEle_]
            # print ([ij for ij in range(ep_nEle) if (ep_eleIsPassLoose[ij])], ep_elePt)
            '''
            -------------------------------------------------------------------------------
            muon VARS
            -------------------------------------------------------------------------------
            '''
            ep_muPt = [getPt(ep_muPx[ij], ep_muPy[ij]) for ij in range(ep_nMu)]
            ep_muEta = [getEta(ep_muPx[ij], ep_muPy[ij], ep_muPz[ij]) for ij in range(ep_nMu)]
            ep_muPhi = [getPhi(ep_muPx[ij], ep_muPy[ij]) for ij in range(ep_nMu)]
            '''

            -------------------------------------------------------------------------------
            photon VARS
            -------------------------------------------------------------------------------
            '''
            ep_phoPt = [getPt(ep_phoPx[ij], ep_phoPy[ij]) for ij in range(ep_nPho)]
            ep_phoEta = [getEta(ep_phoPx[ij], ep_phoPy[ij], ep_phoPz[ij]) for ij in range(ep_nPho)]
            ep_phoPhi = [getPhi(ep_phoPx[ij], ep_phoPy[ij]) for ij in range(ep_nPho)]

            '''
            -------------------------------------------------------------------------------
            THIN JET VARS
            -------------------------------------------------------------------------------
            '''
            ep_THINjetPt = [getPt(ep_THINjetPx[ij], ep_THINjetPy[ij]) for ij in range(ep_THINnJet)]
            ep_THINjetEta = [getEta(ep_THINjetPx[ij], ep_THINjetPy[ij], ep_THINjetPz[ij]) for ij in range(ep_THINnJet)]
            ep_THINjetPhi = [getPhi(ep_THINjetPx[ij], ep_THINjetPy[ij]) for ij in range(ep_THINnJet)]
            ep_THINbjets_index = [ij for ij in range(ep_THINnJet) if (ep_THINjetDeepCSV[ij] > deepCSV_Med and abs(ep_THINjetEta[ij]) < 2.5)]
            nBjets = len(ep_THINbjets_index)

            if len(ep_THINjetPt)==0 : continue

            # WenuPhi = WmunuPhi = ZeePhi = ZmumuPhi = 0.001
            # ep_ZeeRecoil = ep_ZmumuRecoil = ep_WenuRecoil = ep_WmunuRecoil = 200.0
            min_dPhi_jet_MET = min([DeltaPhi(jet_phi,ep_pfMetCorrPhi) for jet_phi in ep_THINjetPhi])

            Jet2Pt  = dummy;Jet2Eta     = dummy
            Jet2Phi = dummy;Jet2deepCSV = dummy
            Jet3Pt  = dummy;Jet3Eta     = dummy
            Jet3Phi = dummy;Jet3deepCSV = dummy

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
            SIGNAL REGION Cutflow
            --------------------------------------------------------------------------------
            '''
            if (ep_pfMetCorrPt > 200.):
               cut_ep_pfMetCorrPt +=1
               if (ep_nEle_index == 0) and (ep_nMu == 0) and (ep_nPho ==0) and (ep_nTau_discBased_TightEleTightMuVeto==0):
                   cut_ep_nLep+=1
                   if (min_dPhi_jet_MET > 0.5):
                       cut_min_dPhi +=1
                       if (ep_THINnJet ==1 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                           cut_ep_THINnJet_1b +=1;
                           if (ep_THINjetDeepCSV[0] > deepCSV_Med):
                               cut_ep_THINjetDeepCSV_1b+=1
                       elif (ep_THINnJet ==3 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                           cut_ep_THINnJet_2b +=1
                           if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med):
                               cut_ep_THINjetDeepCSV_2b+=1
            '''
            --------------------------------------------------------------------------------
            1b SIGNAL REGION
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for 1b SR.
            if (ep_THINnJet ==1 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_nEle_index == 0) and (ep_nMu == 0) and (ep_nPho ==0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_pfMetCorrPt > 200.) and (min_dPhi_jet_MET > 0.5) and mettrigdecision:
                isSR1b=True
                SR1bcount+=1
                if ep_THINnJet==2:
                    Jet2Pt  = ep_THINjetPt[1]; Jet2Eta     = ep_THINjetEta[1]
                    Jet2Phi = ep_THINjetPhi[1];Jet2deepCSV = ep_THINjetDeepCSV[1]
                ## cal function for each of them based on pt and eta
                weightMET=wgt.getMETtrig_First(ep_pfMetCorrPt)
                weightEle=1
                weightMu=1
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)

                weight = weightMET * weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            2b SIGNAL REGION
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for 2b SR.
            if (ep_THINnJet ==3 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and (ep_nEle_index == 0) and (ep_nMu == 0) and (ep_nPho ==0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_pfMetCorrPt > 200.) and (min_dPhi_jet_MET > 0.5) and mettrigdecision:
                isSR2b=True
                SR2bcount+=1
                if ep_THINnJet==3:
                    Jet3Pt  = ep_THINjetPt[2]; Jet3Eta     = ep_THINjetEta[2]
                    Jet3Phi = ep_THINjetPhi[2];Jet3deepCSV = ep_THINjetDeepCSV[2]
                ## cal function for each of them based on pt and eta
                weightMET=wgt.getMETtrig_First(ep_pfMetCorrPt)
                weightEle=1
                weightMu=1
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0

                weightPU=wgt.puweight(ep_pu_nTrueInt)


                weight = weightMET*weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU


            '''
            --------------------------------------------------------------------------------
            ZEE CONTROL REGION Cutflow
            --------------------------------------------------------------------------------
            '''
            if (ep_ZeeRecoil > 200.) :
               Zee_cut_ep_Recoil +=1
               if (ep_nEle_index == 2) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_elePt[0] > 30.) and (ep_eleIsTight[0]) :
                   Zee_cut_ep_nLep+=1
                   if (min_dPhi_jet_MET > 0.5):
                       Zee_cut_min_dPhi +=1
                       if (ep_Zeemass >= 60 and ep_Zeemass <= 120):
                           Zee_cut_ep_Zeemass +=1
                           if (ep_THINnJet ==1 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               Zee_cut_ep_THINnJet_1b +=1;
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med):
                                   Zee_cut_ep_THINjetDeepCSV_1b+=1
                           if (ep_THINnJet ==3 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               Zee_cut_ep_THINnJet_2b +=1
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med):
                                   Zee_cut_ep_THINjetDeepCSV_2b+=1
            '''
            --------------------------------------------------------------------------------
            ZEE CONTROL REGION 1b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Zee SR.
            if (ep_ZeeRecoil > 200.) and (ep_nEle_index == 2) and (ep_elePt[0] > 30.) and (ep_eleIsTight[0]) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Zeemass >= 60 and ep_Zeemass <= 120) and (ep_THINnJet ==1 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and eletrigdecision  :
                ZeeCR1bcount+=1
                is1bCRZee=True
                ## cal function for each of them based on pt and eta
                ele_trig = True
                no_ele_trig = False
                weightEle=wgt.ele_weight(ep_elePt[0],ep_eleEta[0],ele_trig,'T')*wgt.ele_weight(ep_elePt[1],ep_eleEta[1],no_ele_trig,'L')
                weightMu=1
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)


                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            ZEE CONTROL REGION 2b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Zee SR.
            if (ep_ZeeRecoil > 200.) and (ep_nEle_index == 2) and (ep_elePt[0] > 30.) and (ep_eleIsTight[0]) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Zeemass >= 60 and ep_Zeemass <= 120) and (ep_THINnJet ==3 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and eletrigdecision  :
                ZeeCR2bcount+=1
                is2bCRZee=True
                ## cal function for each of them based on pt and eta
                ele_trig = True
                no_ele_trig = False
                weightEle=wgt.ele_weight(ep_elePt[0],ep_eleEta[0],ele_trig,'T')*wgt.ele_weight(ep_elePt[1],ep_eleEta[1],no_ele_trig,'L')
                weightMu=1
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)

                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU

            '''
            --------------------------------------------------------------------------------
            ZMUMU CONTROL REGION Cutflow
            --------------------------------------------------------------------------------
            '''
            if (ep_ZmumuRecoil > 200.) :
                #print 'Zmumu_cut_ep_Recoil'
                Zmumu_cut_ep_Recoil +=1
                if (ep_nEle_index == 0) and (ep_nMu == 2) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_muPt[0] > 30.) and (ep_isTightMuon[0]):
                   #print 'Zmumu_cut_ep_nLep'
                   Zmumu_cut_ep_nLep+=1
                   if (min_dPhi_jet_MET > 0.5):
                       #print 'Zmumu_cut_min_dPhi'
                       Zmumu_cut_min_dPhi +=1
                       if (ep_Zmumumass >= 60 and ep_Zmumumass <= 120):
                           #print 'Zmumu_cut_ep_Zmumumass'
                           Zmumu_cut_ep_Zmumumass +=1
                           if (ep_THINnJet ==1 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               #print 'Zmumu_cut_ep_THINnJet_1b'
                               Zmumu_cut_ep_THINnJet_1b +=1;
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med):
                                   #print 'Zmumu_cut_ep_THINjetDeepCSV_1b'
                                   Zmumu_cut_ep_THINjetDeepCSV_1b+=1
                           if (ep_THINnJet ==3 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               Zmumu_cut_ep_THINnJet_2b +=1
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med):
                                   Zmumu_cut_ep_THINjetDeepCSV_2b+=1
            '''
            --------------------------------------------------------------------------------
            ZMUMU CONTROL REGION 1b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Zmumu SR.
            if (ep_ZmumuRecoil > 200.) and (ep_nMu == 2) and (ep_muPt[0] > 30.) and (ep_isTightMuon[0]) and (ep_nEle_index == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Zmumumass >= 60 and ep_Zmumumass <= 120) and (ep_THINnJet ==1 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and mutrigdecision :
                ZmumuCR1count+=1
                is1bCRZmumu=True
                ## cal function for each of them based on pt and eta
                weightEle=1
                mu_trig = True
                no_mu_trig = False
                weightMu=wgt.mu_weight(ep_muPt[0],ep_muEta[0],mu_trig,'T')*wgt.mu_weight(ep_muPt[1],ep_muEta[1],no_mu_trig,'L')
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightMET=wgt.getMETtrig_First(ep_ZmumuRecoil)

                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            ZMUMU CONTROL REGION 2b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Zmumu SR.
            if (ep_ZmumuRecoil > 200.) and (ep_nMu == 2) and (ep_muPt[0] > 30.) and (ep_isTightMuon[0]) and (ep_nEle_index == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Zmumumass >= 60 and ep_Zmumumass <= 120) and (ep_THINnJet ==3 or ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and mutrigdecision  :
                ZmumuCR2count+=1
                is2bCRZmumu=True
                ## cal function for each of them based on pt and eta
                weightEle=1
                mu_trig = True
                no_mu_trig = False
                weightMu=wgt.mu_weight(ep_muPt[0],ep_muEta[0],mu_trig,'T')*wgt.mu_weight(ep_muPt[1],ep_muEta[1],no_mu_trig,'L')
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightMET=wgt.getMETtrig_First(ep_ZmumuRecoil)

                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            WENU CONTROL REGION Cutflow
            --------------------------------------------------------------------------------
            '''
            if (ep_WenuRecoil > 200.) :
                #print 'Wenu_cut_ep_Recoil'
                Wenu_cut_ep_Recoil +=1
                if (ep_nEle_index == 1) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_elePt[0] > 30.) and (ep_eleIsTight[0]):
                   #print 'Wenu_cut_ep_nLep'
                   Wenu_cut_ep_nLep+=1
                   if (min_dPhi_jet_MET > 0.5):
                       #print 'Wenu_cut_min_dPhi'
                       Wenu_cut_min_dPhi +=1
                       if (ep_Wenumass >= 0 and ep_Wenumass <= 160):
                           #print 'Wenu_cut_ep_Wenumass'
                           Wenu_cut_ep_Wenumass +=1
                           if (ep_THINnJet ==1) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               #print 'Wenu_cut_ep_THINnJet_1b'
                               Wenu_cut_ep_THINnJet_1b +=1;
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med):
                                   #print 'Wenu_cut_ep_THINjetDeepCSV_1b'
                                   Wenu_cut_ep_THINjetDeepCSV_1b+=1
                           elif (ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               Wenu_cut_ep_THINnJet_2b +=1
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med):
                                   Wenu_cut_ep_THINjetDeepCSV_2b+=1
            '''
            --------------------------------------------------------------------------------
            WENU CONTROL REGION 1b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Wenu SR.
            if (ep_WenuRecoil > 200.) and (ep_nEle_index == 1) and (ep_elePt[0] > 30.) and (ep_eleIsTight[0]) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Wenumass >= 0 and ep_Wenumass <= 160) and (ep_THINnJet ==1) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and eletrigdecision  :
                WenuCR1bcount+=1
                is1bCRWenu=True
                ## cal function for each of them based on pt and eta
                ele_trig = True
                weightEle=wgt.ele_weight(ep_elePt[0],ep_eleEta[0],ele_trig,'T')
                weightMu=1
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)


                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            WENU CONTROL REGION 2b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Wenu SR.
            if (ep_WenuRecoil > 200.) and (ep_nEle_index == 1) and (ep_elePt[0] > 30.) and (ep_eleIsTight[0]) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Wenumass >= 0 and ep_Wenumass <= 160) and (ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and eletrigdecision  :
                WenuCR2bcount+=1
                is2bCRWenu=True
                ## cal function for each of them based on pt and eta
                ele_trig = True
                weightEle=wgt.ele_weight(ep_elePt[0],ep_eleEta[0],ele_trig,'T')
                weightMu=1
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)


                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            WMUNU CONTROL REGION Cutflow
            --------------------------------------------------------------------------------
            '''
            if (ep_WmunuRecoil > 200.) :
                #print 'Wmunu_cut_ep_Recoil'
                Wmunu_cut_ep_Recoil +=1
                if (ep_nEle_index == 0) and (ep_nMu == 1) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_muPt[0] > 30.) and (ep_isTightMuon[0]):
                   #print 'Wmunu_cut_ep_nLep'
                   Wmunu_cut_ep_nLep+=1
                   if (min_dPhi_jet_MET > 0.5):
                       #print 'Wmunu_cut_min_dPhi'
                       Wmunu_cut_min_dPhi +=1
                       if (ep_Wmunumass >= 0 and ep_Wmunumass <= 160):
                           #print 'Wmunu_cut_ep_Wmunumass'
                           Wmunu_cut_ep_Wmunumass +=1
                           if (ep_THINnJet ==1) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               #print 'Wmunu_cut_ep_THINnJet_1b'
                               Wmunu_cut_ep_THINnJet_1b +=1;
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med):
                                   #print 'Wmunu_cut_ep_THINjetDeepCSV_1b'
                                   Wmunu_cut_ep_THINjetDeepCSV_1b+=1
                           elif (ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               Wmunu_cut_ep_THINnJet_2b +=1
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med):
                                   Wmunu_cut_ep_THINjetDeepCSV_2b+=1
            '''
            --------------------------------------------------------------------------------
            WMUNU CONTROL REGION 1b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Wmunu SR.
            if (ep_WmunuRecoil > 200.) and (ep_nMu == 1) and (ep_muPt[0] > 30.) and (ep_isTightMuon[0]) and (ep_nEle_index == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Wmunumass >= 0 and ep_Wmunumass <= 160) and (ep_THINnJet ==1) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and mutrigdecision  :
                WmunuCR1bcount+=1
                is1bCRWmunu=True
                ## cal function for each of them based on pt and eta
                weightEle=1
                mu_trig = True
                weightMu=wgt.mu_weight(ep_muPt[0],ep_muEta[0],mu_trig,'T')
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightMET=wgt.getMETtrig_First(ep_WmunuRecoil)

                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            WMUNU CONTROL REGION 2b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Wmunu SR.
            if (ep_WmunuRecoil > 200.) and (ep_nMu == 1) and (ep_muPt[0] > 30.) and (ep_isTightMuon[0]) and (ep_nEle_index == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Wmunumass >= 0 and ep_Wmunumass <= 160) and (ep_THINnJet ==2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and mutrigdecision  :
                WmunuCR2bcount+=1
                is2bCRWmunu=True
                ## cal function for each of them based on pt and eta
                weightEle=1
                mu_trig = True
                weightMu=wgt.mu_weight(ep_muPt[0],ep_muEta[0],mu_trig,'T')
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightMET=wgt.getMETtrig_First(ep_WmunuRecoil)

                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            TOPENU CONTROL REGION Cutflow
            --------------------------------------------------------------------------------
            '''
            if (ep_WenuRecoil > 200.) :
                #print 'Wenu_cut_ep_Recoil'
                Topenu_cut_ep_Recoil +=1
                if (ep_nEle_index == 1) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_elePt[0] > 30.) and (ep_eleIsTight[0]):
                   #print 'Wenu_cut_ep_nLep'
                   Topenu_cut_ep_nLep+=1
                   if (min_dPhi_jet_MET > 0.5):
                       #print 'Topenu_cut_min_dPhi'
                       Topenu_cut_min_dPhi +=1
                       if (ep_Wenumass >= 0 and ep_Wenumass <= 160):
                           #print 'Wenu_cut_ep_Wenumass'
                           Topenu_cut_ep_Topenumass +=1
                           if (ep_THINnJet >1) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               Topenu_cut_ep_THINnJet_1b +=1;
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med):
                                   Topenu_cut_ep_THINjetDeepCSV_1b+=1
                           if (ep_THINnJet > 2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                                Topenu_cut_ep_THINnJet_2b +=1
                                if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med):
                                    Topenu_cut_ep_THINjetDeepCSV_2b+=1
            '''
            --------------------------------------------------------------------------------
            TOPENU CONTROL REGION 1b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Topenu SR.
            if (ep_WenuRecoil > 200.) and (ep_nEle_index == 1) and (ep_elePt[0] > 30.) and (ep_eleIsTight[0]) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Wenumass >= 0 and ep_Wenumass <= 160) and (ep_THINnJet >1) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and eletrigdecision  :
                TopenuCR1bcount+=1
                is1bCRTopenu=True
                ## cal function for each of them based on pt and eta
                ele_trig = True
                weightEle=wgt.ele_weight(ep_elePt[0],ep_eleEta[0],ele_trig,'T')
                weightMu=1
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)


                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            TOPENU CONTROL REGION 2b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Topenu SR.
            if (ep_WenuRecoil > 200.) and (ep_nEle_index == 1) and (ep_elePt[0] > 30.) and (ep_eleIsTight[0]) and (ep_nMu == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Wenumass >= 0 and ep_Wenumass <= 160) and (ep_THINnJet >2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and eletrigdecision  :
                TopenuCR2bcount+=1
                is2bCRTopenu=True
                ## cal function for each of them based on pt and eta
                ele_trig = True
                weightEle=wgt.ele_weight(ep_elePt[0],ep_eleEta[0],ele_trig,'T')
                weightMu=1
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)


                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            TOPMUNU CONTROL REGION Cutflow
            --------------------------------------------------------------------------------
            '''
            if (ep_WmunuRecoil > 200.) :
                #print 'Wenu_cut_ep_Recoil'
                Topmunu_cut_ep_Recoil +=1
                if (ep_nEle_index == 0) and (ep_nMu == 1) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (ep_muPt[0] > 30.) and (ep_isTightMuon[0]):
                   #print 'Wenu_cut_ep_nLep'
                   Topmunu_cut_ep_nLep+=1
                   if (min_dPhi_jet_MET > 0.5):
                       #print 'Topmunu_cut_min_dPhi'
                       Topmunu_cut_min_dPhi +=1
                       if (ep_Wmunumass >= 0 and ep_Wmunumass <= 160):
                           #print 'Wenu_cut_ep_Wenumass'
                           Topmunu_cut_ep_Topmunumass +=1
                           if (ep_THINnJet >1) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               #print 'Wenu_cut_ep_THINnJet_1b'
                               Topmunu_cut_ep_THINnJet_1b +=1;
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med):
                                   #print 'Wenu_cut_ep_THINjetDeepCSV_1b'
                                   Topmunu_cut_ep_THINjetDeepCSV_1b+=1
                           if (ep_THINnJet >2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                               Topmunu_cut_ep_THINnJet_2b +=1
                               if (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med):
                                   Topmunu_cut_ep_THINjetDeepCSV_2b+=1
            '''
            --------------------------------------------------------------------------------
            TOPMUNU CONTROL REGION 1b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Topmunu SR.
            if (ep_WmunuRecoil > 200.) and (ep_nMu == 1) and (ep_muPt[0] > 30.) and (ep_isTightMuon[0]) and (ep_nEle_index == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Wmunumass >= 0 and ep_Wmunumass <= 160) and (ep_THINnJet >1) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and mutrigdecision  :
                TopmunuCR1bcount+=1
                is1bCRTopmunu=True
                ## cal function for each of them based on pt and eta
                weightEle=1
                mu_trig = True
                weightMu=wgt.mu_weight(ep_muPt[0],ep_muEta[0],mu_trig,'T')
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightMET=wgt.getMETtrig_First(ep_WmunuRecoil)

                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU
            '''
            --------------------------------------------------------------------------------
            TOPMUNU CONTROL REGION 2b
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for Topmunu SR.
            if (ep_WmunuRecoil > 200.) and (ep_nMu == 1) and (ep_muPt[0] > 30.) and (ep_isTightMuon[0]) and (ep_nEle_index == 0) and (ep_nTau_discBased_TightEleTightMuVeto==0) and (min_dPhi_jet_MET > 0.5) and (ep_Wmunumass >= 0 and ep_Wmunumass <= 160) and (ep_THINnJet > 2) and (ep_THINjetPt[0] > 50.) and (ep_THINjetDeepCSV[0] > deepCSV_Med) and (ep_THINjetDeepCSV[1] > deepCSV_Med) and mutrigdecision  :
                TopmunuCR2bcount+=1
                is2bCRTopmunu=True
                ## cal function for each of them based on pt and eta
                weightEle=1
                mu_trig = True
                weightMu=wgt.mu_weight(ep_muPt[0],ep_muEta[0],mu_trig,'T')
                weightB=wgt.getBTagSF(ep_THINnJet,ep_THINjetPt,ep_THINjetEta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)
                weightTau=1
                if ep_genParSample==23:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    if len(ep_genParPt)==1: weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                else: weightEWK = 1.0
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightMET=wgt.getMETtrig_First(ep_WmunuRecoil)

                weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU

            if isSR1b:
                df_out_SR_1b = df_out_SR_1b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0],'Jet1Eta':ep_THINjetEta[0],'Jet1Phi':ep_THINjetPhi[0],'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':Jet2Pt,'Jet2Eta':Jet2Eta,'Jet2Phi':Jet2Phi,'Jet2deepCSV':Jet2deepCSV,
                                                    'Jet3Pt':dummy,'Jet3Eta':dummy,'Jet3Phi':dummy,'Jet3deepCSV':dummy,
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('isSR1b')
            if isSR2b:
                df_out_SR_2b = df_out_SR_2b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0], 'Jet1Eta':ep_THINjetEta[0], 'Jet1Phi':ep_THINjetPhi[0], 'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':ep_THINjetPt[1], 'Jet2Eta':ep_THINjetEta[1], 'Jet2Phi':ep_THINjetPhi[1], 'Jet2deepCSV':ep_THINjetDeepCSV[1],
                                                    'Jet3Pt':Jet3Pt, 'Jet3Eta':Jet3Eta, 'Jet3Phi':Jet3Phi, 'Jet3deepCSV':Jet3deepCSV,
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('isSR2b')

            if is1bCRZee:
                df_out_ZeeCR_1b = df_out_ZeeCR_1b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_ZeeRecoil ,'Zmass':ep_Zeemass,'ZpT':ZpT_ee,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0],'Jet1Eta':ep_THINjetEta[0],'Jet1Phi':ep_THINjetPhi[0],'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':Jet2Pt,'Jet2Eta':Jet2Eta,'Jet2Phi':Jet2Phi,'Jet2deepCSV':Jet2deepCSV,
                                                    'Jet3Pt':dummy,'Jet3Eta':dummy,'Jet3Phi':dummy,'Jet3deepCSV':dummy,
                                                    'leadingLepPt':ep_elePt[0],'leadingLepEta':ep_eleEta[0],'leadingLepPhi':ep_elePhi[0],
                                                    'subleadingLepPt':ep_elePt[1],'subleadingLepEta':ep_eleEta[1],'subleadingLepPhi':ep_elePhi[1],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRZee')
            if is2bCRZee:
                df_out_ZeeCR_2b = df_out_ZeeCR_2b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_ZeeRecoil ,'Zmass':ep_Zeemass,'ZpT':ZpT_ee,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0], 'Jet1Eta':ep_THINjetEta[0], 'Jet1Phi':ep_THINjetPhi[0], 'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':ep_THINjetPt[1], 'Jet2Eta':ep_THINjetEta[1], 'Jet2Phi':ep_THINjetPhi[1], 'Jet2deepCSV':ep_THINjetDeepCSV[1],
                                                    'Jet3Pt':Jet3Pt, 'Jet3Eta':Jet3Eta, 'Jet3Phi':Jet3Phi, 'Jet3deepCSV':Jet3deepCSV,
                                                    'leadingLepPt':ep_elePt[0],'leadingLepEta':ep_eleEta[0],'leadingLepPhi':ep_elePhi[0],
                                                    'subleadingLepPt':ep_elePt[1],'subleadingLepEta':ep_eleEta[1],'subleadingLepPhi':ep_elePhi[1],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRZee')
            if is1bCRZmumu:
                df_out_ZmumuCR_1b = df_out_ZmumuCR_1b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_ZmumuRecoil ,'Zmass':ep_Zmumumass,'ZpT':ZpT_mumu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0],'Jet1Eta':ep_THINjetEta[0],'Jet1Phi':ep_THINjetPhi[0],'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':Jet2Pt,'Jet2Eta':Jet2Eta,'Jet2Phi':Jet2Phi,'Jet2deepCSV':Jet2deepCSV,
                                                    'Jet3Pt':dummy,'Jet3Eta':dummy,'Jet3Phi':dummy,'Jet3deepCSV':dummy,
                                                    'leadingLepPt':ep_muPt[0],'leadingLepEta':ep_muEta[0],'leadingLepPhi':ep_muPhi[0],
                                                    'subleadingLepPt':ep_muPt[1],'subleadingLepEta':ep_muEta[1],'subleadingLepPhi':ep_muPhi[1],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRZmumu')
            if is2bCRZmumu:
                df_out_ZmumuCR_2b = df_out_ZmumuCR_2b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_ZmumuRecoil ,'Zmass':ep_Zmumumass,'ZpT':ZpT_mumu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0], 'Jet1Eta':ep_THINjetEta[0], 'Jet1Phi':ep_THINjetPhi[0], 'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':ep_THINjetPt[1], 'Jet2Eta':ep_THINjetEta[1], 'Jet2Phi':ep_THINjetPhi[1], 'Jet2deepCSV':ep_THINjetDeepCSV[1],
                                                    'Jet3Pt':Jet3Pt, 'Jet3Eta':Jet3Eta, 'Jet3Phi':Jet3Phi, 'Jet3deepCSV':Jet3deepCSV,
                                                    'leadingLepPt':ep_muPt[0],'leadingLepEta':ep_muEta[0],'leadingLepPhi':ep_muPhi[0],
                                                    'subleadingLepPt':ep_muPt[1],'subleadingLepEta':ep_muEta[1],'subleadingLepPhi':ep_muPhi[1],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRZmumu')
            if is1bCRWenu:
                df_out_WenuCR_1b = df_out_WenuCR_1b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_WenuRecoil ,'Wmass':ep_Wenumass,'WpT':WpT_enu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0],'Jet1Eta':ep_THINjetEta[0],'Jet1Phi':ep_THINjetPhi[0],'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':Jet2Pt,'Jet2Eta':Jet2Eta,'Jet2Phi':Jet2Phi,'Jet2deepCSV':Jet2deepCSV,
                                                    'Jet3Pt':dummy,'Jet3Eta':dummy,'Jet3Phi':dummy,'Jet3deepCSV':dummy,
                                                    'leadingLepPt':ep_elePt[0],'leadingLepEta':ep_eleEta[0],'leadingLepPhi':ep_elePhi[0],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRWenu')
            if is2bCRWenu:
                df_out_WenuCR_2b = df_out_WenuCR_2b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_WenuRecoil ,'Wmass':ep_Wenumass,'WpT':WpT_enu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0], 'Jet1Eta':ep_THINjetEta[0], 'Jet1Phi':ep_THINjetPhi[0], 'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':ep_THINjetPt[1], 'Jet2Eta':ep_THINjetEta[1], 'Jet2Phi':ep_THINjetPhi[1], 'Jet2deepCSV':ep_THINjetDeepCSV[1],
                                                    'Jet3Pt':Jet3Pt, 'Jet3Eta':Jet3Eta, 'Jet3Phi':Jet3Phi, 'Jet3deepCSV':Jet3deepCSV,
                                                    'leadingLepPt':ep_elePt[0],'leadingLepEta':ep_eleEta[0],'leadingLepPhi':ep_elePhi[0],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRWenu')
            if is1bCRWmunu:
                df_out_WmunuCR_1b = df_out_WmunuCR_1b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_WmunuRecoil ,'Wmass':ep_Wmunumass,'WpT':WpT_munu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0],'Jet1Eta':ep_THINjetEta[0],'Jet1Phi':ep_THINjetPhi[0],'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':Jet2Pt,'Jet2Eta':Jet2Eta,'Jet2Phi':Jet2Phi,'Jet2deepCSV':Jet2deepCSV,
                                                    'Jet3Pt':dummy,'Jet3Eta':dummy,'Jet3Phi':dummy,'Jet3deepCSV':dummy,
                                                    'leadingLepPt':ep_muPt[0],'leadingLepEta':ep_muEta[0],'leadingLepPhi':ep_muPhi[0],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRWmunu')
            if is2bCRWmunu:
                df_out_WmunuCR_2b = df_out_WmunuCR_2b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_WmunuRecoil ,'Wmass':ep_Wmunumass,'WpT':WpT_munu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0], 'Jet1Eta':ep_THINjetEta[0], 'Jet1Phi':ep_THINjetPhi[0], 'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':ep_THINjetPt[1], 'Jet2Eta':ep_THINjetEta[1], 'Jet2Phi':ep_THINjetPhi[1], 'Jet2deepCSV':ep_THINjetDeepCSV[1],
                                                    'Jet3Pt':Jet3Pt, 'Jet3Eta':Jet3Eta, 'Jet3Phi':Jet3Phi, 'Jet3deepCSV':Jet3deepCSV,
                                                    'leadingLepPt':ep_muPt[0],'leadingLepEta':ep_muEta[0],'leadingLepPhi':ep_muPhi[0],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRWmunu')
            if is1bCRTopenu:
                df_out_TopenuCR_1b = df_out_TopenuCR_1b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_WenuRecoil ,'Wmass':ep_Wenumass,'WpT':WpT_enu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0],'Jet1Eta':ep_THINjetEta[0],'Jet1Phi':ep_THINjetPhi[0],'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':Jet2Pt,'Jet2Eta':Jet2Eta,'Jet2Phi':Jet2Phi,'Jet2deepCSV':Jet2deepCSV,
                                                    'Jet3Pt':dummy,'Jet3Eta':dummy,'Jet3Phi':dummy,'Jet3deepCSV':dummy,
                                                    'leadingLepPt':ep_elePt[0],'leadingLepEta':ep_eleEta[0],'leadingLepPhi':ep_elePhi[0],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRTopenu')
            if is2bCRTopenu:
                df_out_TopenuCR_2b = df_out_TopenuCR_2b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_WenuRecoil ,'Wmass':ep_Wenumass,'WpT':WpT_enu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0], 'Jet1Eta':ep_THINjetEta[0], 'Jet1Phi':ep_THINjetPhi[0], 'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':ep_THINjetPt[1], 'Jet2Eta':ep_THINjetEta[1], 'Jet2Phi':ep_THINjetPhi[1], 'Jet2deepCSV':ep_THINjetDeepCSV[1],
                                                    'Jet3Pt':Jet3Pt, 'Jet3Eta':Jet3Eta, 'Jet3Phi':Jet3Phi, 'Jet3deepCSV':Jet3deepCSV,
                                                    'leadingLepPt':ep_elePt[0],'leadingLepEta':ep_eleEta[0],'leadingLepPhi':ep_elePhi[0],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is2bCRTopenu')
            if is1bCRTopmunu:
                df_out_TopmunuCR_1b = df_out_TopmunuCR_1b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_WmunuRecoil ,'Wmass':ep_Wmunumass,'WpT':WpT_munu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0],'Jet1Eta':ep_THINjetEta[0],'Jet1Phi':ep_THINjetPhi[0],'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':Jet2Pt,'Jet2Eta':Jet2Eta,'Jet2Phi':Jet2Phi,'Jet2deepCSV':Jet2deepCSV,
                                                    'Jet3Pt':dummy,'Jet3Eta':dummy,'Jet3Phi':dummy,'Jet3deepCSV':dummy,
                                                    'leadingLepPt':ep_muPt[0],'leadingLepEta':ep_muEta[0],'leadingLepPhi':ep_muPhi[0],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
                                                    },ignore_index=True
                                                   )
                if istest: print ('is1bCRTopmunu')
            if is2bCRTopmunu:
                df_out_TopmunuCR_2b = df_out_TopmunuCR_2b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,'pu_nPUVert':ep_pu_nPUVert,
                                                    'MET':ep_pfMetCorrPt,'Recoil':ep_WmunuRecoil ,'Wmass':ep_Wmunumass,'WpT':WpT_munu,
                                                    'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_TightEleTightMuVeto,'NEle':ep_nEle_index,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0], 'Jet1Eta':ep_THINjetEta[0], 'Jet1Phi':ep_THINjetPhi[0], 'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':ep_THINjetPt[1], 'Jet2Eta':ep_THINjetEta[1], 'Jet2Phi':ep_THINjetPhi[1], 'Jet2deepCSV':ep_THINjetDeepCSV[1],
                                                    'Jet3Pt':Jet3Pt, 'Jet3Eta':Jet3Eta, 'Jet3Phi':Jet3Phi, 'Jet3deepCSV':Jet3deepCSV,
                                                    'leadingLepPt':ep_muPt[0],'leadingLepEta':ep_muEta[0],'leadingLepPhi':ep_muPhi[0],
                                                    'weight':weight,'weightEle':weightEle,'weightMu':weightMu,'weightB':weightB, 'weightTau':weightTau, 'weightEWK':weightEWK, 'weightTop':weightTop, 'weightPU':weightPU
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

    cfsr_list = {1:'MET',2:'nLep',3:'min_dPhi',4:'nJet',5:'nBjets'}
    print ('\n============SR cutflow============')
    print ('cut_ep_pfMetCorrPt,cut_ep_nLep,cut_min_dPhi,cut_ep_THINnJet_1b,cut_ep_THINjetDeepCSV_1b')
    print (cut_ep_pfMetCorrPt,cut_ep_nLep,cut_min_dPhi,cut_ep_THINnJet_1b,cut_ep_THINjetDeepCSV_1b)
    print ('cut_ep_pfMetCorrPt,cut_ep_nLep,cut_min_dPhi,cut_ep_THINnJet_2b,cut_ep_THINjetDeepCSV_2b')
    print (cut_ep_pfMetCorrPt,cut_ep_nLep,cut_min_dPhi,cut_ep_THINnJet_2b,cut_ep_THINjetDeepCSV_2b)
    print ('SR1bcount',SR1bcount,'SR2bcount',SR2bcount)
    sr_1bdict = {1:cut_ep_pfMetCorrPt,2:cut_ep_nLep,3:cut_min_dPhi,4:cut_ep_THINnJet_1b,5:cut_ep_THINjetDeepCSV_1b}
    sr_2bdict = {1:cut_ep_pfMetCorrPt,2:cut_ep_nLep,3:cut_min_dPhi,4:cut_ep_THINnJet_2b,5:cut_ep_THINjetDeepCSV_2b}
    for i in [1,2,3,4,5]:
        h_reg_SR_1b_cutFlow.AddBinContent(i,sr_1bdict[i])
        h_reg_SR_2b_cutFlow.AddBinContent(i,sr_2bdict[i])
        h_reg_SR_1b_cutFlow.GetXaxis().SetBinLabel(i,cfsr_list[i])
        h_reg_SR_2b_cutFlow.GetXaxis().SetBinLabel(i,cfsr_list[i])

    print ('============SR cutflow============')

    cf_list = {1:'Recoil',2:'nLep',3:'min_dPhi',4:'Z/W mass',5:'nJet',6:'nBjets'}
    print ('\n============Zee cutflow============')
    print('Zee_cut_ep_Recoil,Zee_cut_ep_nLep,Zee_cut_min_dPhi,Zee_cut_ep_Zeemass,Zee_cut_ep_THINnJet_1b,Zee_cut_ep_THINjetDeepCSV_1b')
    print(Zee_cut_ep_Recoil,Zee_cut_ep_nLep,Zee_cut_min_dPhi,Zee_cut_ep_Zeemass,Zee_cut_ep_THINnJet_1b,Zee_cut_ep_THINjetDeepCSV_1b)
    print('Zee_cut_ep_Recoil,Zee_cut_ep_nLep,Zee_cut_min_dPhi,Zee_cut_ep_Zeemass,Zee_cut_ep_THINnJet_2b,Zee_cut_ep_THINjetDeepCSV_2b')
    print(Zee_cut_ep_Recoil,Zee_cut_ep_nLep,Zee_cut_min_dPhi,Zee_cut_ep_Zeemass,Zee_cut_ep_THINnJet_2b,Zee_cut_ep_THINjetDeepCSV_2b)
    print('ZeeCR1bcount',ZeeCR1bcount,'ZeeCR2bcount',ZeeCR2bcount)
    Zee_1bdict = {1:Zee_cut_ep_Recoil,2:Zee_cut_ep_nLep,3:Zee_cut_min_dPhi,4:Zee_cut_ep_Zeemass,5:Zee_cut_ep_THINnJet_1b,6:Zee_cut_ep_THINjetDeepCSV_1b}
    Zee_2bdict = {1:Zee_cut_ep_Recoil,2:Zee_cut_ep_nLep,3:Zee_cut_min_dPhi,4:Zee_cut_ep_Zeemass,5:Zee_cut_ep_THINnJet_2b,6:Zee_cut_ep_THINjetDeepCSV_2b}
    print ('============Zee cutflow============')

    print ('\n============Zmumu cutflow============')
    print('Zmumu_cut_ep_Recoil,Zmumu_cut_ep_nLep,Zmumu_cut_min_dPhi,Zmumu_cut_ep_Zmumumass,Zmumu_cut_ep_THINnJet_1b,Zmumu_cut_ep_THINjetDeepCSV_1b')
    print(Zmumu_cut_ep_Recoil,Zmumu_cut_ep_nLep,Zmumu_cut_min_dPhi,Zmumu_cut_ep_Zmumumass,Zmumu_cut_ep_THINnJet_1b,Zmumu_cut_ep_THINjetDeepCSV_1b)
    print('Zmumu_cut_ep_Recoil,Zmumu_cut_ep_nLep,Zmumu_cut_min_dPhi,Zmumu_cut_ep_Zmumumass,Zmumu_cut_ep_THINnJet_2b,Zmumu_cut_ep_THINjetDeepCSV_2b')
    print(Zmumu_cut_ep_Recoil,Zmumu_cut_ep_nLep,Zmumu_cut_min_dPhi,Zmumu_cut_ep_Zmumumass,Zmumu_cut_ep_THINnJet_2b,Zmumu_cut_ep_THINjetDeepCSV_2b)
    print ('ZmumuCR1count', ZmumuCR1count,'ZmumuCR2count', ZmumuCR2count)
    Zmumu_1bdict = {1:Zmumu_cut_ep_Recoil,2:Zmumu_cut_ep_nLep,3:Zmumu_cut_min_dPhi,4:Zmumu_cut_ep_Zmumumass,5:Zmumu_cut_ep_THINnJet_1b,6:Zmumu_cut_ep_THINjetDeepCSV_1b}
    Zmumu_2bdict = {1:Zmumu_cut_ep_Recoil,2:Zmumu_cut_ep_nLep,3:Zmumu_cut_min_dPhi,4:Zmumu_cut_ep_Zmumumass,5:Zmumu_cut_ep_THINnJet_2b,6:Zmumu_cut_ep_THINjetDeepCSV_2b}
    print ('============Zmumu cutflow============')

    print ('\n============Wenu cutflow============')
    print('Wenu_cut_ep_Recoil,Wenu_cut_ep_nLep,Wenu_cut_min_dPhi,Wenu_cut_ep_Wenumass,Wenu_cut_ep_THINnJet_1b,Wenu_cut_ep_THINjetDeepCSV_1b')
    print(Wenu_cut_ep_Recoil,Wenu_cut_ep_nLep,Wenu_cut_min_dPhi,Wenu_cut_ep_Wenumass,Wenu_cut_ep_THINnJet_1b,Wenu_cut_ep_THINjetDeepCSV_1b)
    print('Wenu_cut_ep_Recoil,Wenu_cut_ep_nLep,Wenu_cut_min_dPhi,Wenu_cut_ep_Wenumass,Wenu_cut_ep_THINnJet_2b,Wenu_cut_ep_THINjetDeepCSV_2b')
    print(Wenu_cut_ep_Recoil,Wenu_cut_ep_nLep,Wenu_cut_min_dPhi,Wenu_cut_ep_Wenumass,Wenu_cut_ep_THINnJet_2b,Wenu_cut_ep_THINjetDeepCSV_2b)
    print ('WenuCR1bcount', WenuCR1bcount,'WenuCR2bcount', WenuCR2bcount)
    Wenu_1bdict = {1:Wenu_cut_ep_Recoil,2:Wenu_cut_ep_nLep,3:Wenu_cut_min_dPhi,4:Wenu_cut_ep_Wenumass,5:Wenu_cut_ep_THINnJet_1b,6:Wenu_cut_ep_THINjetDeepCSV_1b}
    Wenu_2bdict = {1:Wenu_cut_ep_Recoil,2:Wenu_cut_ep_nLep,3:Wenu_cut_min_dPhi,4:Wenu_cut_ep_Wenumass,5:Wenu_cut_ep_THINnJet_2b,6:Wenu_cut_ep_THINjetDeepCSV_2b}
    print ('============Wenu cutflow============')

    print ('\n============Wmunu cutflow============')
    print('Wmunu_cut_ep_Recoil,Wmunu_cut_ep_nLep,Wmunu_cut_min_dPhi,Wmunu_cut_ep_Wmunumass,Wmunu_cut_ep_THINnJet_1b,Wmunu_cut_ep_THINjetDeepCSV_1b')
    print(Wmunu_cut_ep_Recoil,Wmunu_cut_ep_nLep,Wmunu_cut_min_dPhi,Wmunu_cut_ep_Wmunumass,Wmunu_cut_ep_THINnJet_1b,Wmunu_cut_ep_THINjetDeepCSV_1b)
    print('Wmunu_cut_ep_Recoil,Wmunu_cut_ep_nLep,Wmunu_cut_min_dPhi,Wmunu_cut_ep_Wmunumass,Wmunu_cut_ep_THINnJet_2b,Wmunu_cut_ep_THINjetDeepCSV_2b')
    print(Wmunu_cut_ep_Recoil,Wmunu_cut_ep_nLep,Wmunu_cut_min_dPhi,Wmunu_cut_ep_Wmunumass,Wmunu_cut_ep_THINnJet_2b,Wmunu_cut_ep_THINjetDeepCSV_2b)
    print ('WmunuCR1bcount', WmunuCR1bcount,'WmunuCR2bcount', WmunuCR2bcount)
    Wmunu_1bdict = {1:Wmunu_cut_ep_Recoil,2:Wmunu_cut_ep_nLep,3:Wmunu_cut_min_dPhi,4:Wmunu_cut_ep_Wmunumass,5:Wmunu_cut_ep_THINnJet_1b,6:Wmunu_cut_ep_THINjetDeepCSV_1b}
    Wmunu_2bdict = {1:Wmunu_cut_ep_Recoil,2:Wmunu_cut_ep_nLep,3:Wmunu_cut_min_dPhi,4:Wmunu_cut_ep_Wmunumass,5:Wmunu_cut_ep_THINnJet_2b,6:Wmunu_cut_ep_THINjetDeepCSV_2b}
    print ('============Wmunu cutflow============')

    print ('\n============Topenu cutflow============')
    print('Topenu_cut_ep_Recoil,Topenu_cut_ep_nLep,Topenu_cut_min_dPhi,Topenu_cut_ep_Topenumass,Topenu_cut_ep_THINnJet_1b,Topenu_cut_ep_THINjetDeepCSV_1b')
    print(Topenu_cut_ep_Recoil,Topenu_cut_ep_nLep,Topenu_cut_min_dPhi,Topenu_cut_ep_Topenumass,Topenu_cut_ep_THINnJet_1b,Topenu_cut_ep_THINjetDeepCSV_1b)
    print('Topenu_cut_ep_Recoil,Topenu_cut_ep_nLep,Topenu_cut_min_dPhi,Topenu_cut_ep_Topenumass,Topenu_cut_ep_THINnJet_2b,Topenu_cut_ep_THINjetDeepCSV_2b')
    print(Topenu_cut_ep_Recoil,Topenu_cut_ep_nLep,Topenu_cut_min_dPhi,Topenu_cut_ep_Topenumass,Topenu_cut_ep_THINnJet_2b,Topenu_cut_ep_THINjetDeepCSV_2b)
    print ('TopenuCR1bcount', TopenuCR1bcount,'TopenuCR2bcount', TopenuCR2bcount)
    Topenu_1bdict = {1:Topenu_cut_ep_Recoil,2:Topenu_cut_ep_nLep,3:Topenu_cut_min_dPhi,4:Topenu_cut_ep_Topenumass,5:Topenu_cut_ep_THINnJet_1b,6:Topenu_cut_ep_THINjetDeepCSV_1b}
    Topenu_2bdict = {1:Topenu_cut_ep_Recoil,2:Topenu_cut_ep_nLep,3:Topenu_cut_min_dPhi,4:Topenu_cut_ep_Topenumass,5:Topenu_cut_ep_THINnJet_2b,6:Topenu_cut_ep_THINjetDeepCSV_2b}
    print ('============Topenu cutflow============')

    print ('\n============Topmunu cutflow============')
    print('Topmunu_cut_ep_Recoil,Topmunu_cut_ep_nLep,Topmunu_cut_min_dPhi,Topmunu_cut_ep_Topmunumass,Topmunu_cut_ep_THINnJet_1b,Topmunu_cut_ep_THINjetDeepCSV_1b')
    print(Topmunu_cut_ep_Recoil,Topmunu_cut_ep_nLep,Topmunu_cut_min_dPhi,Topmunu_cut_ep_Topmunumass,Topmunu_cut_ep_THINnJet_1b,Topmunu_cut_ep_THINjetDeepCSV_1b)
    print('Topmunu_cut_ep_Recoil,Topmunu_cut_ep_nLep,Topmunu_cut_min_dPhi,Topmunu_cut_ep_Topmunumass,Topmunu_cut_ep_THINnJet_2b,Topmunu_cut_ep_THINjetDeepCSV_2b')
    print(Topmunu_cut_ep_Recoil,Topmunu_cut_ep_nLep,Topmunu_cut_min_dPhi,Topmunu_cut_ep_Topmunumass,Topmunu_cut_ep_THINnJet_2b,Topmunu_cut_ep_THINjetDeepCSV_2b)
    print ('TopmunuCR1bcount', TopmunuCR1bcount,'TopmunuCR2bcount', TopmunuCR2bcount)
    Topmunu_1bdict = {1:Topmunu_cut_ep_Recoil,2:Topmunu_cut_ep_nLep,3:Topmunu_cut_min_dPhi,4:Topmunu_cut_ep_Topmunumass,5:Topmunu_cut_ep_THINnJet_1b,6:Topmunu_cut_ep_THINjetDeepCSV_1b}
    Topmunu_2bdict = {1:Topmunu_cut_ep_Recoil,2:Topmunu_cut_ep_nLep,3:Topmunu_cut_min_dPhi,4:Topmunu_cut_ep_Topmunumass,5:Topmunu_cut_ep_THINnJet_2b,6:Topmunu_cut_ep_THINjetDeepCSV_2b}
    print ('============Topmunu cutflow============\n')

    for i in [1,2,3,4,5,6]:
        h_reg_ZeeCR_1b_cutFlow.AddBinContent(i,Zee_1bdict[i] )
        h_reg_ZeeCR_2b_cutFlow.AddBinContent(i,Zee_2bdict[i] )
        h_reg_ZmumuCR_1b_cutFlow.AddBinContent(i,Zmumu_1bdict[i] )
        h_reg_ZmumuCR_2b_cutFlow.AddBinContent(i,Zmumu_2bdict[i] )
        h_reg_WenuCR_1b_cutFlow.AddBinContent(i,Wenu_1bdict[i] )
        h_reg_WenuCR_2b_cutFlow.AddBinContent(i,Wenu_2bdict[i] )
        h_reg_WmunuCR_1b_cutFlow.AddBinContent(i,Wmunu_1bdict[i] )
        h_reg_WmunuCR_2b_cutFlow.AddBinContent(i,Wmunu_2bdict[i] )
        h_reg_TopenuCR_1b_cutFlow.AddBinContent(i,Topenu_1bdict[i] )
        h_reg_TopenuCR_2b_cutFlow.AddBinContent(i,Topenu_2bdict[i] )
        h_reg_TopmunuCR_1b_cutFlow.AddBinContent(i,Topmunu_1bdict[i] )
        h_reg_TopmunuCR_2b_cutFlow.AddBinContent(i,Topmunu_2bdict[i] )

        h_reg_ZeeCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_ZeeCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_ZmumuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_ZmumuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_WenuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_WenuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_WmunuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_WmunuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_TopenuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_TopenuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_TopmunuCR_1b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])
        h_reg_TopmunuCR_2b_cutFlow.GetXaxis().SetBinLabel(i,cf_list[i])

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
        #inputpath= "/afs/cern.ch/work/p/ptiwari/bb+DM_analysis/ntuple_analysis/CMSSW_10_3_0/src/ExoPieSlimmer/SIG_2016_2HDMa_SkimRootFilesALL"
        #inputpath= "/eos/cms/store/group/phys_exotica/bbMET/2017_skimmedFiles/V0/MC_USCM_25Sep"
        #inputpath= "/afs/cern.ch/work/p/ptiwari/bb+DM_analysis/ntuple_analysis/CMSSW_10_3_0/src/ExoPieProducer/ExoPieAnalyzer/CondorJobs_v2/Filelists_v1"
        inputpath="/afs/cern.ch/work/p/ptiwari/bb+DM_analysis/ntuple_analysis/CMSSW_10_3_0/src/ExoPieProducer/ExoPieAnalyzer/test_rootFile/test_root_file"

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
