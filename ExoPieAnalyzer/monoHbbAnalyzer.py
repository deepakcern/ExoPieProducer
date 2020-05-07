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
#import eventSelector
## for parallel threads in interactive running
from multiprocessing import Process
import multiprocessing as mp
from os import sys

isCondor = True
runInteractive = False
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
#import variables as var
import variables as var
import outvars as out
from cutFlow import cutFlow
#import eventSelector as eventSelector_v2


## from analysisutils
if isCondor:sys.path.append('ExoPieUtils/scalefactortools/')
else:sys.path.append('../../ExoPieUtils/scalefactortools/')
year_file= open("Year.py","w")
year_file.write('era="2018"')
year_file.close()

import ana_weight as wgt


import eventSelector as eventSelector_v2
from eventSelector import getJECWeight 

######################################################################################################
## All import are done before this
######################################################################################################

## ----- start of clock
start = time.clock()


year_file= open("Year.py","w")
runOn2017=False
runOn2018=False
## ----- command line argument
usage = "analyzer for bb+DM (debugging) "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--inputfile",  dest="inputfile",default="myfiles.txt")
parser.add_argument("-inDir", "--inputDir",  dest="inputDir",default=".")
parser.add_argument("-runOnTXT", "--runOnTXT",action="store_true", dest="runOnTXT")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="out.root")
parser.add_argument("-D", "--outputdir", dest="outputdir")
parser.add_argument("-F", "--farmout", action="store_true",  dest="farmout")
parser.add_argument("-y", "--year", dest="year", default="Year")

args = parser.parse_args()

if args.farmout==None:
    isfarmout = False
else:
    isfarmout = args.farmout

if args.inputDir and isfarmout:
    dirName=args.inputDir

runOnTxt=False
if args.runOnTXT:
    runOnTxt = True

if args.year=='2017':
    runOn2017=True
    year_file.write('era="2017"')
elif args.year=='2018':
    runOn2018=True 
    year_file.write('era="2018"')

else:
    print 'please provide year'
    sys.exit()

year_file.close()

if isfarmout:
    infile  = args.inputfile

else: print "No file is provided for farmout"


outputdir = '.'
if args.outputdir:
    outputdir = str(args.outputdir)

infilename = "NCUGlobalTuples.root"

debug_ = False

outDir=outputdir

def TextToList(textfile):
    return([iline.rstrip()    for iline in open(textfile)])

## the input file list and key is caught in one variable as a python list,
#### first element is the list of rootfiles
#### second element is the key, user to name output.root

dummy = -9999.0

if runOn2017:
    LWP = 0.1522
    MWP = 0.4941
else:
    LWP = 0.1241
    MWP = 0.4184

def runbbdm(txtfile):
    count = 0


    print "in main function"

    infile_=[]
    outfilename=""
    prefix="Skimmed_"
    ikey_ = ""

    '''
    if  runInteractive:
        print "running for ", txtfile
        infile_  = TextToList(txtfile[0])
        key_=txtfile[1]
        outfilename= prefix+key_+".root"
    '''


    if  runInteractive:
        print "running for ", txtfile
        infile_  = TextToList(txtfile)
        #key_=txtfile[1]

        ''' old
        prefix="Skimmed_"
        outfilename= prefix+infile_.split("/")[-1]
        '''

        outfilename= outDir+'/'+txtfile.split('/')[-1].replace('.txt','.root')#prefix+key_+".root"

    print 'txtfile', txtfile
    if not runInteractive:
	infile_=TextToList(txtfile)
        prefix_ = '' #'/eos/cms/store/group/phys_exotica/bbMET/2017_skimmedFiles/locallygenerated/'
        if outputdir!='.': prefix_ = outputdir+'/'
        print "prefix_", prefix_
        outfilename = prefix_+txtfile.split('/')[-1].replace('.txt','.root')#"SkimmedTree.root"
        print 'outfilename',  outfilename


    ## define global dataframes
    df_out_SR_resolved  = out.df_out_SR_resolved
    df_out_SR_boosted   = out.df_out_SR_boosted

    df_out_SBand_resolved = out.df_out_SBand_resolved
    df_out_Tope_resolved  = out.df_out_Tope_resolved
    df_out_Topmu_resolved = out.df_out_Topmu_resolved
    df_out_We_resolved    = out.df_out_We_resolved
    df_out_Wmu_resolved   = out.df_out_Wmu_resolved
    df_out_Zmumu_resolved = out.df_out_Zmumu_resolved
    df_out_Zee_resolved   = out.df_out_Zee_resolved


    df_out_SBand_boosted= out.df_out_SBand_boosted
    df_out_Tope_boosted = out.df_out_Tope_boosted
    df_out_Topmu_boosted= out.df_out_Topmu_boosted
    df_out_We_boosted   = out.df_out_We_boosted
    df_out_Wmu_boosted  = out.df_out_Wmu_boosted
    df_out_Zmumu_boosted= out.df_out_Zmumu_boosted
    df_out_Zee_boosted  = out.df_out_Zee_boosted
    df_out_TopWmu_boosted = out.df_out_TopWmu_boosted
    df_out_TopWe_boosted  = out.df_out_TopWe_boosted


    #outputfilename = args.outputfile
    Entrees = 0.0
    h_total = TH1F('h_total','h_total',2,0,2)
    h_total_mcweight = TH1F('h_total_mcweight','h_total_mcweight',2,0,2)
    for infl in infile_:
	f_tmp = TFile.Open(infl,'READ')
        h_tmp = f_tmp.Get('h_total')
        Stree = f_tmp.Get('outTree')
        Entrees += Stree.GetEntries()
        h_tmp_weight = f_tmp.Get('h_total_mcweight')
        h_total.Add(h_tmp)
        h_total_mcweight.Add(h_tmp_weight)


    total_events = h_total_mcweight.Integral()
    h_reg_WenuCR_resolved_cutFlow = TH1F("h_reg_WenuCR_resolved_cutFlow","h_reg_WenuCR_resolved_cutFlow",12,0,12)
    h_reg_WmunuCR_resolved_cutFlow = TH1F("h_reg_WmunuCR_resolved_cutFlow","h_reg_WmunuCR_resolved_cutFlow",12,0,12)
    
    h_reg_TopenuCR_resolved_cutFlow = TH1F("h_reg_TopenuCR_resolved_cutFlow","h_reg_TopenuCR_resolved_cutFlow",12,0,12)
    h_reg_TopmunuCR_resolved_cutFlow = TH1F("h_reg_TopmunuCR_resolved_cutFlow","h_reg_TopmunuCR_resolved_cutFlow",12,0,12)

    h_reg_ZeeCR_resolved_cutFlow = TH1F("h_reg_ZeeCR_resolved_cutFlow","h_reg_ZeeCR_resolved_cutFlow",12,0,12)
    h_reg_ZmumuCR_resolved_cutFlow = TH1F("h_reg_ZmumuCR_resolved_cutFlow","h_reg_ZmumuCR_resolved_cutFlow",12,0,12)


    h_reg_WenuCR_boosted_cutFlow = TH1F("h_reg_WenuCR_boosted_cutFlow","h_reg_WenuCR_boosted_cutFlow",12,0,12)
    h_reg_WmunuCR_boosted_cutFlow = TH1F("h_reg_WmunuCR_boosted_cutFlow","h_reg_WmunuCR_boosted_cutFlow",12,0,12)

    h_reg_TopenuCR_boosted_cutFlow = TH1F("h_reg_TopenuCR_boosted_cutFlow","h_reg_TopenuCR_boosted_cutFlow",12,0,12)
    h_reg_TopmunuCR_boosted_cutFlow = TH1F("h_reg_TopmunuCR_boosted_cutFlow","h_reg_TopmunuCR_boosted_cutFlow",12,0,12)

    h_reg_ZeeCR_boosted_cutFlow = TH1F("h_reg_ZeeCR_boosted_cutFlow","h_reg_ZeeCR_boosted_cutFlow",12,0,12)
    h_reg_ZmumuCR_boosted_cutFlow = TH1F("h_reg_ZmumuCR_boosted_cutFlow","h_reg_ZmumuCR_boosted_cutFlow",12,0,12)

    h_reg_SBand_resolved_cutFlow     = TH1F("h_reg_SBand_resolved_cutFlow","h_reg_SBand_resolved_cutFlow",12,0,12)
    h_reg_SBand_boosted_cutFlow     = TH1F("h_reg_SBand_boosted_cutFlow","h_reg_SBand_boosted_cutFlow",12,0,12)

    h_reg_WenuCR_resolved_cutFlow.AddBinContent(1,total_events)
    h_reg_WmunuCR_resolved_cutFlow.AddBinContent(1,total_events)    
    h_reg_WenuCR_resolved_cutFlow.AddBinContent(2,Entrees)
    h_reg_WmunuCR_resolved_cutFlow.AddBinContent(2,Entrees)

    h_reg_ZmumuCR_resolved_cutFlow.AddBinContent(2,total_events)
    h_reg_ZeeCR_resolved_cutFlow.AddBinContent(2,total_events)
    h_reg_ZmumuCR_resolved_cutFlow.AddBinContent(2,Entrees)
    h_reg_ZeeCR_resolved_cutFlow.AddBinContent(2,Entrees)

    h_reg_TopenuCR_resolved_cutFlow.AddBinContent(1,total_events)
    h_reg_TopmunuCR_resolved_cutFlow.AddBinContent(1,total_events)
    h_reg_TopenuCR_resolved_cutFlow.AddBinContent(2,Entrees)
    h_reg_TopmunuCR_resolved_cutFlow.AddBinContent(2,Entrees)

    h_reg_WenuCR_boosted_cutFlow.AddBinContent(1,total_events)
    h_reg_WmunuCR_boosted_cutFlow.AddBinContent(1,total_events)
    h_reg_WenuCR_boosted_cutFlow.AddBinContent(2,Entrees)
    h_reg_WmunuCR_boosted_cutFlow.AddBinContent(2,Entrees)

    h_reg_ZmumuCR_boosted_cutFlow.AddBinContent(2,total_events)
    h_reg_ZeeCR_boosted_cutFlow.AddBinContent(2,total_events)
    h_reg_ZmumuCR_boosted_cutFlow.AddBinContent(2,Entrees)
    h_reg_ZeeCR_boosted_cutFlow.AddBinContent(2,Entrees)

    h_reg_TopenuCR_boosted_cutFlow.AddBinContent(1,total_events)
    h_reg_TopmunuCR_boosted_cutFlow.AddBinContent(1,total_events)
    h_reg_TopenuCR_boosted_cutFlow.AddBinContent(2,Entrees)
    h_reg_TopmunuCR_boosted_cutFlow.AddBinContent(2,Entrees)


    h_reg_SBand_resolved_cutFlow.AddBinContent(1,total_events)
    h_reg_SBand_resolved_cutFlow.AddBinContent(2,Entrees)

    h_reg_SBand_boosted_cutFlow.AddBinContent(1,total_events)
    h_reg_SBand_boosted_cutFlow.AddBinContent(2,Entrees)



    filename = infile_
    ieve = 0;icount = 0
    for df in read_root(filename, 'outTree', columns=var.allvars, chunksize=125000):

        for ep_runId, ep_lumiSection, ep_eventId, \
            ep_pfMetCorrPt, ep_pfMetCorrPhi, ep_pfMetUncJetResUp, ep_pfMetUncJetResDown, ep_pfMetUncJetEnUp, ep_pfMetUncJetEnDown, \
            ep_THINnJet, ep_THINjetPx, ep_THINjetPy, ep_THINjetPz, ep_THINjetEnergy, \
            ep_THINjetDeepCSV, ep_THINjetHadronFlavor, \
            ep_THINjetNHadEF, ep_THINjetCHadEF, ep_THINjetCEmEF, \
            ep_THINjetCorrUnc,ep_RegNNCorr,\
            ep_nfjet, ep_fjetPx, ep_fjetPy, ep_fjetPz, ep_fjetEnergy, \
            ep_fjetDoubleSV, ep_fjetProbQCDb, ep_fjetProbHbb, ep_fjetProbQCDc, ep_fjetProbHcc, ep_fjetProbHbbc, ep_fjetProbbbvsLight, \
            ep_fjetProbccvsLight, ep_fjetProbTvsQCD, ep_fjetProbWvsQCD, ep_fjetProbZHbbvsQCD, \
            ep_fjetSDMass,ep_SDMassCorr, ep_fjetN2b1, ep_fjetN2b2, ep_fjetCHSPRMass, ep_fjetCHSSDMass, \
            ep_nEle, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy, \
            ep_eleIsPasepight, ep_eleIsPassLoose,ep_eleCharge, \
            ep_nPho, ep_phoIsPasepight, ep_phoPx, ep_phoPy, ep_phoPz, ep_phoEnergy, \
            ep_nMu, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy, ep_iepightMuon, ep_muCharge,\
            ep_nTau_DRBased_EleMuVeto,ep_nTau_discBased_looseElelooseMuVeto,ep_nTau_discBased_looseEleTightMuVeto,ep_nTau_discBased_looseEleTightMuVeto,ep_nTau_discBased_mediumElelooseMuVeto,ep_nTau_discBased_TightEleTightMuVeto,\
            ep_pu_nTrueInt, ep_pu_nPUVert, \
            ep_THINjetNPV, \
            ep_mcweight, ep_genParPt,ep_genParSample, isData, eletrigdecision, mutrigdecision, mettrigdecision,\
            in zip(df.st_runId, df.st_lumiSection, df.st_eventId, \
                   df.st_pfMetCorrPt, df.st_pfMetCorrPhi, df.st_pfMetUncJetResUp, df.st_pfMetUncJetResDown, df.st_pfMetUncJetEnUp, df.st_pfMetUncJetEnDown, \
                   df.st_THINnJet, df.st_THINjetPx, df.st_THINjetPy, df.st_THINjetPz, df.st_THINjetEnergy, \
                   df.st_THINjetDeepCSV, df.st_THINjetHadronFlavor, \
                   df.st_THINjetNHadEF, df.st_THINjetCHadEF, df.st_THINjetCEmEF, \
                   df.st_THINjetCorrUnc,df.st_THINbRegNNCorr,\
                   df.st_nfjet, df.st_fjetPx, df.st_fjetPy, df.st_fjetPz, df.st_fjetEnergy, \
                   df.st_fjetDoubleSV, df.st_fjetProbQCDb, df.st_fjetProbHbb, df.st_fjetProbQCDc, df.st_fjetProbHcc, df.st_fjetProbHbbc, df.st_fjetProbbbvsLight, \
                   df.st_fjetProbccvsLight, df.st_fjetProbTvsQCD, df.st_fjetProbWvsQCD, df.st_fjetProbZHbbvsQCD, \
                   df.st_fjetSDMass,df.st_fjetSDMassCorrFact,df.st_fjetN2b1, df.st_fjetN2b2, df.st_fjetCHSPRMass, df.st_fjetCHSSDMass, \
                   df.st_nEle, df.st_elePx, df.st_elePy, df.st_elePz, df.st_eleEnergy, \
                   df.st_eleIsPassTight, df.st_eleIsPassLoose,df.st_eleCharge, \
                   df.st_nPho, df.st_phoIsPassTight, df.st_phoPx, df.st_phoPy, df.st_phoPz, df.st_phoEnergy, \
                   df.st_nMu, df.st_muPx, df.st_muPy, df.st_muPz, df.st_muEnergy, df.st_isTightMuon,df.st_muCharge, \
                   df.st_nTau_DRBased_EleMuVeto,df.st_nTau_discBased_looseElelooseMuVeto,df.st_nTau_discBased_looseEleTightMuVeto,df.st_nTau_discBased_looseEleTightMuVeto,df.st_nTau_discBased_mediumElelooseMuVeto,df.st_nTau_discBased_TightEleTightMuVeto,\
                   df.st_pu_nTrueInt, df.st_pu_nPUVert, \
                   df.st_THINjetNPV, \
                   df.mcweight, df.st_genParPt, df.st_genParSample,df.st_isData,df.st_eletrigdecision,df.st_mutrigdecision,df.st_mettrigdecision, \

            ):


            ieve = ieve + 1
            if ieve%10000==0: print "Processed",ieve,"Events"

            isBoostedSR=False
            isBoostedSBand=False
            isBoostedCRWenu=False
            isBoostedCRWmunu=False
            isBoostedCRZee=False
            isBoostedCRZmumu=False
            isBoostedCRTope=False
            isBoostedCRTopmu=False
            isBoostedSB=False

            isResolvedSR=False
	    isResolvedSBand=False
            isResolvedCRWenu=False
            isResolvedCRWmunu=False
            isResolvedCRZee=False
            isResolvedCRZmumu=False
            isResolvedCRTope=False
            isResolvedCRTopmu=False
            isResolvedSB=False


            '''
            -------------------------------------------------------------------------------
            FAT JET COLLECTION
            -------------------------------------------------------------------------------
            '''


            fatjetpt = [getPt(ep_fjetPx[ij], ep_fjetPy[ij]) for ij in range(ep_nfjet)]
            fatjeteta = [getEta(ep_fjetPx[ij], ep_fjetPy[ij], ep_fjetPz[ij]) for ij in range(ep_nfjet)]
            fatjetphi = [getPhi(ep_fjetPx[ij], ep_fjetPy[ij]) for ij in range(ep_nfjet)]


            pass_nfjetIndex = [index for index in range(ep_nfjet) if ((fatjetpt[index] > 200.0) and (abs(fatjeteta[index])< 2.5) and (ep_fjetSDMass[index] > 100.0) and (ep_fjetSDMass[index] < 150.0) and (ep_fjetProbHbb[index] > 0.86)) ]

            FatJet_SBand_index = [index for index in range(ep_nfjet) if ((fatjetpt[index] > 200.0) and (abs(fatjeteta[index])< 2.5)) and ((ep_fjetSDMass[index] > 50.0) and (ep_fjetSDMass[index] < 100.0) or ((ep_fjetSDMass[index] > 150.0) and (ep_fjetSDMass[index] < 350.0) )) and (ep_fjetProbHbb[index] > 0.86)]

            FatJet_ZCR_index   = [index for index in range(ep_nfjet) if ((fatjetpt[index] > 200.0) and (abs(fatjeteta[index])< 2.5) and (ep_fjetProbHbb[index] > 0.86))]

            nFatJet_SBand  = len(FatJet_SBand_index)
            nFatJet_ZCR    = len(FatJet_ZCR_index)

            sel_nfjets = len(pass_nfjetIndex)
            if (sel_nfjets==1):
		fjetIndex = pass_nfjetIndex[0]


            '''
            ------------------------------------------------------------------------------
            AK4JET COLLECTION FOR BOOSTED CATEGORY
            -----------------------------------------------------------------------------
            '''
 
       
            ak4jetpt  = [getPt(ep_THINjetPx[ij], ep_THINjetPy[ij]) for ij in range(ep_THINnJet)]
            ak4jeteta = [getEta(ep_THINjetPx[ij], ep_THINjetPy[ij], ep_THINjetPz[ij]) for ij in range(ep_THINnJet)]
            ak4jetphi = [getPhi(ep_THINjetPx[ij], ep_THINjetPy[ij]) for ij in range(ep_THINnJet)]

            nBjets_notiso_index = [ij for ij in range(ep_THINnJet) if (ep_THINjetDeepCSV[ij] > MWP and abs(ak4jeteta[ij]) < 2.5)]
            nBjets_notiso = len(nBjets_notiso_index)

            ak4_pt30_eta4p5  = [True for ij in range(ep_THINnJet)] #pt > 30 and eta < 4.5 is already applied at skimmer level
            fjet_pt200_eta_2p5 = [True for ij in range(ep_nfjet)] #pt > 200 and eta < 2.5 is already applied at skimmer level
            pass_ak4jet_index_cleaned = []
            if len(ak4_pt30_eta4p5) > 0 and len(fjet_pt200_eta_2p5) > 0 :
		ak4jet_cleaned_against_fjet = anautil.jetcleaning(ak4_pt30_eta4p5, fjet_pt200_eta_2p5, ak4jeteta, fatjeteta, ak4jetphi, fatjetphi,0.8)
                pass_ak4jet_index_cleaned = boolutil.WhereIsTrue(ak4jet_cleaned_against_fjet)
            #print 'pass_ak4jet_index_cleaned', pass_ak4jet_index_cleaned

	    nJets_cleaned = len(pass_ak4jet_index_cleaned)

            Bjet_index = [ij for ij in pass_ak4jet_index_cleaned if (ep_THINjetDeepCSV[ij] > LWP and abs(ak4jeteta[ij]) < 2.5)]
            nBjets_iso = len(Bjet_index)

            '''
            ------------------------------------------------------------------------------
            #AK4JETS COLLECTION/WORK  FOR RESOLVED CATEGORY
            -----------------------------------------------------------------------------
            '''

            jet1Index_list = []
            jet2Index_list = []

            for firstjet in range(ep_THINnJet):
		for secondjet in range(ep_THINnJet):
		    if (firstjet<secondjet) and (ep_THINjetDeepCSV[firstjet] > MWP and (abs(ak4jeteta[firstjet]) < 2.5) ) and (ep_THINjetDeepCSV[secondjet] > MWP and (abs(ak4jeteta[secondjet]) < 2.5)):
		        jet1Index_list.append(firstjet)
                        jet2Index_list.append(secondjet)
                    else:continue

            
            h_mass= -9999.0

            #nBjets = 0
            if nBjets_notiso==2: 
		jet1Index=jet1Index_list[0]
                jet2Index=jet2Index_list[0]
		h_mass  = InvMass(ep_THINjetPx[jet1Index], ep_THINjetPy[jet1Index], ep_THINjetPz[jet1Index], ep_THINjetEnergy[jet1Index],ep_THINjetPx[jet2Index], ep_THINjetPy[jet2Index], ep_THINjetPz[jet2Index],ep_THINjetEnergy[jet2Index])
                #print 'jet1Index',jet1Index,'jet2Index',jet2Index
                #print 'jet1Index_list',jet1Index_list,'jet2Index_list',jet2Index_list
                #print 'nBjets_notiso_index',nBjets_notiso_index,'jet1Index',jet1Index,'jet2Index',jet2Index

            #print 'h_mass',h_mass

	    '''
            ------------------------------------------------------------------------------
            LEPTON COLLECTIONS
            ------------------------------------------------------------------------------
            '''
            #ep_HPSTau_n = ep_nTau_discBased_looseElelooseMuVeto
            ep_HPSTau_n = ep_nTau_DRBased_EleMuVeto

            elephi = [getPhi(ep_elePx[ij], ep_elePy[ij]) for ij in range(ep_nEle)]
            elept = [getPt(ep_elePx[ij], ep_elePy[ij]) for ij in range(ep_nEle)]
            eleeta = [getEta(ep_elePx[ij], ep_elePy[ij], ep_elePz[ij]) for ij in range(ep_nEle)]

            ele_loose_index = [index for index in range(ep_nEle) if (ep_eleIsPassLoose[index])]
            ele_tight_index = [index for index in range(ep_nEle) if (ep_eleIsPasepight[index]) and getPt(ep_elePx[index], ep_elePy[index]) > 40]
            isTightEles = [(ep_eleIsPasepight[ij]) and (elept[ij] > 40) for ij in range(ep_nEle)]

            '''
            elephi_tight = [getPhi(ep_elePx[ij], ep_elePy[ij]) for ij in range(ep_nEle) if ((ep_eleIsPasepight[ij]) and getPt(ep_elePx[ij], ep_elePy[ij]) > 40) ]
            elept_tight = [getPt(ep_elePx[ij], ep_elePy[ij]) for ij in range(ep_nEle) if (ep_eleIsPasepight[ij] and getPt(ep_elePx[ij], ep_elePy[ij]) > 40)]
            eleeta_tight = [getEta(ep_elePx[ij], ep_elePy[ij], ep_elePz[ij]) for ij in range(ep_nEle) if (ep_eleIsPasepight[ij])]
            '''

            muphi  = [getPhi(ep_muPx[ij],ep_muPy[ij]) for ij in range(ep_nMu)]
            mupt   = [getPt(ep_muPx[ij],ep_muPy[ij]) for ij in range(ep_nMu) ]
            mueta = [getEta(ep_muPx[ij], ep_muPy[ij], ep_muPz[ij]) for ij in range(ep_nMu)]
            muon_tight_index = [index for index in range(ep_nMu) if (ep_iepightMuon[index] and (mupt[index] > 30) )]
            isTightMuons = [(ep_iepightMuon[ij]) and (mupt[ij] > 30) for ij in range(ep_nMu)]

            '''
	    muphi_tight  = [getPhi(ep_muPx[ij],ep_muPy[ij]) for ij in range(ep_nMu) if (ep_iepightMuon[ij]) and (mupt[ij] > 20)]
            mupt_tight   = [getPt(ep_muPx[ij],ep_muPy[ij]) for ij in range(ep_nMu) if (ep_iepightMuon[ij]) and (mupt[ij] > 20)]
            '''

            nEle_loose = len(ele_loose_index)
            nTightEle  = len(ele_tight_index)
            nLooseMu   = ep_nMu
            nTightMu   = len(muon_tight_index)

            #if nEle_loose==2 or nLooseMu==2: print 'nTightEle',nTightEle, 'isTightEles',isTightEles,'nTightMu',nTightMu,'isTightMuons',isTightMuons

            '''
            --------------------------------------------------------------------------
            CLEANED PHOTON COLLECATION
            -------------------------------------------------------------------------
            '''
	    #skip photon veto for now
            pho_pt  = [getPt(ep_phoPx[ij], ep_phoPy[ij]) for ij in range(ep_nPho)]
            pho_phi = [getPhi(ep_phoPx[ij], ep_phoPy[ij]) for ij in range(ep_nPho)]
            pho_eta = [getEta(ep_phoPx[ij], ep_phoPy[ij], ep_phoPz[ij]) for ij in range(ep_nPho)]
            myphotons = [True for ij in range(ep_nPho)]
            myeleBooleans = [True for ij in range(ep_nEle)]
            mymuBooleans = [True for ij in range(ep_nMu)]
            cleanedPho_ag_ele = []; cleanedPho_ag_mu = [];pass_pho_index_cleaned=[]
	    if ep_nPho > 0: #and ep_nEle > 0:
                cleanedPho_ag_ele = anautil.jetcleaning(myphotons, myeleBooleans, pho_eta, eleeta, pho_phi, elephi, 0.4)
		cleanedPho_ag_mu  = anautil.jetcleaning(myphotons, mymuBooleans, pho_eta, mueta, pho_phi, muphi, 0.4) 
                cleanedPhoton     = boolutil.logical_AND_List2(cleanedPho_ag_ele,cleanedPho_ag_mu)
                pass_pho_index_cleaned = boolutil.WhereIsTrue(cleanedPhoton)
                #print 'cleanedPho_ag_ele',cleanedPho_ag_ele, 'cleanedPho_ag_mu', cleanedPho_ag_mu


            nPho = len(pass_pho_index_cleaned)




            '''
            ------------------------------------------------------------------------------
            HADRONIC RECOIL
            ------------------------------------------------------------------------------
            '''

            #======   usage: WRecoil_Phi_Wmass(nEle,elept,elephi,elepx_,elepy_,met_,metphi_) =======
	    Werecoil, WerecoildPhi, WeMass = eventSelector_v2.WRecoil_Phi_Wmass(nEle_loose,elept,elephi,ep_elePx,ep_elePy,ep_pfMetCorrPt, ep_pfMetCorrPhi)
            WerecoilResUp, WerecoildPhiResUp, WeMassResUp = eventSelector_v2.WRecoil_Phi_Wmass(nEle_loose,elept,elephi,ep_elePx,ep_elePy,ep_pfMetUncJetResUp, ep_pfMetCorrPhi)
            WerecoilResDown, WerecoildPhiResDown, WeMassResDown = eventSelector_v2.WRecoil_Phi_Wmass(nEle_loose,elept,elephi,ep_elePx,ep_elePy,ep_pfMetUncJetResDown, ep_pfMetCorrPhi)
            WerecoilEnUp, WerecoildPhiEnUp, WeMassEnUp = eventSelector_v2.WRecoil_Phi_Wmass(nEle_loose,elept,elephi,ep_elePx,ep_elePy,ep_pfMetUncJetEnUp, ep_pfMetCorrPhi)
            WerecoilEnDown, WerecoildPhiEnDown, WeMassEnDown = eventSelector_v2.WRecoil_Phi_Wmass(nEle_loose,elept,elephi,ep_elePx,ep_elePy,ep_pfMetUncJetEnDown, ep_pfMetCorrPhi)

            Wmurecoil, WmurecoildPhi, WmuMass = eventSelector_v2.WRecoil_Phi_Wmass(ep_nMu, mupt, muphi, ep_muPx, ep_muPy,ep_pfMetCorrPt, ep_pfMetCorrPhi)
            WmurecoilResUp, WmurecoildPhiResUp, WmuMassResUp = eventSelector_v2.WRecoil_Phi_Wmass(ep_nMu, mupt, muphi, ep_muPx, ep_muPy,ep_pfMetUncJetResUp, ep_pfMetCorrPhi)
            WmurecoilResDown, WmurecoildPhiResDown, WmuMassResDown = eventSelector_v2.WRecoil_Phi_Wmass(ep_nMu, mupt, muphi, ep_muPx, ep_muPy,ep_pfMetUncJetResDown, ep_pfMetCorrPhi)
            WmurecoilEnUp, WmurecoildPhiEnUp, WmuMassEnUp = eventSelector_v2.WRecoil_Phi_Wmass(ep_nMu, mupt, muphi, ep_muPx, ep_muPy,ep_pfMetUncJetEnUp, ep_pfMetCorrPhi)
            WmurecoilEnDown, WmurecoildPhiEnDown, WmuMassEnDown = eventSelector_v2.WRecoil_Phi_Wmass(ep_nMu, mupt, muphi, ep_muPx, ep_muPy,ep_pfMetUncJetEnDown, ep_pfMetCorrPhi)

            #======   usage: ZRecoil_Phi_Zmass(nEle, eleCharge_, elepx_, elepy_, elepz_, elee_,met_,metphi_)=====
            ZeeRecoil,ZeeRecoil_dPhi,ZeeMass = eventSelector_v2.ZRecoil_Phi_Zmass(nEle_loose, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetCorrPt, ep_pfMetCorrPhi)
            ZeeRecoilResUp,ZeeRecoil_dPhiResUp,ZeeMassResUp = eventSelector_v2.ZRecoil_Phi_Zmass(nEle_loose, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetUncJetResUp, ep_pfMetCorrPhi)
            ZeeRecoilResDown,ZeeRecoil_dPhiResDown,ZeeMassResDown = eventSelector_v2.ZRecoil_Phi_Zmass(nEle_loose, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetUncJetResDown, ep_pfMetCorrPhi)
            ZeeRecoilEnUp,ZeeRecoil_dPhiEnUp,ZeeMassEnUp = eventSelector_v2.ZRecoil_Phi_Zmass(nEle_loose, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetUncJetEnUp, ep_pfMetCorrPhi)
            ZeeRecoilEnDown,ZeeRecoil_dPhiEnDown,ZeeMassEnDown = eventSelector_v2.ZRecoil_Phi_Zmass(nEle_loose, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetUncJetEnDown, ep_pfMetCorrPhi)


            ZmumuRecoil,ZmumuRecoil_dPhi,ZmumuMass = eventSelector_v2.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetCorrPt, ep_pfMetCorrPhi)
            ZmumuRecoilResUp,ZmumuRecoil_dPhiResUp,ZmumuMassResUp = eventSelector_v2.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetUncJetResUp, ep_pfMetCorrPhi)
            ZmumuRecoilResDown,ZmumuRecoil_dPhiResDown,ZmumuMassResDown = eventSelector_v2.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetUncJetResDown, ep_pfMetCorrPhi)
            ZmumuRecoilEnUp,ZmumuRecoil_dPhiEnUp,ZmumuMassEnUp = eventSelector_v2.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetUncJetEnUp, ep_pfMetCorrPhi)
            ZmumuRecoilEnDown,ZmumuRecoil_dPhiEnDown,ZmumuMassEnDown = eventSelector_v2.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetUncJetEnDown, ep_pfMetCorrPhi)

            '''
            -----------------------------------------------------------------------------
            MINIMUM dPhi BETWEEN AK4JETS AND MET/RECOIL FOR BOOSTED CATEGORY
            ----------------------------------------------------------------------------
            '''

            min_ak4jet_MET_dPhi = 3.5; mini_AK8jet_MET_dPhi = 3.5
            minDPhi_ak4jet_Werecoil = minDPhi_ak4jet_Wmurecoil = minDPhi_ak4jet_ZeeRecoil = minDPhi_ak4jet_ZmumuRecoil = -10.0
            minDPhi_ak8jet_Werecoil = minDPhi_ak8jet_Wmurecoil = minDPhi_ak8jet_ZeeRecoil = minDPhi_ak8jet_ZmumuRecoil = -10.0

            if ep_THINnJet>0 and WerecoildPhi > -10.0:
                minDPhi_ak4jet_Werecoil    =  min([DeltaPhi(WerecoildPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])#,ak4jetphi
            if ep_THINnJet>0 and WmurecoildPhi > -10.0:
                minDPhi_ak4jet_Wmurecoil   =  min([DeltaPhi(WmurecoildPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])
            if ep_THINnJet>0 and ZeeRecoil_dPhi > -10.0:
                minDPhi_ak4jet_ZeeRecoil   =  min([DeltaPhi(ZeeRecoil_dPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])
            if ep_THINnJet>0 and ZmumuRecoil_dPhi > -10.0:
                minDPhi_ak4jet_ZmumuRecoil =  min([DeltaPhi(ZmumuRecoil_dPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])

            if ep_THINnJet>0:
                min_ak4jet_MET_dPhi = min([DeltaPhi(ep_pfMetCorrPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])



            '''
            -----------------------------------------------------------------------------
            MINIMUM dPhi BETWEEN AK4JETS AND MET/RECOIL FOR RESOLVED CATEGORY
            ----------------------------------------------------------------------------
            '''

            min_ak4jet_MET_dPhi_R = 3.5
            minDPhi_ak4jet_Werecoil_R = minDPhi_ak4jet_Wmurecoil_R = minDPhi_ak4jet_ZeeRecoil_R = minDPhi_ak4jet_ZmumuRecoil_R = -10.0
            ''' 
            if ep_THINnJet > 0 and WerecoildPhi > -10.0:
                minDPhi_ak4jet_Werecoil_R  = min([DeltaPhi(WerecoildPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])#,ak4jetphi
            if ep_THINnJet > 0 and WmurecoildPhi > -10.0:
                minDPhi_ak4jet_Wmurecoil_R   = min([DeltaPhi(WmurecoildPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])
            if ep_THINnJet > 0 and ZeeRecoil_dPhi > -10.0:
                minDPhi_ak4jet_ZeeRecoil_R = min([DeltaPhi(ZeeRecoil_dPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])
            if ep_THINnJet > 0 and ZeeRecoil_dPhi > -10.0:
                minDPhi_ak4jet_ZmumuRecoil_R = min([DeltaPhi(ZmumuRecoil_dPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])
            '''
            if ep_THINnJet > 0:
                min_ak4jet_MET_dPhi_R  = min([DeltaPhi(ep_pfMetCorrPhi,ak4jetphi[ijet]) for ijet in range(ep_THINnJet)])
            '''
            ------------------------------------------------------------------------------
            MINIMUM dPhi BETWEEN AK4JETS AND AK8JETS
            ------------------------------------------------------------------------------
            '''

            mindPhi_ak4_ak8_sr = mindPhi_ak4_ak8_sb = mindPhi_ak4_ak8_zcr = 3.4
            if ep_THINnJet>0 and sel_nfjets > 0:
                index = pass_nfjetIndex[0]
		mindPhi_ak4_ak8_sr = min([DeltaPhi(fatjetphi[index],ak4jetphi[jj]) for jj in range(ep_THINnJet)])

            if ep_THINnJet>0 and nFatJet_SBand >0:
                index=FatJet_SBand_index[0]
		mindPhi_ak4_ak8_sb = min([DeltaPhi(fatjetphi[index],ak4jetphi[jj]) for jj in range(ep_THINnJet)])

	    if ep_THINnJet>0 and nFatJet_ZCR>0:
                index = FatJet_ZCR_index[0]
		mindPhi_ak4_ak8_zcr = min([DeltaPhi(fatjetphi[index],ak4jetphi[jj]) for jj in range(ep_THINnJet)])

            #print mindPhi_ak4_ak8_sr, mindPhi_ak4_ak8_sb, mindPhi_ak4_ak8_zcr
            '''
            ------------------------------------------------------------------------------
             GET REGION[SR/CR] BASED ON PROPERTIES OF EVENT
            ------------------------------------------------------------------------------
            '''

            # ==== usage: getSel_boosted(nEle,nTightEle,isTightEle,nMu,nTightMu,isTightMuon,nTau,nPho,nBjets,cleaned_ak4jets,nFatJet,pfMet,mini_ak4jet_MET_dPhi,ZeeRecoil,min_ak4jets_ZeeRecoil_dPhi,ZeeMass,ZmumuRecoil,min_ak4jets_ZmumuRecoil_dPhi,ZmumuMass,WenuRecoil,min_ak4jets_WenuRecoil_dPhi,WenuMass,WmunuRecoil,min_ak4jets_WmunuRecoil_dPhi,WmunuMass) ======
            
            region_boosted = eventSelector_v2.getSel_boosted(ep_nEle,nTightEle,isTightEles,ep_nMu,nTightMu,isTightMuons,ep_HPSTau_n,nPho,nBjets_iso,\
                     nJets_cleaned,sel_nfjets,nFatJet_SBand,nFatJet_ZCR,ep_pfMetCorrPt,min_ak4jet_MET_dPhi,mini_AK8jet_MET_dPhi,\
                     ZeeRecoil,minDPhi_ak4jet_ZeeRecoil,ZeeMass,ZmumuRecoil,minDPhi_ak4jet_ZmumuRecoil,ZmumuMass,\
                     Werecoil,minDPhi_ak4jet_Werecoil, WeMass,Wmurecoil, minDPhi_ak4jet_Wmurecoil, WmuMass)

            if region_boosted['boosted_signal'] and mettrigdecision:  isBoostedSR      = True
            if region_boosted['boosted_SBand']  and mettrigdecision  : isBoostedSBand    = True
            if region_boosted['boosted_te']     and eletrigdecision:  isBoostedCRTope  = True
            if region_boosted['boosted_tm']     and mettrigdecision:  isBoostedCRTopmu = True
            if region_boosted['boosted_wen']    and eletrigdecision:  isBoostedCRWenu  = True
            if region_boosted['boosted_wmn']    and mettrigdecision:  isBoostedCRWmunu = True
            if region_boosted['boosted_zee']    and eletrigdecision:  isBoostedCRZee   = True
            if region_boosted['boosted_zmm'] and mettrigdecision:  isBoostedCRZmumu = True

            '''
            getSel_resolved(nEle,nTightEle,isTightEle,nMu,nTightMu,isTightMuon,nTau,nPho,\
                   nBjets,nak4jets,ak4CSV,ak4pt,ak4eta,pfMet,mini_ak4jet_MET_dPhi,h_mass,\
                   ZeeRecoil,min_ak4jets_ZeeRecoil_dPhi,ZeeMass,ZmumuRecoil,min_ak4jets_ZmumuRecoil_dPhi,ZmumuMass,WenuRecoil,\
                   min_ak4jets_WenuRecoil_dPhi,WenuMass,WmunuRecoil,min_ak4jets_WmunuRecoil_dPhi,WmunuMass):
            '''

            region_resolved = eventSelector_v2.getSel_resolved(ep_nEle,nTightEle,isTightEles,ep_nMu,nTightMu,isTightMuons,ep_HPSTau_n,nPho,\
                              nBjets_notiso,ep_THINnJet,ep_THINjetDeepCSV,ak4jetpt,ak4jeteta,ep_pfMetCorrPt,min_ak4jet_MET_dPhi_R,h_mass,\
			      ZeeRecoil,minDPhi_ak4jet_ZeeRecoil_R,ZeeMass,ZmumuRecoil,minDPhi_ak4jet_ZmumuRecoil_R,ZmumuMass,\
			      Werecoil,minDPhi_ak4jet_Werecoil_R,WeMass,Wmurecoil,minDPhi_ak4jet_Wmurecoil_R,WmuMass)



            if region_resolved['resolved_signal'] and mettrigdecision and not isBoostedSR      :  isResolvedSR      = True
            if region_resolved['resolved_SBand']  and mettrigdecision and not isBoostedSBand   :  isResolvedSBand   = True
            if region_resolved['resolved_tm']     and mettrigdecision and not isBoostedCRTopmu :  isResolvedCRTopmu = True
            if region_resolved['resolved_te']     and eletrigdecision and not isBoostedCRTope  :  isResolvedCRTope  = True
            if region_resolved['resolved_wmn']    and mettrigdecision and not isBoostedCRWmunu :  isResolvedCRWmunu = True
            if region_resolved['resolved_wen']    and eletrigdecision and not isBoostedCRWenu  :  isResolvedCRWenu  = True
            if region_resolved['resolved_zmm']    and mettrigdecision and not isBoostedCRZmumu :  isResolvedCRZmumu = True
            if region_resolved['resolved_zee']    and eletrigdecision and not isBoostedCRZee   :  isResolvedCRZee   = True



            '''
            --------------------------------------------------------------------------------
            COMMAN WEIGHT CALCULATION FOR ALL REGIONS
            --------------------------------------------------------------------------------
            '''
            
           

            weight =1.0; JEC_up =1.0; JEC_down =1.0 ;PUweight =1.0; PUweight_up =1.0; PUweight_down =1.0; lepweight =1.0; lepweight_up =1.0; lepweight_down =1.0; recoilweight =1.0; recoil_up =1.0; recoil_down =1.0; recoilweight = 1.0
            btagweight =1.0; btagweight_up =1.0; btagweight_down =1.0; btagweight_B =1.0; btagweight_B_up =1.0 ;btagweight_B_down =1.0 ;ewkweight =1.0; ewkweight_up =1.0; ewkweight_down =1.0; qcdk=1.0
            toppTweight =1.0; toppTweight_up =1.0; toppTweight_down =1.0; METweight =1.0; METweight_up =1.0; METweight_down =1.0; R_weight=1.0;B_weight=1.0

            commanweight=1.0;commanweight_B=1.0;R_weight=1.0;B_weight=1.0;



            if not isData:
                btagweight,btagweight_up,btagweight_down     = wgt.getBTagSF(ep_THINnJet,ak4jetpt,ak4jeteta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV,'MWP',index=False)
                btagweight_B,btagweight_B_up,btagweight_B_down = wgt.getBTagSF(pass_ak4jet_index_cleaned,ak4jetpt,ak4jeteta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV,'LWP',index=True)

                if ep_genParSample   == 23 and len(ep_genParPt) > 0 : 
		    ewkweight = wgt.getEWKZ(ep_genParPt[0])
                    qcdk      = wgt.getQCDZ(ep_genParPt[0])
                if ep_genParSample == 24 and len(ep_genParPt) > 0 : 
                    ewkweight = wgt.getEWKW(ep_genParPt[0])
                    qcdk      = wgt.getQCDW(ep_genParPt[0])
                if ep_genParSample == 6 and len(ep_genParPt) > 0  : toppTweight,toppTweight_up,toppTweight_down = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

	        PUweight, PUweight_up, PUweight_down = wgt.puweight(ep_pu_nTrueInt)
                commanweight = ewkweight*qcdk*toppTweight*PUweight*btagweight
                commanweight_B = ewkweight*qcdk*toppTweight*PUweight*btagweight_B
            
                R_Lepweight =  eventSelector_v2.getweight(ep_pfMetCorrPt,Wmurecoil,ZmumuRecoil,nPho,ep_HPSTau_n,nEle_loose,nTightEle,isTightEles,ele_loose_index,ele_tight_index,elept,eleeta,ep_nMu,nTightMu,isTightMuons,muon_tight_index,mupt,mueta,cat='R')
                B_Lepweight =  eventSelector_v2.getweight(ep_pfMetCorrPt,Wmurecoil,ZmumuRecoil,nPho,ep_HPSTau_n,nEle_loose,nTightEle,isTightEles,ele_loose_index,ele_tight_index,elept,eleeta,ep_nMu,nTightMu,isTightMuons,muon_tight_index,mupt,mueta,cat='B')
               
                R_weight = R_Lepweight*commanweight
                B_weight = B_Lepweight*commanweight_B
	    additional_jets=ep_THINnJet-2
            #cutFlow(nEle_loose,nEle_tight,nMu_loose,nMu_tight,nTau,Werecoil,Wmurecoil,ZeeRecoil,ZmumuRecoil,pfMet,njets,nBjets,h_mass,nAK8JetsSR,nAK8JetsSBand,nAK8JetsZCR,ZeeMass,ZmumuMass)


            '''
            ------------------------------------------------------------
            CUTFLOW PART FOR BOOSTED AND RESOLVED ANALYSIS
            ------------------------------------------------------------
            '''

            cutFlow_bins = cutFlow(nPho,nEle_loose,nTightEle,ep_nMu,nTightMu,ep_HPSTau_n,Werecoil,Wmurecoil,ZeeRecoil,ZmumuRecoil,ep_pfMetCorrPt,ep_THINnJet,nJets_cleaned,nBjets_notiso,nBjets_iso,h_mass,sel_nfjets,nFatJet_SBand,nFatJet_ZCR,ZeeMass,ZmumuMass)
           
            TopmuBins_B = cutFlow_bins.singleLepton_B(B_weight,isEleRegion=False,isTop=True)
            TopeBins_B  = cutFlow_bins.singleLepton_B(B_weight,isEleRegion=True,isTop=True)
            WmuBins_B   = cutFlow_bins.singleLepton_B(B_weight,isEleRegion=False,isTop=False)
            WeBins_B    = cutFlow_bins.singleLepton_B(B_weight,isEleRegion=True,isTop=False)

            ZeeBins_B   = cutFlow_bins.diLepton_B(B_weight,isEleRegion=True)
            ZmumuBins_B = cutFlow_bins.diLepton_B(B_weight,isEleRegion=False)
            SBand_B     = cutFlow_bins.SBand_B(B_weight) 

            TopmuBins_R = cutFlow_bins.singleLepton_R(R_weight,isEleRegion=False,isTop=True)
            TopeBins_R  = cutFlow_bins.singleLepton_R(R_weight,isEleRegion=True,isTop=True)
            WmuBins_R   = cutFlow_bins.singleLepton_R(R_weight,isEleRegion=False,isTop=False)
            WeBins_R    = cutFlow_bins.singleLepton_R(R_weight,isEleRegion=True,isTop=False)

            ZeeBins_R   = cutFlow_bins.diLepton_R(R_weight,isEleRegion=True)
            ZmumuBins_R = cutFlow_bins.diLepton_R(R_weight,isEleRegion=False)
            SBand_R     = cutFlow_bins.SBand_R(R_weight)


            if eletrigdecision:
                h_reg_TopenuCR_resolved_cutFlow.AddBinContent(3,1)
                h_reg_WenuCR_resolved_cutFlow.AddBinContent(3,1)
                h_reg_ZeeCR_resolved_cutFlow.AddBinContent(3,1)
                h_reg_TopenuCR_boosted_cutFlow.AddBinContent(3,1)
                h_reg_WenuCR_boosted_cutFlow.AddBinContent(3,1)
                h_reg_ZeeCR_boosted_cutFlow.AddBinContent(3,1)

                for ibin, binweight in enumerate(TopeBins_R):
                    h_reg_TopenuCR_resolved_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(WeBins_R):
                    h_reg_WenuCR_resolved_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(ZeeBins_R):
                    h_reg_ZeeCR_resolved_cutFlow.AddBinContent(4+ibin,binweight)

                for ibin, binweight in enumerate(TopeBins_B):
                    h_reg_TopenuCR_boosted_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(WeBins_B):
                    h_reg_WenuCR_boosted_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(ZeeBins_B):
                    h_reg_ZeeCR_boosted_cutFlow.AddBinContent(4+ibin,binweight)



            if mettrigdecision:
                h_reg_TopmunuCR_resolved_cutFlow.AddBinContent(3,1)
                h_reg_WmunuCR_resolved_cutFlow.AddBinContent(3,1)
                h_reg_ZmumuCR_resolved_cutFlow.AddBinContent(3,1)
                h_reg_SBand_resolved_cutFlow.AddBinContent(3,1)
                h_reg_TopmunuCR_boosted_cutFlow.AddBinContent(3,1)
                h_reg_WmunuCR_boosted_cutFlow.AddBinContent(3,1)
                h_reg_ZmumuCR_boosted_cutFlow.AddBinContent(3,1)
                h_reg_SBand_boosted_cutFlow.AddBinContent(3,1)
                for ibin, binweight in enumerate(TopmuBins_R):
                    h_reg_TopmunuCR_resolved_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(WmuBins_R):
                    h_reg_WmunuCR_resolved_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(ZmumuBins_R):
                    h_reg_ZmumuCR_resolved_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(SBand_R):
                    h_reg_SBand_resolved_cutFlow.AddBinContent(4+ibin,binweight)

                for ibin, binweight in enumerate(TopmuBins_B):
                    h_reg_TopmunuCR_boosted_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(WmuBins_B):
                    h_reg_WmunuCR_boosted_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(ZmumuBins_B):
                    h_reg_ZmumuCR_boosted_cutFlow.AddBinContent(4+ibin,binweight)
                for ibin, binweight in enumerate(SBand_B):
                    h_reg_SBand_boosted_cutFlow.AddBinContent(4+ibin,binweight)


            '''
            --------------------------------------------------------------------------------
            SIGNAL REGION RESOLVED
            --------------------------------------------------------------------------------
            '''

           
            if  isResolvedSR:
		
                if not isData:
                        ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
                        weightMET,weightMET_up,weightMET_down=wgt.getMETtrig_First(ep_pfMetCorrPt,cat='R')
                        weight = commanweight*weightMET
                        JEC_up,JEC_down = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
               
		#print 'SR_R','R_weight',R_weight,'weight',weight
    
                df_out_SR_resolved = df_out_SR_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                               'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                               'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                               'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet1CSV':ep_THINjetDeepCSV[jet1Index],
                                               'Jet2Pt':ak4jetpt[jet2Index], 'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index], 'Jet2CSV':ep_THINjetDeepCSV[jet2Index],
                                               'Jet3Pt':dummy, 'Jet3Eta':dummy, 'Jet3Phi':dummy, 'Jet3CSV':dummy,
                                               'DiJetMass':h_mass,'nJets':additional_jets,
                                               'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                               'METweight':METweight,'METweight_up':METweight_up,'METweight_down':METweight_down,'METRes_up':ep_pfMetUncJetResUp[0],'METRes_down':ep_pfMetUncJetResDown[0],'METEn_up':ep_pfMetUncJetEnUp[0],'METEn_down':ep_pfMetUncJetEnDown[0],
                                               'btagweight':btagweight,'btagweight_up':btagweight_up,'btagweight_down':btagweight_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                               'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                              ignore_index=True)
          
            if  isResolvedSBand:
                
                if not isData:
                        ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
                        weightMET,weightMET_up,weightMET_down=wgt.getMETtrig_First(ep_pfMetCorrPt,cat='R')
                        weight = commanweight*weightMET
                        JEC_up,JEC_down = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
            
                #print 'SBand_R','R_weight',R_weight,'weight',weight

                df_out_SBand_resolved = df_out_SBand_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                               'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                               'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                               'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet1CSV':ep_THINjetDeepCSV[jet1Index],
                                               'Jet2Pt':ak4jetpt[jet2Index], 'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index], 'Jet2CSV':ep_THINjetDeepCSV[jet2Index],
                                               'Jet3Pt':dummy, 'Jet3Eta':dummy, 'Jet3Phi':dummy, 'Jet3CSV':dummy,
                                               'DiJetMass':h_mass,'nJets':additional_jets,
                                               'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                               'METweight':METweight,'METweight_up':METweight_up,'METweight_down':METweight_down,'METRes_up':ep_pfMetUncJetResUp[0],'METRes_down':ep_pfMetUncJetResDown[0],'METEn_up':ep_pfMetUncJetEnUp[0],'METEn_down':ep_pfMetUncJetEnDown[0],
                                               'btagweight':btagweight,'btagweight_up':btagweight_up,'btagweight_down':btagweight_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                               'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                              ignore_index=True)

	    
	    if isResolvedCRTope:

                ele1_index    = ele_tight_index[0]

                if not isData:
                    lepIdReco,lepIdReco_up,lepIdReco_down = wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],'T')
                    lepTrigSF,lepTrigSF_up,lepTrigSF_down = wgt.eletrig_weight(elept[ele1_index],eleeta[ele1_index])
                    lepweight = lepIdReco*lepTrigSF
                    lepweight_up = lepIdReco_up*lepTrigSF_up
                    lepweight_down = lepIdReco_down*lepTrigSF_down
                    ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
                    weight = lepweight*commanweight
                    JEC_up,JEC_down = getJECWeight(ep_THINnJet,ep_THINjetCorrUnc,index=False)
               
		#print 'Tope','R_weight',R_weight,'weight',weight

                df_out_Tope_resolved  = df_out_Tope_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WerecoildPhi,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],'Wmass':WeMass,
                                                'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':1.0,'recoilweight_up':1.0,'recoilweight_down':1.0,'recoilRes_up':WerecoilResUp,'recoilRes_down':WerecoilResDown,'recoilEn_up':WerecoilEnUp,'recoilEn_down':WerecoilEnDown,
                                                'btagweight':btagweight,'btagweight_up':btagweight_up,'btagweight_down':btagweight_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)

           
            if isResolvedCRTopmu:
                muon1_index          = muon_tight_index[0]

                if not isData:
                    recoilweight,recoil_up,recoil_down=wgt.getMETtrig_First(Wmurecoil,cat='R') 
                    lepweight,lep_up,lep_down=wgt.mu_weight(mupt[0],mueta[0],'T')
                    ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
                    weight = lepweight*commanweight*recoilweight
                    JEC_up,JEC_down = getJECWeight(ep_THINnJet,ep_THINjetCorrUnc,index=False)
                
		#print 'Topmu','R_weight',R_weight,'weight',weight

                df_out_Topmu_resolved  = df_out_Topmu_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WmurecoildPhi,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],'Wmass':WmuMass,
                                                'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':recoilweight,'recoilweight_up':recoil_up,'recoilweight_down':recoil_down,'recoilRes_up':WmurecoilResUp,'recoilRes_down':WmurecoilResDown,'recoilEn_up':WmurecoilEnUp,'recoilEn_down':WmurecoilEnDown,
                                                'btagweight':btagweight,'btagweight_up':btagweight_up,'btagweight_down':btagweight_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)

          
            if isResolvedCRWenu:
                ele1_index    = ele_tight_index[0]
                ele_trig      = True

                if not isData: 
                    lepIdReco,lepIdReco_up,lepIdReco_down = wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],'T')
                    lepTrigSF,lepTrigSF_up,lepTrigSF_down = wgt.eletrig_weight(elept[ele1_index],eleeta[ele1_index])
                    lepweight = lepIdReco*lepTrigSF
                    lepweight_up = lepIdReco_up*lepTrigSF_up
                    lepweight_down = lepIdReco_down*lepTrigSF_down
                    ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
		    weight = lepweight*commanweight
                    JEC_up,JEC_down = getJECWeight(ep_THINnJet,ep_THINjetCorrUnc,index=False)
                
                #print 'We','R_weight',R_weight,'weight',weight

                df_out_We_resolved  = df_out_We_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WerecoildPhi,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],'Wmass':WeMass,
						'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':1.0,'recoilweight_up':1.0,'recoilweight_down':1.0,'recoilRes_up':WerecoilResUp,'recoilRes_down':WerecoilResDown,'recoilEn_up':WerecoilEnUp,'recoilEn_down':WerecoilEnDown,
                                                'btagweight':btagweight,'btagweight_up':btagweight_up,'btagweight_down':btagweight_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)
           
            if isResolvedCRWmunu:
                muon1_index          = muon_tight_index[0]

                if not isData:
                    recoilweight,recoil_up,recoil_down = wgt.getMETtrig_First(Wmurecoil,cat='R') 
                    lepweight,lep_up,lep_down=wgt.mu_weight(mupt[0],mueta[0],'T')
                    ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
                    weight = lepweight*commanweight*recoilweight
                    JEC_up,JEC_down = getJECWeight(ep_THINnJet,ep_THINjetCorrUnc,index=False)
               
                #print 'Wmu','R_weight',R_weight,'weight',weight


                df_out_Wmu_resolved  = df_out_Wmu_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WmurecoildPhi,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],'Wmass':WmuMass,
						'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':recoilweight,'recoilweight_up':recoil_up,'recoilweight_down':recoil_down,'recoilRes_up':WmurecoilResUp,'recoilRes_down':WmurecoilResDown,'recoilEn_up':WmurecoilEnUp,'recoilEn_down':WmurecoilEnDown,
                                                'btagweight':btagweight,'btagweight_up':btagweight_up,'btagweight_down':btagweight_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)


            
            if isResolvedCRZee:
               
                ele1_loose = ele_loose_index[0]
                ele2_loose = ele_loose_index[1]
                if isTightEles[0] and not isData:
		    lepweight1,lep1_up,lep1_down          = (wgt.ele_weight(elept[ele1_loose],eleeta[ele1_loose],'T'))
		    if isTightEles[1]:lepweight2,lep2_up,lep2_down          = (wgt.ele_weight(elept[ele2_loose],eleeta[ele2_loose],'T'))
                    else:lepweight2,lep2_up,lep2_down          = (wgt.ele_weight(elept[ele2_loose],eleeta[ele2_loose],'L'))
                   
                    lepTrigSF,lepTrigSF_up,lepTrigSF_down = wgt.eletrig_weight(elept[ele1_loose],eleeta[ele1_loose])
                    lepweight                             = lepweight1*lepweight2*lepTrigSF
                    lepweight_up                          = lep1_up*lep2_up*lepTrigSF_up
                    lepweight_down                        = lep1_down*lep2_down*lepTrigSF_down

                elif isTightEles[1] and not isData:
                    lepTrigSF,lepTrigSF_up,lepTrigSF_down = wgt.eletrig_weight(elept[ele2_loose],eleeta[ele2_loose])
		    lepweight1,lep1_up,lep1_down=(wgt.ele_weight(elept[ele1_loose],eleeta[ele1_loose],'L'))
		    lepweight2,lep2_up,lep2_down=(wgt.ele_weight(elept[ele2_loose],eleeta[ele2_loose],'T'))
		    lepweight = lepweight1*lepweight2*lepTrigSF
                    lepweight_up                          = lep1_up*lep2_up*lepTrigSF_up
                    lepweight_down                        = lep1_down*lep2_down*lepTrigSF_down

                ZpT = math.sqrt( (ep_elePx[0] + ep_elePx[0])**2 + (ep_elePy[1]+ep_elePy[1])**2 )

                if not isData: 
		    weight = lepweight*commanweight
                    ewkweight_up=ewkweight*1.5; ewkweight_down=ewkweight*0.5
                    JEC_up,JEC_down = getJECWeight(ep_THINnJet,ep_THINjetCorrUnc,index=False)
               
                #print 'Zee','R_weight',R_weight,'weight',weight
               
                df_out_Zee_resolved    = df_out_Zee_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZeeRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':ZeeRecoil_dPhi,
                                                'lep1_pT':elept[ele1_loose],'lep1_eta':eleeta[ele1_loose],'lep1_Phi':elephi[ele1_loose],
                                                'lep2_pT':elept[ele2_loose],'lep2_eta':eleeta[ele2_loose],'lep2_Phi':elephi[ele2_loose],
                                                'Zmass':ZeeMass,'ZpT':ZpT,
						'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':1.0,'recoilweight_up':1.0,'recoilweight_down':1.0,'recoilRes_up':ZeeRecoilResUp,'recoilRes_down':ZeeRecoilResDown,'recoilEn_up':ZeeRecoilEnUp,'recoilEn_down':ZeeRecoilEnDown,
                                                'btagweight':btagweight,'btagweight_up':btagweight_up,'btagweight_down':btagweight_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)


            
            if isResolvedCRZmumu:

                if isTightMuons[0] and not isData:
		    lepweight1,lep1_up,lep1_down = wgt.mu_weight(mupt[0],mueta[0],'T')
                    if isTightMuons[1]:lepweight2,lep2_up,lep2_down = wgt.mu_weight(mupt[1],mueta[1],'T')
                    else:lepweight2,lep2_up,lep2_down = wgt.mu_weight(mupt[1],mueta[1],'L')
                    lepweight = lepweight1*lepweight2

                elif isTightMuons[1] and not isData:
		    lepweight1,lep1_up,lep1_down = wgt.mu_weight(mupt[1],mueta[1],'L')
                    lepweight2,lep2_up,lep2_down = wgt.mu_weight(mupt[0],mueta[0],'T')
                    lepweight = lepweight1*lepweight2

                ZpT  = math.sqrt( (ep_muPx[0] + ep_muPx[1])**2 + (ep_muPy[0]+ep_muPy[1])**2 )
               

                if not isData: 
                    recoilweight,recoil_up,recoil_down = wgt.getMETtrig_First(ZmumuRecoil,cat='R')
                    ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
		    weight = lepweight*commanweight*recoilweight
                    JEC_up,JEC_down = getJECWeight(ep_THINnJet,ep_THINjetCorrUnc,index=False)
             
		#print 'Zmumu','R_weight',R_weight,'weight',weight

                df_out_Zmumu_resolved  = df_out_Zmumu_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZmumuRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
						'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':ZmumuRecoil_dPhi,
                                                'lep1_pT':mupt[0],'lep1_eta':mueta[0],'lep1_Phi':muphi[0],
                                                'lep2_pT':mupt[1],'lep2_eta':mueta[1],'lep2_Phi':muphi[1],
                                                'Zmass':ZmumuMass,'ZpT':ZpT,
						'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':recoilweight,'recoilweight_up':recoil_up,'recoilweight_down':recoil_down,'recoilRes_up':ZmumuRecoilResUp,'recoilRes_down':ZmumuRecoilResDown,'recoilEn_up':ZmumuRecoilEnUp,'recoilEn_down':ZmumuRecoilEnDown,
                                                'btagweight':btagweight,'btagweight_up':btagweight_up,'btagweight_down':btagweight_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)

            '''
            --------------------------------------------------------------------------------
            SIGNAL REGION BOOSTED
            --------------------------------------------------------------------------------
            '''

            jet1pT=jet1Eta=jet1Phi=jet1CSV=jet2CSV= -9999.0
            if nJets_cleaned > 0:
                jet1pT           = ak4jetpt[pass_ak4jet_index_cleaned[0]]
                jet1Eta          = ak4jeteta[pass_ak4jet_index_cleaned[0]]
                jet1Phi          = ak4jetphi[pass_ak4jet_index_cleaned[0]]
                jet1CSV          = ep_THINjetDeepCSV[pass_ak4jet_index_cleaned[0]]
                if nJets_cleaned>1:jet2CSV          = ep_THINjetDeepCSV[pass_ak4jet_index_cleaned[1]]

            if isBoostedSR:
                fjet_index           = pass_nfjetIndex[0]
                min_dPhi_ak4_MET     = min_ak4jet_MET_dPhi
                fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))
		lepweight = lepweight_up = lepweight_down = 1.0
                isAK4jet1EtaMatch=0.0;isAK4jet2EtaMatch=0;M_Jet1AK8Jet=dummy;M_Jet2AK8Jet=dummy
                if nJets_cleaned>0: 
                    M_Jet1AK8Jet = InvMass(ep_THINjetPx[pass_ak4jet_index_cleaned[0]], ep_THINjetPy[pass_ak4jet_index_cleaned[0]], ep_THINjetPz[pass_ak4jet_index_cleaned[0]], ep_THINjetEnergy[pass_ak4jet_index_cleaned[0]],ep_fjetPx[fjet_index], ep_fjetPy[fjet_index], ep_fjetPz[fjet_index],ep_fjetEnergy[fjet_index])
		    if fatjeteta[fjet_index]*ak4jeteta[pass_ak4jet_index_cleaned[0]] > 0  : isAK4jet1EtaMatch=1
                if nJets_cleaned>1:
                    M_Jet2AK8Jet = InvMass(ep_THINjetPx[pass_ak4jet_index_cleaned[1]], ep_THINjetPy[pass_ak4jet_index_cleaned[1]], ep_THINjetPz[pass_ak4jet_index_cleaned[1]], ep_THINjetEnergy[pass_ak4jet_index_cleaned[1]],ep_fjetPx[fjet_index], ep_fjetPy[fjet_index], ep_fjetPz[fjet_index],ep_fjetEnergy[fjet_index])                   
		    if fatjeteta[fjet_index]*ak4jeteta[pass_ak4jet_index_cleaned[1]] > 0  : isAK4jet2EtaMatch=1
                 

                if not isData: 
                        ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
                        weightMET,weightMET_up,weightMET_down=wgt.getMETtrig_First(ep_pfMetCorrPt,cat='B')
			weight = commanweight_B*weightMET
                        JEC_up,JEC_down = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
               
                #print 'SR_B','B_weight',B_weight,'weight',weight
                df_out_SR_boosted = df_out_SR_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,'isAK4jet1EtaMatch':isAK4jet1EtaMatch,'isAK4jet2EtaMatch':isAK4jet2EtaMatch,'M_Jet1AK8Jet':M_Jet1AK8Jet,'M_Jet2AK8Jet':M_Jet2AK8Jet,'Jet1CSV':jet1CSV,'Jet1CSV':jet2CSV,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_dPhi_ak4_MET,'met_Phi':ep_pfMetCorrPhi,'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,'min_dphi_jets':mindPhi_ak4_ak8_sr,
                                               'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':1.0,'lepweight_up':1.0,'lepweight_down':1.0,
                                               'METweight':METweight,'METweight_up':METweight_up,'METweight_down':METweight_down,'METRes_up':ep_pfMetUncJetResUp[0],'METRes_down':ep_pfMetUncJetResDown[0],'METEn_up':ep_pfMetUncJetEnUp[0],'METEn_down':ep_pfMetUncJetEnDown[0],
                                               'btagweight':btagweight_B,'btagweight_up':btagweight_B_up,'btagweight_down':btagweight_B_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                               'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)

            '''
            --------------------------------------------------------------------------------
            SIDEBAND CONTROL REGION BOOSTED
            --------------------------------------------------------------------------------
            '''

            if isBoostedSBand:
                #print 'mini_AK8jet_MET_dPhi', mini_AK8jet_MET_dPhi
                fjet_index           = FatJet_SBand_index[0]
                min_dPhi_ak4_MET     = min_ak4jet_MET_dPhi
                fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))
		lepweight = lepweight_up = lepweight_down = 1.0
            
                isAK4jet1EtaMatch=0.0;isAK4jet2EtaMatch=0;M_Jet1AK8Jet=dummy;M_Jet2AK8Jet=dummy
                if nJets_cleaned>0:
                    M_Jet1AK8Jet = InvMass(ep_THINjetPx[pass_ak4jet_index_cleaned[0]], ep_THINjetPy[pass_ak4jet_index_cleaned[0]], ep_THINjetPz[pass_ak4jet_index_cleaned[0]], ep_THINjetEnergy[pass_ak4jet_index_cleaned[0]],ep_fjetPx[fjet_index], ep_fjetPy[fjet_index], ep_fjetPz[fjet_index],ep_fjetEnergy[fjet_index])
                    if fatjeteta[fjet_index]*ak4jeteta[pass_ak4jet_index_cleaned[0]] > 0  : isAK4jet1EtaMatch=1
                if nJets_cleaned>1:
                    M_Jet2AK8Jet = InvMass(ep_THINjetPx[pass_ak4jet_index_cleaned[1]], ep_THINjetPy[pass_ak4jet_index_cleaned[1]], ep_THINjetPz[pass_ak4jet_index_cleaned[1]], ep_THINjetEnergy[pass_ak4jet_index_cleaned[1]],ep_fjetPx[fjet_index], ep_fjetPy[fjet_index], ep_fjetPz[fjet_index],ep_fjetEnergy[fjet_index])
                    if fatjeteta[fjet_index]*ak4jeteta[pass_ak4jet_index_cleaned[1]] > 0  : isAK4jet2EtaMatch=1
 
                if not isData:
                    weightMET,MET_up,MET_down=wgt.getMETtrig_First(ep_pfMetCorrPt,cat='B')
                    ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5 
		    weight = commanweight_B*weightMET
                    JEC_up,JEC_down = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
             
                #print 'SBand_B','B_weight',B_weight,'weight',weight

                df_out_SBand_boosted = df_out_SBand_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjet_index], 'FJetEta':fatjeteta[fjet_index], 'FJetPhi':fatjetphi[fjet_index], 'FJetCSV':ep_fjetProbHbb[fjet_index],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,'isAK4jet1EtaMatch':isAK4jet1EtaMatch,'isAK4jet2EtaMatch':isAK4jet2EtaMatch,'M_Jet1AK8Jet':M_Jet1AK8Jet,'M_Jet2AK8Jet':M_Jet2AK8Jet,'Jet1CSV':jet1CSV,'Jet1CSV':jet2CSV,
                                                'FJetMass':ep_fjetSDMass[fjet_index], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_dPhi_ak4_MET,'met_Phi':ep_pfMetCorrPhi,'FJetN2b1':ep_fjetN2b1[fjet_index],'FJetN2b2':ep_fjetN2b2[fjet_index],'FJetrho':fatjet_rho,'min_dphi_jets':mindPhi_ak4_ak8_sb,
                                                'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':1.0,'lepweight_up':1.0,'lepweight_down':1.0,
                                                'METweight':METweight,'METweight_up':METweight_up,'METweight_down':METweight_down,'METRes_up':ep_pfMetUncJetResUp[0],'METRes_down':ep_pfMetUncJetResDown[0],'METEn_up':ep_pfMetUncJetEnUp[0],'METEn_down':ep_pfMetUncJetEnDown[0],
                                                'btagweight':btagweight_B,'btagweight_up':btagweight_B_up,'btagweight_down':btagweight_B_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)


            '''
            --------------------------------------------------------------------------------
            TOP EleREGION BOOSTED
            --------------------------------------------------------------------------------
            '''

            if isBoostedCRTope:
                fjet_index           = pass_nfjetIndex[0]
                min_dPhi_ak4_Recoil  = minDPhi_ak4jet_Werecoil
                fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))

                ele1_index           = ele_tight_index[0]

                if not isData:
                    weightele,lep_up,lep_down             = wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],'T')
                    lepTrigSF,lepTrigSF_up,lepTrigSF_down = wgt.eletrig_weight(elept[ele1_index],eleeta[ele1_index])
                    lepweight                             = weightele*lepTrigSF
                    lepweight_up                          = lep_up*lepTrigSF_up
                    lepweight_down                        = lep_down*lepTrigSF_down
                    ewkweight_up=ewkweight*1.5;  ewkweight_down=ewkweight*0.5
                    weight                                = commanweight_B*lepweight
                    JEC_up,JEC_down                       = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
               
                #print 'Tope_B','B_weight',B_weight,'weight',weight 
               
                df_out_Tope_boosted  = df_out_Tope_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WerecoildPhi,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,'min_dphi_jets':mindPhi_ak4_ak8_sr,'Wmass':WeMass,
                                                'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':1.0,'recoilweight_up':1.0,'recoilweight_down':1.0,'recoilRes_up':WerecoilResUp,'recoilRes_down':WerecoilResDown,'recoilEn_up':WerecoilEnUp,'recoilEn_down':WerecoilEnDown,
                                                'btagweight':btagweight_B,'btagweight_up':btagweight_B_up,'btagweight_down':btagweight_B_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)

            '''
            --------------------------------------------------------------------------------
            TOP MuREGION BOOSTED
            --------------------------------------------------------------------------------
            '''


            if isBoostedCRTopmu:
                fjet_index           = pass_nfjetIndex[0]
                min_dPhi_ak4_Recoil  = minDPhi_ak4jet_Wmurecoil
                muon1_index          = muon_tight_index[0]
                fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))

                if not isData: 
                    lepweight,lepweight_up,lepweight_down      =  wgt.mu_weight(mupt[0],mueta[0],'T')
                    recoilweight,recoil_up,recoil_down         =  wgt.getMETtrig_First(Wmurecoil,cat='B')
		    weight                                     = commanweight_B*lepweight*recoilweight
                    JEC_up,JEC_down                            = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
                    ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
             
                #print 'Topmu_B','B_weight',B_weight,'weight',weight
                df_out_Topmu_boosted  = df_out_Topmu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WmurecoildPhi,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,'min_dphi_jets':mindPhi_ak4_ak8_sr,'Wmass':WmuMass,
                                                'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':recoilweight,'recoilweight_up':recoil_up,'recoilweight_down':recoil_down,'recoilRes_up':WmurecoilResUp,'recoilRes_down':WmurecoilResDown,'recoilEn_up':WmurecoilEnUp,'recoilEn_down':WmurecoilEnDown,
                                                'btagweight':btagweight_B,'btagweight_up':btagweight_B_up,'btagweight_down':btagweight_B_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)

            '''
            --------------------------------------------------------------------------------
            W EleREGION BOOSTED
            --------------------------------------------------------------------------------
            '''

            if isBoostedCRWenu:
                fjet_index           = pass_nfjetIndex[0]
                min_dPhi_ak4_Recoil  = minDPhi_ak4jet_Werecoil
                fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))

                ele1_index           = ele_tight_index[0]

                if not isData:
                    weightele,lep_up,lep_down             = wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],'T')
                    lepTrigSF,lepTrigSF_up,lepTrigSF_down = wgt.eletrig_weight(elept[ele1_index],eleeta[ele1_index])
                    lepweight                             = weightele*lepTrigSF
                    lepweight_up                          = lep_up*lepTrigSF_up
                    lepweight_down                        = lep_down*lepTrigSF_down
                    ewkweight_up=ewkweight*1.5;  ewkweight_down=ewkweight*0.5 
		    weight                                = commanweight_B*lepweight
                    JEC_up,JEC_down                       = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
             
                #print 'We_B','B_weight',B_weight,'weight',weight

                df_out_We_boosted  = df_out_We_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WerecoildPhi,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,'min_dphi_jets':mindPhi_ak4_ak8_sr,'Wmass':WeMass,
                                                'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':1.0,'recoilweight_up':1.0,'recoilweight_down':1.0,'recoilRes_up':WerecoilResUp,'recoilRes_down':WerecoilResDown,'recoilEn_up':WerecoilEnUp,'recoilEn_down':WerecoilEnDown,
                                                'btagweight':btagweight_B,'btagweight_up':btagweight_B_up,'btagweight_down':btagweight_B_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)
            

            '''
            --------------------------------------------------------------------------------
            W MuREGION BOOSTED
            --------------------------------------------------------------------------------
            '''


            if isBoostedCRWmunu:
                fjet_index           = pass_nfjetIndex[0]
                min_dPhi_ak4_Recoil  = minDPhi_ak4jet_Wmurecoil
                muon1_index          = muon_tight_index[0]
                fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))


                if not isData: 
                    lepweight,lepweight_up,lepweight_down=wgt.mu_weight(mupt[0],mueta[0],'T')
		    weightRecoil,recoil_up,recoil_down=wgt.getMETtrig_First(Wmurecoil,cat='B')
		    weight = commanweight_B*lepweight*weightRecoil
                    JEC_up,JEC_down = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
             
                #print 'Wmu_B','B_weight',B_weight,'weight',weight

                df_out_Wmu_boosted  = df_out_Wmu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WmurecoildPhi,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,'min_dphi_jets':mindPhi_ak4_ak8_sr,'Wmass':WmuMass,
                                                'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':recoilweight,'recoilweight_up':recoil_up,'recoilweight_down':recoil_down,'recoilRes_up':WmurecoilResUp,'recoilRes_down':WmurecoilResDown,'recoilEn_up':WmurecoilEnUp,'recoilEn_down':WmurecoilEnDown,
                                                'btagweight':btagweight_B,'btagweight_up':btagweight_B_up,'btagweight_down':btagweight_B_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)


            '''
            --------------------------------------------------------------------------------
            ZCR MuREGION BOOSTED
            --------------------------------------------------------------------------------
            '''


            if isBoostedCRZmumu:
                fjet_index           = FatJet_ZCR_index[0]#pass_nfjetIndex[0]
               
                ZpT = math.sqrt( (ep_muPx[0] + ep_muPx[1])**2 + (ep_muPy[0]+ep_muPy[1])**2 )

                #print 'number of muons',len(mupt)
                if isTightMuons[0] and not isData:
		    lepweight1,lep1_up,lep1_down = wgt.mu_weight(mupt[0],mueta[0],'T')
                    if isTightMuons[1]:lepweight2,lep2_up,lep2_down = wgt.mu_weight(mupt[1],mueta[1],'T')
                    else:lepweight2,lep2_up,lep2_down = wgt.mu_weight(mupt[1],mueta[1],'L')
                    lepweight = lepweight1*lepweight2
                    lepweight_up = lep1_up*lep2_up
                    lepweight_down = lep1_down*lep2_down

                elif isTightMuons[1] and not isData:
		    lepweight1,lep1_up,lep1_down = wgt.mu_weight(mupt[1],mueta[1],'L')
                    lepweight2,lep2_up,lep2_down = wgt.mu_weight(mupt[0],mueta[0],'T')
                    lepweight = lepweight1*lepweight2
                    lepweight_up = lep1_up*lep2_up
                    lepweight_down = lep1_down*lep2_down


                if ep_fjetSDMass[fjet_index]==0: fatjet_rho = math.log((0.0001*0.0001)/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))
                else:fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))


                if not isData:
                    ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
		    recoilweight,recoil_up,recoil_down         =  wgt.getMETtrig_First(ZmumuRecoil,cat='B')
	 	    weight                                     = commanweight_B*lepweight*recoilweight
                    JEC_up,JEC_down                            = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
                
                #print 'Zmumu_B','B_weight',B_weight,'weight',weight

                df_out_Zmumu_boosted  = df_out_Zmumu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZmumuRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjet_index], 'FJetEta':fatjeteta[fjet_index], 'FJetPhi':fatjetphi[fjet_index], 'FJetCSV':ep_fjetProbHbb[fjet_index],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjet_index], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':ZmumuRecoil_dPhi,
                                                'lep1_pT':mupt[0],'lep1_eta':mueta[0],'lep1_Phi':muphi[0],
                                                'lep2_pT':mupt[1],'lep2_eta':mueta[1],'lep2_Phi':muphi[1],
                                                'Zmass':ZmumuMass,'ZpT':ZpT,'FJetN2b1':ep_fjetN2b1[fjet_index],'FJetN2b2':ep_fjetN2b2[fjet_index],'FJetrho':fatjet_rho,'min_dphi_jets':mindPhi_ak4_ak8_zcr,
                                                'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':recoilweight,'recoilweight_up':recoil_up,'recoilweight_down':recoil_down,'recoilRes_up':ZmumuRecoilResUp,'recoilRes_down':ZmumuRecoilResDown,'recoilEn_up':ZmumuRecoilEnUp,'recoilEn_down':ZmumuRecoilEnDown,
                                                'btagweight':btagweight_B,'btagweight_up':btagweight_B_up,'btagweight_down':btagweight_B_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)


            '''
            --------------------------------------------------------------------------------
            ZCR EleREGION BOOSTED
            --------------------------------------------------------------------------------
            '''


            if isBoostedCRZee:
                fjet_index           = FatJet_ZCR_index[0]#pass_nfjetIndex[0]
                ele1_index           = ele_loose_index[0]
                ele2_index           = ele_loose_index[1]
                ZpT = math.sqrt( (ep_elePx[ele1_index] + ep_elePx[ele2_index])**2 + (ep_elePy[ele1_index]+ep_elePy[ele2_index])**2 )
                #print 'ep_fjetSDMass',ep_fjetSDMass[fjet_index],'fatjetpt',fatjetpt[fjet_index]
                if ep_fjetSDMass[fjet_index]==0:fatjet_rho = math.log((0.0001*0.0001)/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))
                else:fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))
                #print 'fatjet_rho',fatjet_rho

                ele1_loose = ele_loose_index[0]
                ele2_loose = ele_loose_index[1]
                if isTightEles[0] and not isData:
		    weightele1,ele1_up,ele1_down=(wgt.ele_weight(elept[ele1_loose],eleeta[ele1_loose],'T'))
                    if isTightEles[1]:weightele2,ele2_up,ele2_down=(wgt.ele_weight(elept[ele2_loose],eleeta[ele2_loose],'T'))
                    else:weightele2,ele2_up,ele2_down=(wgt.ele_weight(elept[ele2_loose],eleeta[ele2_loose],'L'))
                    lepTrigSF,lepTrigSF_up,lepTrigSF_down = wgt.eletrig_weight(elept[ele1_loose],eleeta[ele1_loose])
                    lepweight = weightele1*weightele2*lepTrigSF
                    lepweight_up = ele1_up*ele2_up*lepTrigSF_up
                    lepweight_down = ele1_down*ele2_down*lepTrigSF_down

                elif isTightEles[1] and not isData:
                    lepTrigSF,lepTrigSF_up,lepTrigSF_down = wgt.eletrig_weight(elept[ele2_loose],eleeta[ele2_loose])
		    weightele1,ele1_up,ele1_down = (wgt.ele_weight(elept[ele1_loose],eleeta[ele1_loose],'L'))
                    weightele2,ele2_up,ele2_down  = (wgt.ele_weight(elept[ele2_loose],eleeta[ele2_loose],'T'))
                    lepweight = weightele1*weightele2*lepTrigSF
                    lepweight_up = ele1_up*ele2_up*lepTrigSF_up
                    lepweight_down = ele1_down*ele2_down*lepTrigSF_down

                
                if not isData: 
                    ewkweight_up=ewkweight*1.5;ewkweight_down=ewkweight*0.5
		    weight          = commanweight_B*lepweight
                    JEC_up,JEC_down = getJECWeight(pass_ak4jet_index_cleaned,ep_THINjetCorrUnc,index=True)
                
                #print 'Zee_B','B_weight',B_weight,'weight',weight

                df_out_Zee_boosted    = df_out_Zee_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZeeRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjet_index], 'FJetEta':fatjeteta[fjet_index], 'FJetPhi':fatjetphi[fjet_index], 'FJetCSV':ep_fjetProbHbb[fjet_index],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjet_index], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':ZeeRecoil_dPhi,
                                                'lep1_pT':elept[0],'lep1_eta':eleeta[0],'lep1_Phi':elephi[0],
                                                'lep2_pT':elept[1],'lep2_eta':eleeta[1],'lep2_Phi':elephi[1],
                                                'Zmass':ZeeMass,'ZpT':ZpT,'FJetN2b1':ep_fjetN2b1[fjet_index],'FJetN2b2':ep_fjetN2b2[fjet_index],'FJetrho':fatjet_rho,'min_dphi_jets':mindPhi_ak4_ak8_zcr,
						'weight':weight,'puweight':PUweight,'puweight_up':PUweight_up,'puweight_down':PUweight_down,'lepweight':lepweight,'lepweight_up':lepweight_up,'lepweight_down':lepweight_down,
                                                'recoilweight':1.0,'recoilweight_up':1.0,'recoilweight_down':1.0,'recoilRes_up':ZeeRecoilResUp,'recoilRes_down':ZeeRecoilResDown,'recoilEn_up':ZeeRecoilEnUp,'recoilEn_down':ZeeRecoilEnDown,
                                                'btagweight':btagweight_B,'btagweight_up':btagweight_B_up,'btagweight_down':btagweight_B_down,'ewkweight':ewkweight,'ewkweight_up':ewkweight_up,'ewkweight_down':ewkweight_down,
                                                'toppTweight':toppTweight,'toppTweight_up':toppTweight_up,'toppTweight_down':toppTweight_down,'jec':1.0,'jec_up':JEC_up,'jec_down':JEC_down
                                           },
                                                ignore_index=True)


    outfilenameis=outfilename
    #result = df_out_Tope_boosted.empty
    #if not result:df_out_Tope_boosted.fillna(0.0)
    df_out_SR_resolved.to_root(outfilenameis, key='monoHbb_SR_resolved',mode='w')
    df_out_SBand_resolved.to_root(outfilenameis, key='monoHbb_SBand_resolved',mode='a')

    df_out_Tope_resolved.to_root(outfilenameis, key='monoHbb_Tope_resolved',mode='a')
    df_out_Topmu_resolved.to_root(outfilenameis, key='monoHbb_Topmu_resolved',mode='a')
    df_out_We_resolved.to_root(outfilenameis, key='monoHbb_We_resolved',mode='a')
    df_out_Wmu_resolved.to_root(outfilenameis, key='monoHbb_Wmu_resolved',mode='a')
    df_out_Zmumu_resolved.to_root(outfilenameis, key='monoHbb_Zmumu_resolved',mode='a')
    df_out_Zee_resolved.to_root(outfilenameis, key='monoHbb_Zee_resolved',mode='a')


    df_out_SR_boosted.to_root(outfilenameis, key='monoHbb_SR_boosted',mode='a')
    df_out_SBand_boosted.to_root(outfilenameis, key='monoHbb_SBand_boosted',mode='a')
    df_out_Tope_boosted.to_root(outfilenameis, key='monoHbb_Tope_boosted', mode='a')
    df_out_Topmu_boosted.to_root(outfilenameis, key='monoHbb_Topmu_boosted', mode='a')
    df_out_We_boosted.to_root(outfilenameis, key='monoHbb_We_boosted', mode='a')
    df_out_Wmu_boosted.to_root(outfilenameis, key='monoHbb_Wmu_boosted', mode='a')
    df_out_Zmumu_boosted.to_root(outfilenameis, key='monoHbb_Zmumu_boosted', mode='a')
    df_out_Zee_boosted.to_root(outfilenameis, key='monoHbb_Zee_boosted', mode='a')
    df_out_TopWmu_boosted.to_root(outfilenameis, key='monoHbb_TopWmu_boosted',mode='a')
    df_out_TopWe_boosted.to_root(outfilenameis, key='monoHbb_TopWe_boosted',mode='a')

    WCR_CutFlow = {1:"Total",2:"preselection",3:"trigger",4:"lep",5:"lepVeto",6:"Recoil",7:"MET",8:"nJets",9:"nBjets",10:"Mbb",11:"extra",12:"extra"}
    SR_CutFlow = {1:"Total",2:"preselection",3:"trigger",4:"lepVeto",5:"tauVeto",6:"nPho",7:"MET",8:"nJets",9:"nBjets",10:"Mbb",11:"extra",12:"extra"}

    WCR_CutFlow_B = {1:"Total",2:"preselection",3:"trigger",4:"lep",5:"lepVeto",6:"Recoil",7:"MET",8:"nJets",9:"nBjets",10:"Mbb",11:"extra",12:"extra"}

    ZCR_CutFlow_B = {1:"Total",2:"preselection",3:"trigger",4:"lep",5:"lepVeto",6:"Recoil",7:"MET",8:"nJets",9:"nBjets",10:"Zmass",11:"extra",12:"extra"}

    SR_CutFlow_B = {1:"Total",2:"preselection",3:"trigger",4:"lepVeto",5:"tauVeto",6:"nPho",7:"MET",8:"nJets",9:"nBjets",10:"extra",11:"extra",12:"extra"}

    for i in [1,2,3,4,5,6,7,8,9,10,11,12]:
        h_reg_WenuCR_resolved_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow[i])
        h_reg_WmunuCR_resolved_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow[i])
        h_reg_TopenuCR_resolved_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow[i])
        h_reg_TopmunuCR_resolved_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow[i])
        h_reg_ZeeCR_resolved_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow[i])
        h_reg_ZmumuCR_resolved_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow[i])
	h_reg_SBand_resolved_cutFlow.GetXaxis().SetBinLabel(i,SR_CutFlow[i])
        #h_reg_SBand_boosted_cutFlow.GetXaxis().SetBinLabel(i,SR_CutFlow[i])  
      
    for i in [1,2,3,4,5,6,7,8,9,10,11]:
        h_reg_WenuCR_boosted_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow_B[i])
        h_reg_WmunuCR_boosted_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow_B[i])
        h_reg_TopenuCR_boosted_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow_B[i])
        h_reg_TopmunuCR_boosted_cutFlow.GetXaxis().SetBinLabel(i,WCR_CutFlow_B[i])
        h_reg_ZeeCR_boosted_cutFlow.GetXaxis().SetBinLabel(i,ZCR_CutFlow_B[i])
        h_reg_ZmumuCR_boosted_cutFlow.GetXaxis().SetBinLabel(i,ZCR_CutFlow_B[i])
        h_reg_SBand_boosted_cutFlow.GetXaxis().SetBinLabel(i,SR_CutFlow_B[i])


    h_reg_WenuCR_resolved_cutFlow.SetEntries(1)
    h_reg_WmunuCR_resolved_cutFlow.SetEntries(1)
    h_reg_TopenuCR_resolved_cutFlow.SetEntries(1)
    h_reg_TopmunuCR_resolved_cutFlow.SetEntries(1)
    h_reg_ZeeCR_resolved_cutFlow.SetEntries(1)
    h_reg_ZmumuCR_resolved_cutFlow.SetEntries(1)

    h_reg_SBand_resolved_cutFlow.SetEntries(1)

    h_reg_WenuCR_boosted_cutFlow.SetEntries(1)
    h_reg_WmunuCR_boosted_cutFlow.SetEntries(1)
    h_reg_TopenuCR_boosted_cutFlow.SetEntries(1)
    h_reg_ZeeCR_boosted_cutFlow.SetEntries(1)
    h_reg_ZmumuCR_boosted_cutFlow.SetEntries(1)
    h_reg_SBand_boosted_cutFlow.SetEntries(1)




    outfile = TFile(outfilenameis,'UPDATE')
    outfile.cd()

    h_reg_WenuCR_resolved_cutFlow.Write()
    h_reg_WmunuCR_resolved_cutFlow.Write()
    h_reg_TopenuCR_resolved_cutFlow.Write()
    h_reg_TopmunuCR_resolved_cutFlow.Write()
    h_reg_SBand_resolved_cutFlow.Write()
    h_reg_SBand_boosted_cutFlow.Write()

    h_reg_ZeeCR_resolved_cutFlow.Write()
    h_reg_ZmumuCR_resolved_cutFlow.Write()

    h_reg_WenuCR_boosted_cutFlow.Write()
    h_reg_WmunuCR_boosted_cutFlow.Write()
    h_reg_TopenuCR_boosted_cutFlow.Write()
    h_reg_ZeeCR_boosted_cutFlow.Write()
    h_reg_ZmumuCR_boosted_cutFlow.Write()

    h_total_mcweight.Write()
    h_total.Write()
    outfile.Write()

    print "output written to ", outfilename
    end = time.clock()
    print "%.4gs" % (end-start)





if __name__ == '__main__':
    path='/afs/cern.ch/work/d/dekumar/public/monoH/Analyzer/CMSSW_10_3_0/src/ExoPieProducer/ExoPieAnalyzer/Files'
    files=glob.glob(path+'/*txt')
    if not runInteractive:
        txtFile=infile
        #for infile in files:
        runbbdm(infile)

    if runInteractive and runOnTxt:
	filesPath = dirName+'/*txt'
	files     = glob.glob(filesPath)
        n = 8 #submit n txt files at a time, make equal to cores
        final = [files[i * n:(i + 1) * n] for i in range((len(files) + n - 1) // n )]
        method1 =True
        if method1:
            try:
                pool = mp.Pool(mp.cpu_count())
                pool.map(runbbdm,files)
                pool.close()
                pool.join()
            except Exception as e:
		print e
                print 'stoping process because of above error'
        else:
            for i in range(len(final)):
                try:
                    pool = mp.Pool(8)
                    pool.map(runbbdm,final[i])
                    pool.close()
                    pool.join()
	        except Exception as e:
		    print e
		    print "Corrupt file inside set of input txt file is detected! Skipping this txt file:  ", final[i]
		    continue#pass
    '''
    if runInteractive and not runOnTxt:
       #following part is for interactive running. This is still under testing because output file name can't be changed at this moment
        inputpath= "/eos/cms/store/user/khurana/test/"

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
        pool.map(runbbdm, allsample)
        ## this works fine but the output file name get same value becuase it is done via a text file at the moment, need to find a better way,
    '''