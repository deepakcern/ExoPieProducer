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
#import variables as var
import variables_data as var
import outvars_new as out
import eventSelector_v2


## from analysisutils
if isCondor:sys.path.append('ExoPieUtils/scalefactortools/')
else:sys.path.append('../../ExoPieUtils/scalefactortools/')
year_file= open("Year.py","w")
year_file.write('era="2017"')
year_file.close()

import ana_weight as wgt



######################################################################################################
## All import are done before this
######################################################################################################

## ----- start of clock
start = time.clock()





## ----- command line argument
usage = "analyzer for bb+DM (debugging) "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--inputfile",  dest="inputfile",default="myfiles.txt")
parser.add_argument("-inDir", "--inputDir",  dest="inputDir",default=".")
parser.add_argument("-runOnTXT", "--runOnTXT",action="store_true", dest="runOnTXT")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="out.root")
parser.add_argument("-D", "--outputdir", dest="outputdir")
parser.add_argument("-F", "--farmout", action="store_true",  dest="farmout")

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
LWP = 0.1522
MWP = 0.4941

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
    h_total = TH1F('h_total','h_total',2,0,2)
    h_total_mcweight = TH1F('h_total_mcweight','h_total_mcweight',2,0,2)
    for infl in infile_:
	f_tmp = TFile.Open(infl,'READ')
        h_tmp = f_tmp.Get('h_total')
        h_tmp_weight = f_tmp.Get('h_total_mcweight')
        h_total.Add(h_tmp)
        h_total_mcweight.Add(h_tmp_weight)

    filename = infile_
    ieve = 0;icount = 0
    for df in read_root(filename, 'outTree', columns=var.allvars, chunksize=125000):

        for ep_runId, ep_lumiSection, ep_eventId, \
            ep_pfMetCorrPt, ep_pfMetCorrPhi, ep_pfMetUncJetResUp, ep_pfMetUncJetResDown, ep_pfMetUncJetEnUp, ep_pfMetUncJetEnDown, \
            ep_THINnJet, ep_THINjetPx, ep_THINjetPy, ep_THINjetPz, ep_THINjetEnergy, \
            ep_THINjetDeepCSV, ep_THINjetHadronFlavor, \
            ep_THINjetNHadEF, ep_THINjetCHadEF, ep_THINjetCEmEF, ep_THINjetPhoEF, ep_THINjetEleEF, ep_THINjetMuoEF, \
            ep_THINjetCorrUnc, \
            ep_nfjet, ep_fjetPx, ep_fjetPy, ep_fjetPz, ep_fjetEnergy, \
            ep_fjetDoubleSV, ep_fjetProbQCDb, ep_fjetProbHbb, ep_fjetProbQCDc, ep_fjetProbHcc, ep_fjetProbHbbc, ep_fjetProbbbvsLight, \
            ep_fjetProbccvsLight, ep_fjetProbTvsQCD, ep_fjetProbWvsQCD, ep_fjetProbZHbbvsQCD, \
            ep_fjetSDMass, ep_fjetN2b1, ep_fjetN2b2, ep_fjetCHSPRMass, ep_fjetCHSSDMass, \
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
                   df.st_THINjetNHadEF, df.st_THINjetCHadEF, df.st_THINjetCEmEF, df.st_THINjetPhoEF, df.st_THINjetEleEF, df.st_THINjetMuoEF, \
                   df.st_THINjetCorrUnc, \
                   df.st_nfjet, df.st_fjetPx, df.st_fjetPy, df.st_fjetPz, df.st_fjetEnergy, \
                   df.st_fjetDoubleSV, df.st_fjetProbQCDb, df.st_fjetProbHbb, df.st_fjetProbQCDc, df.st_fjetProbHcc, df.st_fjetProbHbbc, df.st_fjetProbbbvsLight, \
                   df.st_fjetProbccvsLight, df.st_fjetProbTvsQCD, df.st_fjetProbWvsQCD, df.st_fjetProbZHbbvsQCD, \
                   df.st_fjetSDMass, df.st_fjetN2b1, df.st_fjetN2b2, df.st_fjetCHSPRMass, df.st_fjetCHSSDMass, \
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
            if ieve%1000==0: print "Processed",ieve,"Events"

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

            FatJet_SBand_index = [index for index in range(ep_nfjet) if ((fatjetpt[index] > 200.0) and (abs(fatjeteta[index])< 2.5)) and ((ep_fjetSDMass[index] > 30.0) and (ep_fjetSDMass[index] < 100.0) or ((ep_fjetSDMass[index] > 150.0) and (ep_fjetSDMass[index] < 350.0) )) and (ep_fjetProbHbb[index] > 0.86)]

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
            #if ep_THINnJet < 2: continue
       
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

            Bjet_index = [ij for ij in pass_ak4jet_index_cleaned if (ep_THINjetDeepCSV[ij] > MWP and abs(ak4jeteta[ij]) < 2.5)]
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

            nPho = ep_nPho #len(pass_pho_index_cleaned)


            '''
            ------------------------------------------------------------------------------
            HADRONIC RECOIL
            ------------------------------------------------------------------------------
            '''

            #======   usage: WRecoil_Phi_Wmass(nEle,elept,elephi,elepx_,elepy_,met_,metphi_) =======
	    Werecoil, WerecoildPhi, WeMass = eventSelector_v2.WRecoil_Phi_Wmass(nEle_loose,elept,elephi,ep_elePx,ep_elePy,ep_pfMetCorrPt, ep_pfMetCorrPhi)
            Wmurecoil, WmurecoildPhi, WmuMass = eventSelector_v2.WRecoil_Phi_Wmass(ep_nMu, mupt, muphi, ep_muPx, ep_muPy,ep_pfMetCorrPt, ep_pfMetCorrPhi)

            #======   usage: ZRecoil_Phi_Zmass(nEle, eleCharge_, elepx_, elepy_, elepz_, elee_,met_,metphi_)=====
            ZeeRecoil,ZeeRecoil_dPhi,ZeeMass = eventSelector_v2.ZRecoil_Phi_Zmass(nEle_loose, ep_eleCharge, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy,ep_pfMetCorrPt, ep_pfMetCorrPhi)
            ZmumuRecoil,ZmumuRecoil_dPhi,ZmumuMass = eventSelector_v2.ZRecoil_Phi_Zmass(ep_nMu, ep_muCharge, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy,ep_pfMetCorrPt, ep_pfMetCorrPhi)

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
            
           
            weight = PUweight = lepweight = recoilweight = recoilweight = btagweight = ewkweight = toppTweight =1.0

            btagweight           = wgt.getBTagSF(ep_THINnJet,ak4jetpt,ak4jeteta,ep_THINjetHadronFlavor,ep_THINjetDeepCSV)

            if ep_genParSample   == 23 and len(ep_genParPt) > 0 : ewkweight = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
            if ep_genParSample == 24 and len(ep_genParPt) > 0 : ewkweight = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
            if ep_genParSample == 6 and len(ep_genParPt) > 0  : toppTweight = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

	    PUweight=wgt.puweight(ep_pu_nTrueInt)
            commanweight = btagweight*ewkweight*toppTweight*PUweight
            commanweight_nobtag = ewkweight*toppTweight*PUweight
            #print 'ewkweight',ewkweight
	    additional_jets=ep_THINnJet-2

            '''
            --------------------------------------------------------------------------------
            SIGNAL REGION RESOLVED
            --------------------------------------------------------------------------------
            '''


            if  isResolvedSR:
                METweight=wgt.getMETtrig_First(ep_pfMetCorrPt)

                if not isData: weight = METweight*commanweight
    
                df_out_SR_resolved = df_out_SR_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                               'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                               'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                               'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet1CSV':ep_THINjetDeepCSV[jet1Index],
                                               'Jet2Pt':ak4jetpt[jet2Index], 'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index], 'Jet2CSV':ep_THINjetDeepCSV[jet2Index],
                                               'Jet3Pt':dummy, 'Jet3Eta':dummy, 'Jet3Phi':dummy, 'Jet3CSV':dummy,
                                               'DiJetMass':h_mass,'nJets':additional_jets,
                                               'weight':weight
                                           },
                                              ignore_index=True)

            if  isResolvedSBand:
                METweight=wgt.getMETtrig_First(ep_pfMetCorrPt)

                if not isData: weight = METweight*commanweight

                df_out_SBand_resolved = df_out_SBand_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                               'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                               'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                               'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet1CSV':ep_THINjetDeepCSV[jet1Index],
                                               'Jet2Pt':ak4jetpt[jet2Index], 'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index], 'Jet2CSV':ep_THINjetDeepCSV[jet2Index],
                                               'Jet3Pt':dummy, 'Jet3Eta':dummy, 'Jet3Phi':dummy, 'Jet3CSV':dummy,
                                               'DiJetMass':h_mass,'nJets':additional_jets,
                                               'weight':weight
                                           },
                                              ignore_index=True)


	    if isResolvedCRTope:

                ele1_index    = ele_tight_index[0]
                ele_trig      = True
                lepweight     =wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],ele_trig,'T')

                if not isData: weight = lepweight*commanweight
                #print 'weight',weight,'PUweight',PUweight,'lepweight',lepweight,'recoilweight',recoilweight,'btagweight',btagweight,'ewkweight',ewkweight,'toppTweight',toppTweight
                df_out_Tope_resolved  = df_out_Tope_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WerecoildPhi,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],
                                                'weight':weight,'puweight':PUweight,'lepweight':lepweight,'recoilweight':recoilweight,'btagweight':btagweight,'ewkweight':ewkweight,'toppTweight':toppTweight
                                           },
                                                ignore_index=True)


            if isResolvedCRTopmu:
                muon1_index          = muon_tight_index[0]
                recoilweight=wgt.getMETtrig_First(Wmurecoil)
                mu_trig = False
		lepweight=wgt.mu_weight(mupt[0],mueta[0],mu_trig,'T')

                if not isData: weight = recoilweight*lepweight*commanweight
                #print 'weight',weight,'PUweight',PUweight,'lepweight',lepweight,'recoilweight',recoilweight,'btagweight',btagweight,'ewkweight',ewkweight,'toppTweight',toppTweight
                df_out_Topmu_resolved  = df_out_Topmu_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WmurecoildPhi,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],
                                                'weight':weight,'puweight':PUweight,'lepweight':lepweight,'recoilweight':recoilweight,'btagweight':btagweight,'ewkweight':ewkweight,'toppTweight':toppTweight
                                           },
                                                ignore_index=True)


            if isResolvedCRWenu:
                ele1_index    = ele_tight_index[0]
                ele_trig      = True
                lepweight     =wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],ele_trig,'T')

                if not isData: weight = lepweight*commanweight
                #print 'weight',weight,'PUweight',PUweight,'lepweight',lepweight,'recoilweight',recoilweight,'btagweight',btagweight,'ewkweight',ewkweight,'toppTweight',toppTweight
                df_out_We_resolved  = df_out_We_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WerecoildPhi,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],
                                                'weight':weight,'puweight':PUweight,'lepweight':lepweight,'recoilweight':recoilweight,'btagweight':btagweight,'ewkweight':ewkweight,'toppTweight':toppTweight
                                           },
                                                ignore_index=True)

            if isResolvedCRWmunu:
                muon1_index          = muon_tight_index[0]
                recoilweight=wgt.getMETtrig_First(Wmurecoil)
                mu_trig = False
                lepweight=wgt.mu_weight(mupt[0],mueta[0],mu_trig,'T')

                if not isData: weight = recoilweight*lepweight*commanweight
                #print 'weight',weight,'PUweight',PUweight,'lepweight',lepweight,'recoilweight',recoilweight,'btagweight',btagweight,'ewkweight',ewkweight,'toppTweight',toppTweight
                df_out_Wmu_resolved  = df_out_Wmu_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WmurecoildPhi,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],
                                                'weight':weight,'puweight':PUweight,'lepweight':lepweight,'recoilweight':recoilweight,'btagweight':btagweight,'ewkweight':ewkweight,'toppTweight':toppTweight
                                           },
                                                ignore_index=True)



            if isResolvedCRZee:
                ele1_trig = True; ele2_trig = False
                ele1_loose = ele_loose_index[0]
                ele2_loose = ele_loose_index[1]
                if isTightEles[0]:lepweight=(wgt.ele_weight(elept[ele1_loose],eleeta[ele1_loose],ele1_trig,'T'))*(wgt.ele_weight(elept[ele2_loose],eleeta[ele2_loose],ele2_trig,'L'))
                if isTightEles[1]:lepweight=(wgt.ele_weight(elept[ele2_loose],eleeta[ele2_loose],ele2_trig,'L'))*(wgt.ele_weight(elept[ele1_loose],eleeta[ele1_loose],ele1_trig,'T'))

                ZpT = math.sqrt( (ep_elePx[0] + ep_elePx[0])**2 + (ep_elePy[1]+ep_elePy[1])**2 )

                if not isData: weight = lepweight*commanweight



                df_out_Zee_resolved    = df_out_Zee_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZeeRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
                                                'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':ZeeRecoil_dPhi,
                                                'lep1_pT':elept[ele1_loose],'lep1_eta':eleeta[ele1_loose],'lep1_Phi':elephi[ele1_loose],
                                                'lep2_pT':elept[ele2_loose],'lep2_eta':eleeta[ele2_loose],'lep2_Phi':elephi[ele2_loose],
                                                'Zmass':ZeeMass,'ZpT':ZpT,
                                                'weight':weight,'puweight':PUweight,'lepweight':lepweight,'recoilweight':recoilweight,'btagweight':btagweight,'ewkweight':ewkweight,'toppTweight':toppTweight
                                           },
                                                ignore_index=True)



            if isResolvedCRZmumu:
                mu1_trig= False; mu2_trig=False
		#print 'number of muons',len(mupt)
                if isTightMuons[0]:lepweight=wgt.mu_weight(mupt[0],mueta[0],mu1_trig,'T')*wgt.mu_weight(mupt[1],mueta[1],mu2_trig,'L')
                if isTightMuons[1]:lepweight=wgt.mu_weight(mupt[1],mueta[1],mu2_trig,'L')*wgt.mu_weight(mupt[0],mueta[0],mu1_trig,'T')

                ZpT          = math.sqrt( (ep_muPx[0] + ep_muPx[1])**2 + (ep_muPy[0]+ep_muPy[1])**2 )
                recoilweight = wgt.getMETtrig_First(ZmumuRecoil)

                if not isData: weight = lepweight*recoilweight*commanweight



                df_out_Zmumu_resolved  = df_out_Zmumu_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZmumuRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_notiso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':dummy, 'FJetEta':dummy, 'FJetPhi':dummy, 'FJetCSV':dummy,
						'Jet1Pt':ak4jetpt[jet1Index], 'Jet1Eta':ak4jeteta[jet1Index], 'Jet1Phi':ak4jetphi[jet1Index], 'Jet2Pt':ak4jetpt[jet2Index],'Jet2Eta':ak4jeteta[jet2Index], 'Jet2Phi':ak4jetphi[jet2Index],
                                                'FJetMass':dummy, 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':additional_jets,'min_dPhi':min_ak4jet_MET_dPhi_R,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':ZmumuRecoil_dPhi,
                                                'lep1_pT':mupt[0],'lep1_eta':mueta[0],'lep1_Phi':muphi[0],
                                                'lep2_pT':mupt[1],'lep2_eta':mueta[1],'lep2_Phi':muphi[1],
                                                'Zmass':ZmumuMass,'ZpT':ZpT,
                                                'weight':weight,'puweight':PUweight,'lepweight':lepweight,'recoilweight':recoilweight,'btagweight':btagweight,'ewkweight':ewkweight,'toppTweight':toppTweight
                                           },
                                                ignore_index=True)

            '''
            --------------------------------------------------------------------------------
            SIGNAL REGION BOOSTED
            --------------------------------------------------------------------------------
            '''

            jet1pT=jet1Eta=jet1Phi= -9999.0
            if nJets_cleaned > 0:
                jet1pT           = ak4jetpt[pass_ak4jet_index_cleaned[0]]
                jet1Eta          = ak4jeteta[pass_ak4jet_index_cleaned[0]]
                jet1Phi          = ak4jetphi[pass_ak4jet_index_cleaned[0]]

            if isBoostedSR:
                fjet_index           = pass_nfjetIndex[0]
                min_dPhi_ak4_MET     = min_ak4jet_MET_dPhi
                fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))

                weightMET=wgt.getMETtrig_First(ep_pfMetCorrPt)
                if not isData: weight = commanweight_nobtag*weightMET
                #print 'weight', weight
                df_out_SR_boosted = df_out_SR_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_dPhi_ak4_MET,'met_Phi':ep_pfMetCorrPhi,'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,
                                                'weight':weight
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

                weightMET=wgt.getMETtrig_First(ep_pfMetCorrPt)
                if not isData: weight = commanweight_nobtag*weightMET
                #print 'weight', weight
                df_out_SBand_boosted = df_out_SBand_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjet_index], 'FJetEta':fatjeteta[fjet_index], 'FJetPhi':fatjetphi[fjet_index], 'FJetCSV':ep_fjetProbHbb[fjet_index],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjet_index], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_dPhi_ak4_MET,'met_Phi':ep_pfMetCorrPhi,'FJetN2b1':ep_fjetN2b1[fjet_index],'FJetN2b2':ep_fjetN2b2[fjet_index],'FJetrho':fatjet_rho,
                                                'weight':weight
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
                ele_trig = True
                weightele=wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],ele_trig,'T')

                if not isData: weight = commanweight_nobtag*weightele
                #print 'weight Tope', weight

                #print 'nPho', nPho
                df_out_Tope_boosted  = df_out_Tope_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WerecoildPhi,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,
                                                'weight':weight
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

                mu_trig = False
                lepweight=wgt.mu_weight(mupt[0],mueta[0],mu_trig,'T')

                weightRecoil=wgt.getMETtrig_First(Wmurecoil)
                if not isData: weight = commanweight_nobtag*weightRecoil*lepweight
                #print 'weight Topmu', weight

                df_out_Topmu_boosted  = df_out_Topmu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WmurecoildPhi,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,
                                                'weight':weight
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
                ele_trig=True
                weightele=wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],ele_trig,'T')

                if not isData: weight = commanweight_nobtag*weightele
                

                df_out_We_boosted  = df_out_We_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WerecoildPhi,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,
                                                'weight':weight
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

                mu_trig = False
                lepweight=wgt.mu_weight(mupt[0],mueta[0],mu_trig,'T')

                weightRecoil=wgt.getMETtrig_First(Wmurecoil)
                if not isData: weight = commanweight_nobtag*weightRecoil*lepweight



                df_out_Wmu_boosted  = df_out_Wmu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':WmurecoildPhi,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],'FJetN2b1':ep_fjetN2b1[fjetIndex],'FJetN2b2':ep_fjetN2b2[fjetIndex],'FJetrho':fatjet_rho,
                                                'weight':weight
                                           },
                                                ignore_index=True)


            '''
            --------------------------------------------------------------------------------
            ZCR MuREGION BOOSTED
            --------------------------------------------------------------------------------
            '''


            if isBoostedCRZmumu:
                fjet_index           = FatJet_ZCR_index[0]#pass_nfjetIndex[0]
                fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))
                ZpT = math.sqrt( (ep_muPx[0] + ep_muPx[1])**2 + (ep_muPy[0]+ep_muPy[1])**2 )

                weightRecoil=wgt.getMETtrig_First(ZmumuRecoil)
                if not isData: weight = commanweight_nobtag*weightRecoil



                df_out_Zmumu_boosted  = df_out_Zmumu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZmumuRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjet_index], 'FJetEta':fatjeteta[fjet_index], 'FJetPhi':fatjetphi[fjet_index], 'FJetCSV':ep_fjetProbHbb[fjet_index],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjet_index], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':ZmumuRecoil_dPhi,
                                                'lep1_pT':mupt[0],'lep1_eta':mueta[0],'lep1_Phi':muphi[0],
                                                'lep2_pT':mupt[1],'lep2_eta':mueta[1],'lep2_Phi':muphi[1],'FJetN2b1':ep_fjetN2b1[fjet_index],'FJetN2b2':ep_fjetN2b2[fjet_index],'FJetrho':fatjet_rho,
                                                'Zmass':ZmumuMass,'ZpT':ZpT,
                                                'weight':weight
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
                fatjet_rho = math.log((ep_fjetSDMass[fjet_index]*ep_fjetSDMass[fjet_index])/(fatjetpt[fjet_index]*fatjetpt[fjet_index]))

                ele_trig = True
                no_ele_trig = False
                weightele=(wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],ele_trig,'T'))*(wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],no_ele_trig,'L'))
                if not isData: weight = commanweight_nobtag*weightele



                df_out_Zee_boosted    = df_out_Zee_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZeeRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets_iso, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjet_index], 'FJetEta':fatjeteta[fjet_index], 'FJetPhi':fatjetphi[fjet_index], 'FJetCSV':ep_fjetProbHbb[fjet_index],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjet_index], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,'min_dPhi':min_ak4jet_MET_dPhi,'met_Phi':ep_pfMetCorrPhi,'RECOIL_Phi':ZeeRecoil_dPhi,
                                                'lep1_pT':elept[0],'lep1_eta':eleeta[0],'lep1_Phi':elephi[0],
                                                'lep2_pT':elept[1],'lep2_eta':eleeta[1],'lep2_Phi':elephi[1],
                                                'Zmass':ZeeMass,'ZpT':ZpT,'FJetN2b1':ep_fjetN2b1[fjet_index],'FJetN2b2':ep_fjetN2b2[fjet_index],'FJetrho':fatjet_rho,
                                                'weight':weight
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

    outfile = TFile(outfilenameis,'UPDATE')
    outfile.cd()
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
        for i in range(len(final)):
            try:
                pool = mp.Pool(8)
                pool.map(runbbdm,files)#final[i])
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
