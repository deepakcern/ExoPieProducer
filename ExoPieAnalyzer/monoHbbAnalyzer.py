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
import outvars as out
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
usage = "analyzer for Hbb+DM (debugging) "
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
            isBoostedSBand= False
            isBoostedCRWenu=False
            isBoostedCRWmunu=False
            isBoostedCRZee=False
            isBoostedCRZmumu=False
            isBoostedCRTope=False
            isBoostedCRTopmu=False
            isBoostedSB=False
            
            isResolvedSR=False
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

            FatJet_SBand_index = [index for index in range(ep_nfjet) if ((fatjetpt[index] > 200.0) and (abs(fatjeteta[index])< 2.5) and ((ep_fjetSDMass[index] > 30.0 and ep_fjetSDMass[index] < 100.0) or (ep_fjetSDMass[index] > 150.0 and ep_fjetSDMass[index] < 350.0)) and (ep_fjetProbHbb[index] > 0.86)) ]

            FatJet_ZCR_index   = [index for index in range(ep_nfjet) if ((fatjetpt[index] > 200.0) and (abs(fatjeteta[index])< 2.5) and (ep_fjetProbHbb[index] > 0.86))]
            
            nFatJet_SBand  = len(FatJet_SBand_index)
            nFatJet_ZCR    = len(FatJet_ZCR_index)

            sel_nfjets = len(pass_nfjetIndex)

            #print 'nFatJet_SBand',nFatJet_SBand,'nFatJet_ZCR',nFatJet_ZCR,'sel_nfjets',sel_nfjets

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


            ak4_pt30_eta4p5  = [True for ij in range(ep_THINnJet)] #pt > 30 and eta < 4.5 is already applied at skimmer level
            fjet_pt200_eta_2p5 = [True for ij in range(ep_nfjet)] #pt > 200 and eta < 2.5 is already applied at skimmer level
            pass_ak4jet_index_cleaned = []
            if ak4_pt30_eta4p5 > 0 and fjet_pt200_eta_2p5 > 0 : 
		ak4jet_cleaned_against_fjet = anautil.jetcleaning(ak4_pt30_eta4p5, fjet_pt200_eta_2p5, ak4jeteta, fatjeteta, ak4jetphi, fatjetphi,0.8)
                pass_ak4jet_index_cleaned = boolutil.WhereIsTrue(ak4jet_cleaned_against_fjet)
            #print 'pass_ak4jet_index_cleaned', pass_ak4jet_index_cleaned
		
	    nJets_cleaned = len(pass_ak4jet_index_cleaned)

            Bjet_index = [ij for ij in pass_ak4jet_index_cleaned if (ep_THINjetDeepCSV[ij] > MWP and abs(ak4jeteta[ij]) < 2.5)]
            nBjets = len(Bjet_index)

            
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

            '''
            --------------------------------------------------------------------------
            CLEANED PHOTON COLLECATION
            -------------------------------------------------------------------------
            '''
            #skip photon part for now
            '''
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
            '''
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

            min_ak4jet_MET_dPhi = 1.0
            minDPhi_ak4jet_Werecoil = minDPhi_ak4jet_Wmurecoil = minDPhi_ak4jet_ZeeRecoil = minDPhi_ak4jet_ZmumuRecoil = -10.0
            
            if nJets_cleaned>0 and WerecoildPhi > -10.0:
 	        minDPhi_ak4jet_Werecoil=min([DeltaPhi(WerecoildPhi,ak4jetphi[ijet]) for ijet in pass_ak4jet_index_cleaned])#,ak4jetphi
            if nJets_cleaned>0 and WmurecoildPhi > -10.0:
                minDPhi_ak4jet_Wmurecoil=min([DeltaPhi(WmurecoildPhi,ak4jetphi[ijet]) for ijet in pass_ak4jet_index_cleaned])
            if nJets_cleaned>0 and ZeeRecoil_dPhi > -10.0:
                minDPhi_ak4jet_ZeeRecoil=min([DeltaPhi(ZeeRecoil_dPhi,ak4jetphi[ijet]) for ijet in pass_ak4jet_index_cleaned])
            if nJets_cleaned>0 and ZeeRecoil_dPhi > -10.0:
                minDPhi_ak4jet_ZmumuRecoil=min([DeltaPhi(ZmumuRecoil_dPhi,ak4jetphi[ijet]) for ijet in pass_ak4jet_index_cleaned])
            if nJets_cleaned>0:
                min_ak4jet_MET_dPhi = min([DeltaPhi(ep_pfMetCorrPhi,ak4jetphi[ijet]) for ijet in pass_ak4jet_index_cleaned])
            #print 'min_ak4jet_MET_dPhi', min_ak4jet_MET_dPhi

            '''
            ------------------------------------------------------------------------------
             GET REGION[SR/CR] BASED ON PROPERTIES OF EVENT
            ------------------------------------------------------------------------------
            '''


            # ==== usage: getSel_boosted(nEle,nTightEle,isTightEle,nMu,nTightMu,isTightMuon,nTau,nPho,nBjets,cleaned_ak4jets,nFatJet,nFatJet_SBand,nFatJet_ZCR,pfMet,mini_ak4jet_MET_dPhi,ZeeRecoil,min_ak4jets_ZeeRecoil_dPhi,ZeeMass,ZmumuRecoil,min_ak4jets_ZmumuRecoil_dPhi,ZmumuMass,WenuRecoil,min_ak4jets_WenuRecoil_dPhi,WenuMass,WmunuRecoil,min_ak4jets_WmunuRecoil_dPhi,WmunuMass) ======


            region_boosted = eventSelector_v2.getSel_boosted(ep_nEle,nTightEle,isTightEles,ep_nMu,nTightMu,isTightMuons,ep_HPSTau_n,nPho,nBjets,\
                     nJets_cleaned,sel_nfjets,nFatJet_SBand,nFatJet_ZCR,ep_pfMetCorrPt,min_ak4jet_MET_dPhi,\
                     ZeeRecoil,minDPhi_ak4jet_ZeeRecoil,ZeeMass,ZmumuRecoil,minDPhi_ak4jet_ZmumuRecoil,ZmumuMass,\
                     Werecoil,minDPhi_ak4jet_Werecoil, WeMass,Wmurecoil, minDPhi_ak4jet_Wmurecoil, WmuMass)


            #print 'mettrigdecision',mettrigdecision,'eletrigdecision',eletrigdecision
            if region_boosted['boosted_signal'] and mettrigdecision :  isBoostedSR      = True 
            if region_boosted['boosted_SBand'] and mettrigdecision  : isBoostedSBand    = True 
            if region_boosted['boosted_te']    and  eletrigdecision :  isBoostedCRTope  = True
            if region_boosted['boosted_tm']    and  mettrigdecision :  isBoostedCRTopmu = True
            if region_boosted['boosted_wen']   and  eletrigdecision :  isBoostedCRWenu  = True
            if region_boosted['boosted_wmn']   and  mettrigdecision :  isBoostedCRWmunu = True
            if region_boosted['boosted_zee']   and eletrigdecision  :  isBoostedCRZee   = True
            if region_boosted['boosted_zmm']   and mettrigdecision  :  isBoostedCRZmumu = True

            ''' 
            --------------------------------------------------------------------------------
            SIGNAL REGION BOOSTED
            --------------------------------------------------------------------------------
            '''
           
            weight = 1.0
	    jet1pT=jet1Eta=jet1Phi= -9999.0
            if nJets_cleaned > 0:
                jet1pT           = ak4jetpt[pass_ak4jet_index_cleaned[0]]
                jet1Eta          = ak4jeteta[pass_ak4jet_index_cleaned[0]]
                jet1Phi          = ak4jetphi[pass_ak4jet_index_cleaned[0]]

            if isBoostedSR:
		fjet_index           = pass_nfjetIndex[0]
                min_dPhi_ak4_MET     = min_ak4jet_MET_dPhi
     		
                weightb=weightEWK=weightTop=weightMET=weightPU=1.0

                if ep_genParSample   == 23 and len(ep_genParPt) > 0 : weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0]) 
                elif ep_genParSample == 24 and len(ep_genParPt) > 0 : weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6 and len(ep_genParPt) > 0 : weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightMET=wgt.getMETtrig_First(ep_pfMetCorrPt)
          	if not isData: weight = weightEWK*weightTop*weightPU*weightMET
		#print 'weight', weight
                df_out_SR_boosted = df_out_SR_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
                                                'weight':weight
                                           },
                                                ignore_index=True)


            '''
            --------------------------------------------------------------------------------
            SIDEBAND CONTROL REGION BOOSTED
            --------------------------------------------------------------------------------
            '''

            if isBoostedSBand:
                fjet_index           = FatJet_SBand_index[0]
                min_dPhi_ak4_MET     = min_ak4jet_MET_dPhi

                weightb=weightEWK=weightTop=weightMET=weightPU=1.0

                if ep_genParSample   == 23 and len(ep_genParPt) > 0 : weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample == 24 and len(ep_genParPt) > 0 : weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6 and len(ep_genParPt) > 0 : weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightMET=wgt.getMETtrig_First(ep_pfMetCorrPt)
                if not isData: weight = weightEWK*weightTop*weightPU*weightMET
                #print 'weight', weight
                df_out_SBand_boosted = df_out_SBand_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjet_index], 'FJetEta':fatjeteta[fjet_index], 'FJetPhi':fatjetphi[fjet_index], 'FJetCSV':ep_fjetProbHbb[fjet_index],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjet_index], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
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
                ele1_index           = ele_tight_index[0]
              
                weightb=weightEWK=weightTop=weightPU=weightele=1.0
		ele_trig = True
                if ep_genParSample   == 23 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample == 24 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6 and len(ep_genParPt) > 0:  weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightele=wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],ele_trig,'T')

                if not isData: weight = weightEWK*weightTop*weightPU*weightele
                #print 'weight Tope', weight

                #print 'nPho', nPho 
                df_out_Tope_boosted  = df_out_Tope_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
						'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],
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


                weightb=weightEWK=weightTop=weightRecoil=weightPU=1.0

                if ep_genParSample   == 23 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample == 24 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6 and len(ep_genParPt) > 0 :  weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightRecoil=wgt.getMETtrig_First(Wmurecoil)
                if not isData: weight = weightEWK*weightTop*weightPU*weightRecoil
                #print 'weight Topmu', weight
	
                df_out_Topmu_boosted  = df_out_Topmu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],
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
                ele1_index           = ele_tight_index[0]
               
                weightb=weightEWK=weightTop=weightPU=weightele=1.0
                ele_trig = True
                if ep_genParSample   == 23 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample == 24 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6 and len(ep_genParPt) > 0:  weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightele=wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],ele_trig,'T')

                if not isData: weight = weightEWK*weightTop*weightPU*weightele
                #print 'weight We', weight
 
                df_out_We_boosted  = df_out_We_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],
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
                
                weightb=weightEWK=weightTop=weightRecoil=weightPU=1.0

                if ep_genParSample   == 23 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample == 24 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6  and len(ep_genParPt) > 0:  weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightRecoil=wgt.getMETtrig_First(Wmurecoil)
                if not isData: weight = weightEWK*weightTop*weightPU*weightRecoil



                df_out_Wmu_boosted  = df_out_Wmu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],
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
               
                ZpT = math.sqrt( (ep_muPx[0] + ep_muPx[1])**2 + (ep_muPy[0]+ep_muPy[1])**2 )

                weightb=weightEWK=weightTop=weightRecoil=weightPU=1.0

                if ep_genParSample   == 23 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample == 24 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6  and len(ep_genParPt) > 0:  weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightRecoil=wgt.getMETtrig_First(ZmumuRecoil)
                if not isData: weight = weightEWK*weightTop*weightPU*weightRecoil



                df_out_Zmumu_boosted  = df_out_Zmumu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZmumuRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjet_index], 'FJetEta':fatjeteta[fjet_index], 'FJetPhi':fatjetphi[fjet_index], 'FJetCSV':ep_fjetProbHbb[fjet_index],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjet_index], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
                                                'lep1_pT':mupt[0],'lep1_eta':mueta[0],'lep1_Phi':muphi[0],
                                                'lep2_pT':mupt[1],'lep2_eta':mueta[1],'lep2_Phi':muphi[1],
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

                weightb=weightEWK=weightTop=weightPU=weightele=1.0

                if ep_genParSample   == 23 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample == 24 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6  and len(ep_genParPt) > 0:  weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                ele_trig = True
                no_ele_trig = False
                weightele=(wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],ele_trig,'T'))*(wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],no_ele_trig,'L'))
                if not isData: weight = weightEWK*weightTop*weightPU*weightele



                df_out_Zee_boosted    = df_out_Zee_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':ZeeRecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjet_index], 'FJetEta':fatjeteta[fjet_index], 'FJetPhi':fatjetphi[fjet_index], 'FJetCSV':ep_fjetProbHbb[fjet_index],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjet_index], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
                                                'lep1_pT':elept[0],'lep1_eta':eleeta[0],'lep1_Phi':elephi[0],
                                                'lep2_pT':elept[1],'lep2_eta':eleeta[1],'lep2_Phi':elephi[1],
                                                'Zmass':ZeeMass,'ZpT':ZpT,
                                                'weight':weight
                                           },
                                                ignore_index=True)



            '''
            --------------------------------------------------------------------------------
            Top+W combined MuREGION BOOSTED
            --------------------------------------------------------------------------------
            '''


            if isBoostedCRWmunu or isBoostedCRTopmu:
                fjet_index           = pass_nfjetIndex[0]
                min_dPhi_ak4_Recoil  = minDPhi_ak4jet_Wmurecoil
                muon1_index          = muon_tight_index[0]

                weightb=weightEWK=weightTop=weightRecoil=weightPU=1.0

                if ep_genParSample   == 23 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample == 24 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6  and len(ep_genParPt) > 0:  weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightRecoil=wgt.getMETtrig_First(Wmurecoil)
                if not isData: weight = weightEWK*weightTop*weightPU*weightRecoil



                df_out_TopWmu_boosted  = df_out_TopWmu_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Wmurecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
                                                'lep1_pT':mupt[muon1_index],'lep1_eta':mueta[muon1_index],'lep1_Phi':muphi[muon1_index],
                                                'weight':weight
                                           },
                                                ignore_index=True)



            '''
            --------------------------------------------------------------------------------
            Top+W COMBINED EleREGION BOOSTED
            --------------------------------------------------------------------------------
            '''

            if isBoostedCRWenu or isBoostedCRTope:
                fjet_index           = pass_nfjetIndex[0]
                min_dPhi_ak4_Recoil  = minDPhi_ak4jet_Werecoil
                ele1_index           = ele_tight_index[0]

                weightb=weightEWK=weightTop=weightPU=weightele=1.0
                ele_trig = True
                if ep_genParSample   == 23 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample == 24 and len(ep_genParPt) > 0: weightEWK = wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                elif ep_genParSample == 6 and len(ep_genParPt) > 0:  weightTop = wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightele=wgt.ele_weight(elept[ele1_index],eleeta[ele1_index],ele_trig,'T')

                if not isData: weight = weightEWK*weightTop*weightPU*weightele
                #print 'weight We', weight

                df_out_TopWe_boosted  = df_out_TopWe_boosted.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                'MET':ep_pfMetCorrPt,'RECOIL':Werecoil ,'Njets_PassID':ep_THINnJet,
                                                'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                                'FJetPt':fatjetpt[fjetIndex], 'FJetEta':fatjeteta[fjetIndex], 'FJetPhi':fatjetphi[fjetIndex], 'FJetCSV':ep_fjetProbHbb[fjetIndex],
                                                'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet2Pt':dummy,'Jet2Eta':dummy, 'Jet2Phi':dummy,
                                                'FJetMass':ep_fjetSDMass[fjetIndex], 'DiJetPt':dummy, 'DiJetEta':dummy,'nJets':nJets_cleaned,
                                                'lep1_pT':elept[ele1_index],'lep1_eta':eleeta[ele1_index],'lep1_Phi':elephi[ele1_index],
                                                'weight':weight
                                           },
                                                ignore_index=True)


            ''' 
            --------------------------------------------------------------------------------
            SIGNAL REGION RESOLVED
            --------------------------------------------------------------------------------
            '''
            
            h_mass=-1.0
            ## Resolved selection will be applied IFF there is no boosted candidate 
            if not isBoostedSR:
                if ep_THINnJet < 2: continue
                h_mass  = InvMass(ep_THINjetPx[0], ep_THINjetPy[0], ep_THINjetPz[0], ep_THINjetEnergy[0],ep_THINjetPx[1], ep_THINjetPy[1], ep_THINjetPz[1],ep_THINjetEnergy[1])
                #print 'h_mass',h_mass 
                if (ep_THINnJet >=2) and (ep_nEle == 0) and (ep_nMu == 0) and (nPho ==0) and (ep_HPSTau_n==0) and (ep_pfMetCorrPt > 200.) and (ep_THINjetDeepCSV[0] > MWP and ep_THINjetDeepCSV[1] > MWP) and (abs(ak4jeteta[0]) < 2.5) and (abs(ak4jeteta[1]) < 2.5) and (h_mass > 100.0 and h_mass < 150.0) :
                    isResolvedSR=True 
                    ## call the proper functions
                    jet1pT=jet2pT=jet3pT=jet1Eta=jet2Eta=jet3Eta=jet1Phi=jet2Phi=jet3Phi=dummy 

                    
                    jet1pT    = ak4jetpt[0]
                    jet2pT    = ak4jetpt[1]
                    if (ep_THINnJet > 2):jet3pT    = ak4jetpt[2]

                    jet1Eta   = ak4jeteta[0]
                    jet2Eta   = ak4jeteta[1]
                    if (ep_THINnJet > 2):jet3Eta   = ak4jeteta[2]

                    jet1Phi = ak4jetphi[0]
                    jet2Phi = ak4jetphi[1]
                    if (ep_THINnJet >2 ):jet3Phi = ak4jetphi[2]
		    #print 'jet1pT', jet1pT,'jet2pT', jet2pT, 'jet3pT',jet3pT
                    additional_jets = ep_THINnJet -2

                    ## cal function for each of them 
                    #h_mass=125.
                    
                    ## cal function for each of them based on pt and eta
                    weightEle=1
                    weightMu=1
                    weightB=1
                    weightTau=1
                    weightEWK=1
                    weightTop=1
                    weightPU=1
                    weightOther=1
                    
                    weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU * weightOther 

                    
            
            if  isResolvedSR: 
                df_out_SR_resolved = df_out_SR_resolved.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId, 
                                               'MET':ep_pfMetCorrPt, 'Njets_PassID':ep_THINnJet, 
                                               'Nbjets_PassID':nBjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':nPho,
                                               'Jet1Pt':jet1pT, 'Jet1Eta':jet1Eta, 'Jet1Phi':jet1Phi, 'Jet1CSV':ep_THINjetDeepCSV[0], 
                                               'Jet2Pt':jet2pT, 'Jet2Eta':jet2Eta, 'Jet2Phi':jet2Phi, 'Jet2CSV':ep_THINjetDeepCSV[1],
                                               'Jet3Pt':jet3pT, 'Jet3Eta':jet3Eta, 'Jet3Phi':jet3Phi, 'Jet3CSV':dummy,
                                               'DiJetMass':h_mass,'nJets':additional_jets,
                                               'weight':weight
                                           },
                                              ignore_index=True)
                
	    

    outfilenameis=outfilename
    #result = df_out_Tope_boosted.empty 
    #if not result:df_out_Tope_boosted.fillna(0.0)
    df_out_SR_resolved.to_root(outfilenameis, key='monoHbb_SR_resolved',mode='w')
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
                pool = mp.Pool()#8)
                pool.map(runbbdm,final[i])
                pool.close()
                pool.join()
	    except Exception as e:
		print e
		print "Corrupt file inside set of input txt file is detected! Skipping this txt file:  ", final[i]
		continue
     
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
    
