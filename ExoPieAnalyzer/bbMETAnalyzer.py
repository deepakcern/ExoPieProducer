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
usage = "analyzer for bb+DM (debugging) "
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

debug_ = False



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
        print "running for ", txtfile[0]
        infile_  = TextToList(txtfile[0])
        #key_=txtfile[1]
        outfilename= txtfile[0].split('/')[-1].replace('.root.txt','.root')#prefix+key_+".root"


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

    for df in read_root(filename, 'outTree', columns=var.allvars_bbDM, chunksize=125000):

        for ep_runId, ep_lumiSection, ep_eventId, \
            ep_pfMetCorrPt, ep_pfMetCorrPhi, ep_pfMetUncJetResUp, ep_pfMetUncJetResDown, ep_pfMetUncJetEnUp, ep_pfMetUncJetEnDown, \
            ep_isData, \
            ep_THINnJet, ep_THINjetPx, ep_THINjetPy, ep_THINjetPz, ep_THINjetEnergy, \
            ep_THINjetDeepCSV, ep_THINjetHadronFlavor, \
            ep_THINjetNHadEF, ep_THINjetCHadEF, ep_THINjetCEmEF, ep_THINjetPhoEF, ep_THINjetEleEF, ep_THINjetMuoEF, \
            ep_THINjetCorrUnc, \
            ep_nEle, ep_elePx, ep_elePy, ep_elePz, ep_eleEnergy, \
            ep_eleIsPasepight, ep_eleIsPassLoose, \
            ep_nPho, ep_phoIsPasepight, ep_phoPx, ep_phoPy, ep_phoPz, ep_phoEnergy, \
            ep_nMu, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy, ep_iepightMuon, \
            ep_nTau_discBased_looseElelooseMuVeto,ep_nTau_discBased_looseEleTightMuVeto,ep_nTau_discBased_looseEleTightMuVeto,ep_nTau_discBased_mediumElelooseMuVeto,ep_nTau_discBased_TightEleTightMuVeto,\
            ep_pu_nTrueInt, ep_pu_nPUVert, \
            ep_THINjetNPV, \
            ep_mcweight, ep_genParPt, ep_genParSample, \
            in zip(df.st_runId, df.st_lumiSection, df.st_eventId, \
                   df.st_pfMetCorrPt, df.st_pfMetCorrPhi, df.st_pfMetUncJetResUp, df.st_pfMetUncJetResDown, df.st_pfMetUncJetEnUp, df.st_pfMetUncJetEnDown, \
                   df.st_isData, \
                   df.st_THINnJet, df.st_THINjetPx, df.st_THINjetPy, df.st_THINjetPz, df.st_THINjetEnergy, \
                   df.st_THINjetDeepCSV, df.st_THINjetHadronFlavor, \
                   df.st_THINjetNHadEF, df.st_THINjetCHadEF, df.st_THINjetCEmEF, df.st_THINjetPhoEF, df.st_THINjetEleEF, df.st_THINjetMuoEF, \
                   df.st_THINjetCorrUnc, \
                   df.st_nEle, df.st_elePx, df.st_elePy, df.st_elePz, df.st_eleEnergy, \
                   df.st_eleIsPassTight, df.st_eleIsPassLoose, \
                   df.st_nPho, df.st_phoIsPassTight, df.st_phoPx, df.st_phoPy, df.st_phoPz, df.st_phoEnergy, \
                   df.st_nMu, df.st_muPx, df.st_muPy, df.st_muPz, df.st_muEnergy, df.st_isTightMuon, \
                   df.st_nTau_discBased_looseElelooseMuVeto,df.st_nTau_discBased_looseEleTightMuVeto,df.st_nTau_discBased_looseEleTightMuVeto,df.st_nTau_discBased_mediumElelooseMuVeto,df.st_nTau_discBased_TightEleTightMuVeto,\
                   df.st_pu_nTrueInt, df.st_pu_nPUVert, \
                   df.st_THINjetNPV, \
                   df.mcweight, df.st_genParPt, df.st_genParSample, \
            ):

            ieve = ieve + 1
            if ieve%1000==0: print "Processed",ieve,"Events"

            isSR1b=False
            is1bCRWenu=False
            is1bCRWmunu=False
            is1bCRZee=False
            is1bCRZmumu=False
            is1bCRTope=False
            is1bCRTopmu=False

            isSR2b=False
            is2bCRWenu=False
            is2bCRWmunu=False
            is2bCRZee=False
            is2bCRZmumu=False
            is2bCRTope=False
            is2bCRTopmu=False

            deepCSV_Med = 0.6321
            '''
            -------------------------------------------------------------------------------
            electron VARS
            -------------------------------------------------------------------------------
            '''
            ep_elePt  = [getPt(ep_elePx[ij], ep_elePy[ij]) for ij in range(ep_nEle)]
            ep_eleEta = [getEta(ep_elePx[ij], ep_elePy[ij], ep_elePz[ij]) for ij in range(ep_nEle)]
            ep_elePhi = [getPhi(ep_elePx[ij], ep_elePy[ij]) for ij in range(ep_nEle)]

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
            ep_THINbjets = [ij for ij in range(ep_THINnJet) if (ep_THINjetDeepCSV[ij] > deepCSV_Med and abs(ep_THINjetEta[ij]) < 2.5)]
            nBjets = len(ep_THINbjets)

            if len(ep_THINjetPt)==0: continue

            min_dPhi_jet_MET = min([DeltaPhi(jet_phi,ep_pfMetCorrPhi) for jet_phi in ep_THINjetPhi])
            # min_dPhi_jet_WenuRecoil = min([DeltaPhi(jet_phi,WenuPhi) for jet_phi in ep_THINjetPhi])
            # min_dPhi_jet_WmunuRecoil = min([DeltaPhi(jet_phi,WmunuPhi) for jet_phi in ep_THINjetPhi])
            # min_dPhi_jet_ZeeRecoil = min([DeltaPhi(jet_phi,ZeePhi) for jet_phi in ep_THINjetPhi])
            # min_dPhi_jet_ZmumuRecoil = min([DeltaPhi(jet_phi,ZmumuPhi) for jet_phi in ep_THINjetPhi])
            # min_dPhi_jet_GammaRecoil = min([DeltaPhi(jet_phi,GammaPhi) for jet_phi in ep_THINjetPhi])

            '''
            --------------------------------------------------------------------------------
            1b SIGNAL REGION
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for 1b SR.
            if (ep_THINnJet ==1 or ep_THINnJet==2 ) and (ep_THINjetPt[0] > 50.) and (nBjets==1) and (ep_nEle == 0) and (ep_nMu == 0) and (ep_nPho ==0) and (ep_nTau_discBased_looseElelooseMuVeto==0) and (ep_pfMetCorrPt > 200.) and (min_dPhi_jet_MET > 0.5) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                isSR1b=True
                ## cal function for each of them based on pt and eta
                weightMET=wgt.getMETtrig_First(ep_pfMetCorrPt)
                weightEle=1
                weightMu=1
                weightB=1
                weightTau=1
                if ep_genParSample==23:
                    weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0
                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightOther=1

                weight = weightMET*weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU * weightOther
            '''
            --------------------------------------------------------------------------------
            2b SIGNAL REGION
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for 2b SR.
            if (ep_THINnJet ==2 or ep_THINnJet==3 ) and (ep_THINjetPt[0] > 50.) and (nBjets==2) and (ep_nEle == 0) and (ep_nMu == 0) and (ep_nPho ==0) and (ep_nTau_discBased_looseElelooseMuVeto==0) and (ep_pfMetCorrPt > 200.) and (min_dPhi_jet_MET > 0.5) and (ep_THINjetCHadEF[0] >0.1) and (ep_THINjetNHadEF[0] < 0.8):
                isSR2b=True
                ## cal function for each of them based on pt and eta
                weightMET=wgt.getMETtrig_First(ep_pfMetCorrPt)
                weightEle=1
                weightMu=1
                weightB=1
                weightTau=1
                if ep_genParSample==23:
                    weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
                elif ep_genParSample==24:
                    weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
                if ep_genParSample==6:
                    weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
                else:
                    weightEWK = 1.0
                    weightTop = 1.0

                weightPU=wgt.puweight(ep_pu_nTrueInt)
                weightOther=1

                weight = weightMET*weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU * weightOther

            # '''
            # --------------------------------------------------------------------------------
            # WENU CONTROL REGION 1b
            # --------------------------------------------------------------------------------
            # '''
            # ## place all the selection for Wenu SR.
            # if (ep_THINnJet ==1) and (ep_THINjetPt[0] > 50.) and (nBjets==1) and (ep_nEle == 1) and (ep_nMu == 0) and (ep_nPho ==0) and (ep_HPSTau_n==0) and (ep_pfMetCorrPt > 50.) and (ep_WenuRecoil > 200.) and (ep_Wenumass <= 160) and (min_dPhi_jet_WenuRecoil > 0.5) and (ep_THINjetCHadEF[0]) >0.1 and (ep_THINjetNHadEF[0] < 0.8) and (ep_elePt[0] > 30.):
            #     is1bCRWenu=True
            #     ## cal function for each of them based on pt and eta
            #     weightEle=wgt.ele_weight()
            #     weightMu=1
            #     weightB=1
            #     weightTau=1
            #     if ep_genParSample==23:
            #         weightEWK=wgt.getEWKZ(ep_genParPt[0])*wgt.getQCDZ(ep_genParPt[0])
            #     elif ep_genParSample==24:
            #         weightEWK=wgt.getEWKW(ep_genParPt[0])*wgt.getQCDW(ep_genParPt[0])
            #     if ep_genParSample==6:
            #         weightTop=wgt.getTopPtReWgt(ep_genParPt[0],ep_genParPt[1])
            #     else:
            #         weightEWK = 1.0
            #         weightTop = 1.0
            #     weightPU=wgt.puweight(ep_pu_nTrueInt)
            #     weightOther=1
            #
            #     weight = weightEle * weightMu * weightB * weightTau * weightEWK * weightTop * weightPU * weightOther
            dummy=-9999.0
            if len(ep_THINjetPt)>1:
                Jet2Pt  = ep_THINjetPt[1]; Jet2Eta     = ep_THINjetEta[1]
                Jet2Phi = ep_THINjetPhi[1];Jet2deepCSV = ep_THINjetDeepCSV[1]
            else:
                Jet2Pt  = dummy;Jet2Eta     = dummy
                Jet2Phi = dummy;Jet2deepCSV = dummy

            if len(ep_THINjetPt)>2:
                Jet3Pt  = ep_THINjetPt[2]; Jet3Eta     = ep_THINjetEta[2]
                Jet3Phi = ep_THINjetPhi[2];Jet3deepCSV = ep_THINjetDeepCSV[2]
            else:
                Jet3Pt  = dummy;Jet3Eta     = dummy
                Jet3Phi = dummy;Jet3deepCSV = dummy

            if isSR1b:
                df_out_SR_1b = df_out_SR_1b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                    'MET':ep_pfMetCorrPt,'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_looseElelooseMuVeto,'NEle':ep_nEle,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0],'Jet1Eta':ep_THINjetEta[0],'Jet1Phi':ep_THINjetPhi[0],'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':Jet2Pt,'Jet2Eta':Jet2Eta,'Jet2Phi':Jet2Phi,'Jet2deepCSV':Jet2deepCSV,
                                                    'Jet3Pt':dummy,'Jet3Eta':dummy,'Jet3Phi':dummy,'Jet3deepCSV':dummy,
                                                    'weight':weight
                                                    },ignore_index=True
                                                   )
            if isSR2b:
                df_out_SR_2b = df_out_SR_2b.append({'run':ep_runId, 'lumi':ep_lumiSection, 'event':ep_eventId,
                                                    'MET':ep_pfMetCorrPt,'dPhi_jetMET':min_dPhi_jet_MET,
                                                    'NTau':ep_nTau_discBased_looseElelooseMuVeto,'NEle':ep_nEle,'NMu':ep_nMu, 'nPho':ep_nPho,
                                                    'Njets_PassID':ep_THINnJet,'Nbjets_PassID':nBjets,
                                                    'Jet1Pt':ep_THINjetPt[0], 'Jet1Eta':ep_THINjetEta[0], 'Jet1Phi':ep_THINjetPhi[0], 'Jet1deepCSV':ep_THINjetDeepCSV[0],
                                                    'Jet2Pt':ep_THINjetPt[1], 'Jet2Eta':ep_THINjetEta[1], 'Jet2Phi':ep_THINjetPhi[1], 'Jet2deepCSV':ep_THINjetDeepCSV[1],
                                                    'Jet3Pt':Jet3Pt, 'Jet3Eta':Jet3Eta, 'Jet3Phi':Jet3Phi, 'Jet3deepCSV':Jet3deepCSV,
                                                    'weight':weight
                                                    },ignore_index=True
                                                   )

    outfilenameis=outfilename
    df_out_SR_1b.to_root(outfilenameis, key='bbDM_SR_1b',mode='a')
    df_out_SR_2b.to_root(outfilenameis, key='bbDM_SR_2b',mode='a')

    outfile = TFile(outfilenameis,'UPDATE')
    outfile.cd()
    h_total_mcweight.Write()
    h_total.Write()
    outfile.Write()

    print "output written to ", outfilename
    end = time.clock()
    print "%.4gs" % (end-start)



if __name__ == '__main__':
    if not runInteractive:
        txtFile='signal_sample.txt'#infile
        runbbdm(txtFile)

    if runInteractive and runOnTxt:
        filesPath = dirName+'/*txt'
        files     = glob.glob(filesPath)
        n = 8 #submit n txt files at a time, make equal to cores
        final = [files[i * n:(i + 1) * n] for i in range((len(files) + n - 1) // n )]
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
        inputpath= "/afs/cern.ch/work/p/ptiwari/bb+DM_analysis/ntuple_analysis/CMSSW_10_3_0/src/ExoPieSlimmer/SIG_2016_2HDMa_SkimRootFiles"

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
