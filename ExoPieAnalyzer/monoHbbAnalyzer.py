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
import outvars as out


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



def TextToList(textfile):
    return([iline.rstrip()    for iline in open(textfile)])

## the input file list and key is caught in one variable as a python list, 
#### first element is the list of rootfiles 
#### second element is the key, user to name output.root

def runbbdm(txtfile):
    

    
    print "in main function"
    
    infile_=[]
    outfilename=""
    prefix="Skimmed_"
    ikey_ = ""

    
    if  runInteractive:
        print "running for ", txtfile
        infile_  = TextToList(txtfile[0])
        key_=txtfile[1]
        outfilename= prefix+key_+".root"
    

    
    if not runInteractive:
        print "running for ", txtfile
        infile_  = TextToList(txtfile[0])
        key_=txtfile[1]
    
        ''' old
        prefix="Skimmed_"
        outfilename= prefix+infile_.split("/")[-1]
        '''
    
        outfilename= prefix+key_+".root"
    
        
    if not runInteractive: 
	infile_=TextToList(txtfile)
        prefix_ = '' #'/eos/cms/store/group/phys_exotica/bbMET/2017_skimmedFiles/locallygenerated/'
        if outputdir!='.': prefix_ = outputdir+'/'
        print "prefix_", prefix_
        outfilename = prefix_+txtfile.split('/')[-1].replace('.txt','.root')#"SkimmedTree.root"
        print 'outfilename',  outfilename 
    
    
    ## define global dataframes 
    df_out_SR_resolved = out.df_out_SR_resolved
    #outputfilename = args.outputfile
    h_total = TH1F('h_total','h_total',2,0,2)
    h_total_mcweight = TH1F('h_total_mcweight','h_total_mcweight',2,0,2)


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
            ep_eleIsPasepight, ep_eleIsPassLoose, \
            ep_nPho, ep_phoIsPasepight, ep_phoPx, ep_phoPy, ep_phoPz, ep_phoEnergy, \
            ep_nMu, ep_muPx, ep_muPy, ep_muPz, ep_muEnergy, ep_iepightMuon, \
            ep_HPSTau_n, ep_nTauTightElectron, ep_nTauTightMuon, ep_nTauTightEleMu, ep_nTauLooseEleMu, \
            ep_pu_nTrueInt, ep_pu_nPUVert, \
            ep_THINjetNPV, \
            ep_mcweight, ep_nGenPar, ep_genParId, ep_genMomParId, ep_genParEp, ep_genParPx, ep_genParPy, ep_genParPz, ep_genParEnergy, \
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
                   df.st_eleIsPassTight, df.st_eleIsPassLoose, \
                   df.st_nPho, df.st_phoIsPassTight, df.st_phoPx, df.st_phoPy, df.st_phoPz, df.st_phoEnergy, \
                   df.st_nMu, df.st_muPx, df.st_muPy, df.st_muPz, df.st_muEnergy, df.st_isTightMuon, \
                   df.st_HPSTau_n, df.st_nTauTightElectron, df.st_nTauTightMuon, df.st_nTauTightEleMu, df.st_nTauLooseEleMu, 
                   df.st_pu_nTrueInt, df.st_pu_nPUVert, \
                   df.st_THINjetNPV, \
                   df.mcweight, df.st_nGenPar, df.st_genParId, df.st_genMomParId, df.st_genParSt, df.st_genParPx, df.st_genParPy, df.st_genParPz, df.st_genParEnergy, \
                   
            ):
            
            #ep_WenuRecoil, ep_Wenumass, ep_WenuPhi, ep_WmunuRecoil, ep_Wmunumass, ep_WmunuPhi, 
            #ep_ZeeRecoil, ep_ZeeMass, ep_ZeePhi, ep_ZmumuRecoil, ep_ZmumuMass, ep_ZmumuPhi,
            #ep_GammaRecoil, ep_GammaPhi, 
            #df.WenuRecoil, df.Wenumass, df.WenuPhi, df.WmunuRecoil, df.Wmunumass, df.WmunuPhi, 
            #df.ZeeRecoil, df.ZeeMass, df.ZeePhi, df.ZmumuRecoil, df.ZmumuMass, df.ZmumuPhi, 
            #df.GammaRecoil, df.GammaPhi 
            
            
            ieve = ieve + 1
            if ieve%1000==0: print "Processed",ieve,"Events"
            
            isBoostedSR=False
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
            --------------------------------------------------------------------------------
            SIGNAL REGION BOOSTED
            --------------------------------------------------------------------------------
            '''
            ## place all the selection for boosted SR. 
            if  ( ep_nfjet == 1) and (ep_nEle == 0) and (ep_nMu == 0) and (ep_nPho ==0) and (ep_HPSTau_n==0) and (ep_pfMetCorrPt > 200.):
                isBoostedSR=True
                
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
            
                
            ''' 
            --------------------------------------------------------------------------------
            SIGNAL REGION RESOLVED
            --------------------------------------------------------------------------------
            '''
            
            
            
            ep_THINjetPt=[]
            ep_THINjetEta=[]
            ep_THINjetPhi=[]
            n_bjets=-1
            h_mass=-1
            ## Resolved selection will be applied IFF there is no boosted candidate 
            if not isBoostedSR:
                if (ep_THINnJet>2):
                    #if  ( ep_nfjet == 1) and (ep_nEle == 0) and (ep_nMu == 0) and (ep_nPho ==0) and (ep_HPSTau_n==0) and (ep_pfMetCorrPt > 170.) :
                    isResolvedSR=True 
                    ## call the proper functions 
                    
                    ep_THINjetPt=[1 for ijet in range(ep_THINnJet)]
                    ep_THINjetEta=[1 for ijet in range(ep_THINnJet)]
                    ep_THINjetPhi=[1 for ijet in range(ep_THINnJet)]
                    
                    ## cal function for each of them 
                    n_bjets=2
                    h_mass=125.
                    
                    
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
                                               'Nbjets_PassID':n_bjets, 'NTauJets':ep_HPSTau_n, 'NEle':ep_nEle, 'NMu':ep_nMu, 'nPho':ep_nPho,
                                               'Jet1Pt':ep_THINjetPt[0], 'Jet1Eta':ep_THINjetEta[0], 'Jet1Phi':ep_THINjetPhi[0], 'Jet1CSV':ep_THINjetDeepCSV[0], 
                                               'Jet2Pt':ep_THINjetPt[1], 'Jet2Eta':ep_THINjetEta[1], 'Jet2Phi':ep_THINjetPhi[1], 'Jet2CSV':ep_THINjetDeepCSV[1],
                                               'Jet3Pt':ep_THINjetPt[2], 'Jet3Eta':ep_THINjetEta[2], 'Jet3Phi':ep_THINjetPhi[2], 'Jet3CSV':ep_THINjetDeepCSV[2],
                                               'DiJetMass':h_mass,
                                               'weight':weight
                                           },
                                              ignore_index=True)
                


    outfilenameis=outfilename
    df_out_SR_resolved.to_root(outfilenameis, key='monoHbb_SR_resolved')
    
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
        txtFile=infile
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
    
