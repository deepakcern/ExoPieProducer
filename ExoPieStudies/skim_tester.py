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


is2016 = True

filename=''


if is2016: filename='/eos/cms/store/group/phys_exotica/bbMET/2016_SkimmedFiles/V0/data/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_191014_204909_0000_0.root'

if not is2016: filename='/eos/cms/store/group/phys_exotica/bbMET/2017_skimmedFiles/V0/MC_merged_V0/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8.root'

'''
if is2016: filename='/tmp/khurana/Merged_TT2016_skim,.root'

if not is2016: filename='/eos/cms/store/group/phys_exotica/bbMET/2017_skimmedFiles/V0/MC_merged_V0/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root'
'''
allvariables=[]

if is2016: allvariables= ['st_THINnJet','st_pfMetCorrPt', 'st_nfjet', 'st_nEle', 'st_nPho', 'st_nMu', 'st_nTau_DRBased_EleMuVeto', 'st_nTau_discBased_looseElelooseMuVeto', 'WenuRecoil', 'ZeeRecoil', 'st_THINjetDeepCSV', 'st_THINjetPx','st_THINjetPy']

if not is2016: allvariables = ['st_THINnJet','st_pfMetCorrPt', 'st_nfjet', 'st_nEle', 'st_nPho', 'st_nMu', 'st_nTau_DRBased_EleMuVeto', 'st_nTau_discBased_looseElelooseMuVeto', 'WenuRecoil', 'ZeeRecoil','st_THINjetDeepCSV','st_THINjetPx','st_THINjetPy']


df = read_root(filename, 'outTree', columns=allvariables)


'''
print "<njet>: ",df.st_THINnJet.mean()
print "<st_pfMetCorrPt>: ",df.st_pfMetCorrPt.mean()
print "<st_nfjet>: ",df.st_nfjet.mean()
print "<st_nEle>: ",df.st_nEle.mean()
print "<st_nPho>: ",df.st_nPho.mean()
print "<st_nMu>: ",df.st_nMu.mean()
print "<st_nTau_DRBased_EleMuVeto>: ",df.st_nTau_DRBased_EleMuVeto.mean()
print "<st_nTau_discBased_looseElelooseMuVeto>: ",df.st_nTau_discBased_looseElelooseMuVeto.mean()
print "<WenuRecoil>: ",df.WenuRecoil.mean()
print "<ZeeRecoil>: ",df.ZeeRecoil.mean()
'''


cutval=1000

if is2016: cutval=0.6321
if not is2016: cutval=0.4941

sel_csv=[]
sel_nb=[]
for df in read_root(filename, 'outTree', columns=allvariables, chunksize=125000):
    
    ## all the code should be inside following for loop. 
    
    ## make a zip of the varianbles which you want to read. This can be same as the allvariables or a subset of it. 
    ## variable in the zip are accessed using the df.VARNAME_YOU_ADDED_IN_THE_LIST
    
    varzip = zip (df.st_THINnJet, df.st_THINjetDeepCSV, df.st_THINjetPx, df.st_THINjetPy)
    
    for njet,  csv, px, py   in varzip:
        
        #print len(csv), njet, len(px) 
        
        ## checking jets
        nbjet=0
        for i in range(len(csv)): 
            if csv[i] > cutval : 
                sel_csv.append(csv[i])
                nbjet = nbjet+1
        sel_nb.append(nbjet)
                

print "<csv>: ", Series(sel_csv).mean()
print "<nbjet> ", Series(sel_nb).mean()

df = DataFrame({ 'nbjets':Series(sel_nb)})

outputfilename=""
if is2016: outputfilename="out2016.root"
if not is2016: outputfilename="out2017.root"

df.to_root(outputfilename, key='testTree',mode='w')





df = read_root(filename, 'outTree', columns=allvariables)

print df.st_THINjetDeepCSV

