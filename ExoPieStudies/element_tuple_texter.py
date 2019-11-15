


''' 

this is just an example about how to access the .root file and then convert it into a dataframe. 

Once converted how to looop over it

And how to apply selection on it 

''' 


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


is2016 =  False

filename=''
'''
if not is2016: filename='/eos/cms/store/group/phys_exotica/bbMET/ExoPieElementTuples_20190821/MC/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/190822_052819/0000/ExoPieElementTuples_1.root'

if is2016: filename='/eos/cms/store/group/phys_exotica/bbMET/ExoPieElementTuples/MC_2016miniaodV3_V2/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/191014_204909/0000/ExoPieElementTuples_1.root'
'''

if not is2016: filename='/tmp/khurana/Merged_TT2017.root'

if is2016: filename='/tmp/khurana/Merged_TT2016.root'


## add the file name 



## make  alist of variables you want to access in this code. 
allvariables=[]

if is2016: allvariables= ['nEle', 'eleIsPassTight','THINnJet', 'THINjetPassIDLoose', 'THINjetPx', 'THINjetPy', 'THINjetDeepCSV_b']

if not is2016: allvariables = ['nEle', 'eleIsPassTight','THINnJet','THINjetPassIDTight', 'THINjetPx', 'THINjetPy', 'THINjetDeepCSV_b']

'''
df = read_root(filename, 'tree/treeMaker', columns=allvariables)


df = df[ (df.eleIsPassTight == True).all()] 

print df.nEle, df.eleIsPassTight

print "mean = ", df.nEle.mean()
'''


cutval=1000

if is2016: cutval=0.6321
if not is2016: cutval=0.4941


## first of all open the rootfile in the chunk of 125000 
## this si just to make sure if the root file is big then we read a part of it and then read the next part of it and so on... 
## variables in the list "allvariables" will be read and saved in the df for 125000 events
ievent = 0
nele_sel=[]
njet_sel=[]
sel_csv=[]
sel_nb=[]

for df in read_root(filename, 'tree/treeMaker', columns=allvariables, chunksize=125000):
    
    ## all the code should be inside following for loop. 
    
    ## make a zip of the varianbles which you want to read. This can be same as the allvariables or a subset of it. 
    ## variable in the zip are accessed using the df.VARNAME_YOU_ADDED_IN_THE_LIST
    
    varzip = zip (df.nEle, df.eleIsPassTight, df.THINnJet)
    if is2016:varzip = zip (df.nEle, df.eleIsPassTight, df.THINnJet, df.THINjetPassIDLoose, df.THINjetPx, df.THINjetPy, df.THINjetDeepCSV_b)
    if not is2016:varzip = zip (df.nEle, df.eleIsPassTight, df.THINnJet, df.THINjetPassIDTight, df.THINjetPx, df.THINjetPy, df.THINjetDeepCSV_b)
    
    
    ## now loop of these events 
    ## note the order of the zip and columns should be same. 
    ## met_, fstatus_, fnames_ are the names given by user, the vars will be refered to these names in rest of the code, 
    for nele, tight_eid, njet, tight_jid, jpx, jpy, csv  in varzip:
        ##print nele, tight_eid
        nsel=0
        for i in range(len(tight_eid)): 
            if (tight_eid[i] == True): nsel = nsel+1
        nele_sel.append(nsel)
        nseljet=0
        
        
        jpt = ((jpx * jpx) + (jpy * jpy))
        #print jpt
        
        nbjet=0
        for i in range(len(tight_jid)): 
            if (tight_jid[i] == True):
                if jpt[i] < 900.0 : continue 
                nseljet= nseljet+1
                if csv[i] >cutval:
                    nbjet = nbjet+1
                    
                    sel_csv.append(csv[i])
        sel_nb.append(nbjet)

                
                
                
        
        njet_sel.append(nseljet)
        
        
        ievent = ievent +1 

print "<ele>: ",Series(nele_sel).mean()
print "<jet>: ",Series(njet_sel).mean()
print "<csv>: ", Series(sel_csv).mean()
print "<nbjet> ", Series(sel_nb).mean()


df = DataFrame({ 'nbjets':Series(sel_nb)})
print df

outputfilename=""
if is2016: outputfilename="eleout2016.root"
if not is2016: outputfilename="eleout2017.root"

df.to_root(outputfilename, key='testTree',mode='w')
