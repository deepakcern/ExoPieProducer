
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
filename = "/afs/cern.ch/work/k/kuchen/public/skimmedfile/ZZ_TuneCP5_13TeV-pythia8_tmp.root"



df = read_root(filename,'outTree')  

#print df
print df.st_THINjetMinDeltaPhiIdx_Recoil.ee


#for recoils in zip(df.st_THINjetMinDeltaPhiIdx_Recoil.ee):
#    print recoils
    
#root [10] outTree->Scan("st_THINjetMinDeltaPhiIdx_Recoil.mumu","st_THINjetMinDeltaPhiIdx_Recoil.mumu>0")

#print df
'''
df = df[df.MET>300]
df = df[df.MET>500]

print df

fout = TFile("optimisation.root","RECREATE")



stop = time.clock()
print "%.4gs" % (stop-start)
'''
