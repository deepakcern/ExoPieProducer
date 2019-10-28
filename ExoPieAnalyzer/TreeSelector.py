
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
filename = "analysis.root"



df = read_root(filename,'monoHbb_SR_resolved')

df = df[df.MET>300]
df = df[df.MET>500]

print df

fout = TFile("optimisation.root","RECREATE")



stop = time.clock()
print "%.4gs" % (stop-start)
