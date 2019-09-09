
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


def SetHist(varname, binning):
    h=TH1F()
    if len(binning) == 3:
        h = TH1F(varname, varname, binning[0], binning[1], binning[2])
    else:
        h = TH1F(varname, varname, binning[0], binning[1], binning[2])  ## make it variable binning histogram
    return h


def VarToHist(df, varname, binning ):
    
    df_var = df[varname]
    h_var  = SetHist(varname, binning)
    
    for i in df_var: h_var.Fill(i)
    
    return h_var

df = read_root(filename,'monoHbb_SR_resolved')

h_list = []


h_list.append(VarToHist(df, "MET", [80,200,800]))
h_list.append(VarToHist(df, "DiJetMass", [80,30,200]))
h_list.append(VarToHist(df, "Jet1CSV", [80,-1,1]))
h_list.append(VarToHist(df, "Jet1Eta", [70,-3.5,3.5]))
h_list.append(VarToHist(df, "Jet1Phi", [70,-3.5,3.5]))
h_list.append(VarToHist(df, "Jet1Pt", [80,200,1000]))

fout = TFile("out_hist.root","RECREATE")

for ih in h_list: ih.Write()


print len(h_list)
stop = time.clock()
print "%.4gs" % (stop-start)
