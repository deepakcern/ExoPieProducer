'''
BuildStatsModel.py: Package to build statistical fitting model for background estimation and limit extraction

Author: Raman Khurana
Date:   26-September-2018
V0:     Simple implementation of the model using transfer factors 

To do:  Write the full package. 
     :  
'''

import os
import sys 


''' import ROOT '''
import ROOT as r


''' import user defined stuff 
    mainly, rootfit, combine etc
'''
from HiggsAnalysis.CombinedLimit.ModelTools import *
ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ");
ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include");
ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libRooFitCore.so")


sys.path.append('/afs/cern.ch/work/k/khurana/RKTools/CommonUtilities/Helpers/')

import  describe as des
import fileutils as fu
#Import = getattr(r.RooWorkspace, 'import')

print des.inputfile, des.outputfile



''' open input file in read mode '''
fin_ = fu.OpenRootFile(des.inputfile)

''' open output root file in recreate mode '''
fout_ = fu.OpenRootFile(des.outputfile, "RECREATE")

'''define the workspace '''
wspace = r.RooWorkspace(des.wsname)


''' A search in a MET tail, define MET as our variable 
    this should be added to the RooArgList once rooreal var is defined
'''

met = r.RooRealVar("met","E_{T}^{miss}",des.met_low, des.met_hi)
vars_ = r.RooArgList(met)

'''
naming convention
Histogram: h_region_process, e.g. h_sr1_data, h_sr2_wjets, h_wenu_wjets
datahist:  use dh_ instead of h_
'''


'''
Get the data histogram from Signal region, this is met in our case. This is a 7 bin histogram in this case. This can be optimised later on.                                         '''
h_sr2_data = fin_.Get("SR_2b_data_obs")
r.SetOwnership(h_sr2_data,False)

'''
convert the histogram into RooDataHist
'''
dh_sr2_data = r.RooDataHist("dh_sr2_data","dh_sr2_data",vars_, h_sr2_data)

'''
Import just created RooDataHist into the workspace.
'''
getattr(wspace, 'import')(vars_)
getattr(wspace,'import')(dh_sr2_data)

wspace.Print()
#wspace.Import(dh_sr2_data)

''' 
----------------------------------------------------------------------------------------------------------------
-----------------------For WJets background---------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------

'''


'''
 copy the signal region WJets histogram
'''
#h_sr2_wjets =  fin_.Get("SR_2b_WJets")

#fu.GetBinContents(h_sr2_wjets)


'''

wjets_SR_content_bin1 = h_sr2_Wjets->GetBinContent(1);
wjets_SR_content_bin2 = h_sr2_Wjets->GetBinContent(2);
wjets_SR_content_bin3 = h_sr2_Wjets->GetBinContent(3);
wjets_SR_content_bin4 = h_sr2_Wjets->GetBinContent(4);
wjets_SR_content_bin5 = h_sr2_Wjets->GetBinContent(5);
wjets_SR_content_bin6 = h_sr2_Wjets->GetBinContent(6);
wjets_SR_content_bin7 = h_sr2_Wjets->GetBinContent(7);
'''




fout_.cd()
wspace.Write()

