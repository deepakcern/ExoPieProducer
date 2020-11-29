import os
import sys
import optparse
import glob

usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)
parser.add_option("-v", "--version", type = "string", dest = "Version")
parser.add_option("-i", "--inputDir", dest = "inputDir", default = ".")
parser.add_option("-y", "--year", dest="year", default="Year")
(options, args) = parser.parse_args()


if options.Version == None:
    print('Please provide which version of histograms are being plotted')
    sys.exit()
else:
    histVersion = options.Version
    
if options.inputDir == None:
    print('Please provide the analysis histogram directory')
    sys.exit()
else:
    inDir = options.inputDir
    

sampList2016 = ['SingleElectron-Run2016B-17Jul2018_ver2-v1',
'SingleElectron-Run2016C-17Jul2018-v1',
'SingleElectron-Run2016D-17Jul2018-v1',
'SingleElectron-Run2016E-17Jul2018-v1',
'SingleElectron-Run2016F-17Jul2018-v1',
'SingleElectron-Run2016G-17Jul2018-v1',
'SingleElectron-Run2016H-17Jul2018-v1',
'MET-Run2016B-17Jul2018_ver2-v1',
'MET-Run2016C-17Jul2018-v1',
'MET-Run2016D-17Jul2018-v1',
'MET-Run2016E-17Jul2018-v1',
'MET-Run2016F-17Jul2018-v1',
'MET-Run2016G-17Jul2018-v1',
'MET-Run2016H-17Jul2018-v2',
'ZJetsToNuNu_HT-800To1200_13TeV-madgraph',
'ZJetsToNuNu_HT-600To800_13TeV-madgraph',
'ZJetsToNuNu_HT-400To600_13TeV-madgraph',
'ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph',
'ZJetsToNuNu_HT-200To400_13TeV-madgraph',
'ZJetsToNuNu_HT-1200To2500_13TeV-madgraph',
'ZJetsToNuNu_HT-100To200_13TeV-madgraph',
'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',
'GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1',
'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1',
'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1',
'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
'ZZ_TuneCUETP8M1_13TeV-pythia8',
'WZ_TuneCUETP8M1_13TeV-pythia8',
'WW_TuneCUETP8M1_13TeV-pythia8',
'ttHTobb_M125_13TeV_powheg_pythia8',
'WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8',
'WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8',
'ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8',
'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8',
'ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8',
'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8',
'bbDM_2HDMa_LO_5f_Ma750_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma700_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma50_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma50_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma500_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma500_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma450_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma400_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma400_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma350_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma350_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma300_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma250_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma250_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma200_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma200_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma150_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma150_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma10_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma100_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma1000_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8']
    
sampList2017 = ['MET-Run2017B-31Mar2018-v1',
'MET-Run2017C-31Mar2018-v1',
'MET-Run2017D-31Mar2018-v1',
'MET-Run2017E-31Mar2018-v1',
'MET-Run2017F-31Mar2018-v1',
'SingleElectron-Run2017B-31Mar2018-v1',
'SingleElectron-Run2017C-31Mar2018-v1',
'SingleElectron-Run2017D-31Mar2018-v1',
'SingleElectron-Run2017E-31Mar2018-v1',
'SingleElectron-Run2017F-31Mar2018-v1',
'DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-40to70_TuneCP5_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8',
'DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
'ZJetsToNuNu_HT-800To1200_13TeV-madgraph',
'ZJetsToNuNu_HT-600To800_13TeV-madgraph',
'ZJetsToNuNu_HT-400To600_13TeV-madgraph',
'ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph',
'ZJetsToNuNu_HT-200To400_13TeV-madgraph',
'ZJetsToNuNu_HT-1200To2500_13TeV-madgraph',
'ZJetsToNuNu_HT-100To200_13TeV-madgraph',
'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
'ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8',
'ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8',
'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8',
'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8',
'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',
'QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8',
'QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8',
'QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8',
'QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8',
'QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8',
'QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8',
'QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8',
'QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8',
'GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
'GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8',
'GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
'GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
'GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8',
'ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8',
'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8',
'ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8',
'WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8',
'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8',
'ttHTobb_ttToSemiLep_M125_TuneCP5_13TeV-powheg-pythia8',
'ttHTobb_ttTo2L2Nu_M125_TuneCP5_13TeV-powheg-pythia8',
'ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8',
'VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix',
'GluGluHToBB_M125_13TeV_powheg_pythia8',
'WW_TuneCP5_13TeV-pythia8',
'WZ_TuneCP5_13TeV-pythia8',
'ZZ_TuneCP5_13TeV-pythia8',
'bbDM_2HDMa_LO_5f_Ma750_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma700_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma50_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma50_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma500_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma500_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma450_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma400_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma400_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma350_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma350_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma300_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma250_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma250_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma200_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma200_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma150_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma150_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma10_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma100_MChi1_MA600_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8',
'bbDM_2HDMa_LO_5f_Ma1000_MChi1_MA1200_tanb35_sint_0p7_MH_600_MHC_600_TuneCP3_13TeV-madgraph-pythia8']



if options.year == '2016':
    print('code is running for 2016')
    sampList = sampList2016
elif options.year == '2017':
    print('code is running for 2017')
    sampList = sampList2017
elif options.year == '2018':
    print('code is running for 2018')
    sampList = sampList2018
else:
    print('please provide year')
    sys.exit()
# In[ ]:


outDir = 'analysis_histo_'+histVersion
for fl in list(set(sampList)):
    if not bool(glob.glob(inDir+'/*'+fl+'*')): continue
    filestr = ""
    filestr=inDir+'/*'+fl+'*.root'
    os.system("mkdir -p "+outDir)
    command = "hadd "+outDir+"/"+fl+".root "+filestr
    os.system(command)
print("+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+\n+=+=+=+=+=+=+ DONE +=+=+=+=+=+=+\n+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+")
