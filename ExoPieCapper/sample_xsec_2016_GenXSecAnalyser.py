import os, sys
#//--------------------------------------------------------------------------------------
def getXsec(samplename):
    samplename = str(samplename)
    #Xsec = 1.0
    zjet_fac = 3.0
    if 'DYJetsToLL_M-50_HT-100to200'   in samplename: Xsec  =	147.3
    if 'DYJetsToLL_M-50_HT-1200to2500'   in samplename: Xsec  =	0.1512
    if 'DYJetsToLL_M-50_HT-200to400'   in samplename: Xsec  =	41.04
    if 'DYJetsToLL_M-50_HT-2500toInf'   in samplename: Xsec  =	0.003659
    if 'DYJetsToLL_M-50_HT-400to600'   in samplename: Xsec  =	5.676
    if 'DYJetsToLL_M-50_HT-600to800'   in samplename: Xsec  =	1.359
    if 'DYJetsToLL_M-50_HT-800to1200'   in samplename: Xsec  =	0.623
    if 'GJets_HT-100To200'   in samplename: Xsec  =	9239
    if 'GJets_HT-200To400'   in samplename: Xsec  =	2303
    if 'GJets_HT-400To600'   in samplename: Xsec  =	274.6
    if 'GJets_HT-40To100'   in samplename: Xsec  =	20810
    if 'GJets_HT-600ToInf'   in samplename: Xsec  =	93.53
    if 'QCD_HT1000to1500'   in samplename: Xsec  =	1208
    if 'QCD_HT1500to2000'   in samplename: Xsec  =	120
    if 'QCD_HT2000toInf'   in samplename: Xsec  =	25.27
    if 'QCD_HT500to700'   in samplename: Xsec  =	32160
    if 'QCD_HT700to1000'   in samplename: Xsec  =	6828
    if 'ST_s-channel_4f_leptonDecays'   in samplename: Xsec  =	10.32
    if 'ST_t-channel_antitop_4f_inclusiveDecays'   in samplename: Xsec  =	80.95
    if 'ST_t-channel_top_4f_inclusiveDecays'   in samplename: Xsec  =	 136.02
    if 'ST_tW_antitop_5f_inclusiveDecays'   in samplename: Xsec  =	38.06
    if 'ST_tW_top_5f_inclusiveDecays'   in samplename: Xsec  =	38.09
    if 'TTTo2L2Nu'   in samplename: Xsec  =	687.1*0.105
    if 'TTToHadronic'   in samplename: Xsec  =	687.1*0.457
    if 'TTToSemiLeptonic'   in samplename: Xsec  =	687.1*0.438
    if 'WJetsToLNu_HT-100To200'   in samplename: Xsec  =	1344
    if 'WJetsToLNu_HT-1200To2500'   in samplename: Xsec  =	1.329
    if 'WJetsToLNu_HT-200To400'   in samplename: Xsec  =	359.3
    if 'WJetsToLNu_HT-2500ToInf'   in samplename: Xsec  =	0.03219
    if 'WJetsToLNu_HT-400To600'   in samplename: Xsec  =	48.86
    if 'WJetsToLNu_HT-600To800'   in samplename: Xsec  =	12.03
    if 'WJetsToLNu_HT-800To1200'   in samplename: Xsec  =	5.482
    if 'WW_TuneCUETP8M1_13TeV'   in samplename: Xsec  =	64.3
    if 'WZ_TuneCUETP8M1_13TeV'   in samplename: Xsec  =	23.43
    if 'ZJetsToNuNu_HT-100To200'   in samplename: Xsec  =	93.53*zjet_fac
    if 'ZJetsToNuNu_HT-1200To2500'   in samplename: Xsec  =	0.09543*zjet_fac
    if 'ZJetsToNuNu_HT-200To400'   in samplename: Xsec  =	25.93*zjet_fac
    if 'ZJetsToNuNu_HT-2500ToInf'   in samplename: Xsec  =	0.002304*zjet_fac
    if 'ZJetsToNuNu_HT-400To600'   in samplename: Xsec  =	3.585*zjet_fac
    if 'ZJetsToNuNu_HT-600To800'   in samplename: Xsec  =	0.8548*zjet_fac
    if 'ZJetsToNuNu_HT-800To1200'   in samplename: Xsec  =	1.499*zjet_fac
    if 'ZZ_TuneCUETP8M1'   in samplename: Xsec  =	10.16

    return Xsec
