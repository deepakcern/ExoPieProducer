import os, sys
#//--------------------------------------------------------------------------------------
def getXsec(samplename):
    samplename = str(samplename)
    if 'DYJetsToLL_M-50_HT-100to200'   in samplename: Xsec  =	161
    elif 'DYJetsToLL_M-50_HT-1200to2500'   in samplename: Xsec  =	0.1927
    elif 'DYJetsToLL_M-50_HT-200to400'   in samplename: Xsec  =	48.58
    elif 'DYJetsToLL_M-50_HT-2500toInf'   in samplename: Xsec  =	0.003478
    elif 'DYJetsToLL_M-50_HT-400to600'   in samplename: Xsec  =	6.982
    elif 'DYJetsToLL_M-50_HT-600to800'   in samplename: Xsec  =	1.747
    elif 'DYJetsToLL_M-50_HT-800to1200'   in samplename: Xsec  =	0.8052
    elif 'GJets_HT-100To200'   in samplename: Xsec  =	8622
    elif 'GJets_HT-200To400'   in samplename: Xsec  =	2193
    elif 'GJets_HT-400To600'   in samplename: Xsec  =	258.5
    elif 'GJets_HT-40To100'   in samplename: Xsec  =	18620
    elif 'GJets_HT-600ToInf'   in samplename: Xsec  =	85.21
    elif 'QCD_HT1000to1500'   in samplename: Xsec  =	1094
    elif 'QCD_HT1500to2000'   in samplename: Xsec  =	98.99
    elif 'QCD_HT2000toInf'   in samplename: Xsec  =	20.23
    elif 'QCD_HT500to700'   in samplename: Xsec  =	29990
    elif 'QCD_HT700to1000'   in samplename: Xsec  =	6351
    elif 'QCD_HT200to300' in samplename: Xsec  =  1735000
    elif 'QCD_HT300to500' in samplename: Xsec  =  366800
    elif 'ST_s-channel_4f_leptonDecays'   in samplename: Xsec  =	10.32
    elif 'ST_t-channel_antitop_4f_inclusiveDecays'   in samplename: Xsec  =	80.95
    elif 'ST_t-channel_top_4f_inclusiveDecays'   in samplename: Xsec  =	 136.02
    elif 'ST_tW_antitop_5f_inclusiveDecays'   in samplename: Xsec  =	38.06
    elif 'ST_tW_top_5f_inclusiveDecays'   in samplename: Xsec  =	38.09
    elif 'TTTo2L2Nu'   in samplename: Xsec  =	687.1*0.105
    elif 'TTToHadronic'   in samplename: Xsec  =	687.1*0.457
    elif 'TTToSemiLeptonic'   in samplename: Xsec  =	687.1*0.438
    elif 'WJetsToLNu_HT-100To200'   in samplename: Xsec  =	1395
    elif 'WJetsToLNu_HT-1200To2500'   in samplename: Xsec  =	1.08
    elif 'WJetsToLNu_HT-200To400'   in samplename: Xsec  =	409.3
    elif 'WJetsToLNu_HT-2500ToInf'   in samplename: Xsec  =	0.008053
    elif 'WJetsToLNu_HT-400To600'   in samplename: Xsec  =	57.91
    elif 'WJetsToLNu_HT-600To800'   in samplename: Xsec  =	12.93
    elif 'WJetsToLNu_HT-800To1200'   in samplename: Xsec  =	5.395
    elif 'WW_TuneCP5_13TeV'   in samplename: Xsec  =	75.9
    elif 'WZ_TuneCP5_13TeV'   in samplename: Xsec  =	27.57
    elif 'ZJetsToNuNu_HT-100To200'   in samplename: Xsec  =	304.5
    elif 'ZJetsToNuNu_HT-1200To2500'   in samplename: Xsec  =	0.343
    elif 'ZJetsToNuNu_HT-200To400'   in samplename: Xsec  =	91.85
    elif 'ZJetsToNuNu_HT-2500ToInf'   in samplename: Xsec  =	0.005146
    elif 'ZJetsToNuNu_HT-400To600'   in samplename: Xsec  =	13.11
    elif 'ZJetsToNuNu_HT-600To800'   in samplename: Xsec  =	3.257
    elif 'ZJetsToNuNu_HT-800To1200'   in samplename: Xsec  =	1.499
    elif 'ZZ_TuneCP5'   in samplename: Xsec  =	12.14

    elif 'ggZH_HToBB_ZToNuNu_M125'  in samplename: Xsec = 0.01222*0.588
    elif 'ggZH_HToBB_ZToLL_M125'    in samplename: Xsec = 0.006185*0.588
    elif 'WminusH_HToBB_WToQQ_M125' in samplename: Xsec = 0.3654*0.588
    elif 'WplusH_HToBB_WToLNu_M125' in samplename: Xsec = 0.2819*0.588
    elif 'ZH_HToBB_ZToLL_M125'      in samplename and 'ggZH_HToBB_ZToLL_M125' not in samplename:   Xsec = 0.07924*0.588
    elif 'ZH_HToBB_ZToNuNu_M125'    in samplename and 'ggZH_HToBB_ZToNuNu_M125' not in samplename: Xsec = 0.1565*0.588
    return Xsec
