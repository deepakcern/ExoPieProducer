import os, sys
#//--------------------------------------------------------------------------------------
def getXsec(samplename):
    if 'DYJetsToLL_M-50_HT-100to200'   in samplename: Xsec  = 161.000000
    if 'DYJetsToLL_M-50_HT-1200to2500'   in samplename: Xsec  = 0.192700
    if 'DYJetsToLL_M-50_HT-200to400'   in samplename: Xsec  = 48.580000
    if 'DYJetsToLL_M-50_HT-2500toInf'   in samplename: Xsec  = 0.003478
    if 'DYJetsToLL_M-50_HT-400to600'   in samplename: Xsec  = 6.982000
    if 'DYJetsToLL_M-50_HT-600to800'   in samplename: Xsec  = 1.747000
    if 'DYJetsToLL_M-50_HT-800to1200'   in samplename: Xsec  = 0.805200

    if 'GJets_HT-40To100'   in samplename: Xsec  = 20790.000000
    if 'GJets_HT-100To200'   in samplename: Xsec  = 9238.000000
    if 'GJets_HT-200To400'   in samplename: Xsec  = 2193.000000
    if 'GJets_HT-400To600'   in samplename: Xsec  = 274.400000
    if 'GJets_HT-600ToInf'   in samplename: Xsec  = 93.460000

    if 'QCD_HT200to300'   in samplename: Xsec  = 1556000.000000
    if 'QCD_HT300to500'   in samplename: Xsec  = 323600.000000
    if 'QCD_HT500to700'   in samplename: Xsec  = 32100.000000
    if 'QCD_HT700to1000'   in samplename: Xsec  = 6831.000000
    if 'QCD_HT1000to1500'   in samplename: Xsec  = 1207.000000
    if 'QCD_HT1500to2000'   in samplename: Xsec  = 119.900000
    if 'QCD_HT2000toInf'   in samplename: Xsec  = 25.240000

    if 'ST_s-channel_4f_leptonDecays'   in samplename: Xsec  = 3.740000
    if 'ST_t-channel_antitop_4f_inclusiveDecays'   in samplename: Xsec  = 67.910000
    if 'ST_t-channel_top_4f_inclusiveDecays'   in samplename: Xsec  = 113.300000
    if 'ST_tW_antitop_5f_inclusiveDecays'   in samplename: Xsec  = 34.970000
    if 'ST_tW_top_5f_inclusiveDecays'   in samplename: Xsec  = 34.910000

    if 'TTTo2L2Nu_TuneCP5_PSweights'   in samplename: Xsec  = 687.100000*0.105
    if 'TTToHadronic_TuneCP5_PSweights'   in samplename: Xsec  = 687.100000*0.457
    if 'TTToSemiLeptonic_TuneCP5'   in samplename: Xsec  = 687.100000*0.438

    if 'WJetsToLNu_HT-100To200'   in samplename: Xsec  = 1395.000000
    if 'WJetsToLNu_HT-200To400'   in samplename: Xsec  = 407.900000
    if 'WJetsToLNu_HT-400To600'   in samplename: Xsec  = 57.480000
    if 'WJetsToLNu_HT-600To800'   in samplename: Xsec  = 12.870000
    if 'WJetsToLNu_HT-800To1200'   in samplename: Xsec  = 5.366000
    if 'WJetsToLNu_HT-1200To2500'   in samplename: Xsec  = 1.074000
    if 'WJetsToLNu_HT-2500ToInf'   in samplename: Xsec  = 0.008001

    if 'WW_TuneCP5_13TeV-'   in samplename: Xsec  = 75.900000
    if 'WZ_TuneCP5_13TeV-'   in samplename: Xsec  = 27.570000

    if 'ZJetsToNuNu_HT-100To200'   in samplename: Xsec  = 280.35
    if 'ZJetsToNuNu_HT-200To400'   in samplename: Xsec  = 77.67
    if 'ZJetsToNuNu_HT-400To600'   in samplename: Xsec  = 10.73
    if 'ZJetsToNuNu_HT-600To800'   in samplename: Xsec  = 2.559
    if 'ZJetsToNuNu_HT-800To1200'   in samplename: Xsec  = 1.1796
    if 'ZJetsToNuNu_HT-1200To2500'   in samplename: Xsec  = 0.28833
    if 'ZJetsToNuNu_HT-2500ToInf'   in samplename: Xsec  = 0.006945

    if 'ZZ_TuneCP5'   in samplename: Xsec  = 12.140000

    if 'ggZH_HToBB_ZToNuNu_M125'  in samplename: Xsec = 0.01222
    if 'ggZH_HToBB_ZToLL_M125'    in samplename: Xsec = 0.006185
    if 'WminusH_HToBB_WToQQ_M125' in samplename: Xsec = 0.3654
    if 'WplusH_HToBB_WToLNu_M125' in samplename: Xsec = 0.2819
    if 'ZH_HToBB_ZToLL_M125'      in samplename and 'ggZH_HToBB_ZToLL_M125' not in samplename:   Xsec = 0.07924
    if 'ZH_HToBB_ZToNuNu_M125'    in samplename and 'ggZH_HToBB_ZToNuNu_M125' not in samplename: Xsec = 0.1565
    return Xsec
