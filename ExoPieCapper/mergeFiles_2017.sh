path='PUperSample_tifrAll/'
outputPath='PUperSample_tifrAll_merged/'
#'2017_TopnJetsLess3_B_merged/'

hadd "$outputPath"combined_data_MET.root  "$path"*MET-Run2017*
hadd "$outputPath"combined_data_SE.root  "$path"*SingleElectron*

hadd "$outputPath"crab_ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8.root "$path"*ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8*
hadd "$outputPath"crab_ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8.root "$path"*ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8*
hadd "$outputPath"crab_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root "$path"*TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8*
hadd "$outputPath"crab_WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8.root "$path"*WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8*
hadd "$outputPath"crab_WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8.root "$path"*WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8*
hadd "$outputPath"crab_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8.root "$path"*ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8*
hadd "$outputPath"crab_ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8.root "$path"*ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8*
hadd "$outputPath"DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8.root  "$path"*GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.root "$path"*QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8*
hadd "$outputPath"QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8.root "$path"*QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8*
hadd "$outputPath"QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root "$path"*QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8*
hadd "$outputPath"QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8.root "$path"*QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8*
hadd "$outputPath"QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8.root "$path"*QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8*
hadd "$outputPath"QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root "$path"*QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8*
hadd "$outputPath"QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8.root "$path"*QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8*
hadd "$outputPath"QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8.root "$path"*QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8*
hadd "$outputPath"ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8.root "$path"*ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8*
hadd "$outputPath"ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root "$path"*ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8*
hadd "$outputPath"ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root "$path"*ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8*
hadd "$outputPath"ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root "$path"*ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8*
hadd "$outputPath"ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root "$path"*ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8*
hadd "$outputPath"TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root "$path"*TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8*
hadd "$outputPath"TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root "$path"*TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8*
hadd "$outputPath"WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8.root "$path"*WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8*
hadd "$outputPath"ZJetsToNuNu_HT-100To200_13TeV-madgraph.root "$path"*ZJetsToNuNu_HT-100To200_13TeV-madgraph*
hadd "$outputPath"ZJetsToNuNu_HT-1200To2500_13TeV-madgraph.root "$path"*ZJetsToNuNu_HT-1200To2500_13TeV-madgraph*
hadd "$outputPath"ZJetsToNuNu_HT-200To400_13TeV-madgraph.root "$path"*ZJetsToNuNu_HT-200To400_13TeV-madgraph*
hadd "$outputPath"ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph.root "$path"*ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph*
hadd "$outputPath"ZJetsToNuNu_HT-400To600_13TeV-madgraph.root "$path"*ZJetsToNuNu_HT-400To600_13TeV-madgraph*
hadd "$outputPath"ZJetsToNuNu_HT-600To800_13TeV-madgraph.root "$path"*ZJetsToNuNu_HT-600To800_13TeV-madgraph*
hadd "$outputPath"ZJetsToNuNu_HT-800To1200_13TeV-madgraph.root "$path"*ZJetsToNuNu_HT-800To1200_13TeV-madgraph*
hadd "$outputPath"ZZ_TuneCP5_13TeV-pythia8.root "$path"*ZZ_TuneCP5_13TeV-pythia8*
hadd "$outputPath"WW_TuneCP5_13TeV-pythia8.root "$path"*WW_TuneCP5_13TeV-pythia8*
hadd "$outputPath"WZ_TuneCP5_13TeV-pythia8.root "$path"*WZ_TuneCP5_13TeV-pythia8*
hadd "$outputPath"ttHTobb_ttToSemiLep_M125_TuneCP5_13TeV-powheg-pythia8.root "$path"*ttHTobb_ttToSemiLep_M125_TuneCP5_13TeV-powheg-pythia8*
