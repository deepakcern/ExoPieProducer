from pandas import  DataFrame

## SR dataframes

df_out_SR_1b = DataFrame(columns=['run','lumi','event','MET',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'weight'])

df_out_SR_2b = DataFrame(columns=['run','lumi','event','MET',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'weight'])


## define more data frames for each region
