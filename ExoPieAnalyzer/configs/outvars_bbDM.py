from pandas import  DataFrame

## SR dataframes

df_out_SR_1b = DataFrame(columns=['run','lumi','event','MET','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'weight'])
df_out_SR_2b = DataFrame(columns=['run','lumi','event','MET','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'weight'])


df_out_ZeeCR_1b = DataFrame(columns=['run','lumi','event','MET','Recoil','Zmass','ZpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'subleadingLepPt','subleadingLepEta','subleadingLepPhi',
                                  'weight'])
df_out_ZeeCR_2b = DataFrame(columns=['run','lumi','event','MET','Recoil','Zmass','ZpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'subleadingLepPt','subleadingLepEta','subleadingLepPhi',
                                  'weight'])
df_out_ZmumuCR_1b = DataFrame(columns=['run','lumi','event','MET','Recoil','Zmass','ZpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'subleadingLepPt','subleadingLepEta','subleadingLepPhi',
                                  'weight'])
df_out_ZmumuCR_2b = DataFrame(columns=['run','lumi','event','MET','Recoil','Zmass','ZpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'subleadingLepPt','subleadingLepEta','subleadingLepPhi',
                                  'weight'])


df_out_WenuCR_1b = DataFrame(columns=['run','lumi','event','MET','Recoil','Wmass','WpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight'])
df_out_WenuCR_2b = DataFrame(columns=['run','lumi','event','MET','Recoil','Wmass','WpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight'])
df_out_WmunuCR_1b = DataFrame(columns=['run','lumi','event','MET','Recoil','Wmass','WpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight'])
df_out_WmunuCR_2b = DataFrame(columns=['run','lumi','event','MET','Recoil','Wmass','WpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight'])


df_out_TopenuCR_1b = DataFrame(columns=['run','lumi','event','MET','Recoil','Wmass','WpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight'])
df_out_TopenuCR_2b = DataFrame(columns=['run','lumi','event','MET','Recoil','Wmass','WpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight'])
df_out_TopmunuCR_1b = DataFrame(columns=['run','lumi','event','MET','Recoil','Wmass','WpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight'])
df_out_TopmunuCR_2b = DataFrame(columns=['run','lumi','event','MET','Recoil','Wmass','WpT','pu_nPUVert',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight'])

## define more data frames for each region
