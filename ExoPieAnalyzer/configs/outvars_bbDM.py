from pandas import  DataFrame

## SR dataframes

df_out_SR_1b = DataFrame(columns=['run','lumi','event','MET','METPhi','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'weight','weightMET','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightMET_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','MET_Res_up','MET_En_up','MET_En_down','MET_Res_down','weightJEC_down','weightMET_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_SR_2b = DataFrame(columns=['run','lumi','event','MET','METPhi','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'weight','weightMET','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightMET_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','MET_Res_up','MET_En_up','MET_En_down','MET_Res_down','weightJEC_down','weightMET_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])


df_out_ZeeCR_1b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Zmass','ZpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'subleadingLepPt','subleadingLepEta','subleadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_ZeeCR_2b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Zmass','ZpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'subleadingLepPt','subleadingLepEta','subleadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_ZmumuCR_1b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Zmass','ZpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'subleadingLepPt','subleadingLepEta','subleadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_ZmumuCR_2b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Zmass','ZpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'subleadingLepPt','subleadingLepEta','subleadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])


df_out_WenuCR_1b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Wmass','WpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_WenuCR_2b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Wmass','WpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_WmunuCR_1b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Wmass','WpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_WmunuCR_2b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Wmass','WpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])


df_out_TopenuCR_1b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Wmass','WpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_TopenuCR_2b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Wmass','WpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_TopmunuCR_1b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Wmass','WpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
df_out_TopmunuCR_2b = DataFrame(columns=['run','lumi','event','MET','METPhi','Recoil','RecoilPhi','Wmass','WpT','nPV',
                                  'dPhi_jetMET','NTau','NEle','NMu','nPho','Njets_PassID','Nbjets_PassID',
                                  'Jet1Pt','Jet1Eta','Jet1Phi','Jet1deepCSV','Jet2Pt','Jet2Eta','Jet2Phi','Jet2deepCSV',
                                  'Jet3Pt','Jet3Eta','Jet3Phi','Jet3deepCSV',
                                  'leadingLepPt','leadingLepEta','leadingLepPhi',
                                  'weight','weightRecoil','weightEle','weightMu','weightB','weightEWK','weightQCD','weightTop','weightPU','weightRecoil_up','weightEle_up','weightMu_up','weightB_up','weightEWK_up','weightQCD_up','weightTop_up','weightPU_up','weightPrefire_up','weightJEC_up','Recoil_Res_up','Recoil_En_up','Recoil_En_down','Recoil_Res_down','weightJEC_down','weightRecoil_down','weightEle_down','weightMu_down','weightB_down','weightEWK_down','weightQCD_down','weightTop_down','weightPU_down','weightPrefire_down'])
## define more data frames for each region
