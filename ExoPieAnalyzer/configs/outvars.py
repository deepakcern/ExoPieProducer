from pandas import  DataFrame

## first two dataframes are just dummy, you can delete from them to get dataframe for a given region,
df_out_boosted = DataFrame(columns=['run', 'lumi', 'event', 'MET',
                                    'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho',
                                    'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV', 'Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                    'FJetMass',
                                    'Mu1Pt', 'Mu1Eta', 'Mu1Phi', 'Mu2Pt', 'Mu2Eta', 'Mu2Phi',
                                    'Ele1Pt', 'Ele1Eta', 'Ele1Phi', 'Ele2Pt', 'Ele2Eta', 'Ele2Phi',
                                    'MT',
                                    'weight'])



df_out_resolved = DataFrame(columns=['run', 'lumi', 'event', 'MET',
                                     'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho',
                                     'Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi', 'Jet3Pt','Jet3Eta','Jet3Phi', 'Jet1CSV', 'Jet2CSV','Jet3CSV',
                                     'DiJetMass', 'DiJetPt', 'DiJetEta',
                                     'Mu1Pt', 'Mu1Eta', 'Mu1Phi', 'Mu2Pt', 'Mu2Eta', 'Mu2Phi',
                                     'Ele1Pt', 'Ele1Eta', 'Ele1Phi', 'Ele2Pt', 'Ele2Eta', 'Ele2Phi',
                                     'MT',
                                     'weight'])


'''
BOOSTED CATEGORY  dataframe
'''

df_out_SR_boosted = DataFrame(columns=['run', 'lumi', 'event', 'pu_nTrueInt','THINjetNPV','MET','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV', 'N2DDT','fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi','isAK4jet1EtaMatch','isAK4jet2EtaMatch','M_Jet1AK8Jet','M_Jet2AK8Jet','Jet1CSV','Jet2CSV',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'METweight','METweight_up','METweight_down','METRes_up','METRes_down','METEn_up','METEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_SBand_boosted = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt', 'THINjetNPV','MET','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV', 'N2DDT','fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi','isAK4jet1EtaMatch','isAK4jet2EtaMatch','M_Jet1AK8Jet','M_Jet2AK8Jet','Jet1CSV','Jet2CSV',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'METweight','METweight_up','METweight_down','METRes_up','METRes_down','METEn_up','METEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])



df_out_Tope_boosted = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt', 'THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV', 'N2DDT','fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi','isAK4jet1EtaMatch','isAK4jet2EtaMatch','M_Jet1AK8Jet','M_Jet2AK8Jet','Jet1CSV','Jet2CSV',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                       'lep1_pT','lep1_eta','lep1_Phi','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets','Wmass',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])


df_out_Topmu_boosted = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt', 'THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV', 'N2DDT','fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi','isAK4jet1EtaMatch','isAK4jet2EtaMatch','M_Jet1AK8Jet','M_Jet2AK8Jet','Jet1CSV','Jet2CSV',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                       'lep1_pT','lep1_eta','lep1_Phi','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets','Wmass',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])


df_out_We_boosted = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt', 'THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV', 'N2DDT','fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi','isAK4jet1EtaMatch','isAK4jet2EtaMatch','M_Jet1AK8Jet','M_Jet2AK8Jet','Jet1CSV','Jet2CSV',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                       'lep1_pT','lep1_eta','lep1_Phi','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets','Wmass',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_Wmu_boosted = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt', 'THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV', 'N2DDT','fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi','isAK4jet1EtaMatch','isAK4jet2EtaMatch','M_Jet1AK8Jet','M_Jet2AK8Jet','Jet1CSV','Jet2CSV',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                       'lep1_pT','lep1_eta','lep1_Phi','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets','Wmass',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])


df_out_Zee_boosted = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt', 'THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV','N2DDT', 'fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                       'lep1_pT','lep1_eta','lep1_Phi',
                                       'lep2_pT','lep2_eta','lep2_Phi',
                                       'Zmass','ZpT','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_Zmumu_boosted = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt', 'THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV','N2DDT', 'fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                       'lep1_pT','lep1_eta','lep1_Phi',
                                       'lep2_pT','lep2_eta','lep2_Phi',
                                       'Zmass','ZpT','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_TopWmu_boosted = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt', 'THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV', 'N2DDT','fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi','isAK4jet1EtaMatch','isAK4jet2EtaMatch','M_Jet1AK8Jet','M_Jet2AK8Jet','Jet1CSV','Jet2CSV',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                       'lep1_pT','lep1_eta','lep1_Phi','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets','Wmass',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_TopWe_boosted = DataFrame(columns=['run', 'lumi', 'event', 'pu_nTrueInt','THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV','N2DDT', 'fjetTau21','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi','isAK4jet1EtaMatch','isAK4jet2EtaMatch','M_Jet1AK8Jet','M_Jet2AK8Jet','Jet1CSV','Jet2CSV',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                       'lep1_pT','lep1_eta','lep1_Phi','FJetN2b1','FJetN2b2','FJetrho','min_dphi_jets','Wmass',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                       'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])
'''
RESOLVED CATEGORY  dataframe
'''


df_out_Tope_resolved  = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt','THINjetNPV','MET','RECOIL' ,'trkMET','trkMETPhi','METSig',
                                         'Njets_PassID','Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                         'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                         'DiJetMass', 'DiJetPt', 'DiJetEta','DiJetPhi','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                         'lep1_pT','lep1_eta','lep1_Phi','Wmass',
                                         'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                         'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                         'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                         'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                         'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_Topmu_resolved  = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt','THINjetNPV','MET','RECOIL' ,'trkMET','trkMETPhi','METSig',
                                         'Njets_PassID','Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                         'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                         'DiJetMass', 'DiJetPt', 'DiJetEta','DiJetPhi','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                         'lep1_pT','lep1_eta','lep1_Phi','Wmass',
                                         'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                         'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                         'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
                                         'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
                                         'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_Wmu_resolved  = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt','THINjetNPV','MET','RECOIL' ,'trkMET','trkMETPhi','METSig',
                                         'Njets_PassID','Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                         'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                         'DiJetMass', 'DiJetPt', 'DiJetEta','DiJetPhi','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                         'lep1_pT','lep1_eta','lep1_Phi','Wmass',
                                         'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                         'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                         'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
					 'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
					 'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_We_resolved  = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt','THINjetNPV','MET','RECOIL' ,'trkMET','trkMETPhi','METSig',
                                         'Njets_PassID','Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                         'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                         'DiJetMass', 'DiJetPt', 'DiJetEta','DiJetPhi','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                         'lep1_pT','lep1_eta','lep1_Phi','Wmass',
                                         'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                         'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                         'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
					 'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
					 'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_Zee_resolved    = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt','THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
					'Njets_PassID','Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                        'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                        'DiJetMass', 'DiJetPt', 'DiJetEta','DiJetPhi','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                        'lep1_pT','lep1_eta','lep1_Phi',
                                        'lep2_pT','lep2_eta','lep2_Phi',
                                        'Zmass','ZpT',
                                        'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                        'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                        'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
					'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
					'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_Zmumu_resolved    = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt','THINjetNPV','MET','RECOIL','trkMET','trkMETPhi','METSig',
                                        'Njets_PassID','Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                        'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV','Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                        'DiJetMass', 'DiJetPt', 'DiJetEta','DiJetPhi','nJets','min_dPhi','met_Phi','RECOIL_Phi',
                                        'lep1_pT','lep1_eta','lep1_Phi',
                                        'lep2_pT','lep2_eta','lep2_Phi',
                                        'Zmass','ZpT',
                                        'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                        'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                        'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
					'recoilweight','recoilweight_up','recoilweight_down','recoilRes_up','recoilRes_down','recoilEn_up','recoilEn_down',
					'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])


df_out_SR_resolved = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt','THINjetNPV','MET','trkMET','trkMETPhi','METSig',
                                       'Njets_PassID','Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                       'Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet1CSV','Jet2Pt', 'Jet2Eta', 'Jet2Phi', 'Jet2CSV',
                                       'Jet3Pt', 'Jet3Eta', 'Jet3Phi', 'Jet3CSV',
                                       'DiJetMass','DiJetPt', 'DiJetEta','DiJetPhi','nJets','met_Phi',
                                       'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                       'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                       'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
				       'METweight','METweight_up','METweight_down','METRes_up','METRes_down','METEn_up','METEn_down',
				       'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])

df_out_SBand_resolved = DataFrame(columns=['run', 'lumi', 'event','pu_nTrueInt', 'THINjetNPV','MET','trkMET','trkMETPhi','METSig',
                                        'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho','st_TopMatching',
                                        'Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi', 'Jet3Pt','Jet3Eta','Jet3Phi','Jet1CSV', 'Jet2CSV','Jet3CSV',
                                        'DiJetMass','DiJetPt', 'DiJetEta','DiJetPhi','nJets','met_Phi',
                                        'isak4JetBasedHemEvent', 'isak8JetBasedHemEvent',
                                        'ismetphiBasedHemEvent1', 'ismetphiBasedHemEvent2',
                                        'weight','puweight','puweight_up','puweight_down','lepweight','lepweight_up','lepweight_down',
					'METweight','METweight_up','METweight_down','METRes_up','METRes_down','METEn_up','METEn_down',
					'btagweight','btagweight_up','btagweight_down','ewkweight','ewkweight_up','ewkweight_down','toppTweight','toppTweight_up','toppTweight_down','jec','jec_up','jec_down','prefiringweight','prefiringweight_up','prefiringweight_down'])
