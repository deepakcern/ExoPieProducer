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


## SR dataframes 

df_out_SR_boosted = DataFrame(columns=['run', 'lumi', 'event', 'MET', 
                                       'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho',
                                       'FJetPt', 'FJetEta', 'FJetPhi', 'FJetCSV', 'Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi',
                                       'FJetMass', 'DiJetPt', 'DiJetEta','nJets',
                                       'weight'])

df_out_SR_resolved = DataFrame(columns=['run', 'lumi', 'event', 'MET', 
                                        'Njets_PassID', 'Nbjets_PassID', 'NTauJets', 'NEle', 'NMu', 'nPho',
                                        'Jet1Pt', 'Jet1Eta', 'Jet1Phi', 'Jet2Pt','Jet2Eta', 'Jet2Phi', 'Jet3Pt','Jet3Eta','Jet3Phi','Jet1CSV', 'Jet2CSV','Jet3CSV',
                                        'DiJetMass','nJets',
                                        'weight'])


## define more data frames for each region
