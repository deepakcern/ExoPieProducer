void copy(){
  
  TFile* f = new TFile("AllMETHistos_180416.root","update");
    TH1F* h_new;
  
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_ZJets_2b");      h_new->SetName("ZJets");    h_new->Write();
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_WJets_2b");	  h_new->SetName("WJets");    h_new->Write();
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_GJets_2b");	  h_new->SetName("GJets");    h_new->Write();
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_TT_2b");	  h_new->SetName("Top");    h_new->Write();
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_STop_2b");	  h_new->SetName("STop");    h_new->Write();
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_GJets_2b");	  h_new->SetName("GJets");    h_new->Write();
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_QCD_2b");	  h_new->SetName("QCD");    h_new->Write();
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_DYJets_2b");	  h_new->SetName("DYJets");    h_new->Write();
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_DIBOSON_2b");    h_new->SetName("DIBOSON");    h_new->Write();
  h_new = (TH1F*) f->Get("bbDM_2016_SignalRegion_ZJets_2b");      h_new->SetName("data_obs");    h_new->Write();
}
