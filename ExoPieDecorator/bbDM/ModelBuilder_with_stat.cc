

#include "ModelBuilder.h"
#include "RooCategory.h"

void ModelBuilder::add_cut(std::string region,std::string ecut){
  std::cout << " Adding cut - " << region << ", " << ecut << std::endl;
  extracuts[region]+=ecut;
}

void ModelBuilder::saveReweightedHist(){
  TFile* freweighted = new TFile("datacards/original.root");
  
  
  TH1F* h_rew = (TH1F*) freweighted->Get("monohiggs_signal_BarZp-500-1_signal");
  h_rew->SetName("reweighted_monohiggs_signal_BarZp-500-1_signal");
  fOut->WriteTObject(h_rew);
  
  RooDataHist tmp_hist(h_rew->GetName(),h_rew->GetTitle(), RooArgList(*(wspace->var(varstring.c_str()))),h_rew);     
  wspace->import(tmp_hist);                                                                                                                                                       
  RooHistPdf tmp_pdf(h_rew->GetName(),h_rew->GetTitle(),RooArgSet(*(wspace->var(varstring.c_str()))),*((RooDataHist*)(wspace->data(tmp_hist.GetName()))));
  wspace->import(tmp_pdf);
    
  
  

  h_rew = (TH1F*) freweighted->Get("monohiggs_signal_BarZp-500-1_signal_btagUp");
  h_rew->SetName("reweighted_monohiggs_signal_BarZp-500-1_signal_btagUp");
  fOut->WriteTObject(h_rew);
  
  RooDataHist tmp_hist1(h_rew->GetName(),h_rew->GetTitle(), RooArgList(*(wspace->var(varstring.c_str()))),h_rew);     
  wspace->import(tmp_hist1);                                                                                                                                                       
  RooHistPdf tmp_pdf1(h_rew->GetName(),h_rew->GetTitle(),RooArgSet(*(wspace->var(varstring.c_str()))),*((RooDataHist*)(wspace->data(tmp_hist1.GetName()))));
  wspace->import(tmp_pdf1);
  

  
  h_rew = (TH1F*) freweighted->Get("monohiggs_signal_BarZp-500-1_signal_btagDown");
  h_rew->SetName("reweighted_monohiggs_signal_BarZp-500-1_signal_btagDown");
  fOut->WriteTObject(h_rew);
  
  RooDataHist tmp_hist2(h_rew->GetName(),h_rew->GetTitle(), RooArgList(*(wspace->var(varstring.c_str()))),h_rew);     
  wspace->import(tmp_hist2);                                                                                                                                                       
  RooHistPdf tmp_pdf2(h_rew->GetName(),h_rew->GetTitle(),RooArgSet(*(wspace->var(varstring.c_str()))),*((RooDataHist*)(wspace->data(tmp_hist2.GetName()))));
  wspace->import(tmp_pdf2);
  
  

  
  h_rew = (TH1F*) freweighted->Get("monohiggs_signal_BarZp-500-1_signal_mistagUp");
  h_rew->SetName("reweighted_monohiggs_signal_BarZp-500-1_signal_mistagUp");
  fOut->WriteTObject(h_rew);
  
  RooDataHist tmp_hist3(h_rew->GetName(),h_rew->GetTitle(), RooArgList(*(wspace->var(varstring.c_str()))),h_rew);     
  wspace->import(tmp_hist3);                                                                                                                                                       
  RooHistPdf tmp_pdf3(h_rew->GetName(),h_rew->GetTitle(),RooArgSet(*(wspace->var(varstring.c_str()))),*((RooDataHist*)(wspace->data(tmp_hist3.GetName()))));
  wspace->import(tmp_pdf3);
  


  h_rew = (TH1F*) freweighted->Get("monohiggs_signal_BarZp-500-1_signal_mistagDown");
  h_rew->SetName("reweighted_monohiggs_signal_BarZp-500-1_signal_mistagDown");
  fOut->WriteTObject(h_rew);
  
  RooDataHist tmp_hist4 (h_rew->GetName(),h_rew->GetTitle(), RooArgList(*(wspace->var(varstring.c_str()))),h_rew);     
  wspace->import(tmp_hist4);                                                                                                                                                       
  RooHistPdf tmp_pdf4 (h_rew->GetName(),h_rew->GetTitle(),RooArgSet(*(wspace->var(varstring.c_str()))),*((RooDataHist*)(wspace->data(tmp_hist4.GetName()))));
  wspace->import(tmp_pdf4);
  
  




}

void ModelBuilder::saveHist(TH1F *histogram){

  fOut->WriteTObject(histogram);
  /*

  for (int i=1; i<=histogram->GetNbinsX();i++){
    string histname = histogram->GetName();
    if (histname.find("data") != std::string::npos)
      continue;
    if (histname.find("Up") != std::string::npos)
      continue;
    if (histname.find("Down") != std::string::npos)
      continue;
    if (!(histname.find("dielectron_fail") != std::string::npos) and !(histname.find("dimuon_fail") != std::string::npos) and !(histname.find("singleelectronw_fail") != std::string::npos) and !(histname.find("singlemuonw_fail") != std::string::npos))
      continue;

    string bin_string = histogram->GetName();

    string chan;
    string proc;
    
    if (bin_string.find("dielectron") != std::string::npos)
      chan = "die";
    if (bin_string.find("dielectron_fail") != std::string::npos)
      chan = "die_f";
    if (bin_string.find("dimuon") != std::string::npos)
      chan = "dim";
    if (bin_string.find("dimuon_fail") != std::string::npos)
      chan = "dim_f";
    if (bin_string.find("singlemuontop") != std::string::npos)
      chan = "tmn";
    if (bin_string.find("singleelectrontop") != std::string::npos)
      chan = "ten";
    if (bin_string.find("singlemuonw_") != std::string::npos)
      chan = "wmn";
    if (bin_string.find("singlemuonw_fail") != std::string::npos)
      chan = "wmn_f";
    if (bin_string.find("singleelectronw_") != std::string::npos)
      chan = "wen";
    if (bin_string.find("singleelectronw_fail") != std::string::npos)
      chan = "wen_f";

    if (bin_string.find("zll") != std::string::npos)
      proc = "zll";
    if (bin_string.find("qcd") != std::string::npos)
      proc = "qcd";
    if (bin_string.find("stop") != std::string::npos)
      proc = "stop";
    if (bin_string.find("dibosons") != std::string::npos)
      proc = "dibosons";
    if (bin_string.find("wjets") != std::string::npos)
      proc = "wjets";
    if (bin_string.find("vh") != std::string::npos)
      proc = "vh";
    if (bin_string.find("zvv") != std::string::npos)
      proc = "zvv";
    if (bin_string.find("ttbar") != std::string::npos)
      proc = "ttbar";


    bin_string = bin_string.append("_");
    bin_string = bin_string.append(chan);
    bin_string = bin_string.append("_");
    bin_string = bin_string.append(proc);
    bin_string = bin_string.append("_bin");
    string bin_id = std::to_string(i);
    bin_string = bin_string.append(bin_id);
    string bin_string_2 = bin_string;
    
    TH1F* err_up = (TH1F*)histogram->Clone((bin_string.append("Up")).c_str());
    TH1F* err_down = (TH1F*)histogram->Clone((bin_string_2.append("Down")).c_str());
    
    err_up->SetBinContent(i,histogram->GetBinContent(i)+histogram->GetBinError(i));
    if (histogram->GetBinContent(i)<=histogram->GetBinError(i))
      err_down->SetBinContent(i,0.001);
    else
      err_down->SetBinContent(i,histogram->GetBinContent(i)-histogram->GetBinError(i));

    fOut->WriteTObject(err_up);
    fOut->WriteTObject(err_down);
    
  }
  */


  
  string histname = histogram->GetName();
  std::cout<<" Raman :: integral of "<<histname<<" is "<<histogram->Integral()<<std::endl;
  //string scale_string = "Scale";
  if (histname.find("Scale") != std::string::npos){
    std::size_t found = histname.find_last_of("_");
    string before = histname.substr(0,found);
    string after = histname.substr(found+1);

    string proc;
    if (histname.find("ttbar") != std::string::npos){
      proc = "_ttbar";
      string final_name = before;
      final_name = final_name.append(proc);
      final_name = final_name.append(after); 
      TH1F* clone = (TH1F*)histogram->Clone(final_name.c_str());
      fOut->WriteTObject(clone);
    }
    if (histname.find("wjets") != std::string::npos){
      proc = "_wjets";
      string final_name = before;
      final_name = final_name.append(proc);
      final_name = final_name.append(after); 
      TH1F* clone = (TH1F*)histogram->Clone(final_name.c_str());
      fOut->WriteTObject(clone);
    }
    if (histname.find("zll") != std::string::npos || histname.find("zvv") != std::string::npos){
      proc = "_zjets";
      string final_name = before;
      final_name = final_name.append(proc);
      final_name = final_name.append(after); 
      TH1F* clone = (TH1F*)histogram->Clone(final_name.c_str());
      fOut->WriteTObject(clone);
    }
    if (histname.find("vh") != std::string::npos || histname.find("zvv") != std::string::npos){
      proc = "_vh";
      string final_name = before;
      final_name = final_name.append(proc);
      final_name = final_name.append(after); 
      TH1F* clone = (TH1F*)histogram->Clone(final_name.c_str());
      fOut->WriteTObject(clone);
    }

  }


  histname = histogram->GetName();
  
  //string scale_string = "Scale";
  if (histname.find("PDF") != std::string::npos){
    std::size_t found = histname.find_last_of("_");
    string before = histname.substr(0,found);
    string after = histname.substr(found+1);

    string proc;
    if (histname.find("ttbar") != std::string::npos){
      proc = "_ttbar";
      string final_name = before;
      final_name = final_name.append(proc);
      final_name = final_name.append(after); 
      TH1F* clone = (TH1F*)histogram->Clone(final_name.c_str());
      fOut->WriteTObject(clone);
    }
    if (histname.find("wjets") != std::string::npos){
      proc = "_wjets";
      string final_name = before;
      final_name = final_name.append(proc);
      final_name = final_name.append(after); 
      TH1F* clone = (TH1F*)histogram->Clone(final_name.c_str());
      fOut->WriteTObject(clone);
    }
    if (histname.find("zll") != std::string::npos || histname.find("zvv") != std::string::npos){
      proc = "_zjets";
      string final_name = before;
      final_name = final_name.append(proc);
      final_name = final_name.append(after); 
      TH1F* clone = (TH1F*)histogram->Clone(final_name.c_str());
      fOut->WriteTObject(clone);
    }
    if (histname.find("vh") != std::string::npos || histname.find("zvv") != std::string::npos){
      proc = "_vh";
      string final_name = before;
      final_name = final_name.append(proc);
      final_name = final_name.append(after); 
      TH1F* clone = (TH1F*)histogram->Clone(final_name.c_str());
      fOut->WriteTObject(clone);
    }

  }
  
   // also make a RooDataHist and RooHistPdf
//   RooDataHist tmp_hist(Form("dhist_%s",histogram->GetName()),histogram->GetTitle(),RooArgList(*(wspace->var(varstring.c_str()))),histogram);
//   wspace->import(tmp_hist);
  
//   RooHistPdf tmp_pdf(Form("hpdf_%s",histogram->GetName()),histogram->GetTitle(),RooArgSet(*(wspace->var(varstring.c_str()))),*((RooDataHist*)(wspace->data(tmp_hist.GetName()))));
//   wspace->import(tmp_pdf);
}
const char * ModelBuilder::turnon(RooWorkspace *ws,RooRealVar &x,std::string ext){
   /*
   RooRealVar m1(Form("plaw_p_%s",ext.c_str()),"p",-4.,-10.,0.);
   RooRealVar eff_mu(Form("eff_mu_%s",ext.c_str()),"eff decay param mu",250,200.,1000.);
   RooRealVar eff_sig(Form("eff_sig_%s",ext.c_str()),"eff decay param sigma ",1,0.,10.);
   //RooGenericPdf *pdf = new RooGenericPdf(Form("powerLaw_%s",ext.c_str()), Form("powerLaw_%s",ext.c_str()),"TMath::Power(@0,@1)",RooArgList(x,m1));
   // Same name as exponential for now
   RooGenericPdf *pdf_pow = new RooGenericPdf(Form("doubleExponential_powerlaw%s",ext.c_str()), Form("doubleExponential_powerlaw%s",ext.c_str()),"TMath::Power(@0,@1)",RooArgList(x,m1));
   RooGenericPdf *eff = new RooGenericPdf(Form("doubleExponential_efficiency%s",ext.c_str()),"Efficiency thing","0.5*(TMath::Erf((@0-@1)/@2)+1)",RooArgList(x,eff_mu,eff_sig));
   //RooEffProd *pdf = new RooEffProd(Form("doubleExponential_%s",ext.c_str()), Form("doubleExponential_%s",ext.c_str()),*pdf_pow,eff);
   RooProdPdf *pdf = new RooProdPdf(Form("doubleExponential_%s",ext.c_str()), Form("doubleExponential_%s",ext.c_str()),RooArgList(*pdf_pow,*eff));
   */
   ext += catname;
   RooRealVar eff_mu(Form("eff_mu_%s",ext.c_str()),"eff decay param mu",240,200.,300.);
   RooRealVar eff_sig(Form("eff_sig_%s",ext.c_str()),"eff decay param sigma ",10,0.,20.);
   RooRealVar frac(Form("frac_%s",ext.c_str()),"Fraction",0.8,0,1);

   RooRealVar m1(Form("plaw_p_%s",ext.c_str()),"p",-4.,-15.,0.);
   RooGenericPdf *pdf_pow = new RooGenericPdf(Form("doubleExponential_powerlaw%s",ext.c_str()), Form("doubleExponential_powerlaw%s",ext.c_str()),"TMath::Power(@0,@1)",RooArgList(x,m1));
   //RooProdPdf *pdf = new RooProdPdf(Form("doubleExponential_%s",ext.c_str()), Form("doubleExponential_%s",ext.c_str()),RooArgList(*pdf_pow,eff));
   RooLognormal *pdf_l = new RooLognormal(Form("doubleExponential_log%s",ext.c_str()), Form("doubleExponential_log%s",ext.c_str()),x,eff_mu,eff_sig);
   RooAddPdf *pdf = new RooAddPdf(Form("doubleExponential_%s",ext.c_str()), Form("doubleExponential_%s",ext.c_str()),RooArgList(*pdf_pow,*pdf_l),RooArgList(frac));
   ws->import(*pdf);
   return pdf->GetName(); 
}

const char * ModelBuilder::powerlaw(RooWorkspace *ws,RooRealVar &x,std::string ext){
   ext += catname;
   RooRealVar p(Form("p_%s",ext.c_str()),"p",-2.,-10.,0.01);
   //RooGenericPdf *pdf = new RooGenericPdf(Form("powerLaw_%s",ext.c_str()), Form("powerLaw_%s",ext.c_str()),"TMath::Power(@0,@1)",RooArgList(x,m1));
   ws->import(p);
   // Same name as exponential for now
   RooGenericPdf *pdf = new RooGenericPdf(Form("doubleExponential_%s",ext.c_str()), Form("doubleExponential_%s",ext.c_str()),"TMath::Power(@0,@1)",RooArgList(x,*(ws->var(p.GetName()))));
   ws->import(*pdf);
   return pdf->GetName(); 
}
/*
const char * ModelBuilder::doubleexp(RooWorkspace *ws,RooRealVar &x,std::string ext){ 
   ext += catname;
   // Double exponential model
  // RooRealVar frac(Form("f_%s",ext.c_str()),"f",0.8,0.,1.);
  // RooRealVar m1(Form("m1_%s",ext.c_str()),"m1",-0.02,-0.3,0.0);
  // RooRealVar m2(Form("m2_%s",ext.c_str()),"m2",-0.01,-0.3,0.0);
   RooRealVar a(Form("a_%s",ext.c_str()),"a",-0.2,-2,2.0);
   RooRealVar b(Form("b_%s",ext.c_str()),"b",0,-2,2.0);

   //RooExponential exp1(Form("exp1_%s",ext.c_str()),Form("exp1_%s",ext.c_str()),x,m1);
   //RooExponential exp2(Form("exp2_%s",ext.c_str()),Form("exp2_%s",ext.c_str()),x,m2);

   //RooAddPdf *sumexp = new RooAddPdf(Form("doubleExponential_%s",ext.c_str()), Form("doubleExponential_%s",ext.c_str()),RooArgList(exp1,exp2),RooArgList(frac));
   RooGenericPdf *sumexp = new RooGenericPdf(Form("doubleExponential_%s",ext.c_str()), Form("doubleExponential_%s",ext.c_str()),"TMath::Exp(@1*@0+@2*@0*@0)",RooArgList(x,a,b));
   ws->import(*sumexp);
   return sumexp->GetName(); 
}
*/

const char * ModelBuilder::doubleexp(RooWorkspace *ws,RooRealVar &x,std::string ext){ 
   ext += catname;
   // Double exponential model
   RooRealVar frac(Form("f_%s",ext.c_str()),"f",0.8,0.,1.);
   RooRealVar m1(Form("m1_%s",ext.c_str()),"m1",-0.002,-0.1,0.0);
   RooRealVar m2(Form("m2_%s",ext.c_str()),"m2",-0.001,-0.1,0.0);

   RooExponential exp1(Form("exp1_%s",ext.c_str()),Form("exp1_%s",ext.c_str()),x,m1);
   RooExponential exp2(Form("exp2_%s",ext.c_str()),Form("exp2_%s",ext.c_str()),x,m2);

   RooAddPdf *sumexp = new RooAddPdf(Form("doubleExponential_%s",ext.c_str()), Form("doubleExponential_%s",ext.c_str()),RooArgList(exp1,exp2),RooArgList(frac));
   ws->import(*sumexp);
   return sumexp->GetName(); 
}

const char * ModelBuilder::singleexp(RooWorkspace *ws,RooRealVar &x,std::string ext){ 
   ext += catname;
   RooRealVar m1(Form("m3_%s",ext.c_str()),"m3",-0.02,-0.3,0.05);
   RooExponential* exp1 =new RooExponential(Form("doubleExponential_%s",ext.c_str()),Form("exp3_%s",ext.c_str()),x,m1);
   ws->import(*exp1);
   return exp1->GetName(); 
}
void ModelBuilder::save(){
  std::map<std::string, TH1F*>::iterator it = save_hists.begin();
  for (;it!=save_hists.end();it++){
    //fOut->WriteTObject((*it).second);
    saveHist((*it).second);
  }

  std::map<std::string, RooDataSet*>::iterator itd = save_datas.begin();
  for (;itd!=save_datas.end();itd++){
    if (!wspace->data((*itd).second->GetName())) wspace->import(*((*itd).second));
    // should also make the basic dataset into a histogram for each additional var
    for (std::map<std::string,TH1F*>::iterator additional_var = additional_vars.begin()
        ; additional_var!=additional_vars.end()
	; additional_var++){
    	std::string var = additional_var->first;
    	TH1F *tmph = additional_var->second;
    	TH1F *hist_c = (TH1F*)generateTemplate(tmph,(RooFormulaVar*) 0
   	  , *(wspace->var(varstring.c_str()))
	  , (*itd).second
	  ,0,1, var, var); // extension is var name
     //fOut->WriteTObject(hist_corrected);
    	saveHist(hist_c);
    }
  }
  
  //saveReweightedHist();
  fOut->cd();
  wspace->Write();
}

bool ModelBuilder::has_process(ControlRegion &cr,std::string proc){
   bool ret = false;
   std::vector<std::pair<std::string, int> >::iterator it = cr.procs.begin();
   for (;it!=cr.procs.end();it++){
	if ((*it).first==proc) { 
		ret=true;
		break;
	}
   }
   return ret;
}
int ModelBuilder::process_type(ControlRegion &cr,std::string proc){
   std::vector<std::pair<std::string, int> >::iterator it = cr.procs.begin();
   for (;it!=cr.procs.end();it++){
	if ((*it).first==proc) { 
		return (*it).second;
	}
   }
   std::cout<< "Error - No Process Found -- " << proc  << std::endl;
   return 0;
}
void ModelBuilder::apply_corrections(std::string correction, std::string region, std::string process, bool run_systematics){
   // For now just run one, but of course can correct ALL processes if we need;
   // Find and import the dataset region_process
   std::map<std::string,RooDataSet*>::iterator it = save_datas.find(region+std::string("_")+process);

   wspace->import(*(*it).second);
   TH1F *hist_corrected = (TH1F*)generateTemplate(lTmp, (RooFormulaVar*)wspace->function(correction.c_str())
   	, *(wspace->var(varstring.c_str()))
	, (RooDataSet*) wspace->data(Form("%s_%s",region.c_str(),process.c_str()))
	,1,1);
     //fOut->WriteTObject(hist_corrected);
   saveHist(hist_corrected);


   for (std::map<std::string,TH1F*>::iterator additional_var = additional_vars.begin()
        ; additional_var!=additional_vars.end()
	; additional_var++){
     std::string var = additional_var->first;
     TH1F *tmph = additional_var->second;
     TH1F *hist_c = (TH1F*)generateTemplate(tmph,(RooFormulaVar*)wspace->function(correction.c_str())
   	, *(wspace->var(varstring.c_str()))
	, (RooDataSet*) wspace->data(Form("%s_%s",region.c_str(),process.c_str()))
	,1,1,var, var); // extension is same name as var
     //fOut->WriteTObject(hist_corrected);
     saveHist(hist_c);
   }

   //save_datas.erase(it);

   if (! run_systematics) return;
   // Systematics
   std::vector<TH1F> v_th1f_;
   std::cout << " Generating Systematics " << std::endl;
   generateVariations(lTmp,(RooFitResult*)fOut->Get(Form("fitResult_%s",correction.c_str()))
   	,(RooFormulaVar*)wspace->function(correction.c_str())
	,wspace->var(varstring.c_str()),v_th1f_,wspace,region+std::string("_")+process,fOut);

   // Book saving of histograms
   for (std::vector<TH1F>::iterator its = v_th1f_.begin();its!=v_th1f_.end();its++){
      // rename to include category name so that combine doesn't correlate them!
      (*its).SetName(Form("%s_%s",catname.c_str(),its->GetName()));
      saveHist(&(*its));
   }


}

void ModelBuilder::run_corrections(std::string correction_name,std::string region){


   // First find the region some variables the region
   ControlRegion cr = v_samples[region];
   // get the format from the first sample
   if (!cr.procs.size()) { 
   	std::cout << "Error, no processes for Region " << region << std::endl;
	assert(0);
   }

   fOut->cd();
   RooRealVar *x = wspace->var(varstring.c_str());
   RooRealVar *w = wspace->var(weightname.c_str());

   RooAbsPdf *pdf 	       ;
   RooAbsPdf *pdf_mc 	       ;
   RooAbsPdf *pdf_background_mc;

   // Build and fit the model for the background
   if (_pdfmodel==0){  // Double exponential
    std::cout << "Using Double Exponential for fitting functions -- " << _pdfmodel << std::endl;
     pdf 	       = wspace->pdf(doubleexp(wspace,*x,Form("%s_data",cr.name.c_str())));
     pdf_mc 	       = wspace->pdf(doubleexp(wspace,*x,Form("%s_mc",cr.name.c_str())));
     pdf_background_mc = wspace->pdf(doubleexp(wspace,*x,Form("%s_bkg_mc",cr.name.c_str())));
   } else if (_pdfmodel==1) {
     // Why not use one power law for the resolved (lower stat) case?
     std::cout << "Using PowerLaw for fitting functions -- " <<  _pdfmodel <<std::endl;
     pdf 		= wspace->pdf(powerlaw(wspace,*x,Form("%s_data",cr.name.c_str())));
     pdf_mc 	   	= wspace->pdf(powerlaw(wspace,*x,Form("%s_mc",cr.name.c_str())));
     pdf_background_mc = wspace->pdf(powerlaw(wspace,*x,Form("%s_bkg_mc",cr.name.c_str())));
   
   } else if (_pdfmodel==2) {
    std::cout << "Using PowerLaw + Log-normal for fitting functions -- " << _pdfmodel <<std::endl;
     pdf 		= wspace->pdf(turnon(wspace,*x,Form("%s_data",cr.name.c_str())));
     pdf_mc 	   	= wspace->pdf(turnon(wspace,*x,Form("%s_mc",cr.name.c_str())));
     pdf_background_mc  = wspace->pdf(turnon(wspace,*x,Form("%s_bkg_mc",cr.name.c_str())));
   } else if (_pdfmodel==3) {
     std::cout << "Using Single Exponential for fitting functions -- " << _pdfmodel << std::endl;
     pdf 	       = wspace->pdf(singleexp(wspace,*x,Form("%s_data"  ,cr.name.c_str())));
     pdf_mc 	       = wspace->pdf(singleexp(wspace,*x,Form("%s_mc"    ,cr.name.c_str())));
     pdf_background_mc = wspace->pdf(singleexp(wspace,*x,Form("%s_bkg_mc",cr.name.c_str())));
   } else {
	std::cout << "! No Pdf Model for type " << _pdfmodel << std::endl;
	assert(0);
   }

   // Build a dataset which is all data in control region
   RooArgSet *treevars = (RooArgSet*) wspace->genobj("treevars");
   RooDataSet all_data(Form("%s_all_data",region.c_str()),Form("%s_all_data",region.c_str()),RooArgSet(*treevars));
   RooDataSet all_sig_mc(Form("%s_all_signal",region.c_str()),Form("%s_all_signal",region.c_str()),RooArgSet(*treevars),w->GetName());
   RooDataSet all_bkg_mc(Form("%s_all_background",region.c_str()),Form("%s_all_background",region.c_str()),RooArgSet(*treevars),w->GetName());

   for (std::vector<std::pair<std::string,int> >::iterator it_p = cr.procs.begin();
   	it_p!=cr.procs.end();it_p++){
	int type = (*it_p).second;
	std::string prc = (*it_p).first;
	std::cout << "Adding DataSet -- " << cr.name+std::string("_")+prc << " For fitting of the region - " << region << " type " << type << ", sum-weights = " << (*save_datas[cr.name+std::string("_")+prc]).sumEntries() << std::endl; 
	if (type==0)	 all_data.append(*save_datas[cr.name+std::string("_")+prc]);
	else if (type>0) all_bkg_mc.append(*save_datas[cr.name+std::string("_")+prc]);
	else if (type<0) all_sig_mc.append(*save_datas[cr.name+std::string("_")+prc]);
   }
   wspace->import(all_data);
   wspace->import(all_sig_mc);
   wspace->import(all_bkg_mc);

   /*
   // Fit control region MC
   RooFitResult *fit_res_control_mc  = pdf_mc->fitTo(all_sig_mc,RooFit::Save(1),RooFit::SumW2Error(false));
   fit_res_control_mc->SetName(Form("fitResult_%s_all_signal",region.c_str())); fit_res_control_mc->Write();

   // Fit background MC and then fix it
   pdf_background_mc->fitTo(all_bkg_mc,RooFit::SumW2Error(true)); // we dont use the errors so who cares about sumW2?
   freezeParameters(pdf_background_mc->getParameters(RooArgSet(*x)));
   
   // Now fit the Zvv Data 
   double nbkgcont = all_bkg_mc.sumEntries();
   double ncont    = all_data.sumEntries()-nbkgcont;

   RooRealVar num_contamination(Form("num_contamination_%s",region.c_str()),Form("num_contamination_%s",region.c_str()),nbkgcont,0,10E10);
   num_contamination.setConstant();
   RooRealVar num(Form("num_%s",region.c_str()),Form("num_%s",region.c_str()),ncont,0,10E10);
   //num.setConstant(true);// Float in the fit -> additional parameter

   RooAddPdf model(Form("model_%s",region.c_str()),Form("model_%s",region.c_str()),RooArgList(*pdf_background_mc,*pdf),RooArgList(num_contamination,num));

   RooFitResult *fit_res_control = model.fitTo(all_data,RooFit::Save(1),RooFit::Strategy(2));
   fit_res_control->SetName(Form("fitResult_%s",correction_name.c_str())); fit_res_control->Write();

   // Find the scale of ndata/nmc to normalize the yields
   double nmontecarlo = all_sig_mc.sumEntries();
   double ndata = num.getVal();
   RooRealVar num_MC(Form("num_MC_sigproc_%s",region.c_str()),Form("num_MC_sigproc_%s",region.c_str()),nmontecarlo); num_MC.setConstant();
   //RooRealVar scalef(Form("scalef_%s",region.c_str()),"scalef",ndata/nmontecarlo);
   RooFormulaVar scalef(Form("scalef_%s",region.c_str()),"@0/@1",RooArgList(num,num_MC));
   // uncomment make this ONLY a shape correction!
   // scalef.setVal(1);
   //scalef.setConstant(true);

   std::cout << region << " N Control Data (all) == " << all_data.sumEntries()       << std::endl;
   std::cout << region << " N Control Bkg  == " << nbkgcont    << std::endl;
   std::cout << region << " N Control MC   == " << nmontecarlo << std::endl;

   // Plot the fits
   TCanvas can_datafit(Form("%s_datafit",region.c_str()),"Data Fit",800,600); 
   RooPlot *pl = x->frame();
   all_bkg_mc.plotOn(pl,RooFit::MarkerStyle(kOpenCircle));
   all_data.plotOn(pl);
   model.plotOn(pl);
   model.plotOn(pl,RooFit::Components(pdf_background_mc->GetName()),RooFit::LineStyle(kDashed),RooFit::LineColor(1));
   model.plotOn(pl,RooFit::VisualizeError(*fit_res_control,2),RooFit::FillColor(kYellow));
   model.plotOn(pl,RooFit::VisualizeError(*fit_res_control,1),RooFit::FillColor(kGreen) );
   model.plotOn(pl);
   all_bkg_mc.plotOn(pl,RooFit::MarkerStyle(kOpenCircle));
   all_data.plotOn(pl);
   //model.paramOn(pl);
   pl->SetMinimum(0.001);
   pl->SetTitle("");
   pl->GetXaxis()->SetTitle("fake MET (GeV)");
   pl->Draw();
   can_datafit.SetLogy();
   can_datafit.Write();

   TCanvas can_mcfit(Form("%s_mcfit",region.c_str()),"MC Fit",800,600); 
   RooPlot *plmc = x->frame();
   all_sig_mc.plotOn(plmc,RooFit::MarkerColor(kBlack));
   pdf_mc->plotOn(plmc,RooFit::LineStyle(1),RooFit::LineColor(2));
   //pdf_mc->paramOn(plmc);
   plmc->SetMinimum(0.001);
   plmc->SetTitle("");
   plmc->GetXaxis()->SetTitle("fake MET (GeV)");
   plmc->Draw();
   can_mcfit.SetLogy();
   can_mcfit.Write();

   TCanvas can_mcdatafit(Form("%s_mcdatafit",region.c_str()),"MC and Data Fits",800,600); 
   RooPlot *plmcdata = x->frame();
   pdf_mc->plotOn(plmcdata,RooFit::LineColor(2),RooFit::Normalization(nmontecarlo));
   pdf->plotOn(plmcdata,RooFit::Normalization(ndata));
   plmcdata->SetMinimum(0.001);
   plmcdata->GetXaxis()->SetTitle("fake MET (GeV)");
   plmcdata->Draw();
   can_mcdatafit.Write();

   // Import the correction and the models 
   RooFormulaVar ratio(Form("%s",correction_name.c_str()),"@0*@1/@2",RooArgList(scalef,*pdf,*pdf_mc));
   wspace->import(ratio);
   wspace->import(num);

   TCanvas can_ra(Form("%s_ratio",region.c_str()),"MC Fit",800,600); 
   RooPlot *plra = x->frame();
   ratio.plotOn(plra,RooFit::LineStyle(1));
   plra->SetTitle("");
   plra->GetYaxis()->SetTitle("r(MET) Correction Data/MC");
   plra->Draw();
   can_ra.Write();
   */
}

void ModelBuilder::addSample(std::string name, std::string region, std::string process, bool is_mc, bool is_signal, bool saveDataset){
 
   TH1F *tmp_hist;

//   std::cout << "VAR - " << varstring << std::endl;
//   std::cout << "CUT - " << cutstring << std::endl;
//   std::cout << "W   - " << weightname << std::endl;

   if (! fIn->Get(name.c_str())){
   	std::cout << Form("Error, no tree %s found in %s", name.c_str(), fIn->GetName()) << std::endl;
	assert(0);
   }

   // This is needed to make all things into a Dataset (for the fitting)
   std::cout << "Appending Sample " << name << std::endl;

   RooDataSet *tmp_data;
   RooRealVar *var    = wspace->var(varstring.c_str());
   RooRealVar *weight = wspace->var(weightname.c_str());

   RooArgSet treevariables(*var,*weight);
   // List all branches in the tree and add them as variables 
   TObjArray *branches = (TObjArray*) ((TTree*)fIn->Get(name.c_str()))->GetListOfBranches();
//   branches->Print();
   TIter next(branches); TBranch *br;
   while ( (br = (TBranch*)next()) ) {
	const char *brname = br->GetName();
	if ( std::strcmp(brname,weightname.c_str())!=0 && std::strcmp(brname,varstring.c_str())!=0 ){
	  RooRealVar *vartmp = new RooRealVar(brname,brname,0,1); vartmp->removeRange();
	  std::cout << "Adding variable" << vartmp->GetName() <<  std::endl;
	  treevariables.add(*vartmp);
	}
   }
   /**************************************************************************/
   if (! (wspace->genobj("treevars"))) {
   	//treevariables.SetName("treevars");
     wspace->import(treevariables);
   }
   
   std::string lcutstring = cutstring;
   std::map<std::string,std::string>::iterator it_ecut = extracuts.find(region);
   if ( it_ecut != extracuts.end() ) {
   	lcutstring+=" && "+(*it_ecut).second;
   	// Also allow only a specific process to be cut 
   }
   it_ecut = extracuts.find(process);
   if ( it_ecut != extracuts.end() ) {	
   	   lcutstring+=" && "+(*it_ecut).second;
   }

   std::cout << " CUT STRING FOR " << process <<  ", in " << region  << " : " << lcutstring.c_str() <<"  "<<saveDataset<< std::endl;
   if (is_mc) {
     std::cout<<" generating template for"<<name.c_str()<<"  "<<varstring<<"  "<<weightname<<"  "<<lcutstring<<"  "<<Form("_tmphist%s",catname.c_str())<<std::endl;
     tmp_hist = (TH1F*)generateTemplate(lTmp,(TTree *)fIn->Get(name.c_str()),varstring,weightname,lcutstring,Form("_tmphist%s",catname.c_str()));
        if (saveDataset)  tmp_data = new RooDataSet("tmpdata","dataset",treevariables,RooFit::Import(*(TTree*)fIn->Get(name.c_str())),RooFit::Cut(lcutstring.c_str()),RooFit::WeightVar(weightname.c_str()));
   } else {
   	tmp_hist = (TH1F*)generateTemplate(lTmp,(TTree *)fIn->Get(name.c_str()),varstring,"",lcutstring,Form("_tmphist%s",catname.c_str()));
   	if (saveDataset) tmp_data = new RooDataSet("tmpdata","dataset",treevariables,RooFit::Import(*(TTree*)fIn->Get(name.c_str())),RooFit::Cut(lcutstring.c_str()));
   }

   std::map<std::string,ControlRegion>::iterator it_sample = v_samples.find(region);
   
   int type;
   if (!is_mc) type=0;
   else if (is_signal) type=-1;
   else type=1;

   TString tRegion = region;
   //tRegion.ReplaceAll("_Met","");
   std::string pRegion = tRegion.Data();
   if (it_sample!=v_samples.end()) {

     bool proc_exists = has_process((*it_sample).second,process);
     std::cout<<" proc_exists for = "<<pRegion+std::string("_")+process<<"  "<<proc_exists<<std::endl;
     if (proc_exists){
       std::cout<<" proc_exists for = "<<pRegion+std::string("_")+process<<std::endl;
        save_hists[pRegion+std::string("_")+process]->Add(tmp_hist);


	if (saveDataset) {
	  save_datas[pRegion+std::string("_")+process]->append(*tmp_data);
	  //std::cout << "Adding to existing sample -- " <<pRegion+std::string("_")+process << " -> " <<  name.c_str() << ", dataset=" << tmp_data->sumEntries()  << " histogram=" <<tmp_hist->Integral() << std::endl; 
	}
     } else {
     	//std::cout << "CREATE NEW PROCESS" << process << ", sample" << name << std::endl;
   	tmp_hist->SetName(Form("%s_%s",pRegion.c_str(),process.c_str()));
	if (saveDataset) tmp_data->SetName(Form("%s_%s",pRegion.c_str(),process.c_str()));
	//fOut->WriteTObject(tmp_hist);
	save_hists[pRegion+std::string("_")+process] = tmp_hist;



	if (saveDataset) save_datas[pRegion+std::string("_")+process] = tmp_data;

	((*it_sample).second).procs.push_back(std::pair<std::string,int> (process,type));

     }
   }
   else {
	ControlRegion cregion;
	cregion.name = region;
	cregion.procs.push_back(std::pair<std::string,int> (process,type));
        v_samples[region] = cregion;

   	tmp_hist->SetName(Form("%s_%s",pRegion.c_str(),process.c_str()));
	save_hists.insert(std::pair<std::string,TH1F*> (pRegion+std::string("_")+process,tmp_hist));


	
	if (saveDataset){
   	 tmp_data->SetName(Form("%s_%s",pRegion.c_str(),process.c_str()));
	 save_datas.insert(std::pair<std::string,RooDataSet*> (pRegion+std::string("_")+process,tmp_data));
	 //fOut->WriteTObject(tmp_hist);
	}
   } 

}


/*


void ModelBuilder::buildAndFitModels(std::string region){

   // Loop though processes for this region
   // Make 3 datasets, signal, bkg and data
   

}
*/
