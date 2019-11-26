#include <iostream>
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"

using namespace RooFit ;

/*
BuildStatsModel.py: Package to build statistical fitting model for background estimation and limit extraction
Author: Raman Khurana                                                                                                                                                               Date:   26-September-2018
V0:     Simple implementation of the model using transfer factors
*/


/*

naming convention: h_region_ process,

e.g.

bincontent: bincontents_sr1_wjets

Histogram:  h_sr1_data, h_sr2_wjets, h_wenu_wjets

Histogram Transfer Factor: htf_

datahist:  use dh_wenu_wjets

RooRealVar: rrv_wenu_wjets
RooRealVar: rrvbc_wenu_wjets (for boncontent)

RooArgList: ral_wenu_wjets
RooArgList: ralbc_wenu_wjets (for the bin content)

RooParamHist: rph_wenu_wjets

RooAddition: rph_norm_wenu_wjets

RooFormulaVar: rfv_

RooFormulaVar:
*/

//gInterpreter->GenerateDictionary("vector<RooFormulaVar>","RooFormulaVar.h;vector");

TFile* OpenRootFile(TString filename, TString mode="READ"){
  TFile* f = new TFile (filename, mode);
  return f;
}


std::vector<float> GetBinContents(TH1F* h){
  std::vector<float> bcs;
  bcs.clear();

  int nbinsx = h->GetNbinsX();
  for (auto ibin=1; ibin<=nbinsx; ibin++){
    bcs.push_back(h->GetBinContent(ibin));
  }

  return bcs;

}




void addTemplate(RooWorkspace& ws,  RooArgList& vars, TH1F* hist) {
  std::cout<<" name = "<<hist->GetName()<<std::endl;
  RooDataHist rhist(hist->GetName(), hist->GetName(),  vars, hist);
  std::cout<<" integral of the histogram for "<<hist->GetName()<<" is "<<rhist.sumEntries()<<"  "<<hist->Integral()<<std::endl;
  ws.import(rhist);
}


// following two functions do essentially same things, but instead of returning RooArgList, one return vector, becuase I couldn't find any way to get roorealvar back from the rooarglist.

RooArgList GetRooArgList(std::vector<float>  bcs, TString name){
  RooArgList ral_ ;
  TString postfix;

  for (auto i=0; i<bcs.size(); i++){
    postfix.Form("%d", i+1);
    RooRealVar rrv_(name+postfix,"Background yield in signal region, bin 1", bcs[i], bcs[i]/2, bcs[i]*2);
    ral_.addClone(rrv_);
  }
  return ral_;
}


// This function is overloaded
std::vector <RooRealVar> GetRooRealVar(std::vector<float>  bcs, TString name){
  std::vector<RooRealVar> rrvV_ ;
  rrvV_.clear();
  TString postfix;

  for (auto i=0; i<bcs.size(); i++){
    postfix.Form("%d", i+1);
    // fix the naming here using some automation  and also in the next function

    rrvV_.push_back(RooRealVar(name+postfix,"Background yield in signal region, bin 1", bcs[i], bcs[i]/5, bcs[i]*5));

  }

  return rrvV_;
}



std::vector <RooRealVar> GetRooRealVar(std::vector<float>  bcs, bool setconstant, TString name){
  std::vector<RooRealVar> rrvV_ ;
  rrvV_.clear();
  TString postfix;
  for (auto i=0; i<bcs.size(); i++){
    std::cout<<" ---------- inside GetRooRealVar "<<std::endl;
    // fix the naming here using some automation
    postfix.Form("%d", i+1);
    RooRealVar rrv_(name+postfix, "Background yield in signal region, bin 1", bcs[i]);
    rrv_.setConstant(1);
    rrvV_.push_back(rrv_);
  }

  return rrvV_;
}



// I still don't know why I can't return simple RooFormulaVar from a function but the std::vector<RooFormulaVar> works fine.
std::vector<RooFormulaVar> GetRooFormulaVar(std::vector<RooRealVar> nuisances, std::vector <RooRealVar> rrv_htf_wenu_2b_wjets, TString name){
  std::vector<RooFormulaVar> test_;
  test_.clear();
  TString postfix;

  for (auto i=0; i<rrv_htf_wenu_2b_wjets.size(); i++){
    postfix.Form("%d", i+1);
    // write correct naming here
    RooFormulaVar tf(name+postfix,"Trasnfer factor","@2*TMath::Power(1.01,@0)*TMath::Power(1.02,@1)",RooArgList(nuisances[0], nuisances[1], rrv_htf_wenu_2b_wjets[i]) );
    test_.push_back(tf);
  }
  return test_;

}






/*
1: roorealvar, here it is mer
2: background histogram in signal region
3: background histogram in CR
4: data histogram in signal region
5: workspace
6: string to save names etc, as per convention, region_proc for the CR
7: string to save names etc, as per convention, region_proc for the SR
8: output file
*/
void createRegion(RooRealVar met, TH1F* h_sr2_wjets , TH1F* h_wenu_2b_wjets, TH1F* h_sr2_data, RooWorkspace& wspace, TString region_proc_cr, TString region_proc_sr, TFile* fOut){

  RooArgList vars(met);

  /* Get the bin content of each bin of the histogram in a vector which can be used later */
  std::vector<float> bincontents_sr2_wjets = GetBinContents(h_sr2_wjets);

  // This will create the RooRealVar with a/5 to a*5 range.
  RooArgList ralbc_sr2_wjets  = GetRooArgList(bincontents_sr2_wjets, "ralbc_"+region_proc_sr);

  // create a vector of RooRealVar, this is needed because I didn't find  way to retrive the RooRealVar back from the RooArgList
  std::vector<RooRealVar> rrvbc_sr2_wjets = GetRooRealVar(bincontents_sr2_wjets, "rrvbc_"+region_proc_sr);

  // Create a RooParametericHist which contains those yields, last argument is just for the binning, we can use the data TH1 for that
  RooParametricHist rph_sr2_wjets("rph_"+region_proc_sr, "wjets PDF in signal region",met,ralbc_sr2_wjets, *h_sr2_data);

  // Always include a _norm term which should be the sum of the yields (thats how combine likes to play with pdfs)
  RooAddition rph_sr2_wjets_norm("rph_"+region_proc_sr+"_norm","Total Number of events from background in signal region",ralbc_sr2_wjets);


  wspace.import(rph_sr2_wjets);
  wspace.import(rph_sr2_wjets_norm,RooFit::RecycleConflictNodes());


  /*
    For the control region, the background process will be dependent on the yields of the background in the signal region using a transfer factor.
    The transfer factor TF must account for acceptance/efficiency etc differences in the signal to control regions.
    In this case we define the transfer factor as: ratio of the WJets (electron) yield in the WJets control region and
    WJets in the Signal region.

    For each bin a transfer factor is calculated and the nuisance parameters are associated with this.

    We could imagine that the transfer factor could be associated with some uncertainty - lets say a 1% uncertainty due to efficiency and 2% due to acceptance.
    We need to make nuisance parameters ourselves to model this and give them a nominal value of 0.
  */

  RooRealVar efficiency("efficiency", "efficiency nuisance parameter",0);
  RooRealVar acceptance("acceptance", "acceptance nuisance parameter",0);


  /*
  We need to make the transfer factor a function of these parameters since variations in these uncertainties will lead to variations of the transfer factor. Here we've assumed Log-normal effects (i.e the same as putting lnN in the CR datacard) but we could use any function which could be used to parameterise the effect - eg if the systematic is due to some alternate template, we could use polynomials for example.
  */


  // --------------------------------------------------------------
  // ------------------------Wjets (electrron) Control region ------
  // --------------------------------------------------------------

  // create roodatahist of the background histogram in CR.
  RooDataHist dh_wenu_2b_wjets("dh_"+region_proc_cr,"dh_"+region_proc_cr, vars, h_wenu_2b_wjets);

  // another copy fo the wjets in wenu CR for division and saving thr TFs central value.
  TH1F* htf_wenu_2b_wjets = (TH1F*) h_wenu_2b_wjets->Clone();
  htf_wenu_2b_wjets->Divide(h_sr2_wjets);

  // Get bin content of each bin of this ratio histogram and save it in the RooRealVar which will be used later for the Actual Transfer Factor with effect of Nuisance parameters included
  // idelaly each of these rooreal var in following vector should be setConstat(1) otherwise it may be treated as free parameter however it should be fixed.
  std::vector <float> bincontents_htf_wenu_2b_wjets =  GetBinContents(htf_wenu_2b_wjets);

  // the rooreal vars are set to constant values instead of free values.
  // this is an overloaded function
  std::vector <RooRealVar> rrv_htf_wenu_2b_wjets = GetRooRealVar(bincontents_htf_wenu_2b_wjets, true);

  // the nuisance part of the code has to be updated after some more discussion.
  std::vector<RooRealVar> nuisances;
  nuisances.clear();
  nuisances.push_back(efficiency);
  nuisances.push_back(acceptance);




  // It crashes if I use the element of vector
  std::vector<RooFormulaVar> rfv_tf_wenu_2b_wjets = GetRooFormulaVar(nuisances, rrv_htf_wenu_2b_wjets, "rfv_tf_wenu_2b_wjets");

  //However following work fine if I do it one by one
  // this is the only part which is not automatic becuase vector of RooFormulaVar is not working, Once this is done i should use the binning of TF also automatic


  /*
  RooFormulaVar TF[7];
  TString istr;
  for (auto i=0; i<bincontents_htf_wenu_2b_wjets.size(); i++){
    istr.Form("%d", i+1);
    TF[i] = RooFormulaVar("TF"+istr+region_proc_cr,"Trasnfer factor","@2*TMath::Power(1.01,@0)*TMath::Power(1.02,@1)",RooArgList(efficiency,acceptance,rrv_htf_wenu_2b_wjets[i]));
  }
  */

  RooFormulaVar TF1("TF1"+region_proc_cr,"Trasnfer factor","@2*TMath::Power(1.01,@0)*TMath::Power(1.02,@1)",RooArgList(efficiency,acceptance,rrv_htf_wenu_2b_wjets[0]) );
  RooFormulaVar TF2("TF2"+region_proc_cr,"Trasnfer factor","@2*TMath::Power(1.01,@0)*TMath::Power(1.02,@1)",RooArgList(efficiency,acceptance,rrv_htf_wenu_2b_wjets[1]) );
  RooFormulaVar TF3("TF3"+region_proc_cr,"Trasnfer factor","@2*TMath::Power(1.01,@0)*TMath::Power(1.02,@1)",RooArgList(efficiency,acceptance,rrv_htf_wenu_2b_wjets[2]) );
  RooFormulaVar TF4("TF4"+region_proc_cr,"Trasnfer factor","@2*TMath::Power(1.01,@0)*TMath::Power(1.02,@1)",RooArgList(efficiency,acceptance,rrv_htf_wenu_2b_wjets[3]) );
  RooFormulaVar TF5("TF5"+region_proc_cr,"Trasnfer factor","@2*TMath::Power(1.01,@0)*TMath::Power(1.02,@1)",RooArgList(efficiency,acceptance,rrv_htf_wenu_2b_wjets[4]) );
  RooFormulaVar TF6("TF6"+region_proc_cr,"Trasnfer factor","@2*TMath::Power(1.01,@0)*TMath::Power(1.02,@1)",RooArgList(efficiency,acceptance,rrv_htf_wenu_2b_wjets[5]) );
  //RooFormulaVar TF7("TF7"+region_proc_cr,"Trasnfer factor","@2*TMath::Power(1.01,@0)*TMath::Power(1.02,@1)",RooArgList(efficiency,acceptance,rrv_htf_wenu_2b_wjets[6]) );


  std::cout<<" size of rfv_tf_wenu_2b_wjets = "<<rfv_tf_wenu_2b_wjets.size()<<std::endl;
  /*
    Then need to make each bin of the background in the control region a function of the background in the signal and the transfer factor -
    i.e NCR=NSR x TF
  */




  RooFormulaVar rfv_wenu_2b_wjets1("rfv_"+region_proc_cr+"1","Background yield in control region, bin 1","@0*@1",RooArgList(TF1, rrvbc_sr2_wjets.at(0)));
  RooFormulaVar rfv_wenu_2b_wjets2("rfv_"+region_proc_cr+"2","Background yield in control region, bin 2","@0*@1",RooArgList(TF2,rrvbc_sr2_wjets.at(1)));
  RooFormulaVar rfv_wenu_2b_wjets3("rfv_"+region_proc_cr+"3","Background yield in control region, bin 3","@0*@1",RooArgList(TF3,rrvbc_sr2_wjets.at(2)));
  RooFormulaVar rfv_wenu_2b_wjets4("rfv_"+region_proc_cr+"4","Background yield in control region, bin 4","@0*@1",RooArgList(TF4,rrvbc_sr2_wjets.at(3)));
  RooFormulaVar rfv_wenu_2b_wjets5("rfv_"+region_proc_cr+"5","Background yield in control region, bin 5","@0*@1",RooArgList(TF5,rrvbc_sr2_wjets.at(4)));
  RooFormulaVar rfv_wenu_2b_wjets6("rfv_"+region_proc_cr+"6","Background yield in control region, bin 6","@0*@1",RooArgList(TF6,rrvbc_sr2_wjets.at(5)));
  //RooFormulaVar rfv_wenu_2b_wjets7("rfv_"+region_proc_cr+"7","Background yield in control region, bin 7","@0*@1",RooArgList(TF7,rrvbc_sr2_wjets.at(6)));


  /*
  RooFormulaVar CRbin1("wjets_CR_bin1","Background yield in control region, bin 1","@0*@1",RooArgList((rfv_tf_wenu_2b_wjets[0]), rrvbc_sr2_wjets.at(0)));
  RooFormulaVar CRbin2("wjets_CR_bin2","Background yield in control region, bin 2","@0*@1",RooArgList((rfv_tf_wenu_2b_wjets[1]),rrvbc_sr2_wjets.at(1)));
  RooFormulaVar CRbin3("wjets_CR_bin3","Background yield in control region, bin 3","@0*@1",RooArgList((rfv_tf_wenu_2b_wjets[2]),rrvbc_sr2_wjets.at(2)));
  RooFormulaVar CRbin4("wjets_CR_bin4","Background yield in control region, bin 4","@0*@1",RooArgList((rfv_tf_wenu_2b_wjets[3]),rrvbc_sr2_wjets.at(3)));
  RooFormulaVar CRbin5("wjets_CR_bin5","Background yield in control region, bin 5","@0*@1",RooArgList((rfv_tf_wenu_2b_wjets[4]),rrvbc_sr2_wjets.at(4)));
  RooFormulaVar CRbin6("wjets_CR_bin6","Background yield in control region, bin 6","@0*@1",RooArgList((rfv_tf_wenu_2b_wjets[5]),rrvbc_sr2_wjets.at(5)));
  RooFormulaVar CRbin7("wjets_CR_bin7","Background yield in control region, bin 7","@0*@1",RooArgList((rfv_tf_wenu_2b_wjets[6]),rrvbc_sr2_wjets.at(6)));
  */


  std::cout<<" aftercrbin1 "<<std::endl;


  // --------------------------------------------------------------
  // ------------------------WJets (muon ) Control region ---------
  // --------------------------------------------------------------

  RooArgList ral_wenu_2b_wjets;
  ral_wenu_2b_wjets.add(rfv_wenu_2b_wjets1);
  ral_wenu_2b_wjets.add(rfv_wenu_2b_wjets2);
  ral_wenu_2b_wjets.add(rfv_wenu_2b_wjets3);
  ral_wenu_2b_wjets.add(rfv_wenu_2b_wjets4);
  ral_wenu_2b_wjets.add(rfv_wenu_2b_wjets5);
  ral_wenu_2b_wjets.add(rfv_wenu_2b_wjets6);
  //ral_wenu_2b_wjets.add(rfv_wenu_2b_wjets7);


  std::cout<<" after RooArgList wjets_CR_bins"<<std::endl;
  RooParametricHist rph_wenu_2b_wjets("rph_"+region_proc_cr, "Background PDF in control region",met,ral_wenu_2b_wjets, *h_sr2_data);

  std::cout<<" after RooParametricHist"<<std::endl;
  RooAddition rph_wenu_2b_wjets_norm("rph_"+region_proc_cr+"_norm","Total Number of events from background in control region", ral_wenu_2b_wjets);

  std::cout<<" after RooAddition"<<std::endl;


  wspace.import(rph_wenu_2b_wjets);

  std::cout<<" imported p_CR_wjets "<<std::endl;
  wspace.import(rph_wenu_2b_wjets_norm ,RooFit::RecycleConflictNodes());

  std::cout<<" norm of imported p_CR_wjets "<<std::endl;


  //fOut->cd();
  //wspace.Write();



}







void PrepareWS(){
  ///afs/cern.ch/work/p/ptiwari/public/bbDM/WCr_Split/AllMETHistos.root
  TString inputfile   = "AllMETHistos.root";
  TString outputfile  = "bbDM_WS.root";
  TString year        = "_2016";
  TString version     = "_V0";

  int met_low = 200;
  int met_hi = 1000;

  Double_t bins[]={200,250,350,500,1000};
  Int_t  binnum = sizeof(bins)/sizeof(Double_t) - 1;

  // As usual, load the combine library to get access to the RooParametricHist
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  // Output file and workspace
  TFile* fOut = OpenRootFile(outputfile,"RECREATE");
  RooWorkspace wspace("ws_bbDM","ws_bbDM");

  // A search in a MET tail, define MET as our variable
  RooRealVar met("met","p_{T}^{miss}",met_low, met_hi);
  RooArgList vars(met);

  std::cout<<" debug 2" <<std::endl;


  // Open input file with all the histograms.
  TFile* fin = OpenRootFile(inputfile);

  /*
  // --------------------------------------------------------------
  // ------------------------Signal region ------------------------
  // --------------------------------------------------------------
  // Get the data histogram from Signal region, this is met in our case. This is a 7 bin histogram in this case. This can be optimised later on.
  TH1F* h_sr2_data = (TH1F*) fin->Get("SR_2b_data_obs");

  std::cout<<" debug 3a" <<std::endl;
  // convert the histogram into RooDataHist
  RooDataHist dh_sr2_data("dh_sr2_data","dh_sr2_data",vars, h_sr2_data);

    std::cout<<" debug 3b" <<std::endl;

  // Import just created RooDataHist into the workspace.
  wspace.import(dh_sr2_data);
  std::cout<<" debug 3c" <<std::endl;
  */


  //the following lines create a freely floating parameter for each of our bins (in this example, there are only 7 bins, defined for our observable met.
  // In this case we vary the normalisation in each bin of the background from N/3 to 3*N, e.g. if actual content in the histogram is 55 then we initialize
  // it with 55 and vary it from 55/3 to 55*3. which is very close to freely floating. This can be checked if this works for the cases when bin content is very low,
  // specially in the tails.

  /*
 KEY: TH1FWmunu_2b_DIBOSON;1Wmunu_2b_DIBOSON;
 KEY: TH1FWmunu_2b_ZJets;1Wmunu_2b_ZJets;
 KEY: TH1FWmunu_2b_GJets;1Wmunu_2b_GJets;
 KEY: TH1FWmunu_2b_STop;1Wmunu_2b_STop;
 KEY: TH1FWmunu_2b_Top;1Wmunu_2b_Top;
 KEY: TH1FWmunu_2b_WJets;1Wmunu_2b_WJets;
 KEY: TH1FWmunu_2b_DYJets;1Wmunu_2b_DYJets;
 KEY: TH1FWmunu_2b_QCD;1Wmunu_2b_QCD;
 KEY: TH1FWmunu_2b_data_obs;1Wmunu_2b_data_obs;
  */



  //-------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------- W enu CR -----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------
  /*

  std::cout<<" calling function for Wenu"<<std::endl;
  // Get the wjets histogram in signal region
  TH1F* h_sr2_wjets = (TH1F*) fin->Get("SR_2b_WJets");
  // Get the wjets hostogram in the Wenu CR
  TH1F* h_wenu_2b_wjets = (TH1F*) fin->Get("Wenu_2b_WJets");
  // Create all the inputs needed for this CR
  createRegion(met, h_sr2_wjets, h_wenu_2b_wjets, h_sr2_data, wspace, "wenu_2b_wjets", "sr2_wjets",  fOut);



  //-------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------- W munu CR -----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------


  std::cout<<" calling function for Wmunu"<<std::endl;
  // Get the wjets hostogram in the Wmunu CR
  TH1F* h_wmunu_2b_wjets = (TH1F*) fin->Get("Wmunu_2b_WJets");
  // Create all the inputs needed for this CR
  createRegion(met, h_sr2_wjets, h_wmunu_2b_wjets, h_sr2_data, wspace, "wmunu_2b_wjets", "sr2_wjets",  fOut);


  //-------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------- Top mu CR -----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------

  std::cout<<" calling function for Top mu"<<std::endl;
  TH1F* h_sr2_top = (TH1F*) fin->Get("SR_2b_Top");
  // Get the top hostogram in the Top mu CR
  TH1F* h_topmu_2b_top = (TH1F*) fin->Get("Topmunu_2b_Top");
  // Create all the inputs needed for this CR
  createRegion(met, h_sr2_top, h_topmu_2b_top, h_sr2_data, wspace, "topmu_2b_top", "sr2_top",  fOut);



  //-------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------- Top e CR -----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------

  std::cout<<" calling function for Top e"<<std::endl;
  //TH1F* h_sr2_top = (TH1F*) fin->Get("SR2_2b_Top");
  // Get the top hostogram in the Top mu CR
  TH1F* h_tope_2b_top = (TH1F*) fin->Get("Topenu_2b_Top");
  // Create all the inputs needed for this CR
  createRegion(met, h_sr2_top, h_tope_2b_top, h_sr2_data, wspace, "tope_2b_top", "sr2_top",  fOut);



  //-------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------- Zmumu CR -----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------

  std::cout<<" calling function for Zmumu"<<std::endl;
  TH1F* h_sr2_Z = (TH1F*) fin->Get("SR_2b_ZJets");
  // Get the top hostogram in the Top mu CR
  TH1F* h_Zmumu_2b_Z = (TH1F*) fin->Get("Zmumu_2b_ZJets");
  // Create all the inputs needed for this CR
  createRegion(met, h_sr2_Z, h_Zmumu_2b_Z, h_sr2_data, wspace, "Zmumu_2b_Z", "sr2_Z",  fOut);


  //-------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------- Zee CR -----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------

  // Get the top hostogram in the Top mu CR
  TH1F* h_Zee_2b_Z = (TH1F*) fin->Get("Zee_2b_ZJets");
  // Create all the inputs needed for this CR
  createRegion(met, h_sr2_Z, h_Zee_2b_Z, h_sr2_data, wspace, "Zee_2b_Z", "sr2_Z",  fOut);

  */

  //-------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------- Signal -------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------

  int signalpoint[]={50,100,150,250,300,350,400};
  Int_t  nsig = sizeof(signalpoint)/sizeof(int);
  std::vector<TString> category;
  category.push_back("2b");
  category.push_back("1b");

  TString mps;
  for (auto is=0; is<nsig; is++){
    mps.Form("%d",signalpoint[is]);
    //bbDM2016_2b_SR_2HDMa_Ma750_MChi1_MA1200_tb35_st_0p7
    //addTemplate(wspace, vars, (TH1F*) fin->Get("bbNLO_pseudo_2b_Mchi_1_Mphi_"+mps ) );
    addTemplate(wspace, vars, (TH1F*) fin->Get("bbDM2016_2b_SR_2HDMa_Ma"+mps+"_MChi1_MA600_tb35_st_0p7" ) );
    //addTemplate(wspace, vars, (TH1F*) fin->Get("bbDM2016_2b_SR_2HDMa_Ma"+mps+"_MChi1_MA1200_tb35_st_0p7" ) );

  }

  /*
  addTemplate(wspace, vars, (TH1F*) fin->Get("SR_2b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("SR_1b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Wenu_1b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Wenu_2b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Topenu_1b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Topenu_2b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Wmunu_1b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Wmunu_2b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Topmunu_1b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Topmunu_2b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Zee_1b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Zee_2b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Zmumu_1b_data_obs" ) );
  addTemplate(wspace, vars, (TH1F*) fin->Get("Zmumu_2b_data_obs" ) );
  */

  // all other histograms
  std::vector<TString> regions;
  //regions.push_back("SR");
  regions.push_back("TOPENUCR");
  regions.push_back("TOPMUNUCR");
  regions.push_back("WENUCR");
  regions.push_back("WMUNUCR");
  regions.push_back("ZEECR");
  regions.push_back("ZMUMUCR");

  std::vector<TString> process;
  process.push_back("qcd");
  process.push_back("zjets");
  process.push_back("gjets");
  process.push_back("wjets");
  process.push_back("dyjets");
  process.push_back("tt");
  process.push_back("singlet");


  TString tempname;
  for (auto ir=0; ir<regions.size(); ir++){
    for (auto ip=0; ip<process.size(); ip++){
      for (auto ic=0; ic<category.size(); ic++){
        tempname = "bbDM2016_"+category[ic] + "_" + regions[ir] + "_" + process[ip];
        std::cout<<" saving "<<tempname<<std::endl;
        addTemplate(wspace, vars, (TH1F*) fin->Get(tempname)  );
      }
    }
  }



  // write the workspace at the very end, once everthing has been imported to the workspace
  fOut->cd();
  wspace.Write();

}
