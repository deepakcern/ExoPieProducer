# Configuration for a simple monojet topology. Use this as a template for your own Run-2 mono-X analysis
# First provide ouput file name in out_file_name field 
out_file_name = 'mono-x.root'

# can define any thing useful here which may be common to several categories, eg binning in MET 
#bins = range(200,1200,200)
bins = [250,300,350,400,500,1000]
#bins = [250,300,350,400,500,600,1000]
# will expect samples with sample_sys_Up/Down but will skip if not found 
#systematics=["Met","FP","btag","mistag","wjethf","zjethf","gjethf"]
systematics=["Met","FP","btag","mistag"]
# Define each of the categories in a dictionary of the following form .. 
#	'name' : the category name 
#	'in_file_name' : input ntuple file for this category 
#	'cutstring': add simple cutrstring, applicable to ALL regions in this category (eg mvamet > 200)
#	'varstring': the main variable to be fit in this category (eg mvamet), must be named as the branch in the ntuples
#	'weightname': name of the weight variable 
#	'bins': binning given as a python list
#	'additionalvars': list additional variables to be histogrammed by the first stage, give as a list of lists, each list element 
#			  as ['variablename',nbins,min,max]
#	'pdfmodel': integer --> N/A  redudant for now unless we move back to parameteric fitting estimates
# 	'samples' : define tree->region/process map given as a dictionary with each entry as follows 
#		TreeName : ['region','process',isMC,isSignal] --> Note isSignal means DM/Higgs etc for signal region but Z-jets/W-jets for the di/single-muon regions !!!

#  OPTIONAL --> 'extra_cuts': additional cuts maybe specific to this control region (eg ptphoton cuts) if this key is missing, the code will not complain   
 
monotop_category = {
	    'name':"monotop"
	   #,'in_file_name':"/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_7_4_12_patch1/src/MonoX/../BaconAnalyzer/MJSelection/skim/monojet-combo-electron.root"  # Without recoil corrections
	   #,'in_file_name':"/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_7_4_12_patch1/src/MonoX/monojet-combo-electron.root_recoil"
       #,'in_file_name':"files/monotop-boosted-combo-weight.root"
#       ,'in_file_name':"files/monotop-boosted-combo-mar9.root"
       ,'in_file_name':"files/monotop-boosted-combo-apr6.root"
	   ,"cutstring":"met>250 && met<1000"
	   ,"varstring":["met",250,1000]
	   ,"weightname":"weight"
	   ,"bins":bins[:]
	   #,"bins":[200.0 , 210.0 , 220.0 , 230.0 , 240.0 , 250.0 , 260.0 , 270.0 , 280.0 , 290.0 , 300.0 , 310.0 , 320.0 , 330.0,340,360,380,420,510,1000]
  	   #,"additionalvars":[['jet1pt',25,150,1000]]
	   ,"additionalvars":[]
           ,"pdfmodel":0
	   #,"extra_cuts":[["singleelectron","rmet>40"],["photon","ptpho>200"]]
	   #,"extra_cuts":[["singleelectron","rmet>40"]]
	   ,"samples":
	   	{  
		  # Signal Region
		   "ZnunuLO_signal"  	           :['signal','zjets',1,0]
       ,"ZllLO_signal"	           :['signal','zll',1,0]
 		  ,"Wjets_signal"  	           :['signal','wjets',1,0]
		  ,"WW_signal"  	           :['signal','dibosons',1,0]
		  ,"WZ_signal"  	           :['signal','dibosons',1,0]
		  ,"ZZ_signal"  	           :['signal','dibosons',1,0]
		  ,"ttbar_signal"   	           :['signal','ttbar',1,0]
		  ,"SingleTop_signal"              :['signal','stop',1,0]
		  ,"QCD_signal"		   	   :['signal','qcd',1,0]
		  ,"data_signal"	           :['signal','data',0,0]

		  # some signals 
		  ,"Mres-900_Mchi-100_signal"		   :['signal','Mres900_Mchi100',1,1]
		  ,"Mres-1100_Mchi-100_signal"		   :['signal','Mres1100_Mchi100',1,1]
		  ,"Mres-1300_Mchi-100_signal"		   :['signal','Mres1300_Mchi100',1,1]
		  ,"Mres-1500_Mchi-100_signal"		   :['signal','Mres1500_Mchi100',1,1]
		  ,"Mres-1700_Mchi-100_signal"		   :['signal','Mres1700_Mchi100',1,1]
		  ,"Mres-1900_Mchi-100_signal"		   :['signal','Mres1900_Mchi100',1,1]
		  ,"Mres-2100_Mchi-100_signal"		   :['signal','Mres2100_Mchi100',1,1]
		  ,"Mchi-300_signal"		   :['signal','Mchi300',1,1]
      ,"Mchi-500_signal"		   :['signal','Mchi500',1,1]
      ,"Mchi-700_signal"		   :['signal','Mchi700',1,1]
      ,"Mchi-900_signal"		   :['signal','Mchi900',1,1]
      ,"Mchi-1100_signal"		   :['signal','Mchi1100',1,1]
      ,"Mchi-1300_signal"		   :['signal','Mchi1300',1,1]
      ,"Mchi-1500_signal"		   :['signal','Mchi1500',1,1]
#      ,"tHq_signal"            :['signal','tHq_signal',1,1]

		  #,"A_2000_150_g025_signal"		   :['signal','axialvector_2000150',1,1]
		  #,"A_2000_100_g025_signal"		   :['signal','axialvector_2000100',1,1]
		  #,"A_2000_10_g025_signal"		   :['signal','axialvector_2000010',1,1]
		  #,"A_2000_1_g025_signal"		   :['signal','axialvector_2000001',1,1]

		  # Di muon-Control
		  ,"ZllLO_di_muon_control"	   :['dimuon','zll',1,1]
		  ,"Wjets_di_muon_control"  	   :['dimuon','wjets',1,0]
		  ,"WW_di_muon_control"  	   :['dimuon','dibosons',1,0]
		  ,"WZ_di_muon_control"  	   :['dimuon','dibosons',1,0]
		  ,"ZZ_di_muon_control"  	   :['dimuon','dibosons',1,0]
		  ,"ttbar_di_muon_control"         :['dimuon','ttbar',1,0]
		  ,"SingleTop_di_muon_control"     :['dimuon','stop',1,0]
		  ,"data_di_muon_control"	   :['dimuon','data',0,0]

		  # Single muon (top) control
		  ,"ZllLO_single_muon_top_control"       :['singlemuontop','zll',1,0]
		  ,"Wjets_single_muon_top_control"     :['singlemuontop','wjets',1,0]
		  ,"ZZ_single_muon_top_control"        :['singlemuontop','dibosons',1,0]
		  ,"WW_single_muon_top_control"        :['singlemuontop','dibosons',1,0]
		  ,"WZ_single_muon_top_control"        :['singlemuontop','dibosons',1,0]
		  ,"SingleTop_single_muon_top_control" :['singlemuontop','stop',1,0]
		  ,"ttbar_single_muon_top_control"     :['singlemuontop','ttbar',1,1]
		  ,"QCD_single_muon_top_control"       :['singlemuontop','qcd',1,0]
		  ,"data_single_muon_top_control"      :['singlemuontop','data',0,0]

                   # Single muon (w) control
		  ,"ZllLO_single_muon_w_control"	   :['singlemuonw','zll',1,0]
		  ,"Wjets_single_muon_w_control"     :['singlemuonw','wjets',1,1]
		  ,"ZZ_single_muon_w_control"        :['singlemuonw','dibosons',1,0]
		  ,"WW_single_muon_w_control"        :['singlemuonw','dibosons',1,0]
		  ,"WZ_single_muon_w_control"        :['singlemuonw','dibosons',1,0]
		  ,"SingleTop_single_muon_w_control" :['singlemuonw','stop',1,0]
		  ,"ttbar_single_muon_w_control"     :['singlemuonw','ttbar',1,0]
		  ,"QCD_single_muon_w_control"	   :['singlemuonw','qcd',1,0]
		  ,"data_single_muon_w_control"	   :['singlemuonw','data',0,0]

		  # Photon control region
		  ,"data_photon_control"	   :['photon','data',0,0]
		  ,"Photon_photon_control"	   :['photon','gjets',1,1]
                  ,"QCD_photon_control"	   	   :['photon','qcd',1,0]

		  # Di electron-Control
		  ,"ZllLO_di_electron_control"	   :['dielectron','zll',1,1]
		  ,"Wjets_di_electron_control"     :['dielectron','wjets',1,0]
		  ,"WW_di_electron_control"  	   :['dielectron','dibosons',1,0]
		  ,"WZ_di_electron_control"  	   :['dielectron','dibosons',1,0]
		  ,"ZZ_di_electron_control"  	   :['dielectron','dibosons',1,0]
		  ,"ttbar_di_electron_control"      :['dielectron','ttbar',1,0]
		  ,"SingleTop_di_electron_control"  :['dielectron','stop',1,0]
		  ,"data_di_electron_control"	   :['dielectron','data',0,0]

		  # Single electron (top) control
		  ,"ZllLO_single_electron_top_control"       :['singleelectrontop','zll',1,0]
		  ,"Wjets_single_electron_top_control"     :['singleelectrontop','wjets',1,0]
		  ,"ZZ_single_electron_top_control"        :['singleelectrontop','dibosons',1,0]
		  ,"WW_single_electron_top_control"        :['singleelectrontop','dibosons',1,0]
		  ,"WZ_single_electron_top_control"        :['singleelectrontop','dibosons',1,0]
		  ,"SingleTop_single_electron_top_control" :['singleelectrontop','stop',1,0]
		  ,"ttbar_single_electron_top_control"     :['singleelectrontop','ttbar',1,1]
		  ,"QCD_single_electron_top_control"       :['singleelectrontop','qcd',1,0]
		  ,"data_single_electron_top_control"      :['singleelectrontop','data',0,0]

                   # Single electron (w) control
		  ,"ZllLO_single_electron_w_control"       :['singleelectronw','zll',1,0]
		  ,"Wjets_single_electron_w_control"     :['singleelectronw','wjets',1,1]
		  ,"ZZ_single_electron_w_control"        :['singleelectronw','dibosons',1,0]
		  ,"WW_single_electron_w_control"        :['singleelectronw','dibosons',1,0]
		  ,"WZ_single_electron_w_control"        :['singleelectronw','dibosons',1,0]
		  ,"SingleTop_single_electron_w_control" :['singleelectronw','stop',1,0]
		  ,"ttbar_single_electron_w_control"     :['singleelectronw','ttbar',1,0]
		  ,"QCD_single_electron_w_control"       :['singleelectronw','qcd',1,0]
		  ,"data_single_electron_w_control"      :['singleelectronw','data',0,0]
	   	},

}
categories = [monotop_category]
