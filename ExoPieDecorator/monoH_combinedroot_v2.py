import sys,os,array

from array import array
from glob import glob

from ROOT import TFile, gROOT, kBlack,TH1F

gROOT.SetBatch(True)

#CRSRPath = '/Users/dekumar/MEGA/Fullwork/2017_Plotting/22102019/monoHROOT'
#CRSRPath = '/Users/dekumar/MEGA/Fullwork/2017_Plotting/LimitFiles/bkg'

#oldPath=12052020_fixedEra_v2_B
boosted_CRSRPath = '/afs/cern.ch/work/d/dekumar/public/monoH/monoHbbPlottingFiles/CMSSW_10_3_0/src/2017/PlotFiles/28102020_setup4_NoJER_B/monoHROOT'
Resolved_CRSRPath = '/afs/cern.ch/work/d/dekumar/public/monoH/monoHbbPlottingFiles/CMSSW_10_3_0/src/2017/PlotFiles/28102020_setup4_NoJER_R/monoHROOT'
#SignalPath = '/Users/dekumar/MEGA/Fullwork/2017_Plotting/rootFiles_Signal'
SignalPath = '/afs/cern.ch/work/d/dekumar/public/monoH/monoHbbPlottingFiles/CMSSW_10_3_0/src/2017/HistFiles/monohbb.v07.02.00.2017_NoJER_Signal/'#2017_fixedEra_v2_signal'

CRSRFiles_boosted = [boosted_CRSRPath+'/'+fl for fl in os.listdir(boosted_CRSRPath) if 'Recoil' in fl or 'MET' in fl ]
CRSRFiles_resolved = [Resolved_CRSRPath+'/'+fl for fl in os.listdir(Resolved_CRSRPath) if 'Recoil' in fl or 'MET' in fl]
SignalFiles = [SignalPath+'/'+fl for fl in os.listdir(SignalPath) if '.root' in fl]

dirName = 'DataCardRootFiles_monohbb.v07.02.00_2017'#'DataCardRootFiles_monohbb.v06.00.05_2017'

os.system('rm -rf '+dirName)
os.system('mkdir '+dirName)

def SetCanvas():

    # CMS inputs
    # -------------
    H_ref = 1000;
    W_ref = 1000;
    W = W_ref
    H  = H_ref

    T = 0.08*H_ref
    B = 0.21*H_ref
    L = 0.12*W_ref
    R = 0.08*W_ref
    # --------------

    c1 = TCanvas("c2","c2",0,0,2000,1500)
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetFrameFillStyle(0)
    c1.SetFrameBorderMode(0)
    c1.SetLeftMargin( L/W )
    c1.SetRightMargin( R/W )
    c1.SetTopMargin( T/H )
    c1.SetBottomMargin( B/H )
    c1.SetTickx(0)
    c1.SetTicky(0)
    c1.SetTickx(1)
    c1.SetTicky(1)
    c1.SetGridy()
    c1.SetGridx()
    #c1.SetLogy(1)
    return c1

def getLegend():
    legend=TLegend(.10,.79,.47,.89)
    legend.SetTextSize(0.038)
    legend.SetFillStyle(0)

    return legend


def getLatex():
    latex =  TLatex()
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(31);
    latex.SetTextAlign(11);
    return latex



def setHistStyle(h_temp,newname):

    #h_temp=h_temp2.Rebin(len(bins)-1,"h_temp",array.array('d',bins))
    h_temp.SetName(newname)
    h_temp.SetTitle(newname)
    h_temp.SetLineWidth(1)
    h_temp.SetBinContent(h_temp.GetXaxis().GetNbins(),h_temp.GetBinContent(h_temp.GetXaxis().GetNbins())+h_temp.GetBinContent(h_temp.GetXaxis().GetNbins()+1)) #Add overflow bin content to last bin
    h_temp.SetMarkerColor(kBlack);
    h_temp.SetMarkerStyle(2);
    return h_temp



def fillEmptyBin(hist,nBins,value):
    trueBins = hist.GetXaxis().GetNbins()
    if nBins>=trueBins:
	for i in range(nBins+1):
            if hist.GetBinError(i+1)<0:print hist.GetBinError(i+1)
    	    if hist.GetBinContent(i+1)==0.0:
	        hist.SetBinContent(i+1,value)
    else:
        for i in range(trueBins+1):
            if hist.GetBinContent(i+1)==0.0:
	        hist.SetBinContent(i+1,value)


CSList = {'ma_150_mA_300':1.606,'ma_150_mA_400':0.987,'ma_150_mA_500':0.5074,'ma_150_mA_600':0.2984,'ma_150_mA_1000':0.0419,'ma_150_mA_1200':0.0106,'ma_150_mA_1600':0.07525}


print ('CSList',CSList)

SRCRhistos=['bkgSum','DIBOSON','ZJets','GJets','QCD','SMH','STop','Top','WJets','DYJets','data_obs']

bins= [200,270,345,480,1000]

f=TFile(dirName+'/'+"AllMETHistos.root","RECREATE")

for infile in CRSRFiles_boosted:
    #print ('checking code for ',infile)
    fin       =   TFile(infile,"READ")
    rootFile  = infile.split('/')[-1]
    reg       = rootFile.split('_')[2]
    #cat       = rootFile.split('_')[-1].replace('.root','')
    #print ('cat',cat)
    print ('running code for ',infile)
    syst = ''
    if '_up.root' in infile or '_down.root' in infile:
	laststr = infile.split('/')[-1]
        syst    = '_'+laststr.split("_")[-2]+'_'+laststr.split("_")[-1].replace('.root','')

    if ('MET' in infile and (not 'SR'  in infile)):continue# or ('Recoil' not in infile): continue
    #if 'TopWmu' in infile or 'TopWe' in infile:continue

    reg = reg.replace('Zmumu','ZMUMU').replace('Zee','ZEE').replace('Wmu','WMU').replace('We','WE').replace('Topmu','TOPMU').replace('Tope','TOPE')

    for hist in SRCRhistos:
        temp   = fin.Get(hist)
        hist=hist.replace('DIBOSON','diboson').replace('ZJets','zjets').replace('GJets','gjets').replace('QCD','qcd').replace('SMH','smh').replace('STop','singlet').replace('Top','tt').replace('WJets','wjets').replace('DYJets','dyjets')
        #if 'boosted' in cat:
        newName   = 'monoHbb2017_B_'+reg+'_'+str(hist)+syst
        #if 'resolved' in cat:
        #    newName   = 'monoHbb2017_R_'+reg+'_'+str(hist)
        #print (temp.GetXaxis().GetNbins())
        if not syst=='' and hist=='data_obs':continue
        #print 'newName',newName
        fillEmptyBin(temp,4,0.00001)
        '''
        if temp.Integral() == 0.0:
            HISTNAME=newName
            temp = TH1F(newName, newName, 4, array('d',bins))
            # print ('=================',hist)
            # print ('=================',temp.GetXaxis().GetNbins())
            for bin in range(4):
                temp.SetBinContent(bin+1,0.00001)
        '''
        myHist = setHistStyle(temp,newName)
        f.cd()
        myHist.Write()




for infile in CRSRFiles_resolved:
    #print ('checking code for ',infile)
    fin       =   TFile(infile,"READ")
    rootFile  = infile.split('/')[-1]
    reg       = rootFile.split('_')[2]
    #cat       = rootFile.split('_')[-1].replace('.root','')
    #print ('cat',cat)

    if ('MET' in infile and (not 'SR' in infile)):continue# or ('Recoil' not in infile): continue
    #if 'TopWmu' in infile or 'TopWe' in infile:continue
    print ('running code for ',infile)
    syst = ''
    if '_up.root' in infile or '_down.root' in infile:
        laststr = infile.split('/')[-1]
        syst    = '_'+laststr.split("_")[-2]+'_'+laststr.split("_")[-1].replace('.root','')

    reg = reg.replace('Zmumu','ZMUMU').replace('Zee','ZEE').replace('Wmu','WMU').replace('We','WE').replace('Topmu','TOPMU').replace('Tope','TOPE')

    for hist in SRCRhistos:
        temp   = fin.Get(hist)
        hist=hist.replace('DIBOSON','diboson').replace('ZJets','zjets').replace('GJets','gjets').replace('QCD','qcd').replace('SMH','smh').replace('STop','singlet').replace('Top','tt').replace('WJets','wjets').replace('DYJets','dyjets')
        # if 'boosted' in cat:
        #     newName   = 'monoHbb2017_B_'+reg+'_'+str(hist)
        # if 'resolved' in cat:
        newName   = 'monoHbb2017_R_'+reg+'_'+str(hist)+syst
        #print (temp.GetXaxis().GetNbins())
        if not syst=='' and hist=='data_obs':continue
        fillEmptyBin(temp,4,0.00001)
        '''
        if temp.Integral() == 0.0:
            HISTNAME=newName
            temp = TH1F(newName, newName, 4, array('d',bins))
            # print ('=================',hist)
            # print ('=================',temp.GetXaxis().GetNbins())
            for bin in range(4):
                temp.SetBinContent(bin+1,0.00001)
        '''
        myHist = setHistStyle(temp,newName)
        f.cd()
        myHist.Write()


lumi = 41.5*1000

BR = 0.588

for infile in SignalFiles:
    print ('infile',infile)
    fin       =   TFile(infile,"READ")
    rootFile = infile.split('/')[-1]
    print rootFile.split('_')
    ma=rootFile.split('_')[12]
    mA=rootFile.split('_')[10]

    if mA=='1400': continue

    sampStr = 'ma_'+ma+'_mA_'+mA
    CS = CSList[sampStr]
    for hst in ["boosted","resolved"]:
        for syst in ["h_reg_SR_MET","h_reg_SR_MET_btagweight_up","h_reg_SR_MET_btagweight_down","h_reg_SR_MET_ewkweight_up","h_reg_SR_MET_ewkweight_down","h_reg_SR_MET_toppTweight_up","h_reg_SR_MET_toppTweight_down","h_reg_SR_MET_metTrigweight_up","h_reg_SR_MET_metTrigweight_down","h_reg_SR_MET_puweight_down","h_reg_SR_MET_puweight_up","h_reg_SR_MET_jec_up","h_reg_SR_MET_jec_down","h_reg_SR_MET_Res_up","h_reg_SR_MET_Res_down","h_reg_SR_MET_En_up","h_reg_SR_MET_En_down"]:
            temp = fin.Get(syst+"_"+hst)
            print 'hist',syst+"_"+hst

            fillEmptyBin(temp,4,0.00001)
            '''
            if  temp.Integral() == 0.0:
                for bin in range(temp.GetXaxis().GetNbins()):
                    temp.SetBinContent(bin,0.00001)
            '''
            h_total = fin.Get('h_total_mcweight')
            totalEvents = h_total.Integral()
            temp.Scale((lumi*CS*BR)/(totalEvents))
            if hst=="boosted":
                cat="B"
            if hst=="resolved":
                cat="R"
            if syst=="h_reg_SR_MET":
                samp = 'monoHbb2017_'+cat+'_SR_ggF_sp_0p35_tb_1p0_mXd_10_mA_'+mA+'_ma_'+ma
            else:samp = 'monoHbb2017_'+cat+'_SR_ggF_sp_0p35_tb_1p0_mXd_10_mA_'+mA+'_ma_'+ma+"_"+syst.split("_")[-2]+"_"+syst.split("_")[-1]
            myHist = setHistStyle(temp,samp)
            f.cd()
            myHist.Write()


f.Close()
