import sys,os,array

from array import array
from glob import glob

from ROOT import TFile, gROOT, kBlack,TH1F

gROOT.SetBatch(True)

#CRSRPath = '/Users/dekumar/MEGA/Fullwork/2017_Plotting/22102019/monoHROOT'
#CRSRPath = '/Users/dekumar/MEGA/Fullwork/2017_Plotting/LimitFiles/bkg'

#oldPath=12052020_fixedEra_v2_B
boosted_CRSRPath = '/afs/cern.ch/work/d/dekumar/public/monoH/monoHbbPlottingFiles/CMSSW_10_3_0/src/2017/PlotFiles/08072021_v12.07.07_noDDBOnZCR_B/monoHROOT'
#boosted_CRSRPath='/Users/deepakkumar/monoHbb/fitModel/rootFiles/bkg/'
Resolved_CRSRPath = '/afs/cern.ch/work/d/dekumar/public/monoH/monoHbbPlottingFiles/CMSSW_10_3_0/src/2017/PlotFiles/09072021_v12.07.07_R_HP/monoHROOT'#09072021_v12.07.07_R_LP/monoHROOT'
#Resolved_CRSRPath = '/afs/cern.ch/work/d/dekumar/public/monoH/monoHbbPlottingFiles/CMSSW_10_3_0/src/2017/PlotFiles/02072021_change_nBcond_R/monoHROOT'
#Resolved_CRSRPath ='/Users/deepakkumar/monoHbb/fitModel/rootFiles/bkg/'
#SignalPath = '/Users/deepakkumar/monoHbb/fitModel/rootFiles/signal/'
SignalPath = '/afs/cern.ch/work/d/dekumar/public/monoH/monoHbbPlottingFiles/CMSSW_10_3_0/src/2017/HistFiles/monohbb.v12.07.07.2017_signal_HP'#monohbb.v12.07.07.2017_signal_HP/'

CRSRFiles_boosted = [boosted_CRSRPath+'/'+fl for fl in os.listdir(boosted_CRSRPath) if 'Recoil' in fl or 'MET' in fl ]
CRSRFiles_resolved = [Resolved_CRSRPath+'/'+fl for fl in os.listdir(Resolved_CRSRPath) if 'Recoil' in fl or 'MET' in fl]
THDMaSignalFiles = [SignalPath+'/'+fl for fl in os.listdir(SignalPath) if 'ggTomonoH' in fl]
Zp2HDMSignalFiles = [SignalPath+'/'+fl for fl in os.listdir(SignalPath) if 'ZprimeToA0hToA0chichihbb' in fl]
ZpBaryonicSignalFiles = [SignalPath+'/'+fl for fl in os.listdir(SignalPath) if 'ZpBaryonic' in fl]

# print SignalFiles
dirName = 'DataCardRootFiles_monohbb.v12.07.07_NoddbOnZCR_renorm_PDF_scale_HP'#

os.system('rm -rf '+dirName)
os.system('mkdir '+dirName)

def getXsec(ma='',mA='',tb='',st=''):
    fs = open('genLevel2HDMa_XSec.txt')
    xsec=0.0
    isXsecPresent=False
    for line in fs:
        line = line.rstrip()
        # print 'checking line %s  for provided inputs : '%line , 'ma,mA,tb,st  ', ma,mA,tb,st
        ma_= line.split()[5]
        mA_= line.split()[4]
        tb_= line.split()[2]
        st_= line.split()[1]
        Xsec = line.split()[6]
        proc = line.split()[0]

        if proc=='gg' and int(ma)==int(ma_) and int(mA)==int(mA_) and float(tb)==float(tb_) and float(st)==float(st_):
            isXsecPresent = True
            break

    if isXsecPresent:
        xsec = float(Xsec)
    return xsec
        

def getZpBxSec(mZp, mChi):
    fs = open('ZpBxsection.txt','r')  
    xsec=0.0
    isXsecPresent=False
    for line in fs:
        zp = line.split()[0].split('_')[1].replace('MZp','')
        chi= line.split()[0].split('_')[2].replace('MChi','')
        #print "checking masses  zp, chi  ", zp, chi
        if mZp == zp and mChi == chi:
            Xsec = line.split()[4] 
            isXsecPresent=True
    if isXsecPresent:
        xsec = float(Xsec)
    return xsec


def getZp2HDMxSec(mZp, mA0):
    fs = open('Zprime_xsecs.txt','r')
    xsec=0.0
    isXsecPresent=False
    lines = fs.readlines()
    for i, line in enumerate(lines):
        if (not 'name' in line) or 'Zinv' in line:continue
    #	print line,  'line  ', line.split()
        zp = line.split()[-1].split('_')[3]
        a0 = line.split()[-1].split('_')[4]
        #print "checking masses  zp, chi  ", zp, chi
        if mZp == zp and mA0 == a0:
            Xsec = lines[i+1].split()[2]
            isXsecPresent=True
    if isXsecPresent:
        xsec = float(Xsec)
    return xsec


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
                if hist.GetBinContent(i+1)==0.0:
                    hist.SetBinContent(i+1,value)
    else:
        for i in range(trueBins+1):
            if hist.GetBinContent(i+1)==0.0:
                hist.SetBinContent(i+1,value)

def renameSyst(syst,year="2017"):
    syst = syst.replace('_MET','').replace('_Recoil','').replace('btagweight','CMS2017_eff_b').replace('metTrigweight','CMS2017_trig_met').replace('pdfWeight','CMS2017_pdf').replace('puweight','CMS2017_PU').replace('scaleWeight','CMS2017_mu_scale')
    syst = syst.replace('toppTweight','CMS2017_eff_t').replace('year',year).replace('ewkweight','ewk')
    syst = syst.replace('eleID','CMS2017_EleID')
    syst = syst.replace('eleReco','CMS2017_EleRECO')
    syst = syst.replace('muID','CMS2017_MuID').replace('muIso','CMS2017_MuISO')
    syst = syst.replace('_up','Up').replace('_down','Down')
    syst = syst.replace('_Absolute','_JECAbsolute').replace('_BBEC1','_JECBBEC1').replace('_EC2','_JECEC2').replace('_FlavorQCD','_JECFlavorQCD').replace('_HF','_JECHF').replace('_Relative','_JECRelative')

    return syst


def getExtraHists(hist,pefix):
    print ('setting extra unc using this hist', hist)
    metbins = [250, 300, 370,500,1000]
    nbins = len(metbins) -1 
    h_CMS2017_prefireUp = TH1F(pefix+'_CMS2017_prefireUp',pefix+'_CMS2017_prefireUp',nbins,array('d',metbins))
    h_CMS2017_prefireDown = TH1F(pefix+'_CMS2017_prefireDown',pefix+'_CMS2017_prefireDown',nbins,array('d',metbins))
    h_CMS2017_fake_bUp = TH1F(pefix+'_CMS2017_fake_bUp',pefix+'_CMS2017_fake_bUp',nbins,array('d',metbins))
    h_CMS2017_fake_bDown = TH1F(pefix+'_CMS2017_fake_bDown',pefix+'_CMS2017_fake_bDown',nbins,array('d',metbins))

    hists = [h_CMS2017_prefireUp,h_CMS2017_prefireDown,h_CMS2017_fake_bUp,h_CMS2017_fake_bDown]

    for ij in range(nbins):
        value = hist.GetBinContent(ij+1)
        value_prefire_up = value + value * 0.01
        value_prefire_down = value - value * 0.01
        h_CMS2017_prefireUp.SetBinContent(ij+1,value_prefire_up)
        h_CMS2017_prefireDown.SetBinContent(ij+1,value_prefire_down)


        value_fakeb_up = value + value * 0.05
        value_fakeb_down = value - value * 0.05

        h_CMS2017_fake_bUp.SetBinContent(ij+1,value_fakeb_up)
        h_CMS2017_fake_bDown.SetBinContent(ij+1,value_fakeb_down)

        #print 'value ', value, 'value_prefire_up', value_prefire_up,'value_prefire_down',value_prefire_down
        #print 'value',value,'value_fakeb_up',value_fakeb_up,'value_fakeb_down',value_fakeb_down

    return hists


def renormHist(h_center, h_up_down):
    print (h_center, h_center.Integral())
    print (h_up_down, h_up_down.Integral())
    h_up_down.Scale(h_center.Integral()/h_up_down.Integral())
    print (h_up_down, h_up_down.Integral())
    return h_up_down


CSList = {'ma_150_mA_300':1.606,'ma_150_mA_400':0.987,'ma_150_mA_500':0.5074,'ma_150_mA_600':0.2984,'ma_150_mA_1000':0.0419,'ma_150_mA_1200':0.0106,'ma_150_mA_1600':0.07525}


# print ('CSList',CSList)

SRCRhistos=['bkgSum','DIBOSON','ZJets','GJets','QCD','SMH','STop','Top','WJets','DYJets','data_obs']

bins= [200,270,345,480,1000]

f=TFile(dirName+'/'+"AllMETHistos.root","RECREATE")


for infile in CRSRFiles_boosted:
    print "running code for boosted"
    #print ('checking code for ',infile)
    fin       =   TFile(infile,"READ")
    rootFile  = infile.split('/')[-1]
    reg       = rootFile.split('_')[2]
    #cat       = rootFile.split('_')[-1].replace('.root','')
    #print ('cat',cat)
    print ('running code for ',infile)
    syst = ''
    if '_up.root' in infile or '_down.root' in infile or 'Up' in infile or 'Down' in infile:
	laststr = infile.split('/')[-1]
        syst    = '_'+laststr.split("_")[-2]+'_'+laststr.split("_")[-1].replace('.root','').replace('Up','_up').replace('Down','_down')

    if ('MET' in infile and (not 'SR'  in infile)):continue# or ('Recoil' not in infile): continue
    #if 'TopWmu' in infile or 'TopWe' in infile:continue

    reg = reg.replace('Zmumu','ZMUMU').replace('Zee','ZEE').replace('Wmu','WMU').replace('We','WE').replace('Topmu','TOPMU').replace('Tope','TOPE')

    for hist in SRCRhistos:
        temp   = fin.Get(hist)
        hist=hist.replace('DIBOSON','diboson').replace('ZJets','zjets').replace('GJets','gjets').replace('QCD','qcd').replace('SMH','smh').replace('STop','singlet').replace('Top','tt').replace('WJets','wjets').replace('DYJets','dyjets')
        #if 'boosted' in cat:
	syst = renameSyst(syst)
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
	temp_unc_hists=[]
        if syst=='' and hist!='data_obs':
            temp_unc_hists=getExtraHists(temp,newName)
        f.cd()
        myHist.Write()
        for ihst in temp_unc_hists:
            ihst.Write()




for infile in CRSRFiles_resolved:
    print ('checking code for ',infile)
    fin       =   TFile(infile,"READ")
    rootFile  = infile.split('/')[-1]
    reg       = rootFile.split('_')[2]
    if 'bdtscore' in rootFile:
        reg =reg+"_bdtscore"
    #cat       = rootFile.split('_')[-1].replace('.root','')
    #print ('cat',cat)

    if ('MET' in infile and (not 'SR' in infile)):continue# or ('Recoil' not in infile): continue
    #if 'TopWmu' in infile or 'TopWe' in infile:continue
    if 'bdtscore' in rootFile:print ('running code for ',infile)
    syst = ''
    if '_up.root' in infile or '_down.root' in infile or 'Up' in infile or 'Down' in infile:
        laststr = infile.split('/')[-1]
        syst    = '_'+laststr.split("_")[-2]+'_'+laststr.split("_")[-1].replace('.root','').replace('Up','_up').replace('Down','_down')
	# print 'systematics ', syst

    reg = reg.replace('Zmumu','ZMUMU').replace('Zee','ZEE').replace('Wmu','WMU').replace('We','WE').replace('Topmu','TOPMU').replace('Tope','TOPE')

    for hist in SRCRhistos:
        temp   = fin.Get(hist)
        hist=hist.replace('DIBOSON','diboson').replace('ZJets','zjets').replace('GJets','gjets').replace('QCD','qcd').replace('SMH','smh').replace('STop','singlet').replace('Top','tt').replace('WJets','wjets').replace('DYJets','dyjets')
        # if 'boosted' in cat:
        #     newName   = 'monoHbb2017_B_'+reg+'_'+str(hist)
        # if 'resolved' in cat:
        syst = renameSyst(syst)
        newName   = 'monoHbb2017_R_'+reg+'_'+str(hist)+syst
        # print newName
        #print (temp.GetXaxis().GetNbins())
        if not syst=='' and hist=='data_obs':continue
        fillEmptyBin(temp,4,0.00001)

        myHist = setHistStyle(temp,newName)
	temp_unc_hists =[]
        if syst=='' and hist!='data_obs':
	    # print 'temp',temp,'newName',newName
            temp_unc_hists=getExtraHists(temp,newName)
        f.cd()
        myHist.Write()
        for ihst in temp_unc_hists:
	    # print 'ihst',ihst
            ihst.Write()


lumi = 41.5*1000

BR = 0.588



for infile in THDMaSignalFiles:
    print ('infile',infile)
    fin       =   TFile(infile,"READ")
    fin2       =   TFile(infile,"READ")
    rootFile = infile.split('/')[-1]
    # print rootFile.split('_')
    ma=rootFile.split('_')[-8]
    mA=rootFile.split('_')[-4]
    tb=rootFile.split('_')[6]
    st=rootFile.split('_')[4]

    # if mA=='1400': continue

    sampStr = 'ma_'+ma+'_mA_'+mA
    CS = getXsec(ma=ma,mA=mA,tb=tb.split('p')[0]+'.'+tb.split('p')[1],st=st.split('p')[0]+'.'+st.split('p')[1]) #CSList[sampStr]
    # print ("CS  %s"%CS)
    if CS==0.0:
        print ("Xsec is not found")
        exit()

    h_boost_tmp = fin2.Get("h_reg_SR_MET_boosted")
    h_resol_tmp = fin2.Get("h_reg_SR_MET_resolved")
    h_total_ = fin2.Get('h_total_mcweight')
    h_boost_tmp.Scale((lumi*CS*BR)/(h_total_.Integral()))
    h_resol_tmp.Scale((lumi*CS*BR)/(h_total_.Integral()))


    for hst in ["boosted","resolved"]:
        for syst in ["h_reg_SR_MET","h_reg_SR_MET_btagweight_up","h_reg_SR_MET_btagweight_down","h_reg_SR_MET_ewkweight_up","h_reg_SR_MET_ewkweight_down","h_reg_SR_MET_toppTweight_up","h_reg_SR_MET_toppTweight_down","h_reg_SR_MET_metTrigweight_up","h_reg_SR_MET_metTrigweight_down","h_reg_SR_MET_puweight_down","h_reg_SR_MET_puweight_up","h_reg_SR_MET_jec_up","h_reg_SR_MET_jec_down","h_reg_SR_MET_Res_up","h_reg_SR_MET_Res_down","h_reg_SR_MET_En_up","h_reg_SR_MET_En_down",
"h_reg_SR_MET_AbsoluteUp", "h_reg_SR_MET_Absolute_yearUp", "h_reg_SR_MET_BBEC1Up", "h_reg_SR_MET_BBEC1_yearUp", "h_reg_SR_MET_EC2Up", "h_reg_SR_MET_EC2_yearUp","h_reg_SR_MET_FlavorQCDUp", "h_reg_SR_MET_HFUp", "h_reg_SR_MET_HF_yearUp", "h_reg_SR_MET_RelativeBalUp", "h_reg_SR_MET_RelativeSample_yearUp","h_reg_SR_MET_AbsoluteDown", "h_reg_SR_MET_Absolute_yearDown", "h_reg_SR_MET_BBEC1Down", "h_reg_SR_MET_BBEC1_yearDown", "h_reg_SR_MET_EC2Down", "h_reg_SR_MET_EC2_yearDown","h_reg_SR_MET_FlavorQCDDown", "h_reg_SR_MET_HFDown", "h_reg_SR_MET_HF_yearDown", "h_reg_SR_MET_RelativeBalDown", "h_reg_SR_MET_RelativeSample_yearDown","h_reg_SR_MET_scaleWeightUp", "h_reg_SR_MET_scaleWeightDown", "h_reg_SR_MET_pdfWeightUp", "h_reg_SR_MET_pdfWeightDown"]:

            temp = fin.Get(syst+"_"+hst)
            # print ('hist',syst+"_"+hst)

            fillEmptyBin(temp,4,0.00001)

            h_total = fin.Get('h_total_mcweight')
            totalEvents = h_total.Integral()
            temp.Scale((lumi*CS*BR)/(totalEvents))
        
            temp_unc_hists=[]
            if hst=="boosted":
                cat="B"
            if hst=="resolved":
                cat="R"
            if syst=="h_reg_SR_MET":
                samp = 'monoHbb2017_'+cat+'_SR_ggF_sp_'+st+'_tb_'+tb+'_mXd_10_mA_'+mA+'_ma_'+ma
                temp_unc_hists=getExtraHists(temp,samp)
		# print "samp  ",samp,'temp_unc_hists',temp_unc_hists
            else:
                samp = 'monoHbb2017_'+cat+'_SR_ggF_sp_'+st+'_tb_'+tb+'_mXd_10_mA_'+mA+'_ma_'+ma+"_"+syst.split("_")[-2]+"_"+syst.split("_")[-1].replace('Down','_down').replace('Up','_up')
                samp=renameSyst(samp)
                if "pdfUp" in samp or "pdfDown" in samp or 'mu_scale' in samp:
                    print ("renormalizing pdf hist")
                    print ("center  :   ",'monoHbb2017_'+cat+'_SR_ggF_sp_'+st+'_tb_'+tb+'_mXd_10_mA_'+mA+'_ma_'+ma)
                    print ("up or down  :  ",temp)
                    if cat=="R":
                        temp = renormHist(h_resol_tmp,temp)
                    if cat=="B":
                        temp = renormHist(h_boost_tmp,temp)

            myHist = setHistStyle(temp,samp)
            f.cd()
            myHist.Write()
            for ihist in temp_unc_hists:
                ihist.Write()

    # for fittype in ["bdtscore"]:
        
    #     temp = fin.Get("h_reg_SR_bdtscore_resolved")
    #     fillEmptyBin(temp,4,0.00001)
    #     h_total = fin.Get('h_total_mcweight')
    #     totalEvents = h_total.Integral()
	# # if 'MH4_150' in infile.split('/')[-1] and '0p35_tn_1p0' in infile.split('/')[-1]:
	# # 	print ("===============")
	# # 	print ("ifile   :  ",infile.split('/')[-1])
	# # 	print ("sel events  ",temp.Integral()," total events  ",totalEvents , " eff     ",temp.Integral()/totalEvents)
	# # 	print ("===================")
    #     temp.Scale((lumi*CS*BR)/(totalEvents))
    #     samp = 'monoHbb2017_R_SR_bdtscore_ggF_sp_'+st+'_tb_'+tb+'_mXd_10_mA_'+mA+'_ma_'+ma
    #     myHist = setHistStyle(temp,samp)
    #     f.cd()
    #     myHist.Write()


for infile in Zp2HDMSignalFiles:
    print ('infile  ',infile)
    fin = TFile(infile,"READ")
    fin2 = TFile(infile,"READ")
    rootFile = infile.split('/')[-1]
    
    MZp = rootFile.split('_')[3]
    MA0 = rootFile.split('_')[4]
    CS = getZp2HDMxSec(MZp, MA0)#1.0

    # print ('CS  ', CS,  ' MZp ',MZp,' MA0   ',MA0)
    h_boost_tmp = fin2.Get("h_reg_SR_MET_boosted")
    h_resol_tmp = fin2.Get("h_reg_SR_MET_resolved")
    h_total_ = fin2.Get('h_total_mcweight')
    h_boost_tmp.Scale((lumi*CS*BR)/(h_total_.Integral()))
    h_resol_tmp.Scale((lumi*CS*BR)/(h_total_.Integral()))


    for hst in ["boosted","resolved"]:
        for syst in ["h_reg_SR_MET","h_reg_SR_MET_btagweight_up","h_reg_SR_MET_btagweight_down","h_reg_SR_MET_ewkweight_up","h_reg_SR_MET_ewkweight_down","h_reg_SR_MET_toppTweight_up","h_reg_SR_MET_toppTweight_down","h_reg_SR_MET_metTrigweight_up","h_reg_SR_MET_metTrigweight_down","h_reg_SR_MET_puweight_down","h_reg_SR_MET_puweight_up","h_reg_SR_MET_jec_up","h_reg_SR_MET_jec_down","h_reg_SR_MET_Res_up","h_reg_SR_MET_Res_down","h_reg_SR_MET_En_up","h_reg_SR_MET_En_down",
"h_reg_SR_MET_AbsoluteUp", "h_reg_SR_MET_Absolute_yearUp", "h_reg_SR_MET_BBEC1Up", "h_reg_SR_MET_BBEC1_yearUp", "h_reg_SR_MET_EC2Up", "h_reg_SR_MET_EC2_yearUp","h_reg_SR_MET_FlavorQCDUp", "h_reg_SR_MET_HFUp", "h_reg_SR_MET_HF_yearUp", "h_reg_SR_MET_RelativeBalUp", "h_reg_SR_MET_RelativeSample_yearUp","h_reg_SR_MET_AbsoluteDown", "h_reg_SR_MET_Absolute_yearDown", "h_reg_SR_MET_BBEC1Down", "h_reg_SR_MET_BBEC1_yearDown", "h_reg_SR_MET_EC2Down", "h_reg_SR_MET_EC2_yearDown","h_reg_SR_MET_FlavorQCDDown", "h_reg_SR_MET_HFDown", "h_reg_SR_MET_HF_yearDown", "h_reg_SR_MET_RelativeBalDown", "h_reg_SR_MET_RelativeSample_yearDown","h_reg_SR_MET_scaleWeightUp", "h_reg_SR_MET_scaleWeightDown", "h_reg_SR_MET_pdfWeightUp", "h_reg_SR_MET_pdfWeightDown"]:
            temp = fin.Get(syst+"_"+hst)


            fillEmptyBin(temp,4,0.00001)

            h_total = fin.Get('h_total_mcweight')
            totalEvents = h_total.Integral()
            temp.Scale((lumi*CS*BR)/(totalEvents))
	    temp_unc_hists=[]
            if hst=="boosted":
                cat="B"
            if hst=="resolved":
                cat="R"
            if syst=="h_reg_SR_MET":
                samp = 'monoHbb2017_'+cat+'_SR_'+MZp+'_'+MA0
                temp_unc_hists=getExtraHists(temp,samp)
            else:
                samp = 'monoHbb2017_'+cat+'_SR_'+MZp+'_'+MA0+"_"+syst.split("_")[-2]+"_"+syst.split("_")[-1].replace('Down','_down').replace('Up','_up')
                samp=renameSyst(samp)

                if "pdfUp" in samp or "pdfDown" in samp or 'mu_scale' in samp:
                    print ("renormalizing pdf hist")
                    print ("center  :   ",'monoHbb2017_'+cat+'_SR_ggF_sp_'+st+'_tb_'+tb+'_mXd_10_mA_'+mA+'_ma_'+ma)
                    print ("up or down  :  ",temp)
                    if cat=="R":
                        temp = renormHist(h_resol_tmp,temp)
                    if cat=="B":
                        temp = renormHist(h_boost_tmp,temp)


            myHist = setHistStyle(temp,samp)
            f.cd()
            myHist.Write()
            for ihist in temp_unc_hists:
                ihist.Write()


for infile in ZpBaryonicSignalFiles:
    print ('infile  ',infile)
    fin = TFile(infile,"READ")
    fin2 = TFile(infile,"READ")
    rootFile = infile.split('/')[-1]
    
    MZp = rootFile.split('_')[-3]
    Mchi = rootFile.split('_')[-1].replace('.root','')
    CS = getZpBxSec(MZp, Mchi)#1.0
    # print ('ZpB Xsec', CS, 'MZp', MZp, 'Mchi', Mchi)
    h_boost_tmp = fin2.Get("h_reg_SR_MET_boosted")
    h_resol_tmp = fin2.Get("h_reg_SR_MET_resolved")
    h_total_ = fin2.Get('h_total_mcweight')
    h_boost_tmp.Scale((lumi*CS*BR)/(h_total_.Integral()))
    h_resol_tmp.Scale((lumi*CS*BR)/(h_total_.Integral()))


    for hst in ["boosted","resolved"]:
        for syst in ["h_reg_SR_MET","h_reg_SR_MET_btagweight_up","h_reg_SR_MET_btagweight_down","h_reg_SR_MET_ewkweight_up","h_reg_SR_MET_ewkweight_down","h_reg_SR_MET_toppTweight_up","h_reg_SR_MET_toppTweight_down","h_reg_SR_MET_metTrigweight_up","h_reg_SR_MET_metTrigweight_down","h_reg_SR_MET_puweight_down","h_reg_SR_MET_puweight_up","h_reg_SR_MET_jec_up","h_reg_SR_MET_jec_down","h_reg_SR_MET_Res_up","h_reg_SR_MET_Res_down","h_reg_SR_MET_En_up","h_reg_SR_MET_En_down",
"h_reg_SR_MET_AbsoluteUp", "h_reg_SR_MET_Absolute_yearUp", "h_reg_SR_MET_BBEC1Up", "h_reg_SR_MET_BBEC1_yearUp", "h_reg_SR_MET_EC2Up", "h_reg_SR_MET_EC2_yearUp","h_reg_SR_MET_FlavorQCDUp", "h_reg_SR_MET_HFUp", "h_reg_SR_MET_HF_yearUp", "h_reg_SR_MET_RelativeBalUp", "h_reg_SR_MET_RelativeSample_yearUp","h_reg_SR_MET_AbsoluteDown", "h_reg_SR_MET_Absolute_yearDown", "h_reg_SR_MET_BBEC1Down", "h_reg_SR_MET_BBEC1_yearDown", "h_reg_SR_MET_EC2Down", "h_reg_SR_MET_EC2_yearDown","h_reg_SR_MET_FlavorQCDDown", "h_reg_SR_MET_HFDown", "h_reg_SR_MET_HF_yearDown", "h_reg_SR_MET_RelativeBalDown", "h_reg_SR_MET_RelativeSample_yearDown","h_reg_SR_MET_scaleWeightUp", "h_reg_SR_MET_scaleWeightDown", "h_reg_SR_MET_pdfWeightUp", "h_reg_SR_MET_pdfWeightDown"]:
            temp = fin.Get(syst+"_"+hst)
            # print 'hist',syst+"_"+hst

            fillEmptyBin(temp,4,0.00001)

            h_total = fin.Get('h_total_mcweight')
            totalEvents = h_total.Integral()
            temp.Scale((lumi*CS*BR)/(totalEvents))
	    temp_unc_hists=[]
            if hst=="boosted":
                cat="B"
            if hst=="resolved":
                cat="R"
            if syst=="h_reg_SR_MET":
                samp = 'monoHbb2017_'+cat+'_SR_MZp_'+MZp+'_Mchi_'+Mchi
                temp_unc_hists=getExtraHists(temp,samp)
            else:
                samp = 'monoHbb2017_'+cat+'_SR_MZp_'+MZp+'_Mchi_'+Mchi+"_"+syst.split("_")[-2]+"_"+syst.split("_")[-1].replace('Down','_down').replace('Up','_up')
                samp=renameSyst(samp)

                if "pdfUp" in samp or "pdfDown" in samp or 'mu_scale' in samp:
                    print ("renormalizing pdf hist")
                    print ("center  :   ",'monoHbb2017_'+cat+'_SR_ggF_sp_'+st+'_tb_'+tb+'_mXd_10_mA_'+mA+'_ma_'+ma)
                    print ("up or down  :  ",temp)
                    if cat=="R":
                        temp = renormHist(h_resol_tmp,temp)
                    if cat=="B":
                        temp = renormHist(h_boost_tmp,temp)


            myHist = setHistStyle(temp,samp)
            f.cd()
            myHist.Write()
            for ihist in temp_unc_hists:
                ihist.Write()

f.Close()


