import os
import sys
import datetime
import sys, optparse
import ROOT as ROOT
import array

usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)

parser.add_option("-d", "--data", dest="datasetname")
parser.add_option("-s", "--sr", action="store_true", dest="plotSRs")
parser.add_option("-m", "--mu", action="store_true", dest="plotMuRegs")
parser.add_option("-e", "--ele", action="store_true", dest="plotEleRegs")
parser.add_option("-p", "--pho", action="store_true", dest="plotPhoRegs")
parser.add_option("-q", "--qcd", action="store_true", dest="plotQCDRegs")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose")
parser.add_option("-y", "--year", dest="year", default="Year")

(options, args) = parser.parse_args()

if options.plotSRs==None:
    makeSRplots = False
else:
    makeSRplots = options.plotSRs

if options.plotMuRegs==None:
    makeMuCRplots = False
else:
    makeMuCRplots = options.plotMuRegs

if options.plotEleRegs==None:
    makeEleCRplots = False
else:
    makeEleCRplots = options.plotEleRegs

if options.plotPhoRegs==None:
    makePhoCRplots = False
else:
    makePhoCRplots = options.plotPhoRegs

if options.plotQCDRegs==None:
    makeQCDCRplots = False
else:
    makeQCDCRplots = options.plotQCDRegs

if options.verbose==None:
    verbose = False
else:
    verbose = options.verbose

if options.datasetname.upper()=="SE":
    dtset="SE"
elif options.datasetname.upper()=="SP":
    dtset="SP"
elif options.datasetname.upper()=="SM":
    dtset="SM"
else:
    dtset="MET"
print ("Using dataset "+dtset)

runOn2016 = False
runOn2017 = False
if options.year=='2016':
    runOn2016=True
elif options.year=='2017':
    runOn2017=True
else:
    print('Please provide on which year you want to run?')

if runOn2016:
    import sample_xsec_2016 as sample_xsec
    luminosity = 35.82 * 1000
    luminosity_ = 35.82
elif runOn2017:
    import sample_xsec_2017 as sample_xsec
    luminosity = 41.0 * 1000
    luminosity_ = 41.0


datestr = str(datetime.date.today().strftime("%d%m%Y"))

#path='/Users/dekumar/MEGA/Fullwork/2017_Plotting/rootFiles_Oct5'
path='hadd_outputs_2017_14112019'
#path='hadd_outputs_13112019'
#path='/home/deepak/MEGA/Fullwork/2017_Plotting/rootFiles'
#path='/Users/dekumar/Desktop/test/bkg_data'
#os.system("ls "+path+" | cat > samplelist.txt")

def set_overflow(hist):
    bin_num = hist.GetXaxis().GetNbins()
    #print (bin_num)
    hist.SetBinContent(bin_num,hist.GetBinContent(bin_num+1)+hist.GetBinContent(bin_num)) #Add overflow bin content to last bin
    hist.SetBinContent(bin_num+1,0.)
    return hist
# bins=[200,270,345,480,1000]

def setHistStyle(h_temp2,hist):
    dovarbin=False
    h_temp_=h_temp2
    # elif 'hadrecoil' in hist or 'met' in hist:
    #     # xv = 680.0/15.0;
    #     bins=[200,270,345,480,1000]
    #     dovarbin= True
    #
    # if dovarbin:
    #     h_temp_=h_temp2.Rebin(len(bins)-1,"h_temp",array.array('d',bins))
    #     h_temp_.SetBinContent(len(bins)-1,h_temp_.GetBinContent(len(bins)-1)+h_temp_.GetBinContent(len(bins))) #Add overflow bin content to last bin
    #     h_temp_.SetBinContent(len(bins),0.)
    #     # h_temp_.GetXaxis().SetRangeUser(200,1000)
    #     # h_temp_.SetMarkerColor(kBlack);
    #     # h_temp_.SetMarkerStyle(2);
    #
    # if not dovarbin:
    #     h_temp_=h_temp2
    return h_temp_



def makeplot(loc,hist,titleX,XMIN,XMAX,Rebin,ISLOG,NORATIOPLOT,reg,varBin):
    # try:

    print ('plotting histogram:   ',hist)
    isrebin=False #bool(varBin)
    if runOn2016:
        files=open("samplelist_2016.txt","r")
    elif runOn2017:
        files=open("samplelist_2017.txt","r")


    ROOT.gStyle.SetOptStat(0);
    ROOT.gStyle.SetOptTitle(0);
    ROOT.gStyle.SetFrameLineWidth(3);
    #gStyle->SetErrorX(0);
    ROOT.gStyle.SetLineWidth(1);

    if '_SR_1b' in hist:
        histolabel="SR_1b"
    elif '_SR_2b' in hist:
        histolabel="SR_2b"
    elif 'ZmumuCR_1b' in hist:
        histolabel="Z CR (#mu#mu) 1b"
    elif 'ZeeCR_1b' in hist:
        histolabel="Z CR (ee) 1b"
    elif 'WmunuCR_1b' in hist :
        histolabel="W CR (#mu) 1b"
    elif 'WenuCR_1b' in hist:
        histolabel="W CR (e) 1b"
    elif 'TopmunuCR_1b' in hist:
        histolabel="Top CR (#mu) 1b"
    elif 'TopenuCR_1b' in hist:
        histolabel="Top CR (e) 1b"

    elif 'ZmumuCR_2b' in hist:
        histolabel="Z CR (#mu#mu) 2b"
    elif 'ZeeCR_2b' in hist:
        histolabel="Z CR (ee) 2b"
    elif 'WmunuCR_2b' in hist :
        histolabel="W CR (#mu) 2b"
    elif 'WenuCR_2b' in hist:
        histolabel="W CR (e) 2b"
    elif 'TopmunuCR_2b' in hist:
        histolabel="Top CR (#mu) 2b"
    elif 'TopenuCR_2b' in hist:
        histolabel="Top CR (e) 2b"

    else:
        histolabel="testing"

    xsec=1.0
    norm = 1.0
    BLINDFACTOR = 1.0
    r_fold = 'rootFiles/'
    DIBOSON = ROOT.TH1F()
    Top = ROOT.TH1F()
    WJets = ROOT.TH1F()
    DYJets = ROOT.TH1F()
    ZJets = ROOT.TH1F()
    STop = ROOT.TH1F()
    GJets = ROOT.TH1F()
    QCD = ROOT.TH1F()
    SMH = ROOT.TH1F()

    DYJets_Hits   = []; ZJets_Hits   = []; WJets_Hists   = []; GJets_Hists  = []; DIBOSON_Hists = []; STop_Hists   = []; Top_Hists     = []; QCD_Hists    = [];
    MET_Hist      = []; SE_Hist      = []

    count=0
    for file in files.readlines()[:]:
        myFile=path+'/'+file.rstrip()
        print ('running for file',myFile)
        print ('histName',hist)
        Str=str(count)
        exec("f"+Str+"=ROOT.TFile(myFile,'READ')",locals(), globals())
        exec("h_temp=f"+Str+".Get("+"\'"+str(hist)+"\'"+")",locals(), globals())
        exec("h_total_weight=f"+Str+".Get('h_total_mcweight')",locals(), globals())
        total_events = h_total_weight.Integral()
        #print ('selected events',h_temp.Integral())

        if 'WJetsToLNu_HT' in file:
            xsec = sample_xsec.getXsec(file)
            #print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0): normlisation=(xsec*luminosity)/(total_events)
            else: normlisation=0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                WJets_Hists.append(h_temp2)
            else:WJets_Hists.append(h_temp)

        elif 'DYJetsToLL_M-50' in file:
            xsec = sample_xsec.getXsec(file)
            #print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0): normlisation=(xsec*luminosity)/(total_events)
            else: normlisation=0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                DYJets_Hits.append(h_temp2)
            else:DYJets_Hits.append(h_temp)

        elif 'ZJetsToNuNu' in file:
            xsec = sample_xsec.getXsec(file)
            #print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0): normlisation=(xsec*luminosity)/(total_events)
            else: normlisation=0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                ZJets_Hits.append(h_temp2)
            else:ZJets_Hits.append(h_temp)

        elif 'GJets_HT' in file:
            xsec = sample_xsec.getXsec(file)
            #print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0): normlisation=(xsec*luminosity)/(total_events)
            else: normlisation=0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                GJets_Hists.append(h_temp2)
            else:GJets_Hists.append(h_temp)

        elif ('WWTo' in file) or ('WZTo' in file) or ('ZZTo' in file) or ('WW_' in file) or ('ZZ_' in file) or ('WZ_' in file) :
            xsec = sample_xsec.getXsec(file)
            #print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0): normlisation=(xsec*luminosity)/(total_events)
            else: normlisation=0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                DIBOSON_Hists.append(h_temp2)
            else:DIBOSON_Hists.append(h_temp)


        elif ('ST_t' in file) or ('ST_s' in file):
            xsec = sample_xsec.getXsec(file)
            #print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0): normlisation=(xsec*luminosity)/(total_events)
            else: normlisation=0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                STop_Hists.append(h_temp2)
            else:STop_Hists.append(h_temp)

        elif 'TTTo' in file:
            xsec = sample_xsec.getXsec(file)
            #print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0): normlisation=(xsec*luminosity)/(total_events)
            else: normlisation=0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                Top_Hists.append(h_temp2)

            else:Top_Hists.append(h_temp)

        elif 'QCD_HT' in file:
            xsec = sample_xsec.getXsec(file)
            #print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0): normlisation=(xsec*luminosity)/(total_events)
            else: normlisation=0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                QCD_Hists.append(h_temp2)
            else:QCD_Hists.append(h_temp)

        elif 'combined_data_MET' in file:
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                MET_Hist.append(h_temp2)
            else:MET_Hist.append(h_temp)

        elif 'combined_data_SE' in file:
            if isrebin:
                h_temp2=setHistStyle(h_temp,hist)
                SE_Hist.append(h_temp2)
            else:SE_Hist.append(h_temp)

        count+=1

###==========================================================add all the histograms regional based ======================================

    for i in range(len(WJets_Hists)):
        if i==0:
            WJets=WJets_Hists[i]
        else:WJets.Add(WJets_Hists[i])
    WJets.Sumw2()

    for i in range(len(DYJets_Hits)):
        if i==0:
            DYJets=DYJets_Hits[i]
        else:DYJets.Add(DYJets_Hits[i])
    DYJets.Sumw2()

    for i in range(len(ZJets_Hits)):
        if i==0:
            ZJets=ZJets_Hits[i]
        else:ZJets.Add(ZJets_Hits[i])
    ZJets.Sumw2()


    for i in range(len(GJets_Hists)):
        if i==0:
            GJets=GJets_Hists[i]
        else:GJets.Add(GJets_Hists[i])
    GJets.Sumw2()


    for i in range(len(DIBOSON_Hists)):
        if i==0:
            DIBOSON=DIBOSON_Hists[i]
        else:DIBOSON.Add(DIBOSON_Hists[i])
    DIBOSON.Sumw2()

    for i in range(len(STop_Hists)):
        if i==0:
            STop=STop_Hists[i]
        else:STop.Add(STop_Hists[i])
    STop.Sumw2()

    for i in range(len(Top_Hists)):
        if i==0:
            Top=Top_Hists[i]
        else:Top.Add(Top_Hists[i])
    Top.Sumw2()

    for i in range(len(QCD_Hists)):
        if i==0:
            QCD=QCD_Hists[i]
        else:QCD.Add(QCD_Hists[i])
    QCD.Sumw2()


##=================================================================

    ZJetsCount    =   ZJets.Integral();
    DYJetsCount   =   DYJets.Integral();
    WJetsCount    =   WJets.Integral();
    STopCount     =   STop.Integral();
    GJetsCount    =   GJets.Integral();
    TTCount       =   Top.Integral();
    VVCount       =   DIBOSON.Integral();
    QCDCount      =   QCD.Integral();


    mcsum = ZJetsCount + DYJetsCount + WJetsCount + STopCount + GJetsCount + TTCount + VVCount + QCDCount
    total_hists = WJets_Hists + DYJets_Hits + ZJets_Hits + GJets_Hists + DIBOSON_Hists + STop_Hists + Top_Hists + QCD_Hists

    if '_cutFlow' not in str(hist):
        for histo in total_hists:
            histo = set_overflow(histo)

    ROOT.gStyle.SetHistTopMargin(0.)

#============== CANVAS DECLARATION ===================
    c12 = ROOT.TCanvas("Hist", "Hist", 0,0,1000,1000);

#==================Stack==============================
    hs = ROOT.THStack("hs"," ");

#============Colors for Histos
    DYJets.SetFillColor(ROOT.kGreen+1);
    DYJets.SetLineWidth(0);
    ZJets.SetFillColor(ROOT.kAzure-4);
    ZJets.SetLineWidth(0);
    DIBOSON.SetFillColor(ROOT.kBlue+1);
    DIBOSON.SetLineWidth(0);
    Top.SetFillColor(ROOT.kOrange-1);
    Top.SetLineWidth(0);
    WJets.SetFillColor(ROOT.kViolet-2);
    WJets.SetLineWidth(0);
    STop.SetFillColor(ROOT.kOrange+2);
    STop.SetLineWidth(0);
    GJets.SetFillColor(ROOT.kCyan-8);
    GJets.SetLineWidth(0);
    QCD.SetFillColor(ROOT.kGray+2);
    QCD.SetLineWidth(0);
    SMH.SetFillColor(ROOT.kRed-1);
    SMH.SetLineWidth(0);

#=====================Stack all the histogram =========================

    ZJetsCount    =   ZJets.Integral();
    DYJetsCount   =   DYJets.Integral();
    WJetsCount    =   WJets.Integral();
    STopCount     =   STop.Integral();
    GJetsCount    =   GJets.Integral();
    TTCount       =   Top.Integral();
    VVCount       =   DIBOSON.Integral();
    QCDCount      =   QCD.Integral();



    print('#============================= Yeid ===========================')
    print ('ZJetsCount',ZJetsCount   )
    print ('DYJetsCount',DYJetsCount)
    print ('WJetsCount',  WJetsCount)
    print ('STopCount',STopCount)
    print ('GJetsCount',GJetsCount)
    print ('TTCount',TTCount)
    print ('VVCount',VVCount)
    print ('QCDCount',QCDCount)
    print('#========================================================')


    if (QCDCount > 0):     hs.Add(QCD,"hist");
    if (DYJetsCount > 0):  hs.Add(DYJets,"hist");
    if (ZJetsCount > 0):   hs.Add(ZJets,"hist");
    if (GJetsCount > 0):   hs.Add(GJets,"hist");
    if (VVCount > 0):      hs.Add(DIBOSON,"hist");
    if (WJetsCount > 0):   hs.Add(WJets,"hist");
    if (STopCount > 0):    hs.Add(STop,"hist");
    if (TTCount > 0):      hs.Add(Top,"hist");

    hasNoEvents=False
    Stackhist = hs.GetStack().Last()
    print('hs_Integral',Stackhist.Integral())
    print('Stackhist.GetEntries()',Stackhist.GetEntries())
    maxi = Stackhist.GetMaximum()
    Stackhist.SetLineWidth(2)
    if (Stackhist.Integral()==0):
        hasNoEvents=True
        print ('No events found! for '+hist+'\n')

# =====================histogram for systematic/ statistical uncertainty ========================

    h_err = total_hists[0].Clone("h_err");
    h_err.Reset()
    for i in range(len(total_hists)):
        if i==0: continue
        else:
            if (total_hists[i].Integral()>0):
                h_err.Add(total_hists[i])
    h_err.Sumw2()
    h_err.SetFillColor(ROOT.kGray+3)
    h_err.SetLineColor(ROOT.kGray+3)
    h_err.SetMarkerSize(0)
    h_err.SetFillStyle(3013)

    if(NORATIOPLOT):
        c1_2 = ROOT.TPad("c1_2","newpad",0,0.05,1,1);   #0.993);
        c1_2.SetRightMargin(0.04);

    else:
        c1_2 =  ROOT.TPad("c1_2","newpad",0,0.28,1,1);

    c1_2.SetBottomMargin(0.09);
    c1_2.SetTopMargin(0.06);
    c1_2.SetLogy(ISLOG);
    #if(VARIABLEBINS){ c1_2->SetLogx(0);}
    c1_2.Draw();
    c1_2.cd();
    hs.Draw()
    if ('MET' in hist) and ('SR' in hist):
        sig_leg1b = ROOT.TLegend(0.23, 0.62, 0.60,0.90,'',"brNDC");
        sig_leg1b.SetTextSize(0.030);sig_leg1b.SetBorderSize(0)
        sig_leg1b.SetFillStyle(0);sig_leg1b.SetTextFont(42)
        sig_leg1b.SetHeader("2HDM+a model")
        sig_leg2b = ROOT.TLegend(0.23, 0.62, 0.60,0.90,'',"brNDC");
        sig_leg2b.SetTextSize(0.030);sig_leg2b.SetBorderSize(0)
        sig_leg2b.SetFillStyle(0);sig_leg2b.SetTextFont(42)
        sig_leg2b.SetHeader("2HDM+a model")
        mass_points = [50,250,500]
        signal_files_name = [name for name in os.listdir(sig_path) if (('Ma50' in name) or ('Ma500' in name)) ]
        signal_files = [ROOT.TFile(sig_path+'/'+filename,'READ') for filename in os.listdir(sig_path) if (('Ma50' in filename) or ('Ma500' in filename))]
        total = [fname.Get('h_total_mcweight') for fname in signal_files]
        sig_hist1b = [fname.Get('h_reg_SR_1b_MET') for fname in signal_files]
        sig_hist2b = [fname.Get('h_reg_SR_2b_MET') for fname in signal_files]
        sig_hist1b_list = [i.Scale(luminosity*sig_sample_xsec.getSigXsec(j)/k.Integral()) for i,j,k in zip(sig_hist1b,signal_files_name,total)]
        sig_hist2b_list = [i.Scale(luminosity*sig_sample_xsec.getSigXsec(j)/k.Integral()) for i,j,k in zip(sig_hist2b,signal_files_name,total)]
        LineStyle = [[i.SetLineStyle(2), j.SetLineStyle(2)] for i,j in zip(sig_hist1b,sig_hist2b)]
        LineWidth = [[i.SetLineWidth(2), j.SetLineWidth(2)] for i,j in zip(sig_hist1b,sig_hist2b)]
        LineColor = [[i.SetLineColor(n), j.SetLineColor(n)] for i,j,n in zip(sig_hist1b,sig_hist2b,range(2,len(sig_hist2b)+2))]
        MarkerColor = [[i.SetMarkerColor(n), j.SetMarkerColor(n)] for i,j,n in zip(sig_hist1b,sig_hist2b,range(2,len(sig_hist2b)+2))]
        MarkerStyle = [[i.SetMarkerStyle(n), j.SetMarkerStyle(n)] for i,j,n in zip(sig_hist1b,sig_hist2b,range(len(sig_hist2b)))]
        MarkerSize = [[i.SetMarkerSize(1.5), j.SetMarkerSize(1.5)] for i,j in zip(sig_hist1b,sig_hist2b)]
        sig_leg1b_list = [sig_leg1b.AddEntry(his_list,"ma = "+filename.split('_')[6].strip('Ma')+" GeV, mA = "+filename.split('_')[8].strip('MA')+" GeV","l") for his_list,filename in zip(sig_hist1b,signal_files_name)]
        sig_leg2b_list = [sig_leg2b.AddEntry(his_list,"ma = "+filename.split('_')[6].strip('Ma')+" GeV, mA = "+filename.split('_')[8].strip('MA')+" GeV","l") for his_list,filename in zip(sig_hist2b,signal_files_name)]
        if ('1b' in hist):
            draw_hist1b = [i.Draw("same") for i in sig_hist1b]
            sig_leg1b.Draw()
        if ('2b' in hist):
            draw_hist1b = [i.Draw("same") for i in sig_hist2b]
            sig_leg2b.Draw()
#####================================= data section =========================
    if 'SR' in reg:
        h_data=hs.GetStack().Last()
    else:
        if dtset=="SE":
            h_data=SE_Hist[0]
        elif dtset=="MET":
            h_data=MET_Hist[0]
    h_data.Sumw2()
    h_data.SetLineColor(1)
    h_data.SetLineWidth(2)
    h_data.SetMarkerSize(1.5)
    h_data.SetMarkerStyle(20)

    if(not NORATIOPLOT):
        h_data.Draw("same p e1");
    if (ISLOG):
        if '_cutFlow' in str(hist):
            hs.SetMaximum(10000000)
            hs.SetMinimum(120)
        else:
            hs.SetMaximum(maxi * 50)
            hs.SetMinimum(1)
    else:
        hs.SetMaximum(maxi * 1.35)
        hs.SetMinimum(0)
    print ('Data Integral',h_data.Integral())
##=============================== hs setting section =====================
#
    if (not hasNoEvents):
        hs.GetXaxis().SetNdivisions(508)
        if(NORATIOPLOT):
            hs.GetXaxis().SetTitleOffset(1.05)
            hs.GetXaxis().SetTitleFont(42)
            hs.GetXaxis().SetLabelFont(42)
            hs.GetXaxis().SetLabelSize(.03)
            hs.GetYaxis().SetTitle("Events")
            hs.GetYaxis().SetTitleSize(0.12)
            hs.GetYaxis().SetTitleOffset(1.5)
            hs.GetYaxis().SetTitleFont(42)
            hs.GetYaxis().SetLabelFont(42)
            hs.GetYaxis().SetLabelSize(0.05)
            hs.GetXaxis().SetTitle(str(titleX))
            hs.GetXaxis().SetTitleFont(42)
            hs.GetXaxis().SetLabelFont(42);
            hs.GetXaxis().SetLabelOffset(.01);
            hs.GetXaxis().SetLabelSize(0.03);
            hs.GetYaxis().SetTitle("Events");
            hs.GetYaxis().SetTitleSize(0.05);
            hs.GetYaxis().SetTitleFont(42);
            hs.GetYaxis().SetLabelFont(42);
            hs.GetYaxis().SetLabelSize(.03);
        else:
            #hs.GetXaxis().SetTitle(str(titleX))
            # hs.GetXaxis().SetTitleSize(0.00);
            hs.GetXaxis().SetTitleOffset(0.00);
            hs.GetXaxis().SetTitleFont(42);
            hs.GetXaxis().SetLabelFont(42);
            hs.GetXaxis().SetLabelOffset(.01);
            hs.GetXaxis().SetLabelSize(0.04);
            hs.GetYaxis().SetTitle("Events");
            hs.GetYaxis().SetTitleSize(0.045);
            hs.GetYaxis().SetTitleOffset(1);
            hs.GetYaxis().SetTitleFont(42);
            hs.GetYaxis().SetLabelFont(42);
            hs.GetYaxis().SetLabelSize(.03);

        if not isrebin: hs.GetXaxis().SetRangeUser(XMIN,XMAX);
        hs.GetXaxis().SetNdivisions(508)

#=============================  legend section =========================================
    DYLegend    =   "Z(ll) + jets "
    WLegend     =   "W(l#nu) + jets "
    GLegend     =   "G jets "
    ZLegend     =   "Z(#nu#nu) + jets "
    STLegend    =   "Single t "
    TTLegend    =   "Top "
    VVLegend    =   "DIBOSON "
    QCDLegend   =   "QCD "

    x1_l = 0.96; dx_l = 0.40
    y1_l = 0.94; dy_l = 0.21
    x0_l = x1_l-dx_l
    y0_l = y1_l-dy_l

    legend = ROOT.TLegend(x0_l,y0_l,x1_l, y1_l,"", "brNDC")
    legend.SetNColumns(2)
    legend.SetTextSize(0.0265);
    legend.SetBorderSize(0);
    legend.SetLineColor(1);
    legend.SetLineStyle(1);
    legend.SetLineWidth(1);
    legend.SetFillColor(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);

    if(not NORATIOPLOT):
        if 'SR' in reg:
            legend.AddEntry(h_data,"bkgSum","PEL")
        else:
            legend.AddEntry(h_data,"Data","PEL")
    legend.AddEntry(Top,TTLegend,"f");
    legend.AddEntry(STop,STLegend,"f");
    legend.AddEntry(WJets,WLegend,"f");
    legend.AddEntry(DIBOSON,VVLegend,"f");
    if GJetsCount > 0:legend.AddEntry(GJets,GLegend,"f");
    if ZJetsCount > 0:legend.AddEntry(ZJets,ZLegend,"f");
    legend.AddEntry(DYJets,DYLegend,"f");
    legend.AddEntry(QCD,QCDLegend,"f");

    legend.Draw('same')

#=================================================latex section =====================
    latexPreCMSname= "#bf{CMS} #it{Preliminary}"
    t2c =  ROOT.TLatex(0.10,0.97,latexPreCMSname);
    t2c.SetTextSize(0.045)

    t2a =  ROOT.TLatex(0.7,0.97,str(luminosity_)+' fb^{-1} (13TeV )');
    t2a.SetTextSize(0.040);

    t2b = ROOT.TLatex(0.22,0.88,'');
    t2b.SetTextSize(0.03);

    t2d = ROOT.TLatex(0.30,0.85,str(histolabel));
    t2d.SetTextSize(0.045);

    t2a.SetTextAlign(12);
    t2a.SetNDC(ROOT.kTRUE);
    t2a.SetTextFont(42);
    t2b.SetTextAlign(12);
    t2b.SetNDC(ROOT.kTRUE);
    t2b.SetTextFont(61);
    t2c.SetTextAlign(12);
    t2c.SetNDC(ROOT.kTRUE);
    t2c.SetTextFont(42);
    t2d.SetTextAlign(12);
    t2d.SetNDC(ROOT.kTRUE);
    t2d.SetTextFont(42);
    t2a.Draw("same");
    t2b.Draw("same");
    t2c.Draw("same");
    t2d.Draw("same");
#======================================== ratio log ================

    ratioleg =  ROOT.TLegend(0.6, 0.88, 0.89, 0.98);
    #//ratioleg->SetFillColor(0);
    ratioleg.SetLineColor(0);
    ratioleg.SetShadowColor(0);
    ratioleg.SetTextFont(42);
    ratioleg.SetTextSize(0.09);
    ratioleg.SetBorderSize(1);
    ratioleg.SetNColumns(2);

#============================================= statistical error section ======================

    ratiostaterr = h_err.Clone("ratiostaterr")
    ratiostaterr.Sumw2()
    ratiostaterr.SetStats(0);
    ratiostaterr.SetMinimum(0);
    ratiostaterr.SetMarkerSize(0);
    ratiostaterr.SetFillColor(ROOT.kBlack);
    ratiostaterr.SetFillStyle(3013);
    for i in range(h_err.GetNbinsX()+2):
        ratiostaterr.SetBinContent(i, 0.0)

        if (h_err.GetBinContent(i) > 1e-6 ):
            binerror = h_err.GetBinError(i)/h_err.GetBinContent(i)
            ratiostaterr.SetBinError(i, binerror)
        else:ratiostaterr.SetBinError(i, 999.)

    ratioleg.AddEntry(ratiostaterr, "stat", "f")

 #============================================= Lower Tpad Decalaration ====================================
    if(not NORATIOPLOT):
        c12.cd()
        DataMC    = h_data.Clone()
        DataMC.Add(Stackhist,-1)   # remove for data/mc
        DataMCPre = h_data.Clone();
        DataMC.Divide(Stackhist);
        DataMC.GetYaxis().SetTitle("#frac{Data-Pred}{Pred}.");
        DataMC.GetYaxis().SetTitleSize(0.1);
        DataMC.GetYaxis().SetTitleOffset(0.42);
        DataMC.GetYaxis().SetTitleFont(42);
        DataMC.GetYaxis().SetLabelSize(0.08);
        DataMC.GetYaxis().CenterTitle();
        DataMC.GetXaxis().SetTitle(str(titleX))
        DataMC.GetXaxis().SetLabelSize(0.1);
        DataMC.GetXaxis().SetTitleSize(0.1);
        DataMC.GetXaxis().SetTitleOffset(1);
        DataMC.GetXaxis().SetTitleFont(42);
        DataMC.GetXaxis().SetTickLength(0.07);
        DataMC.GetXaxis().SetLabelFont(42);
        DataMC.GetYaxis().SetLabelFont(42);


    c1_1 = ROOT.TPad("c1_1", "newpad",0,0.00,1,0.3);
    if (not NORATIOPLOT): c1_1.Draw();
    c1_1.cd();
    c1_1.Range(-7.862408,-629.6193,53.07125,486.5489);
    c1_1.SetFillColor(0);
    c1_1.SetTicky(1);
    c1_1.SetLeftMargin(0.1);
    c1_1.SetRightMargin(0.1);
    c1_1.SetTopMargin(0.0);
    c1_1.SetBottomMargin(0.32);
    c1_1.SetFrameFillStyle(0);
    c1_1.SetFrameBorderMode(0);
    c1_1.SetFrameFillStyle(0);
    c1_1.SetFrameBorderMode(0);
    c1_1.SetLogy(0);

    if(not NORATIOPLOT):
        if (0): # if(VARIABLEBINS)
            c1_1.SetLogx(0)
            DataMC.GetXaxis().SetMoreLogLabels()
            DataMC.GetXaxis().SetNoExponent()
            DataMC.GetXaxis().SetNdivisions(508)
        if not isrebin: DataMC.GetXaxis().SetRangeUser(XMIN,XMAX)
        DataMC.SetMarkerSize(1.5)
        DataMC.SetMarkerStyle(20)
        DataMC.SetMarkerColor(1)
        DataMC.SetMinimum(-1.08)
        DataMC.SetMaximum(1.08)
        DataMC.GetXaxis().SetNdivisions(508)
        DataMC.GetYaxis().SetNdivisions(505)
        DataMC.Draw("P e1")
        ratiostaterr.Draw("e2 same")
        DataMC.Draw("P e1 same")
        line1=  ROOT.TLine(XMIN,0.2,XMAX,0.2)
        line2=  ROOT.TLine(XMIN,-0.2,XMAX,-0.2)
        line1.SetLineStyle(2)
        line1.SetLineColor(2)
        line1.SetLineWidth(2)
        line2.SetLineStyle(2)
        line2.SetLineColor(2)
        line2.SetLineWidth(2)
        line1.Draw("same")
        line2.Draw("same")
        ratioleg.Draw("same")
    c12.Draw()

    plot=str(hist)

    if not os.path.exists('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPng/'+reg):
        os.makedirs('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPng/'+reg)
    if not os.path.exists('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPdf/'+reg):
        os.makedirs('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPdf/'+reg)
    if not os.path.exists('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMRoot/'+reg):
        os.makedirs('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMRoot/'+reg)
    if (ISLOG == 0):
        c12.SaveAs('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPdf/'+reg+'/'+plot+'.pdf')
        c12.SaveAs('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPng/'+reg+'/'+plot+'.png')
        print("Saved. \n")
    if (ISLOG == 1):
        c12.SaveAs('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPdf/'+reg+'/'+plot+'_log.pdf')
        c12.SaveAs('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPng/'+reg+'/'+plot+'_log.png')
        print("Saved. \n")

    fshape = ROOT.TFile('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMRoot/'+reg+'/'+plot+'.root', "RECREATE");

    fshape.cd()
    Stackhist.SetNameTitle("bkgSum","bkgSum")
    Stackhist.Write()
    DIBOSON.SetNameTitle("DIBOSON","DIBOSON");
    DIBOSON.Write()
    ZJets.SetNameTitle("ZJets","ZJets");
    ZJets.Write()
    GJets.SetNameTitle("GJets","GJets");
    GJets.Write()
    QCD.SetNameTitle("QCD","QCD");
    QCD.Write()
    SMH.SetNameTitle("SMH","SMH");
    SMH.Write();
    STop.SetNameTitle("STop","STop");
    STop.Write();
    Top.SetNameTitle("Top","Top");
    Top.Write();
    WJets.SetNameTitle("WJets","WJets");
    WJets.Write();
    DYJets.SetNameTitle("DYJets","DYJets");
    DYJets.Write();
    data_obs=h_data
    data_obs.SetNameTitle("data_obs","data_obs");
    data_obs.Write();
    fshape.Write();
    fshape.Close();
 #=======================================================================


######################################################################

regions=[]
PUreg=[]
if makeMuCRplots:
    regions+=['SR_1b','SR_2b','ZmumuCR_1b','ZmumuCR_2b','TopmunuCR_1b','WmunuCR_1b','TopmunuCR_2b','WmunuCR_2b']
    PUreg+=['mu_']
if makeEleCRplots:
    regions+=['ZeeCR_1b','ZeeCR_2b','TopenuCR_1b','WenuCR_1b','TopenuCR_2b','WenuCR_2b']
'''
    PUreg+=['ele_']
if makePhoCRplots:
    regions+=['1gamma2b']
    PUreg+=['pho_']
if makeQCDCRplots:
    regions+=['QCD2b']
    PUreg+=[]
'''
# makeplot("reg_WenuCR_1b_cutFlow",'h_reg_WenuCR_1b_cutFlow','Cutflow',0,6,1,1,0,'WenuCR_1b',varBin=False)
# makeplot("reg_ZeeCR_1b_pu_nPUVert",'h_reg_ZeeCR_1b_pu_nPUVert','pu_nPUVert',0,70,1,0,0,'ZeeCR_1b',varBin=False)
# makeplot("reg_ZeeCR_1b_nopu_nPUVert",'h_reg_ZeeCR_1b_nopu_nPUVert','nopu_nPUVert',0,70,1,0,0,'ZeeCR_1b',varBin=False)
# makeplot("reg_TopmunuCR_1b_cutFlow",'h_reg_TopmunuCR_1b_cutFlow','Cutflow',0,6,1,1,0,'TopmunuCR_1b',varBin=False)

for reg in regions:
    try:
        if 'SR_' in reg:
            makeplot("reg_"+reg+"_MET",'h_reg_'+reg+'_MET','Real MET (GeV)',200.,1000.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_min_dPhi",'h_reg_'+reg+'_min_dPhi','min_dPhi',0,4,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet1Pt",'h_reg_'+reg+'_Jet1Pt','JET1 p_{T} (GeV)',30.,800.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet1Eta",'h_reg_'+reg+'_Jet1Eta','JET1 #eta',-2.5,2.5,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet1Phi",'h_reg_'+reg+'_Jet1Phi','JET1 #phi',-3.14,3.14,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet1deepCSV",'h_reg_'+reg+'_Jet1deepCSV','JET1 deepCSV',0,1.2,1,0,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet2Pt",'h_reg_'+reg+'_Jet2Pt','JET2 p_{T} (GeV)',30.,800.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet2Eta",'h_reg_'+reg+'_Jet2Eta','JET2 #eta',-2.5,2.5,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet2Phi",'h_reg_'+reg+'_Jet2Phi','JET2 #phi',-3.14,3.14,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet2deepCSV",'h_reg_'+reg+'_Jet2deepCSV','JET2 deepCSV',0,1.2,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_nJets",'h_reg_'+reg+'_nJets','nJets',0.,10.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_cutFlow",'h_reg_'+reg+'_cutFlow','CutFlow',0,7,1,1,0,reg,varBin=False)
        else:
            makeplot("reg_"+reg+"_MET",'h_reg_'+reg+'_MET','Real MET (GeV)',0.,700.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Recoil",'h_reg_'+reg+'_Recoil','Hadronic Recoil (GeV)',200.,1000.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_min_dPhi",'h_reg_'+reg+'_min_dPhi','min_dPhi',0,4,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet1Pt",'h_reg_'+reg+'_Jet1Pt','JET1 p_{T} (GeV)',30.,800.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet1Eta",'h_reg_'+reg+'_Jet1Eta','JET1 #eta',-2.5,2.5,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet1Phi",'h_reg_'+reg+'_Jet1Phi','JET1 #phi',-3.14,3.14,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet1deepCSV",'h_reg_'+reg+'_Jet1deepCSV','JET1 deepCSV',0,1.2,1,1,0,reg,varBin=False)
            if 'W' in reg :
                makeplot("reg_"+reg+"_Wmass",'h_reg_'+reg+'_Wmass','W candidate mass (GeV)',0.,165.,1,1,0,reg,varBin=False)
                makeplot("reg_"+reg+"_WpT",'h_reg_'+reg+'_WpT','W candidate p_{T} (GeV)',0.,700.,1,1,0,reg,varBin=False)
            if 'Z' in reg:
                makeplot("reg_"+reg+"_Zmass",'h_reg_'+reg+'_Zmass','Z candidate mass (GeV)',70.,110.,1,0,0,reg,varBin=False)
                makeplot("reg_"+reg+"_ZpT",'h_reg_'+reg+'_ZpT','Z candidate p_{T} (GeV)',0.,700.,1,1,0,reg,varBin=False)
                makeplot("reg_"+reg+"_lep2_pT",'h_reg_'+reg+'_lep2_pT','lepton2 p_{T}',0,500,1,1,0,reg,varBin=False)

            makeplot("reg_"+reg+"_lep1_pT",'h_reg_'+reg+'_lep1_pT','lepton1 p_{T}',0,500,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_cutFlow",'h_reg_'+reg+'_cutFlow','CutFlow',0,9,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet2Pt",'h_reg_'+reg+'_Jet2Pt','JET2 p_{T} (GeV)',30.,800.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet2Eta",'h_reg_'+reg+'_Jet2Eta','JET2 #eta',-2.5,2.5,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet2Phi",'h_reg_'+reg+'_Jet2Phi','JET2 #phi',-3.14,3.14,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_Jet2deepCSV",'h_reg_'+reg+'_Jet2deepCSV','JET2 deepCSV',0,1.2,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_nJets",'h_reg_'+reg+'_nJets','nJets',0.,10.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_nPV",'h_reg_'+reg+'_nPV','After PU reweighting',0.,70.,1,1,0,reg,varBin=False)
            makeplot("reg_"+reg+"_PUnPV",'h_reg_'+reg+'_PUnPV','After PU reweighting',0.,70.,1,1,0,reg,varBin=False)
    except Exception as e:
        print (e)
        print ("Cannot Plot")
        pass
