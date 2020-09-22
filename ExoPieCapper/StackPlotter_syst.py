import os
import sys
import datetime
import sys
import optparse
import ROOT as ROOT
import array
import string

#command  python StackPlotter_syst.py -d <DATASET_NAME> -m -y <Year> -D <histo_DIR> -S <signal_Dir>
usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)

parser.add_option("-d", "--data", dest="datasetname")
parser.add_option("-D", "--pDir", type="string",
                  dest="rootFileDir", help="directory containing histogram")
parser.add_option("-S", "--sigDir", type="string", dest="SIGrootFileDir",
                  help="directory containing signal histogram")
parser.add_option("-s", "--sr", action="store_true", dest="plotSIGNAL")
parser.add_option("-m", "--mu", action="store_true", dest="plotMuChannels")
parser.add_option("-e", "--ele", action="store_true", dest="plotEleChannels")
parser.add_option("-p", "--pho", action="store_true", dest="plotPhoChannels")
parser.add_option("-q", "--qcd", action="store_true", dest="plotQCDChannels")
parser.add_option("-y", "--year", dest="year", default="Year")

(options, args) = parser.parse_args()

if options.plotSIGNAL == None:
    makeSIGplots = False
else:
    makeSIGplots = options.plotSIGNAL

if options.plotMuChannels == None:
    makeMuCHplots = False
else:
    makeMuCHplots = options.plotMuChannels

if options.plotEleChannels == None:
    makeEleCHplots = False
else:
    makeEleCHplots = options.plotEleChannels

if options.plotPhoChannels == None:
    makePhoCRplots = False
else:
    makePhoCRplots = options.plotPhoChannels

if options.plotQCDChannels == None:
    makeQCDCRplots = False
else:
    makeQCDCRplots = options.plotQCDChannels

if options.datasetname.upper() == "SE":
    dtset = "SE"
elif options.datasetname.upper() == "SP":
    dtset = "SP"
elif options.datasetname.upper() == "SM":
    dtset = "SM"
else:
    dtset = "MET"

print("Using dataset "+dtset)

runOn2016 = False
runOn2017 = False
runOn2018 = False
if options.year == '2016':
    runOn2016 = True
elif options.year == '2017':
    runOn2017 = True
elif options.year == '2018':
    runOn2018 = True
else:
    print('Please provide on which year you want to run?')

if runOn2016:
    #import sample_xsec_2016 as sample_xsec
    import sample_xsec_2016_GenXSecAnalyser as sample_xsec
    import sig_sample_xsec_2016 as sig_sample_xsec
    luminosity = 35.82 * 1000
    luminosity_ = '{0:.2f}'.format(35.82)
elif runOn2017:
    #import sample_xsec_2017 as sample_xsec
    import sample_xsec_2017_GenXSecAnalyser as sample_xsec
    import sig_sample_xsec_2017 as sig_sample_xsec
    luminosity = 41.5 * 1000
    luminosity_ = '{0:.2f}'.format(41.50)
elif runOn2018:
    import sample_xsec_2018 as sample_xsec
    import sig_sample_xsec_2018 as sig_sample_xsec
    luminosity = 59.64 * 1000
    luminosity_ = '{0:.2f}'.format(59.64)


datestr = str(datetime.date.today().strftime("%d%m%Y"))

if options.rootFileDir == None:
    print('Please provide histogram directory name')
    sys.exit()
else:
    path = options.rootFileDir

sig_path = options.SIGrootFileDir

print("sig_path sig_path sig_path", sig_path)
if makeMuCHplots:
    yield_outfile_binwise = open(
        'bbDM'+str(options.year)+'_Mu_binwise_yield.txt', 'w')
    yield_outfile = open('bbDM'+str(options.year)+'_Mu_yield.txt', 'w')
if makeEleCHplots:
    yield_outfile_binwise = open(
        'bbDM'+str(options.year) + '_Ele_binwise_yield.txt', 'w')
    yield_outfile = open('bbDM'+str(options.year)+'_Ele_yield.txt', 'w')

alpha_list = list(string.ascii_uppercase)


def set_overflow(hist):
    bin_num = hist.GetXaxis().GetNbins()
    #print (bin_num)
    hist.SetBinContent(bin_num, hist.GetBinContent(
        bin_num+1)+hist.GetBinContent(bin_num))  # Add overflow bin content to last bin
    hist.SetBinContent(bin_num+1, 0.)
    return hist

def setHistStyle(h_temp2, hist, rebin):
    dovarbin = False
    bins = h_temp2.GetNbinsX()
    if rebin > 1:
        if bins%rebin == 0:
            h_temp_ = h_temp2.Rebin(rebin)
        elif bins%(rebin+1) == 0:
            h_temp_ = h_temp2.Rebin(rebin+1)
        elif bins%(rebin-1) == 0:
            h_temp_ = h_temp2.Rebin(rebin-1)
    else:
        h_temp_ = h_temp2
    return h_temp_


def SetCMSAxis(h, xoffset=1., yoffset=1.):
    h.GetXaxis().SetTitleSize(0.047)
    h.GetYaxis().SetTitleSize(0.047)

    print(type(h))
    if type(h) is ((not ROOT.TGraphAsymmErrors) or (not ROOT.TGraph)):
        h.GetZaxis().SetTitleSize(0.047)

    h.GetXaxis().SetLabelSize(0.047)
    h.GetYaxis().SetLabelSize(0.047)
    if type(h) is ((not ROOT.TGraphAsymmErrors) or (not ROOT.TGraph)):
        h.GetZaxis().SetLabelSize(0.047)

    h.GetXaxis().SetTitleOffset(xoffset)
    h.GetYaxis().SetTitleOffset(yoffset)
    return h


def ExtraText(text_, x_, y_):
    if not text_:
        print("nothing provided as text to ExtraText, function crashing")
    ltx = ROOT.TLatex(x_, y_, text_)

    if len(text_) > 0:
        ltx.SetTextFont(42)
        ltx.SetTextSize(0.049)
        #ltx.Draw(x_,y_,text_)
        ltx.Draw('same')
    return ltx


def myCanvas1D():
    c = ROOT.TCanvas("myCanvasName", "The Canvas Title", 650, 600)
    c.SetBottomMargin(0.050)
    c.SetRightMargin(0.050)
    c.SetLeftMargin(0.050)
    c.SetTopMargin(0.050)
    return c

#def SetLegend(coordinate_=[.38,.7,.890,.87],ncol=2):


def SetLegend(coordinate_=[.50, .65, .90, .90], ncol=2):
    c_ = coordinate_
    legend = ROOT.TLegend(c_[0], c_[1], c_[2], c_[3])
    legend.SetBorderSize(0)
    legend.SetNColumns(ncol)
    legend.SetLineColor(1)
    legend.SetLineStyle(1)
    legend.SetLineWidth(1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.045)
    return legend


def drawenergy1D(is2017, text_="Work in progress 2018", data=True):
    #pt = ROOT.TPaveText(0.0877181,0.9,0.9580537,0.96,"brNDC")
    pt = ROOT.TPaveText(0.0877181, 0.95, 0.9580537, 0.96, "brNDC")
    pt.SetBorderSize(0)
    pt.SetTextAlign(12)
    pt.SetFillStyle(0)
    pt.SetTextFont(52)

    cmstextSize = 0.07
    preliminarytextfize = cmstextSize * 0.7
    lumitextsize = cmstextSize * 0.7
    pt.SetTextSize(cmstextSize)
    text = pt.AddText(0.03, 0.57, "#font[61]{CMS}")

    #pt1 = ROOT.TPaveText(0.0877181,0.9,0.9580537,0.96,"brNDC")
    pt1 = ROOT.TPaveText(0.0877181, 0.95, 0.9580537, 0.96, "brNDC")
    pt1.SetBorderSize(0)
    pt1.SetTextAlign(12)
    pt1.SetFillStyle(0)
    pt1.SetTextFont(52)

    pt1.SetTextSize(preliminarytextfize)
    #text1 = pt1.AddText(0.215,0.4,text_)
    text1 = pt1.AddText(0.15, 0.4, text_)

    #pt2 = ROOT.TPaveText(0.0877181,0.9,0.9580537,0.96,"brNDC")
    pt2 = ROOT.TPaveText(0.0877181, 0.95, 0.9580537, 0.96, "brNDC")
    pt2.SetBorderSize(0)
    pt2.SetTextAlign(12)
    pt2.SetFillStyle(0)
    pt2.SetTextFont(52)
    pt2.SetTextFont(42)
    pt2.SetTextSize(lumitextsize)

    pavetext = ''
    if is2017 and data:
        pavetext = str(luminosity_)+' fb^{-1}'+" (13 TeV)"
    if (not is2017) and data:
        pavetext = str(luminosity_)+' fb^{-1}'+"(13 TeV)"

    if is2017 and not data:
        pavetext = "13 TeV"
    if (not is2017) and not data:
        pavetext = "13 TeV"

    if data:
        text3 = pt2.AddText(0.68, 0.5, pavetext)
    if not data:
        text3 = pt2.AddText(0.68, 0.5, pavetext)

    return [pt, pt1, pt2]


def makeplot(loc, hist, titleX, XMIN, XMAX, Rebin, ISLOG, NORATIOPLOT, reg, varBin, row=2):
    # try:

    print('plotting histogram:   ', hist)
    isrebin = True  # bool(varBin)
    if runOn2016:
        files = open("samplelist_2016.txt", "r")
    elif runOn2017:
        files = open("samplelist_2017.txt", "r")
    elif runOn2018:
        files = open("samplelist_2018.txt", "r")

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetFrameLineWidth(3)
    #gStyle->SetErrorX(0)
    ROOT.gStyle.SetLineWidth(1)

    if 'preselR' in hist:
        histolabel = "Pre Selection"
    elif '_SR_1b' in hist:
        histolabel = "SR(1b)"
    elif '_SR_2b' in hist:
        histolabel = "SR(2b)"
    elif 'ZmumuCR_1b' in hist:
        histolabel = "Z(#mu#mu)+1b CR"
    elif 'ZeeCR_1b' in hist:
        histolabel = "Z(ee)+1b CR"
    elif 'WmunuCR_1b' in hist:
        histolabel = "W(#mu#nu)+1b CR"
    elif 'WenuCR_1b' in hist:
        histolabel = "W(e#nu)+1b CR"
    elif 'TopmunuCR_1b' in hist:
        histolabel = "t#bar{t}(#mu#nu)+1b CR"
    elif 'TopenuCR_1b' in hist:
        histolabel = "t#bar{t}(e#nu)+1b CR"
    elif 'ZmumuCR_2b' in hist:
        histolabel = "Z(#mu#mu)+2b CR"
    elif 'ZeeCR_2b' in hist:
        histolabel = "Z(ee)+2b CR"
    elif 'WmunuCR_2b' in hist:
        histolabel = "W(#mu#nu)+2b CR"
    elif 'WenuCR_2b' in hist:
        histolabel = "W(e#nu)+2b CR"
    elif 'TopmunuCR_2b' in hist:
        histolabel = "t#bar{t}(#mu#nu)+2b CR"
    elif 'TopenuCR_2b' in hist:
        histolabel = "t#bar{t}(e#nu)+2b CR"
    else:
        histolabel = "testing"

    xsec = 1.0
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

    DYJets_Hits = []
    ZJets_Hits = []
    WJets_Hists = []
    GJets_Hists = []
    DIBOSON_Hists = []
    STop_Hists = []
    Top_Hists = []
    QCD_Hists = []
    SMH_Hists = []
    MET_Hist = []
    SE_Hist = []

    count = 0
    for file in files.readlines()[:]:
        myFile = path+'/'+file.rstrip()
        Str = str(count)
        exec("f"+Str+"=ROOT.TFile(myFile,'READ')", locals(), globals())
        exec("h_temp=f"+Str+".Get("+"\'"+str(hist)+"\'"+")", locals(), globals())
        exec("h_total_weight=f"+Str+".Get('h_total_mcweight')", locals(), globals())
        total_events = h_total_weight.Integral()

        if 'WJetsToLNu_HT' in file:
            xsec = sample_xsec.getXsec(file)
            # print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0):
                normlisation = (xsec*luminosity)/(total_events)
            else:
                normlisation = 0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                WJets_Hists.append(h_temp2)
            else:
                WJets_Hists.append(h_temp)

        elif 'DYJetsToLL_M-50' in file:
            xsec = sample_xsec.getXsec(file)
            # print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0):
                normlisation = (xsec*luminosity)/(total_events)
            else:
                normlisation = 0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                DYJets_Hits.append(h_temp2)
            else:
                DYJets_Hits.append(h_temp)

        elif 'ZJetsToNuNu' in file:
            xsec = sample_xsec.getXsec(file)
            # print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0):
                normlisation = (xsec*luminosity)/(total_events)
            else:
                normlisation = 0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                ZJets_Hits.append(h_temp2)
            else:
                ZJets_Hits.append(h_temp)

        elif 'GJets_HT' in file:
            xsec = sample_xsec.getXsec(file)
            # print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0):
                normlisation = (xsec*xsec)/(total_events)
            else:
                normlisation = 0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                GJets_Hists.append(h_temp2)
            else:
                GJets_Hists.append(h_temp)

        elif ('WWTo' in file) or ('WZTo' in file) or ('ZZTo' in file) or ('WW_' in file) or ('ZZ_' in file) or ('WZ_' in file):
            xsec = sample_xsec.getXsec(file)
            # print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0):
                normlisation = (xsec*luminosity)/(total_events)
            else:
                normlisation = 0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                DIBOSON_Hists.append(h_temp2)
            else:
                DIBOSON_Hists.append(h_temp)

        elif ('ST_t' in file) or ('ST_s' in file):
            xsec = sample_xsec.getXsec(file)
            # print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0):
                normlisation = (xsec*luminosity)/(total_events)
            else:
                normlisation = 0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                STop_Hists.append(h_temp2)
            else:
                STop_Hists.append(h_temp)

        elif ('TTTo' in file) or ('TT_TuneCUETP8M2T4' in file):
            xsec = sample_xsec.getXsec(file)
            # print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0):
                normlisation = (xsec*luminosity)/(total_events)
            else:
                normlisation = 0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                Top_Hists.append(h_temp2)

            else:
                Top_Hists.append(h_temp)

        elif 'QCD_HT' in file:
            xsec = sample_xsec.getXsec(file)
            # print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0):
                normlisation = (xsec*luminosity)/(total_events)
            else:
                normlisation = 0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                QCD_Hists.append(h_temp2)
            else:
                QCD_Hists.append(h_temp)

        elif 'HToBB' in file:
            xsec = sample_xsec.getXsec(file)
            # print ('file', file ,'xsec', xsec,'\n')
            if (total_events > 0):
                normlisation = (xsec*luminosity)/(total_events)
            else:
                normlisation = 0
            h_temp.Scale(normlisation)
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                SMH_Hists.append(h_temp2)
            else:
                SMH_Hists.append(h_temp)

        elif 'combined_data_MET' in file:
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                MET_Hist.append(h_temp2)
            else:
                MET_Hist.append(h_temp)

        elif 'combined_data_SE' in file:
            if isrebin:
                h_temp2 = setHistStyle(h_temp, hist, Rebin)
                SE_Hist.append(h_temp2)
            else:
                SE_Hist.append(h_temp)

        count += 1

###==========================================================add all the histograms regional based ======================================

    for i in range(len(WJets_Hists)):
        if i == 0:
            WJets = WJets_Hists[i]
        else:
            WJets.Add(WJets_Hists[i])
    WJets.Sumw2()

    for i in range(len(DYJets_Hits)):
        if i == 0:
            DYJets = DYJets_Hits[i]
        else:
            DYJets.Add(DYJets_Hits[i])
    DYJets.Sumw2()

    for i in range(len(ZJets_Hits)):
        if i == 0:
            ZJets = ZJets_Hits[i]
        else:
            ZJets.Add(ZJets_Hits[i])
    ZJets.Sumw2()

    for i in range(len(GJets_Hists)):
        if i == 0:
            GJets = GJets_Hists[i]
        else:
            GJets.Add(GJets_Hists[i])
    GJets.Sumw2()

    for i in range(len(DIBOSON_Hists)):
        if i == 0:
            DIBOSON = DIBOSON_Hists[i]
        else:
            DIBOSON.Add(DIBOSON_Hists[i])
    DIBOSON.Sumw2()

    for i in range(len(STop_Hists)):
        if i == 0:
            STop = STop_Hists[i]
        else:
            STop.Add(STop_Hists[i])
    STop.Sumw2()

    for i in range(len(Top_Hists)):
        if i == 0:
            Top = Top_Hists[i]
        else:
            Top.Add(Top_Hists[i])
    Top.Sumw2()

    for i in range(len(QCD_Hists)):
        if i == 0:
            QCD = QCD_Hists[i]
        else:
            QCD.Add(QCD_Hists[i])
    QCD.Sumw2()

    for i in range(len(SMH_Hists)):
        if i == 0:
            SMH = SMH_Hists[i]
        else:
            SMH.Add(SMH_Hists[i])
    SMH.Sumw2()

##=================================================================

    ZJetsCount = ZJets.Integral()
    DYJetsCount = DYJets.Integral()
    WJetsCount = WJets.Integral()
    STopCount = STop.Integral()
    GJetsCount = GJets.Integral()
    TTCount = Top.Integral()
    VVCount = DIBOSON.Integral()
    QCDCount = QCD.Integral()
    SMHCount = SMH.Integral()

    mcsum = ZJetsCount + DYJetsCount + WJetsCount + \
        STopCount + GJetsCount + TTCount + VVCount + QCDCount
    total_hists = WJets_Hists + DYJets_Hits + ZJets_Hits + \
        GJets_Hists + DIBOSON_Hists + STop_Hists + Top_Hists + QCD_Hists

    if '_cutFlow' not in str(hist):
        for histo in total_hists:
            histo = set_overflow(histo)

    ROOT.gStyle.SetHistTopMargin(0.1)

#============== CANVAS DECLARATION ===================
    #c12 = ROOT.TCanvas("Hist", "Hist", 0,0,1000,1000)
    c12 = myCanvas1D()

#==================Stack==============================
    hs = ROOT.THStack("hs", " ")

#============Colors for Histos
    DYJets.SetFillColor(ROOT.kGreen+1)
    DYJets.SetLineWidth(0)
    ZJets.SetFillColor(ROOT.kAzure-4)
    ZJets.SetLineWidth(0)
    DIBOSON.SetFillColor(ROOT.kBlue+1)
    DIBOSON.SetLineWidth(0)
    Top.SetFillColor(ROOT.kOrange-1)
    Top.SetLineWidth(0)
    WJets.SetFillColor(ROOT.kViolet-2)
    WJets.SetLineWidth(0)
    STop.SetFillColor(ROOT.kOrange+2)
    STop.SetLineWidth(0)
    GJets.SetFillColor(ROOT.kCyan-8)
    GJets.SetLineWidth(0)
    QCD.SetFillColor(ROOT.kGray+2)
    QCD.SetLineWidth(0)
    SMH.SetFillColor(ROOT.kRed-1)
    SMH.SetLineWidth(0)

#=====================Stack all the histogram =========================

    ZJetsCount = ZJets.Integral()
    DYJetsCount = DYJets.Integral()
    WJetsCount = WJets.Integral()
    STopCount = STop.Integral()
    GJetsCount = GJets.Integral()
    TTCount = Top.Integral()
    VVCount = DIBOSON.Integral()
    QCDCount = QCD.Integral()
    SMHCount = SMH.Integral()

    if (SMHCount > 0):
        hs.Add(SMH, "hist")
    if (QCDCount > 0):
        hs.Add(QCD, "hist")
    if (DYJetsCount > 0):
        hs.Add(DYJets, "hist")
    if (ZJetsCount > 0):
        hs.Add(ZJets, "hist")
    if (GJetsCount > 0):
        hs.Add(GJets, "hist")
    if (VVCount > 0):
        hs.Add(DIBOSON, "hist")
    if (WJetsCount > 0):
        hs.Add(WJets, "hist")
    if (STopCount > 0):
        hs.Add(STop, "hist")
    if (TTCount > 0):
        hs.Add(Top, "hist")

    hasNoEvents = False
    Stackhist = hs.GetStack().Last()
    maxi = Stackhist.GetMaximum()
    Stackhist.SetLineWidth(2)
    if (Stackhist.Integral() == 0):
        hasNoEvents = True
        print('No events found! for '+hist+'\n')

# =====================histogram for systematic/ statistical uncertainty ========================

    h_err = total_hists[0].Clone("h_err")
    h_err.Reset()
    for i in range(len(total_hists)):
        if i == 0:
            continue
        else:
            if (total_hists[i].Integral() > 0):
                h_err.Add(total_hists[i])
    h_err.Sumw2()
    h_err.SetFillColor(ROOT.kGray+3)
    h_err.SetLineColor(ROOT.kGray+3)
    h_err.SetMarkerSize(0)
    h_err.SetFillStyle(3013)

    if(NORATIOPLOT):
        c1_2 = ROOT.TPad("c1_2", "newpad", 0, 0.05, 1, 1)  # 0.993)
        c1_2.SetRightMargin(0.06)
    else:
        c1_2 = ROOT.TPad("c1_2", "newpad", 0, 0.20, 1, 1)

    c1_2.SetBottomMargin(0.09)
    c1_2.SetTopMargin(0.08)
    c1_2.SetLeftMargin(0.12)
    c1_2.SetRightMargin(0.06)
    c1_2.SetLogy(ISLOG)
    c1_2.Draw()
    c1_2.cd()
    for h in hs:
        h = SetCMSAxis(h)
    hs.Draw()
    noYieldHisto = bool(('weight' in hist) or ('_up' in hist)
                        or ('_down' in hist) or ('_Recoil' in hist))
    if makeSIGplots:
        if ('MET' in hist) and ('SR' in hist) and not noYieldHisto:
            ma_points = [50, 200]  # how many signal points you want to include
            sig_leg = SetLegend([.50, .38, .60, .58], ncol=1)
            sig_leg.SetHeader("2HDM+a model")
            signal_files_name = [name for name in os.listdir(
                sig_path) for mapoint in ma_points if 'Ma'+str(mapoint)+'_' in name]

            signal_files = [ROOT.TFile(sig_path+'/'+name, 'READ') for name in os.listdir(
                sig_path) for mapoint in ma_points if 'Ma'+str(mapoint)+'_' in name]
            total = [fname.Get('h_total_mcweight') for fname in signal_files]
            sig_hist = [fname.Get(hist) for fname in signal_files]

            sig_hist_list = [i.Scale(luminosity*sig_sample_xsec.getSigXsec(j)/k.Integral())
                             for i, j, k in zip(sig_hist, signal_files_name, total)]

            LineStyling = [(i.SetLineStyle(n), i.SetLineWidth(6), i.SetLineColor(
                n)) for i, n in zip(sig_hist, range(2, len(sig_hist)+2))]

            MarkerStyling = [(i.SetMarkerColor(n), i.SetMarkerStyle(n), i.SetMarkerSize(
                1.5)) for i, n in zip(sig_hist, range(2, len(sig_hist)+2))]

            sig_leg_list = [sig_leg.AddEntry(his_list, "ma = "+filename.split('_')[6].strip('Ma')+" GeV, mA = "+filename.split(
                '_')[8].strip('MA')+" GeV", "lp") for his_list, filename in zip(sig_hist, signal_files_name)]

            draw_hist = [i.Draw(" same Ehist") for i in sig_hist]
            sig_leg.Draw('same')
#####================================= data section =========================
    if 'SR' in reg:
        h_data = hs.GetStack().Last()
    else:
        if dtset == "SE":
            h_data = SE_Hist[0]
        elif dtset == "MET":
            h_data = MET_Hist[0]
    h_data.Sumw2()
    h_data.SetLineColor(1)
    h_data.SetLineWidth(2)
    h_data.SetMarkerSize(1.3)
    h_data.SetMarkerStyle(20)
    h_data = SetCMSAxis(h_data)
    if '_cutFlow' not in str(hist):
        h_data = set_overflow(h_data)
    if(not NORATIOPLOT):
        h_data.Draw("same p e1")
    if (ISLOG):
        if '_cutFlow' in str(hist):
            hs.SetMaximum(1000000000)
            hs.SetMinimum(100)
        else:
            hs.SetMaximum(maxi * 50)
            hs.SetMinimum(1)
    else:
        hs.SetMaximum(maxi * 1.35)
        hs.SetMinimum(0)
    #print ('Data Integral',h_data.Integral())
##=============================== hs setting section =====================
#
    if (not hasNoEvents):
        hs.GetXaxis().SetNdivisions(508)
        if(NORATIOPLOT):
            hs.GetXaxis().SetTitleOffset(1.05)
            hs.GetXaxis().SetTitleFont(42)
            hs.GetXaxis().SetLabelFont(42)
            hs.GetXaxis().SetLabelSize(.03)
            hs.GetXaxis().SetTitle(str(titleX))
            hs.GetXaxis().SetTitleFont(42)
            hs.GetXaxis().SetLabelOffset(.01)
            hs.GetYaxis().SetTitleOffset(0.7)
            hs.GetYaxis().SetTitle("Events/bin")
            hs.GetYaxis().SetTitleSize(0.08)
            hs.GetYaxis().SetTitleFont(42)
            hs.GetYaxis().SetLabelFont(42)
            hs.GetYaxis().SetLabelSize(.04)
        else:
            hs.GetXaxis().SetTitle(str(titleX))
            hs.GetXaxis().SetTitleOffset(0.00)
            hs.GetXaxis().SetTitleFont(42)
            hs.GetXaxis().SetTitleSize(0.05)
            hs.GetXaxis().SetLabelFont(42)
            hs.GetXaxis().SetLabelOffset(.01)
            hs.GetXaxis().SetLabelSize(0.04)
            hs.GetYaxis().SetTitle("Events/bin")
            hs.GetYaxis().SetTitleSize(0.08)
            hs.GetYaxis().SetTitleOffset(0.7)
            hs.GetYaxis().SetTitleFont(42)
            hs.GetYaxis().SetLabelFont(42)
            hs.GetYaxis().SetLabelSize(.05)

        # if not isrebin:
        if True:
            hs.GetXaxis().SetRangeUser(XMIN, XMAX)
        hs.GetXaxis().SetNdivisions(508)

#=============================  legend section =========================================
    DYLegend = "Z(ll)+jets "
    WLegend = "W(l#nu)+jets "
    GLegend = "#gamma+jets "
    ZLegend = "Z(#nu#nu)+jets "
    STLegend = "Single t "
    TTLegend = "t#bar{t} "
    VVLegend = "WW/WZ/ZZ "
    QCDLegend = "QCD "
    SMHLegend = "SMH "

    legend = SetLegend([.50, .58, .93, .92], ncol=2)

    if(not NORATIOPLOT):
        if 'SR' in reg:
            legend.AddEntry(h_data, "bkgSum", "PEL")
        else:
            legend.AddEntry(h_data, "Data", "PEL")
    legend.AddEntry(Top, TTLegend, "f")
    legend.AddEntry(STop, STLegend, "f")
    legend.AddEntry(WJets, WLegend, "f")
    legend.AddEntry(DIBOSON, VVLegend, "f")
    if GJetsCount > 0:
        legend.AddEntry(GJets, GLegend, "f")
    if ZJetsCount > 0:
        legend.AddEntry(ZJets, ZLegend, "f")
    legend.AddEntry(DYJets, DYLegend, "f")
    legend.AddEntry(QCD, QCDLegend, "f")
    legend.AddEntry(SMH, SMHLegend, "f")

    legend.Draw('same')

#=================================================latex section =====================
    t2d = ExtraText(str(histolabel), 0.20, 0.80)
    t2d.SetTextSize(0.06)

    t2d.SetTextAlign(12)
    t2d.SetNDC(ROOT.kTRUE)
    t2d.SetTextFont(42)
    t2d.Draw("same")

    pt = drawenergy1D(True, text_="Internal", data=True)
    for ipt in pt:
        ipt.Draw()
#======================================== ratio log ================

    ratioleg = SetLegend([.72, .80, .90, .90], 1)
    ratioleg.SetTextSize(0.15)

#============================================= statistical error section ======================

    ratiostaterr = h_err.Clone("ratiostaterr")
    ratiostaterr.Sumw2()
    ratiostaterr.SetStats(0)
    ratiostaterr.SetMinimum(0)
    ratiostaterr.SetMarkerSize()
    ratiostaterr.SetFillColor(ROOT.kBlack)
    ratiostaterr.SetFillStyle(3013)
    for i in range(h_err.GetNbinsX()+2):
        ratiostaterr.SetBinContent(i, 0.0)

        if (h_err.GetBinContent(i) > 1e-6):
            binerror = h_err.GetBinError(i)/h_err.GetBinContent(i)
            ratiostaterr.SetBinError(i, binerror)
        else:
            ratiostaterr.SetBinError(i, 999.)

    ratioleg.AddEntry(ratiostaterr, "stat", "f")

 #============================================= Lower Tpad Decalaration ====================================
    if(not NORATIOPLOT):
        c12.cd()
        DataMC = h_data.Clone()
        DataMC.Add(Stackhist, -1)   # remove for data/mc
        #DataMCPre = h_data.Clone()
        DataMC.Divide(Stackhist)
        DataMC.GetYaxis().SetTitle("#frac{Data-Pred}{Pred}")
        DataMC.GetYaxis().SetTitleSize(0.12)
        DataMC.GetYaxis().SetTitleOffset(0.42)
        DataMC.GetYaxis().SetTitleFont(42)
        DataMC.GetYaxis().SetLabelSize(0.12)
        DataMC.GetYaxis().CenterTitle()
        DataMC.GetXaxis().SetTitle(str(titleX))
        DataMC.GetXaxis().SetLabelSize(0.14)
        DataMC.GetXaxis().SetTitleSize(0.16)
        DataMC.GetXaxis().SetTitleOffset(1)
        DataMC.GetXaxis().SetTitleFont(42)
        DataMC.GetXaxis().SetTickLength(0.07)
        DataMC.GetXaxis().SetLabelFont(42)
        DataMC.GetYaxis().SetLabelFont(42)

    c1_1 = ROOT.TPad("c1_1", "newpad", 0, 0.00, 1, 0.3)
    if (not NORATIOPLOT):
        c1_1.Draw()
    c1_1.cd()
    c1_1.Range(-7.862408, -629.6193, 53.07125, 486.5489)
    c1_1.SetFillColor(0)
    c1_1.SetTicky(1)
    c1_1.SetLeftMargin(0.12)
    c1_1.SetRightMargin(0.06)
    c1_1.SetTopMargin(0.00)
    c1_1.SetBottomMargin(0.42)
    c1_1.SetFrameFillStyle(0)
    c1_1.SetFrameBorderMode(0)
    c1_1.SetFrameFillStyle(0)
    c1_1.SetFrameBorderMode(0)
    c1_1.SetLogy(0)

    if(not NORATIOPLOT):
        if (0):  # if(VARIABLEBINS)
            c1_1.SetLogx(0)
            DataMC.GetXaxis().SetMoreLogLabels()
            DataMC.GetXaxis().SetNoExponent()
            DataMC.GetXaxis().SetNdivisions(508)
        # if not isrebin:
        if True
            DataMC.GetXaxis().SetRangeUser(XMIN, XMAX)
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
        line1 = ROOT.TLine(XMIN, 0.2, XMAX, 0.2)
        line2 = ROOT.TLine(XMIN, -0.2, XMAX, -0.2)
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
    plot = str(hist)
    noPdfPng = True
    if ('_up' in str(hist) or '_down' in str(hist)):
        noPdfPng = False
    if not os.path.exists('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPng/'+reg):
        os.makedirs('plots_norm/'+datestr+'_' +
                    str(options.year)+'/bbDMPng/'+reg)
    if not os.path.exists('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMPdf/'+reg):
        os.makedirs('plots_norm/'+datestr+'_' +
                    str(options.year)+'/bbDMPdf/'+reg)
    if not os.path.exists('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMRoot/'):
        os.makedirs('plots_norm/'+datestr+'_'+str(options.year)+'/bbDMRoot/')
    if (ISLOG == 0) and noPdfPng:
        c12.SaveAs('plots_norm/'+datestr+'_'+str(options.year) +
                   '/bbDMPdf/'+reg+'/'+plot+'.pdf')
        c12.SaveAs('plots_norm/'+datestr+'_'+str(options.year) +
                   '/bbDMPng/'+reg+'/'+plot+'.png')
        print("Saved. \n")
    if (ISLOG == 1) and noPdfPng:
        c12.SaveAs('plots_norm/'+datestr+'_'+str(options.year) +
                   '/bbDMPdf/'+reg+'/'+plot+'_log.pdf')
        c12.SaveAs('plots_norm/'+datestr+'_'+str(options.year) +
                   '/bbDMPng/'+reg+'/'+plot+'_log.png')
    fshape = ROOT.TFile('plots_norm/'+datestr+'_' +
                        str(options.year)+'/bbDMRoot/'+plot+'.root', "RECREATE")
    fshape.cd()
    Stackhist.SetNameTitle("bkgSum", "bkgSum")
    Stackhist.Write()
    DIBOSON.SetNameTitle("DIBOSON", "DIBOSON")
    DIBOSON.Write()
    ZJets.SetNameTitle("ZJets", "ZJets")
    ZJets.Write()
    GJets.SetNameTitle("GJets", "GJets")
    GJets.Write()
    QCD.SetNameTitle("QCD", "QCD")
    QCD.Write()
    SMH.SetNameTitle("SMH", "SMH")
    SMH.Write()
    STop.SetNameTitle("STop", "STop")
    STop.Write()
    Top.SetNameTitle("Top", "Top")
    Top.Write()
    WJets.SetNameTitle("WJets", "WJets")
    WJets.Write()
    DYJets.SetNameTitle("DYJets", "DYJets")
    DYJets.Write()
    data_obs = h_data
    data_obs.SetNameTitle("data_obs", "data_obs")
    data_obs.Write()
    if makeSIGplots and ('MET' in hist) and ('SR' in hist) and not noYieldHisto:
        [hist.SetNameTitle(sname.split('5f_')[-1].split('_tanb35_')[0], sname.split('5f_')
                           [-1].split('_tanb35_')[0])for (sname, hist) in zip(signal_files_name, sig_hist)]
        [hist.Write() for hist in sig_hist]
    fshape.Write()
    fshape.Close()
    c12.Close()
    print('\n')
    #noYieldHisto = bool(('weight' in hist) or ('_up' in hist) or ('_down' in hist))
    if (('MET' in hist and 'SR' in hist) or ('Recoil' in hist)) and not noYieldHisto:
        bkg_list = {'ZJets': ZJets, 'DYJets': DYJets, 'WJets': WJets, 'STop': STop, 'GJets': GJets,
                    'Top': Top, 'DIBOSON': DIBOSON, 'QCD': QCD, 'SMH': SMH, 'bkgSum': Stackhist, 'data_obs': h_data}
        yield_outfile.write('region '+str(hist)+'\n')
        yield_outfile_binwise.write('region '+str(hist)+'\n')
        if makeSIGplots:
            for (sname, hist) in zip(signal_files_name, sig_hist):
                bkg_list.update(
                    {sname.split('5f_')[-1].split('_tanb35_')[0]: hist})
        yield_outfile_binwise.write('       Bin1   Bin2   Bin3   Bin4\n')
        for key in bkg_list:
            yield_outfile_binwise.write(str(key)+'   '+str.format('{0:.1f}', bkg_list[key].GetBinContent(1))+'\xb1'+str.format('{0:.1f}', bkg_list[key].GetBinError(1))+'   '+str.format('{0:.1f}', bkg_list[key].GetBinContent(2))+'\xb1'+str.format('{0:.1f}', bkg_list[key].GetBinError(
                2))+'   '+str.format('{0:.1f}', bkg_list[key].GetBinContent(3))+'\xb1'+str.format('{0:.1f}', bkg_list[key].GetBinError(3))+'   '+str.format('{0:.1f}', bkg_list[key].GetBinContent(4))+'\xb1'+str.format('{0:.1f}', bkg_list[key].GetBinError(4))+'\n')
        for key in bkg_list:
            binerror = 0.00
            bkg_list[key].Rebin(bkg_list[key].GetNbinsX())
            binerror = (bkg_list[key].GetBinError(1))
            print(str(key)+' '+str.format('{0:.3f}', bkg_list[key].Integral(
            ))+' +/- ' + str.format('{0:.3f}', binerror)+'\n')
            yield_outfile.write(str(key)+' '+str.format('{0:.1f}', bkg_list[key].GetBinContent(
                1))+'\xb1'+str.format('{0:.1f}', bkg_list[key].GetBinError(1))+'\n')
        yield_outfile_binwise.write('\n')
        yield_outfile.write('\n')

 #=======================================================================


######################################################################

regions = []
PUreg = []
if makeMuCHplots:
    regions += ['preselR', 'SR_1b', 'SR_2b', 'ZmumuCR_1b', 'ZmumuCR_2b',
                'TopmunuCR_1b', 'TopmunuCR_2b', 'WmunuCR_1b', 'WmunuCR_2b']
if makeEleCHplots:
    regions += ['ZeeCR_1b', 'ZeeCR_2b', 'TopenuCR_1b',
                'TopenuCR_2b', 'WenuCR_1b', 'WenuCR_2b']

# makeplot("reg_WenuCR_1b_Recoil",'h_reg_WenuCR_1b_Recoil','Recoil (GeV)',200.,1000.,1,1,0,'WenuCR_1b',varBin=False)
# makeplot("reg_SR_2b_MET",'h_reg_SR_2b_MET','p_{T}^{miss} (GeV)',200.,1000.,rebin,1,0,'SR_2b',varBin=False)
for reg in regions:
    if '_2b' in reg:
        rebin = 3
    else:
        rebin = 1
    try:
    # if True:
        if 'SR_' in reg or 'preselR' in reg:
            makeplot("reg_"+reg+"_cutFlow", 'h_reg_'+reg+'_cutFlow',
                     'CutFlow', 0, 7, 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET", 'h_reg_'+reg+'_MET',
                     ' p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightB_up", 'h_reg_'+reg+'_MET_weightB_up',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightB_down", 'h_reg_'+reg+'_MET_weightB_down',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightEWK_up", 'h_reg_'+reg+'_MET_weightEWK_up',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightEWK_down", 'h_reg_'+reg+'_MET_weightEWK_down',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightTop_up", 'h_reg_'+reg+'_MET_weightTop_up',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightTop_down", 'h_reg_'+reg+'_MET_weightTop_down',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightMET_up", 'h_reg_'+reg+'_MET_weightMET_up',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightMET_down", 'h_reg_'+reg+'_MET_weightMET_down',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightPU_up", 'h_reg_'+reg+'_MET_weightPU_up',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightPU_down", 'h_reg_'+reg+'_MET_weightPU_down',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightJEC_up", 'h_reg_'+reg+'_MET_weightJEC_up',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_weightJEC_down", 'h_reg_'+reg+'_MET_weightJEC_down',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_Res_up", 'h_reg_'+reg+'_MET_Res_up',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_Res_down", 'h_reg_'+reg+'_MET_Res_down',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_En_up", 'h_reg_'+reg+'_MET_En_up',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET_En_down", 'h_reg_'+reg+'_MET_En_down',
                     'p_{T}^{miss} (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_min_dPhi", 'h_reg_'+reg+'_min_dPhi',
                     'min_dPhi', 0, 4, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet1Pt", 'h_reg_'+reg+'_Jet1Pt',
                     'JET1 p_{T} (GeV)', 30., 800., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet1Eta", 'h_reg_'+reg+'_Jet1Eta',
                     'JET1 #eta', -2.5, 2.5, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet1Phi", 'h_reg_'+reg+'_Jet1Phi',
                     'JET1 #phi', -3.14, 3.14, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet1deepCSV", 'h_reg_'+reg+'_Jet1deepCSV',
                     'JET1 deepCSV', 0, 1.2, rebin, 0, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet2Pt", 'h_reg_'+reg+'_Jet2Pt',
                     'JET2 p_{T} (GeV)', 30., 800., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet2Eta", 'h_reg_'+reg+'_Jet2Eta',
                     'JET2 #eta', -2.5, 2.5, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet2Phi", 'h_reg_'+reg+'_Jet2Phi',
                     'JET2 #phi', -3.14, 3.14, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet2deepCSV", 'h_reg_'+reg+'_Jet2deepCSV',
                     'JET2 deepCSV', 0, 1.2, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_nJets", 'h_reg_'+reg+'_nJets',
                     'nJets', 0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_NEle", 'h_reg_'+reg+'_NEle',
                     'NEle', 0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_NMu", 'h_reg_'+reg+'_NMu', 'NMu',
                     0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_nPho", 'h_reg_'+reg+'_nPho',
                     'nPho', 0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_NTau", 'h_reg_'+reg+'_NTau',
                     'NTau', 0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_ratioPtJet21", 'h_reg_'+reg+'_ratioPtJet21',
                     'JET2 p_{T}/JET1 p_{T} ', 0., 1.0, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_dPhiJet12", 'h_reg_'+reg+'_dPhiJet12',
                     'JET1#eta - JET2#eta', -7.5, 7.5, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_dEtaJet12", 'h_reg_'+reg+'_dEtaJet12',
                     'JET1#phi - JET2#phi', -7.5, 7.5, rebin, 1, 0, reg, varBin=False)
            if ('SR_1b' in reg):
                makeplot("reg_"+reg+"_isjet1EtaMatch", 'h_reg_'+reg+'_isjet1EtaMatch',
                         'JET1#eta X JET2#eta', 0, 1, rebin, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_M_Jet1Jet2", 'h_reg_'+reg+'_M_Jet1Jet2',
                         'Inv Mass(Jet1, Jet2)', 0, 2000, 5, 1, 0, reg, varBin=False)
            elif ('SR_2b' in reg):
                makeplot("reg_"+reg+"_isjet2EtaMatch", 'h_reg_'+reg+'_isjet2EtaMatch',
                         'JET1#eta X JET3#eta', -1, 1, 1, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_M_Jet1Jet3", 'h_reg_'+reg+'_M_Jet1Jet3',
                         'Inv Mass(Jet1, Jet3)', 0, 2000, 10, 1, 0, reg, varBin=False)
            elif ('preselR' in reg):
                makeplot("reg_"+reg+"_isjet1EtaMatch", 'h_reg_'+reg+'_isjet1EtaMatch',
                         'JET1#eta X JET3#eta', -1, 1, 1, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_M_Jet1Jet3", 'h_reg_'+reg+'_M_Jet1Jet3',
                         'Inv Mass(Jet1, Jet2)', 0, 2000, 5, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_isjet2EtaMatch", 'h_reg_'+reg+'_isjet2EtaMatch',
                         'JET1#eta X JET3#eta', -1, 1, 1, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_M_Jet1Jet3", 'h_reg_'+reg+'_M_Jet1Jet3',
                         'Inv Mass(Jet1, Jet3)', 0, 2000, 5, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_rJet1PtMET", 'h_reg_'+reg+'_rJet1PtMET',
                     'Jet1 p_{T}/MET', 0, 20, 5, 1, 0, reg, varBin=False)
        else:
            makeplot("reg_"+reg+"_cutFlow", 'h_reg_'+reg+'_cutFlow',
                     'CutFlow', 0, 11, 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_MET", 'h_reg_'+reg+'_MET',
                     'Real p_{T}^{miss} (GeV)', 0., 700., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil", 'h_reg_'+reg+'_Recoil',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightB_up", 'h_reg_'+reg+'_Recoil_weightB_up',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightB_down", 'h_reg_'+reg+'_Recoil_weightB_down',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightEWK_up", 'h_reg_'+reg+'_Recoil_weightEWK_up',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightEWK_down", 'h_reg_'+reg+'_Recoil_weightEWK_down',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightTop_up", 'h_reg_'+reg+'_Recoil_weightTop_up',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightTop_down", 'h_reg_'+reg+'_Recoil_weightTop_down',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightRecoil_up", 'h_reg_'+reg +
                     '_Recoil_weightRecoil_up', 'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightRecoil_down", 'h_reg_'+reg +
                     '_Recoil_weightRecoil_down', 'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightPU_up", 'h_reg_'+reg+'_Recoil_weightPU_up',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightPU_down", 'h_reg_'+reg+'_Recoil_weightPU_down',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightEle_up", 'h_reg_'+reg+'_Recoil_weightEle_up',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightEle_down", 'h_reg_'+reg+'_Recoil_weightEle_down',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightMu_up", 'h_reg_'+reg+'_Recoil_weightMu_up',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightMu_down", 'h_reg_'+reg+'_Recoil_weightMu_down',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightJEC_up", 'h_reg_'+reg+'_Recoil_weightJEC_up',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_weightJEC_down", 'h_reg_'+reg+'_Recoil_weightJEC_down',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_Res_up", 'h_reg_'+reg+'_Recoil_Res_up',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_Res_down", 'h_reg_'+reg+'_Recoil_Res_down',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_En_up", 'h_reg_'+reg+'_Recoil_En_up',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Recoil_En_down", 'h_reg_'+reg+'_Recoil_En_down',
                     'Recoil (GeV)', 200., 1000., 1, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_min_dPhi", 'h_reg_'+reg+'_min_dPhi',
                     'min_dPhi', 0, 4, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet1Pt", 'h_reg_'+reg+'_Jet1Pt',
                     'JET1 p_{T} (GeV)', 30., 800., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet1Eta", 'h_reg_'+reg+'_Jet1Eta',
                     'JET1 #eta', -2.5, 2.5, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet1Phi", 'h_reg_'+reg+'_Jet1Phi',
                     'JET1 #phi', -3.14, 3.14, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_Jet1deepCSV", 'h_reg_'+reg+'_Jet1deepCSV',
                     'JET1 deepCSV', 0, 1.2, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_lep1_pT", 'h_reg_'+reg+'_lep1_pT',
                     'lepton1 p_{T}', 0, 500, rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_nJets", 'h_reg_'+reg+'_nJets',
                     'nJets', 0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_NEle", 'h_reg_'+reg+'_NEle',
                     'NEle', 0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_NMu", 'h_reg_'+reg+'_NMu', 'NMu',
                     0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_nPho", 'h_reg_'+reg+'_nPho',
                     'nPho', 0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_NTau", 'h_reg_'+reg+'_NTau',
                     'NTau', 0., 10., rebin, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_nPV", 'h_reg_'+reg+'_nPV',
                     'Before PU reweighting', 0., 70., 1, 0, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_PUnPV", 'h_reg_'+reg+'_PUnPV',
                     'After PU reweighting', 0., 70., 1, 0, 0, reg, varBin=False)
            if ('W' in reg) or ('Top' in reg):
                makeplot("reg_"+reg+"_Wmass", 'h_reg_'+reg+'_Wmass',
                         'W candidate mass (GeV)', 0., 165., rebin, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_WpT", 'h_reg_'+reg+'_WpT',
                         'W candidate p_{T} (GeV)', 0., 700., rebin, 1, 0, reg, varBin=False)
            if ('Z' in reg):
                makeplot("reg_"+reg+"_Zmass", 'h_reg_'+reg+'_Zmass',
                         'Z candidate mass (GeV)', 70., 110., 1, 0, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_ZpT", 'h_reg_'+reg+'_ZpT',
                         'Z candidate p_{T} (GeV)', 0., 700., rebin, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_lep2_pT", 'h_reg_'+reg+'_lep2_pT',
                         'lepton2 p_{T}', 0, 500, rebin, 1, 0, reg, varBin=False)
            if ('WmunuCR_1b' not in reg) and ('WenuCR_1b' not in reg):
                makeplot("reg_"+reg+"_Jet2Pt", 'h_reg_'+reg+'_Jet2Pt',
                        'JET2 p_{T} (GeV)', 30., 800., rebin, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_Jet2Eta", 'h_reg_'+reg+'_Jet2Eta',
                        'JET2 #eta', -2.5, 2.5, rebin, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_Jet2Phi", 'h_reg_'+reg+'_Jet2Phi',
                        'JET2 #phi', -3.14, 3.14, rebin, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_Jet2deepCSV", 'h_reg_'+reg+'_Jet2deepCSV',
                        'JET2 deepCSV', 0, 1.2, rebin, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_ratioPtJet21", 'h_reg_'+reg+'_ratioPtJet21',
                         'JET2 p_{T}/JET1 p_{T} ', 0., 1.0, rebin, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_dPhiJet12", 'h_reg_'+reg+'_dPhiJet12',
                         'JET1#eta - JET2#eta', -7.5, 7.5, rebin, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_dEtaJet12", 'h_reg_'+reg+'_dEtaJet12',
                         'JET1#phi - JET2#phi', -7.5, 7.5, rebin, 1, 0, reg, varBin=False)
            if ('1b' in reg) and ('WmunuCR_1b' not in reg) and ('WenuCR_1b' not in reg):
                makeplot("reg_"+reg+"_isjet1EtaMatch", 'h_reg_'+reg+'_isjet1EtaMatch',
                         'JET1#eta X JET2#eta', -1, 1, 1, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_M_Jet1Jet2", 'h_reg_'+reg+'_M_Jet1Jet2',
                         'Inv Mass(Jet1, Jet2)', 0, 2000, rebin, 1, 0, reg, varBin=False)
            elif ('WmunuCR_2b' in reg or 'WenuCR_2b' in reg):
                makeplot("reg_"+reg+"_isjet1EtaMatch", 'h_reg_'+reg+'_isjet1EtaMatch',
                         'JET1#eta X JET2#eta', -1, 1, 1, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_M_Jet1Jet2", 'h_reg_'+reg+'_M_Jet1Jet2',
                         'Inv Mass(Jet1, Jet2)', 0, 2000, 5, 1, 0, reg, varBin=False)
            elif ('2b' in reg):
                makeplot("reg_"+reg+"_isjet2EtaMatch", 'h_reg_'+reg+'_isjet2EtaMatch',
                         'JET1#eta X JET3#eta', -1, 1, 1, 1, 0, reg, varBin=False)
                makeplot("reg_"+reg+"_M_Jet1Jet3", 'h_reg_'+reg+'_M_Jet1Jet3',
                         'Inv Mass(Jet1, Jet3)', 0, 2000, 5, 1, 0, reg, varBin=False)
            makeplot("reg_"+reg+"_rJet1PtMET", 'h_reg_'+reg+'_rJet1PtMET',
                     'Jet1 p_{T}/MET', 0, 20, 5, 1, 0, reg, varBin=False)
    except Exception as e:
        print(e)
        print("Cannot Plot")
        pass
yield_outfile.close()
yield_outfile_binwise.close()