#!/usr/bin/env python
import os
import sys
import datetime
import sys, optparse
import ROOT as ROOT
import array
import string

datestr = str(datetime.date.today().strftime("%d%m%Y"))

#command  python StackPlotter_syst.py  -y <Year>
usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)

usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)
parser.add_option("-y", "--year", dest="year", default="Year")
(options, args) = parser.parse_args()
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
    luminosity_ = '{0:.2f}'.format(35.81)
elif runOn2017:
    luminosity_ = '{0:.2f}'.format(41.50)
elif runOn2018:
    luminosity_ = '{0:.2f}'.format(59.64)

def SetCMSAxis(h, xoffset=1., yoffset=1.):
    h.GetXaxis().SetTitleSize(0.047)
    h.GetYaxis().SetTitleSize(0.047)

    print (type(h))
    if type(h) is ( (not ROOT.TGraphAsymmErrors) or (not ROOT.TGraph)):  h.GetZaxis().SetTitleSize(0.047)

    h.GetXaxis().SetLabelSize(0.047)
    h.GetYaxis().SetLabelSize(0.047)
    if type(h) is ( (not ROOT.TGraphAsymmErrors) or (not ROOT.TGraph)): h.GetZaxis().SetLabelSize(0.047)

    h.GetXaxis().SetTitleOffset(xoffset)
    h.GetYaxis().SetTitleOffset(yoffset)
    return h

def ExtraText(text_,x_, y_):
    if not text_: print ("nothing provided as text to ExtraText, function crashing")
    ltx = ROOT.TLatex(x_,y_,text_)

    if len(text_)>0:
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
    c.SetTicky(1)
    c.SetGrid(1)
    return c

def SetLegend(coordinate_=[.50,.65,.90,.90],ncol=2):
    c_=coordinate_
    legend=ROOT.TLegend(c_[0], c_[1],c_[2],c_[3])
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
    pt = ROOT.TPaveText(0.0877181,0.95,0.9580537,0.96,"brNDC")
    pt.SetBorderSize(0)
    pt.SetTextAlign(12)
    pt.SetFillStyle(0)
    pt.SetTextFont(52)
    cmstextSize = 0.07
    preliminarytextfize = cmstextSize * 0.7
    lumitextsize = cmstextSize *0.7
    pt.SetTextSize(cmstextSize)
    text = pt.AddText(0.03,0.57,"#font[61]{CMS}")

    #pt1 = ROOT.TPaveText(0.0877181,0.9,0.9580537,0.96,"brNDC")
    pt1 = ROOT.TPaveText(0.0877181,0.95,0.9580537,0.96,"brNDC")
    pt1.SetBorderSize(0)
    pt1.SetTextAlign(12)
    pt1.SetFillStyle(0)
    pt1.SetTextFont(52)

    pt1.SetTextSize(preliminarytextfize)
    #text1 = pt1.AddText(0.215,0.4,text_)
    text1 = pt1.AddText(0.15,0.4,text_)

    #pt2 = ROOT.TPaveText(0.0877181,0.9,0.9580537,0.96,"brNDC")
    pt2 = ROOT.TPaveText(0.0877181,0.95,0.9580537,0.96,"brNDC")
    pt2.SetBorderSize(0)
    pt2.SetTextAlign(12)
    pt2.SetFillStyle(0)
    pt2.SetTextFont(52)
    pt2.SetTextFont(42)
    pt2.SetTextSize(lumitextsize)

    pavetext = ''
    if is2017 and data: pavetext = str(luminosity_)+' fb^{-1}'+" (13 TeV)"
    if (not is2017) and data: pavetext = str(luminosity_)+' fb^{-1}'+"(13 TeV)"

    if is2017 and not data: pavetext = "13 TeV"
    if (not is2017) and not data: pavetext = "13 TeV"

    if data: text3 = pt2.AddText(0.68,0.5,pavetext)
    if not data: text3 = pt2.AddText(0.68,0.5,pavetext)

    return [pt,pt1,pt2]


def CustomiseHistogram(h, titleX, titleY, color, lineStyle,title):
    h.SetMarkerColor(color)
    h.SetMarkerSize(1.0)
    h.SetLineColor(color)
    h.SetLineWidth(2)
    h.SetLineStyle(lineStyle)
    h = SetCMSAxis(h)
    h.GetYaxis().SetTitle(titleY)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetLabelSize(0.055)
    h.SetTitle(title)
    return

def CustomiseRatio(h1,h2,h3,titleX):
    h1.SetMarkerSize(0.7)
    h1.SetMarkerStyle(20)
    h1.SetMarkerColor(ROOT.kRed)
    h1.SetLineColor(ROOT.kRed)
    h2.SetMarkerSize(0.7)
    h2.SetMarkerStyle(20)
    h2.SetMarkerColor(ROOT.kBlue)
    h2.SetLineColor(ROOT.kBlue)
    h3.SetMarkerSize(0.7)
    h3.SetMarkerStyle(20)
    h3.SetMarkerColor(ROOT.kBlack)
    h3.SetLineColor(ROOT.kBlack)
    h1 = SetCMSAxis(h1)
    h1.Draw("P e1")
    h2.Draw("P e1 same")
    h3.Draw("P e1 same")
    h1.SetMinimum(-0.2)
    h1.SetMaximum(2.2)
    h1.GetXaxis().SetTitle(titleX)
    h1.GetXaxis().SetTitleSize(0.16)
    h1.GetXaxis().SetLabelSize(0.14)
    h1.GetYaxis().SetTitle()
    h1.GetYaxis().SetLabelSize(0.12)
    h1.GetYaxis().SetNdivisions(505)
    return

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetGridStyle(2)
uncertfile=open("uncert.txt","w")
CRSRPath = '/Users/ptiwari/cernBox/Documents/ExoPieCapper/plots_norm/analysis_plots_v18_06-04-00-HemAtSkim/bbDMRoot'
for jetprop in ['weightB','weightEWK','weightTop','weightMET','weightEle','weightMu','weightPU','weightJEC','Res','En']:
    for reg in ['SR_1b','SR_2b','ZmumuCR_1b','ZmumuCR_2b','TopmunuCR_1b','TopmunuCR_2b','WmunuCR_1b','WmunuCR_2b','ZeeCR_1b','ZeeCR_2b','TopenuCR_1b','TopenuCR_2b','WenuCR_1b','WenuCR_2b']:
        if 'SR' in reg:
            miss_En = 'MET'
            if jetprop=='weightEle' or jetprop=='weightMu': continue
        else:
            miss_En = 'Recoil'
            if jetprop == 'weightMET': jetprop='weightRecoil'
        #try:
        for syst in ['up','down']:
            exec("systematics_"+reg+"_"+miss_En+"_"+jetprop+"_"+syst+"_file = ROOT.TFile('"+CRSRPath+"/h_reg_"+reg+"_"+miss_En+"_"+jetprop+"_"+syst+".root')")
            exec(jetprop+"_syst_"+reg+"_"+miss_En+"_"+syst+" = systematics_"+reg+"_"+miss_En+"_"+jetprop+"_"+syst+"_file.Get('bkgSum')")

        exec("central_"+reg+"_"+jetprop+"_file = ROOT.TFile('"+CRSRPath+"/h_reg_"+reg+"_"+miss_En+".root')")
        exec(jetprop+"_syst_"+reg+"_"+miss_En+"_central = central_"+reg+"_"+jetprop+"_file.Get('bkgSum')")
        colors = {
                "up"     : ROOT.kRed,
                "down"     : ROOT.kBlue,
                "central" : ROOT.kBlack,
                }
        if 'SR' in reg:
            titleX = "E_{T}^{miss} (GeV)"
        else:
            titleX = "Hadronic Recoil (GeV)"
        titleY = "#Events"
        # Set Canvas
        c1 = myCanvas1D()
        c1_1 =  ROOT.TPad("c1_1","newpad",0,0.20,1,1);
        c1_1.SetBottomMargin(0.0);
        c1_1.SetTopMargin(0.08);
        c1_1.SetLeftMargin(0.12);
        c1_1.SetRightMargin(0.06);
        c1_1.SetLogy(1);
        c1_1.SetGrid(1)
        c1_1.Draw();
        c1_1.cd();
        exec("CustomiseHistogram("+jetprop+"_syst_"+reg+"_"+miss_En+"_up, titleX, titleY, colors['up'], 1,'Up')")
        exec("CustomiseHistogram("+jetprop+"_syst_"+reg+"_"+miss_En+"_down, titleX, titleY, colors['down'], 1,'Down')")
        exec("CustomiseHistogram("+jetprop+"_syst_"+reg+"_"+miss_En+"_central, titleX, titleY, colors['central'], 1,'Central')")

        exec("upUnc = "+jetprop+"_syst_"+reg+"_"+miss_En+"_up.Integral()")
        exec("downUnc = "+jetprop+"_syst_"+reg+"_"+miss_En+"_down.Integral()")
        exec("centralUnc = "+jetprop+"_syst_"+reg+"_"+miss_En+"_central.Integral()")

        uncertfile.write(jetprop+"_syst_"+reg+"_"+miss_En+": ")
        uncertfile.write(str((max(abs(upUnc-centralUnc),abs(centralUnc-downUnc))/centralUnc)*100)+"\n")

        # exec(jetprop+"_syst_"+reg+"_"+miss_En+"_up.Rebin(1)")
        # exec(jetprop+"_syst_"+reg+"_"+miss_En+"_down.Rebin(1)")
        # exec(jetprop+"_syst_"+reg+"_"+miss_En+"_central.Rebin(1)")

        exec(jetprop+"_syst_"+reg+"_"+miss_En+"_up.Draw()")
        exec(jetprop+"_syst_"+reg+"_"+miss_En+"_down.Draw('same')")
        exec(jetprop+"_syst_"+reg+"_"+miss_En+"_central.Draw('same')")

        leg  = SetLegend([0.7, 0.75, 0.95, 0.9],ncol=1)
        exec("leg.AddEntry("+jetprop+"_syst_"+reg+"_"+miss_En+"_up, 'Up' , 'l')")
        exec("leg.AddEntry("+jetprop+"_syst_"+reg+"_"+miss_En+"_central, 'Central' , 'l')")
        exec("leg.AddEntry("+jetprop+"_syst_"+reg+"_"+miss_En+"_down, 'Down' , 'l')")
        leg.Draw("same")

        #texcms.Draw("same")
        #texCat.Draw("same")
        latex1 = ROOT.TLatex();
        latex1.SetNDC();
        latex1.SetTextSize(0.06);
        latex1.SetTextAlign(31);
        latex1.SetTextAlign(11);
        latex1.DrawLatex(0.30, .80, jetprop+'_'+reg)
        
        pt = drawenergy1D(True, text_="Internal", data=True)
        for ipt in pt:
            ipt.Draw()

        exec(jetprop+"_syst_"+reg+"_"+miss_En+"_up.GetXaxis().SetRangeUser(200, 1000)")
        exec(jetprop+"_syst_"+reg+"_"+miss_En+"_down.GetXaxis().SetRangeUser(200, 1000)")
        exec(jetprop+"_syst_"+reg+"_"+miss_En+"_central.GetXaxis().SetRangeUser(200, 1000)")

        exec("ratioUp = "+jetprop+"_syst_"+reg+"_"+miss_En+"_up.Clone()")
        exec("ratioDown = "+jetprop+"_syst_"+reg+"_"+miss_En+"_down.Clone()")
        exec("ratioCentral = "+jetprop+"_syst_"+reg+"_"+miss_En+"_central.Clone()")

        exec("ratioUp.Divide("+jetprop+"_syst_"+reg+"_"+miss_En+"_central)")
        exec("ratioCentral.Divide("+jetprop+"_syst_"+reg+"_"+miss_En+"_central)")
        exec("ratioDown.Divide("+jetprop+"_syst_"+reg+"_"+miss_En+"_central)")

        c1.cd()

        c1_2 = ROOT.TPad("c1_2", "newpad",0,0.00,1,0.3)
        c1_2.Draw()
        c1_2.cd()
        #c1_2.Range(-7.862408,-629.6193,53.07125,486.5489);
        c1_2.SetFillColor(0);
        c1_2.SetLeftMargin(0.12);
        c1_2.SetRightMargin(0.06);
        c1_2.SetTopMargin(0.00);
        c1_2.SetBottomMargin(0.42);
        c1_2.SetFrameFillStyle(0);
        c1_2.SetFrameBorderMode(0);
        c1_2.SetFrameFillStyle(0);
        c1_2.SetFrameBorderMode(0);
        c1_2.SetLogy(0);
        exec("CustomiseRatio(ratioUp,ratioDown,ratioCentral, titleX)")
        c1_2.SetGrid(1)
        c1_2.Update()
        c1.Update()
        c1.Draw()
        if not os.path.exists('syst_plots/'+datestr+'_'+str(options.year)+'/bbDMPng/'+reg):
            os.makedirs('syst_plots/'+datestr+'_'+str(options.year)+'/bbDMPng/'+reg)
        if not os.path.exists('syst_plots/'+datestr+'_'+str(options.year)+'/bbDMPdf/'+reg):
            os.makedirs('syst_plots/'+datestr+'_'+str(options.year)+'/bbDMPdf/'+reg)

        exec("c1.SaveAs('syst_plots/"+datestr+'_'+str(options.year)+"/bbDMPdf/"+reg+"/"+jetprop+"_syst_"+reg+"_"+miss_En+".pdf')")
        exec("c1.SaveAs('syst_plots/"+datestr+'_'+str(options.year)+"/bbDMPng/"+reg+"/"+jetprop+"_syst_"+reg+"_"+miss_En+".png')")
        c1.Close()
        exec("central_"+reg+"_"+jetprop+"_file.Close()")
        exec("systematics_"+reg+"_"+miss_En+"_"+jetprop+"_up_file.Close()")
        exec("systematics_"+reg+"_"+miss_En+"_"+jetprop+"_down_file.Close()")
