#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import glob
import math
import operator
import optparse
import ROOT as ROOT
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TFile, TGraphAsymmErrors, TLatex, TLine, gStyle, TLegend, gROOT, TGraph, kBlack, kBlue, kRed
from array import array
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import text
from matplotlib.colors import LogNorm

def SetCanvas():
    c = TCanvas("myCanvasName","The Canvas Title",650,600)
    c.SetBottomMargin(0.050)
    c.SetRightMargin(0.050)
    c.SetLeftMargin(0.050)
    c.SetTopMargin(0.050)
    return c


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
    pt = ROOT.TPaveText(0.0997181,0.95,0.9580537,0.96,"brNDC")
    pt.SetBorderSize(0)
    pt.SetTextAlign(12)
    pt.SetFillStyle(0)
    pt.SetTextFont(52)

    cmstextSize = 0.07
    preliminarytextfize = cmstextSize * 0.7
    lumitextsize = cmstextSize *0.7
    pt.SetTextSize(cmstextSize)
    text = pt.AddText(0.03,0.57,"#font[60]{CMS}")

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
    pt2 = ROOT.TPaveText(0.0997181,0.95,0.9580537,0.96,"brNDC")
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

def getLatex():
    latex =  TLatex()
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(31);
    latex.SetTextAlign(11);
    latex.SetTextColor(1);
    return latex


def getGraph(n,x,y,lc,mc,ms):
    gr =TGraph(n,x,y)
    gr.SetFillColor(4)
    #gr.SetFillStyle(3004)
    gr.SetLineColor(4)
    gr.SetLineWidth(2)
    gr.SetMarkerStyle(ms)
    gr.SetMarkerSize(1.5)
    gr.SetLineColor(lc)
    gr.SetLineWidth(1)
    gr.SetMarkerColor(mc)
    gr.GetYaxis().SetTitle("Signal Efficiency")
    gr.GetXaxis().SetTitle("M_{a} (GeV)")
#     gr.SetTitle("") 
    return gr
    

# In[2]:
usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)

parser.add_option("-S", "--sigDir", type="string", dest="SIGrootFileDir",
                  help="directory containing signal histogram")
parser.add_option("-t", "--tag", type="string",
                  dest="plot_tag", help="version of histogram")
parser.add_option("-e", "--era", type="string",
                  dest="era_year", help="year of histogram")
(options, args) = parser.parse_args()

runOn2016 = False
runOn2017 = False
runOn2018 = False

if options.era_year == '2016':
    runOn2016 = True
elif options.era_year == '2017':
    runOn2017 = True
elif options.era_year == '2018':
    runOn2018 = True
else:
    print('Please provide on which year you want to run?')

if runOn2016:
    luminosity_ = '{0:.2f}'.format(35.82)
    era_name = 'bbDM2016_'
elif runOn2017:
    luminosity_ = '{0:.2f}'.format(41.50)
    era_name = 'bbDM2017_'
elif runOn2018:
    luminosity_ = '{0:.2f}'.format(59.64)
    era_name = 'bbDM2018_'

if options.SIGrootFileDir == None:
    print('Please provide signal histogram directory name')
    sys.exit()
else:
    SignalPath = options.SIGrootFileDir

if options.plot_tag == None:
    print('Please provide histogram directory name')
    sys.exit()
else:
    plot_tag = options.plot_tag

# In[3]:


gStyle.SetErrorX(0.5)
gStyle.SetFrameLineWidth(3)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetLegendBorderSize(0)
gStyle.SetFillColor(2)
gStyle.SetLineWidth(1)
gStyle.SetHistFillStyle(2)
gROOT.SetBatch(True)

sig_list = [SignalPath+'/' +fl for fl in os.listdir(SignalPath) if '.root' in fl and 'bbDM_2HDMa' in fl]

if not os.path.exists('Signal_Efficiency_Plots'):
    os.makedirs('Signal_Efficiency_Plots')
    
for cat in ('1b','2b'):
    sig_eff_ma_600 = {}
    sig_eff_ma_1200 = {}
    for iFile in sig_list:
        fin = TFile(iFile,"READ")
        rootFile = iFile.split('/')[-1]
        #print ('rootFile', rootFile.split('_'))
        mA_ = int(rootFile.split('_')[8].strip('MA'))
        #print (mA_)
        if mA_ == 600:
            #print (mA_)
            hist = fin.Get("h_reg_SR_"+cat+"_MET")
            hist_total = fin.Get("h_total_mcweight")
            hist_eff = hist.Integral()/hist_total.Integral()
            sig_eff_ma_600.update({int(rootFile.split('_')[6].strip('Ma')):hist_eff})
        if mA_ == 1200:
            #print (mA_)
            hist = fin.Get("h_reg_SR_"+cat+"_MET")
            hist_total = fin.Get("h_total_mcweight")
            hist_eff = hist.Integral()/hist_total.Integral()
            sig_eff_ma_1200.update({int(rootFile.split('_')[6].strip('Ma')):hist_eff})        
    
    c1 = SetCanvas()
    c1.SetTickx()
    c1.SetTicky()
    c1.SetGridx()
    c1.SetGridy()

    c1_2 = TPad("c1_2","newpad",0,0.05,1,1)
    c1_2.SetBottomMargin(0.10)
    c1_2.SetTopMargin(0.08)
    c1_2.SetLeftMargin(0.15)
    c1_2.SetRightMargin(0.01)
    c1_2.Draw()
    c1_2.cd()

    sig_eff_ma_1200_sorted = sorted(sig_eff_ma_1200.items(), key=operator.itemgetter(0))
    x12, y12 = zip(*sig_eff_ma_1200_sorted)
    x12 = array('d',x12)
    y12 = array('d',y12)
    gr12 = getGraph(len(x12),x12,y12,4,4,20)
    gr12 = SetCMSAxis(gr12,1,1.6)
    gr12.Draw()

    sig_eff_ma_600_sorted = sorted(sig_eff_ma_600.items(), key=operator.itemgetter(0))
    x6, y6 = zip(*sig_eff_ma_600_sorted)
    x6 = array('d',x6)
    y6 = array('d',y6)
    gr6=getGraph(len(x6),x6,y6,2,2,22)
    gr6.Draw('pl same')
    
    if cat == '1b': sr = 'SR1'
    elif cat == '2b': sr = 'SR2'
        
    latex=getLatex()
    latex.DrawLatex(0.22, 0.74,'#splitline{2HDM+a model, '+sr+'}{tan#beta = 35, sin#theta = 0.7}')

    legend = SetLegend([.55,.30,.90,.40],ncol=1)
    legend.AddEntry(gr12,"M_{A} = 1200 GeV","PEL")
    legend.AddEntry(gr6,"M_{A} = 600 GeV","PEL")
    legend.Draw('p same')

    pt = drawenergy1D(True,text_="   Internal",data=True)
    for ipt in pt: ipt.Draw()

    c1.Update()
    c1.SaveAs('Signal_Efficiency_Plots/Sig_Eff_'+sr+'_'+plot_tag+'.pdf')
    c1.Close()


# In[ ]:




