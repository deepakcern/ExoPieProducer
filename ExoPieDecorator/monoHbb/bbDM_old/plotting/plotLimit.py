import os
import sys 
from ROOT import * 

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gROOT.SetBatch(1)

c = TCanvas("c","c",1500, 950)

c.SetLogy(1)
c.SetLogx(1)
c.SetGridx(1)
c.SetGridy(1)
rootfilepath=""
if len(sys.argv)<2: 
    print ("tell me which mediator you need, s or ps? ")
    sys.exit()
if sys.argv[1]=="s":
    rootfilename = "scalar.root"
    
if sys.argv[1]=="ps":
    rootfilename = "pseudo.root"


f = TFile(rootfilepath + rootfilename,"read")


exp2s =  f.Get("exp2")
exp2s.SetMarkerStyle(20)
exp2s.SetMarkerSize(1.1)
exp2s.SetLineWidth(2)
exp2s.SetFillColor(kYellow);
exp2s.SetLineColor(kYellow)
exp2s.GetXaxis().SetRangeUser(0,1200)
#exp2s.GetYaxis().SetRangeUser(10,100000)
exp2s.GetXaxis().SetTitle("M_{#phi} [GeV]");
exp2s.GetXaxis().SetTitleOffset(1.4)
exp2s.GetYaxis().SetTitle("#sigma/#sigma_{th}"); 
exp2s.GetYaxis().SetTitleOffset(1.2)
exp2s.GetYaxis().SetNdivisions(505);
exp2s.GetXaxis().SetNdivisions(505);
exp2s.GetYaxis().SetMoreLogLabels()
exp2s.GetXaxis().SetMoreLogLabels()
exp2s.Draw("A 3")

exp1s =  f.Get("exp1")
exp1s.SetMarkerStyle(20)
exp1s.SetMarkerSize(1.1)
exp1s.SetLineWidth(2)
exp1s.SetFillColor(kGreen);
exp1s.SetLineColor(kGreen)
exp1s.Draw("3 same")

exp =  f.Get("expmed")
exp.SetMarkerStyle(1)
exp.SetMarkerSize(1.1)
exp.SetLineWidth(3)
exp.Draw("L same")

obs =  f.Get("obs")
obs.SetMarkerStyle(20)
#obs.SetMarkerColor(4)
obs.SetMarkerSize(1.1)
#obs.SetLineColor(2)
obs.SetLineWidth(3)
obs.SetLineStyle(9)
obs.Draw("P same")

#hr = c.DrawFrame(fr_left, fr_down, fr_right, fr_up, "");
#hr.SetXTitle("M_{Z'} [GeV]");
#hr.SetYTitle("95% CLs on #sigma(Z`#rightarrow#chi#bar{#chi}H)#timesBR(H#rightarrowb#bar{b})[pb]");
#hr.SetMinimum(0.001);
#hr.SetMaximum(1000);
#hr.Draw(same)


leg = TLegend(.30, .65, .55, .890);
leg.SetFillColor(0);
leg.SetShadowColor(0);
leg.SetTextFont(42);
leg.SetTextSize(0.03);
leg.AddEntry(obs, "CL_{S} Observed", "LP");
leg.AddEntry(exp1s, "CL_{S}  Expected #pm 1#sigma", "LF");
leg.AddEntry(exp2s, " CL_{S}  Expected #pm 2#sigma", "LF");
leg.AddEntry(exp, " CL_{S}  Expected ", "L");

leg.Draw("same")
  

latex =  TLatex();
latex.SetNDC();
latex.SetTextSize(0.04);
latex.SetTextAlign(31);
latex.SetTextAlign(11);
model_ = rootfilename.replace(".root","")
if model_ == "pseudo":
    latex.DrawLatex(0.18, 0.93, "bb+DM pseudoscalar");

if model_ == "scalar":
    latex.DrawLatex(0.18, 0.93, "bb+DM scalar");


name = "limit_"+model_
c.SaveAs(name+".png")
c.SaveAs(name+".pdf")
