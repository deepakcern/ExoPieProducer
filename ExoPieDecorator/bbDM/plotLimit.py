import os
import sys
import ROOT as rt

rt.gStyle.SetOptTitle(0)
rt.gStyle.SetOptStat(0)
rt.gROOT.SetBatch(1)

c = rt.TCanvas("c","c",1500, 950)
c.SetLogy(1)
#c.SetLogx(0)
c.SetGrid(1,1);
leg = rt.TLegend(.15, .65, .35, .890);

rootfilename = 'limit_bbDM2016_2HDMa.root'


f = rt.TFile(rootfilename,"read")


exp2s =  f.Get("exp2")
exp2s.SetMarkerStyle(20)
exp2s.SetMarkerSize(1.1)
exp2s.SetLineWidth(2)
exp2s.SetFillColor(rt.kYellow);
exp2s.SetLineColor(rt.kYellow)
exp2s.GetXaxis().SetTitle("Ma [GeV]");
exp2s.GetYaxis().SetRangeUser(.1,100)
exp2s.GetXaxis().SetTitleOffset(1.4)
exp2s.GetYaxis().SetTitle("95% C.L. asymptotic limit on #mu=#sigma/#sigma_{theory}");
exp2s.GetYaxis().SetTitleOffset(1.2)
exp2s.GetYaxis().SetNdivisions(20,5,0);
#exp2s.GetXaxis().SetNdivisions(505);
exp2s.GetYaxis().SetMoreLogLabels()
#exp2s.GetXaxis().SetMoreLogLabels()
#exp2s.GetXaxis().SetRangeUser(30,500)
exp2s.Draw("A 3")

exp1s =  f.Get("exp1")
exp1s.SetMarkerStyle(20)
exp1s.SetMarkerSize(1.1)
exp1s.SetLineWidth(2)
exp1s.SetFillColor(rt.kGreen);
exp1s.SetLineColor(rt.kGreen)
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

leg = rt.TLegend(.15, .65, .40, .890);
leg.SetFillColor(0);
leg.SetShadowColor(0);
leg.SetTextFont(42);
leg.SetTextSize(0.03);
leg.AddEntry(exp, " CL_{S}  Expected ", "LP");
leg.AddEntry(exp1s, "CL_{S}  Expected #pm 1#sigma", "LF");
leg.AddEntry(exp2s, " CL_{S}  Expected #pm 2#sigma", "LF");
# leg.AddEntry(obs, "CL_{S} Observed", "LP");

leg.Draw("same")

# line = TLine(max(xmin,exp2s.GetXaxis().GetXmin()),1,exp2s.GetXaxis().GetXmax(),1)
c.Update()
line = rt.TLine(c.GetUxmin(),1.0,c.GetUxmax(),1.0);
line.SetLineColor(rt.kRed)
line.SetLineWidth(2)
line.Draw('same ')

latex =  rt.TLatex();
latex.SetNDC();
latex.SetTextSize(0.04);
latex.SetTextAlign(31);
latex.SetTextAlign(11);
model_ = '2HDM+a'
latex.DrawLatex(0.11, 0.91, "bb+DM");
latex.DrawLatex(0.55, 0.91, "M_{A}=600 GeV, tan#beta = 35, sin#theta = .7");

name = "limit_"+model_

c.SaveAs(name+".png")
c.SaveAs(name+".pdf")
