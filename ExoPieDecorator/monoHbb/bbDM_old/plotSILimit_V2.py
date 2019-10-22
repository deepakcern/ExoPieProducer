from __future__ import print_function
import sys
import math
import pandas as pd
import numpy as np
import sklearn as skl

#Global variables
g_q = 0.25
g_DM = 1.0
m_n = 0.939
convert_to_cm2 = 0.3894e-27

mediator_type = "V"

luxdm = [3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,12.0,14.0,17.0,21.0,27.0,33.0,40.0,50.0,60.0,70.0,100.0,300.0,1000.0,4000.0,20000.0]
luxxsec = [3.05386E-039,7.32711E-040,2.7257E-040,1.1551E-040,6.1778E-041,3.34041E-041,1.23152E-041,
           3.85766E-042,9.19548E-043,2.9655E-043,1.18273E-043,6.06967E-044,3.35316E-044,1.35406E-044,
           7.01193E-045,4.3259E-045,2.07382E-045,1.40966E-045,9.79264E-046,7.32674E-046,6.21113E-046,
           5.62354E-046,6.06104E-046,6.17893E-046,6.8994E-046,7.44363E-046,9.92359E-046,2.67048E-045,
           8.43005E-045,3.40367E-044,1.66915E-043]
cdmslitedm = [1.429,1.473,1.574,1.654,1.731,1.838,1.986,2.169,2.384,2.632,2.857,3.093,3.415,3.827,4.207,4.625,5.106,5.613,6.157,6.768,7.504,8.161,9.009,9.968,11.100,12.307,13.645,14.871]
cdmslitexsec = [9.880e-38,3.162e-38,1.075e-38,4.334e-39,1.924e-39,8.338e-40,3.279e-40,1.566e-40,8.338e-41,5.012e-41,3.527e-41,2.701e-41,2.120e-41,1.789e-41,1.566e-41,1.492e-41,1.438e-41,1.321e-41,1.143e-41,9.527e-42,8.040e-42,7.122e-42,6.387e-42,5.796e-42,5.261e-42,5.073e-42,4.775e-42,4.892e-42]
mp_2012_mass = [ 1. , 10. , 100. , 200., 300. , 500. , 1000.  ]
mp_2012_xsec = [8.22306e-40,2.65427e-39,3.1704e-39,3.41978e-39,3.73833e-39,4.95775e-39,3.07741e-38]
lux_new_dm = [5,7,10,12,14,17,21,33,50,100,200,1000,4000,20000,100000]
lux_new_xsec = [3400.17e-45,35.8824e-45,4.2465e-45,1.83477e-45,1.09486e-45,0.696387e-45,0.475726e-45,0.277026e-45,0.21929e-45,
                0.264416e-45,0.428959e-45,2.03858e-45,8.01987e-45,41.0509e-45,211.864e-45]
pandaXdm = [990.3314945554625,925.2193898972866,872.8271965619264,823.4018043468851,776.7752128626581,732.788935040558,684.6096838574658,645.8424430196985,609.2704661368674,574.7694424835353,542.2221006501592,511.517809929314,482.5522042741279,455.22682755073095,425.2966708288167,397.334355690884,367.62145015755084,336.8419402964968,311.652694492943,288.34711585860623,261.6504698748821,239.74349678010785,219.67070907932353,203.24359943267606,186.2268057438061,170.63476180478258,154.8365256585499,140.5009707543591,128.73737259389227,114.57021549500398,102.95755673125127,91.62739011886737,78.43587900951637,68.46096838574658,61.52187115995493,52.15543573159424,46.4158883361278,38.9688239742875,34.01304938279253,30.86394787194334,28.00640623313109,25.41343036702638,23.060525425632015,21.129757595204445,19.36064542292897,17.91284454622005,16.573311102378653,15.48365256585499,14.606863203649896,13.646494265700854,12.873737259389234,12.144738992807904,11.568875283162827,11.020317142808507,10.497769831163431,10,9.525832782420359,9.07414901986344,8.643882620598271,8.23401804346886,7.920164050192553,7.692649574879144,7.399430787258665,7.117388546363303,6.8460968385746614,6.714353762575055,6.458424430196983,6.272899858196248,6.092704661368679,5.917685748188826,5.747694424835354,5.5825862688626975,5.475157597543134,5.3697962331787465,5.215543573159428,5.1151780992931455,5.016744011523652,4.9682395947343885]
pandaXxsec = [3.357062045166831e-45,3.0891712650846654e-45,2.842657948119843e-45,2.670770489307242e-45,2.4576452088760942e-45,2.3090383777601755e-45,2.169417379983613e-45,2.0382388677057253e-45,1.875589088548534e-45,1.7621775391848543e-45,1.655623664355352e-45,1.5555127999426033e-45,1.4614553553916502e-45,1.3730852975827298e-45,1.2900587263800697e-45,1.2120525363130744e-45,1.1153317444167957e-45,1.0263291918746627e-45,9.6426997830899e-46,8.873220300067351e-46,8.16514464461493e-46,7.513572841979244e-46,6.9139959313469694e-46,6.362264763256088e-46,5.977556668786939e-46,5.6161107810093825e-46,5.1679497765625874e-46,4.75555165033135e-46,4.467996938035354e-46,4.111454470535231e-46,3.7833637976298236e-46,3.4814544895963734e-46,3.137718623969464e-46,2.8873312678395163e-46,2.712742540163987e-46,2.548710697377543e-46,2.4962679247207412e-46,2.4449042248699515e-46,2.548710697377543e-46,2.656924616043907e-46,2.8279209717161206e-46,3.137718623969464e-46,3.4814544895963734e-46,4.026856374600323e-46,4.6577004841269566e-46,5.5005525808547175e-46,6.632395697431267e-46,7.832586357813745e-46,9.444289695572828e-46,1.16268683214822e-45,1.4313841625215726e-45,1.7621775391848543e-45,2.169417379983613e-45,2.6158161903617592e-45,3.220332219668606e-45,3.9645521131213423e-45,5.1948812059322726e-45,6.666958706748279e-45,8.556179946276695e-45,1.1211440064157123e-44,1.4092375238430638e-44,1.8085742497312438e-44,2.419620070785238e-44,3.041375039853739e-44,4.0689355664047476e-44,4.906196304291031e-44,6.563806297053605e-44,8.423797158552623e-44,1.0810855067477154e-43,1.3874335420260504e-43,1.8561925481719616e-43,2.432229292982666e-43,3.121451986840807e-43,3.923552410357784e-43,5.141157965763318e-43,6.5980118705069475e-43,8.122814703356869e-43,9.59270986458227e-43]
cresstdm = [0.502,0.535,0.566,0.603,0.638,0.675,0.720,0.767,0.811,0.879,0.942,1.004,1.106,1.242,1.405,1.517,1.727,2.000,2.159,2.339,2.483,2.992,3.285,3.523,4.273,5.166,6.184,7.703,10.393,13.654,17.466,23.330,26.476,29.748]
cresstxsec = [1.929e-36,1.164e-36,6.032e-37,2.755e-37,1.428e-37,7.402e-38,4.036e-38,2.257e-38,1.468e-38,9.317e-39,6.708e-39,4.953e-39,3.390e-39,2.206e-39,1.629e-39,1.331e-39,1.087e-39,8.442e-40,7.073e-40,5.223e-40,3.760e-40,1.237e-40,6.914e-41,4.978e-41,3.079e-41,2.274e-41,1.637e-41,1.209e-41,1.065e-41,1.149e-41,1.304e-41,1.637e-41,1.480e-41,9.628e-42]

interpolations_on_line = []
interpolations_by_eye = [
        ((10., 150.), (100., 50.), (90., 60.)),
        ((175., 50.), (300., 125.), (250., 110.)),
        ((100., 50.), (200., 100.), (120., 60.)),
        ((125., 50.), (200., 100.), (150., 70.)),
        ((325., 100.), (400., 150.), (380., 140.)),
        ((400., 100.), (400., 150.), (400., 125.)),
]
if mediator_type == "V":
  interpolations_by_eye.append(((325.,150.),(400.,200.),(360.,175.)))
  interpolations_by_eye.append(((200.,100.),(300.,150.),(250.,125.)))
  interpolations_by_eye.append(((80.,25.),(125.,25.),(100.,25.)))
  interpolations_by_eye.append(((300.,200.),(525.,300.),(325.,220.)))
  interpolations_by_eye.append(((175.,50.),(300.,100.),(200.,60.)))
  interpolations_by_eye.append(((175.,50.),(300.,125.),(200.,70.)))
  
  interpolations_by_eye.append(((50.,25.),(100.,50.),(90.,45.)))
  interpolations_by_eye.append(((50.,25.),(100.,50.),(70,35.)))
  interpolations_by_eye.append(((20.,10.),(50.,25.),(40.,20.)))
  interpolations_by_eye.append(((20.,10.),(50.,25.),(30.,15.)))
  interpolations_by_eye.append(((10.,5.),(20.,10.),(15.,7.)))
  interpolations_by_eye.append(((325.,150.),(525.,200.),(425.,175.)))
  interpolations_by_eye.append(((10.,5.),(20.,10.),(18.,9.)))
  interpolations_by_eye.append(((20.,10.),(40.,25.),(28.,16.)))
  interpolations_by_eye.append(((20.,10.),(40.,25.),(32.,19.)))
  interpolations_by_eye.append(((40.,10.),(80.,25.),(56.,16.)))
  interpolations_by_eye.append(((40.,10.),(80.,25.),(64.,19.)))
  interpolations_by_eye.append(((10.,10.),(50.,50.),(18.,18.)))
  interpolations_by_eye.append(((10.,10.),(50.,50.),(15.,15.)))
  interpolations_by_eye.append(((10.,10.),(50.,50.),(12.,12.)))
  

for p in interpolations_by_eye:
  p0x = 99999.9
  if ((p[1][0]-p[0][0])**2 + (p[1][1]-p[0][1])**2) != 0.0:
    p0x = (p[2][0]*(p[1][0]-p[0][0])**2 + p[0][0]*(p[1][1]-p[2][1])*(p[1][1]-p[0][1]) + p[1][0]*(p[2][1]-p[0][1])*(p[1][1]-p[0][1]))/((p[1][0]-p[0][0])**2 + (p[1][1]-p[0][1])**2)
  if p0x == 99999.9:
    print("Error: p0x == 99999.9")
    exit(1)
  p0y = 99999.9
  if p[0][1]-p[1][1] != 0.0:
    p0y = ((p[2][0]-p0x)*(p[0][0]-p[1][0]))/(p[0][1]-p[1][1]) + p[2][1]
  elif p[1][0]-p[0][0] != 0.0:
    p0y = (p[1][1]*(p0x-p[0][0])+p[0][1]*(p[1][0]-p0x))/(p[1][0]-p[0][0])
  if p0y == 99999.9:
    print("Error: p0y == 99999.9")
    exit(1)
  interpolations_on_line.append(((p[0][0],p[0][1]),(p[1][0],p[1][1]),(p0x,p0y)))

preface = """#include <fstream>
#include <vector>
#include <iomanip>
#include "TFile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"

void plot()
{
  TCanvas *c = new TCanvas("c", "canvas",700,640);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  c->SetLeftMargin(0.15);
  //c->SetLogx();
  //c->SetLogy();
  c->SetLogz();
  c->cd();
  TGraph* lux_limit = new TGraph();
  TGraph* cdmslite_limit = new TGraph();
  TGraph* lux_limit_new = new TGraph();
  TGraph* pandaX_limit = new TGraph();
  TGraph* cresst_limit = new TGraph();
  //TGraph* mp_limit = new TGraph();
  //TGraph2D* r_limit_histo = new TGraph2D();
  TGraph2D* r_limit_histo_copy = new TGraph2D();
  TGraph2D* r_limit_histo_expected50 = new TGraph2D();
  //TGraph2D* r_limit_histo_expected50_plus = new TGraph2D();
  //TGraph2D* r_limit_histo_expected50_minus = new TGraph2D();
  TGraph* masspoints_histo = new TGraph();
  TGraph* masspoints_remaining_histo = new TGraph();

  TGraph* range_histo = new TGraph();
  range_histo->SetPoint(0,log10(1.0),-50.0);
  range_histo->SetPoint(1,log10(900.0),-33.0);
"""

# masspoints_histo = interpolation points
# masspoints_remaining_histo = sample points

print(preface)

for i in range(0,len(luxdm)):
  print("lux_limit->SetPoint("+str(i)+","+str(math.log10(luxdm[i]))+","+str(math.log10(luxxsec[i]))+");")
for i in range(0,len(cdmslitedm)):
  print("cdmslite_limit->SetPoint("+str(i)+","+str(math.log10(cdmslitedm[i]))+","+str(math.log10(cdmslitexsec[i]))+");")
for i in range(0,len(lux_new_dm)):
  print("lux_limit_new->SetPoint("+str(i)+","+str(math.log10(lux_new_dm[i]))+","+str(math.log10(lux_new_xsec[i]))+");")
for i in range(0,len(pandaXdm)):
  print("pandaX_limit->SetPoint("+str(i)+","+str(math.log10(pandaXdm[i]))+","+str(math.log10(pandaXxsec[i]))+");")
for i in range(0,len(cresstdm)):
  print("cresst_limit->SetPoint("+str(i)+","+str(math.log10(cresstdm[i]))+","+str(math.log10(cresstxsec[i]))+");")
#for i in range(0,len(mp_2012_mass)):
#  print("mp_limit->SetPoint("+str(i)+","+str(math.log10(mp_2012_mass[i]))+","+str(math.log10(mp_2012_xsec[i]))+");")


already_printed = [(0,0)]
printed_limits = {}
i = 0
for line in open("limits_barzp_cleaned_NoDuplicate.txt").readlines():
#for line in open("limits_benedikt.txt").readlines():
#for line in open("limits.txt").readlines():
  line = line.rstrip()
  label = line.split()
  mx = int(label[1])
  mv = int(label[0])
  r_limit = float(label[4]) ## expected limit
  r_limit_obs = float(label[7])
  
  
  if (mv,mx) not in already_printed and 1 <= mv <= 2000 and 1 <= mx <= 700:
    xsec = 0.0
    xsec = convert_to_cm2*9*g_q**2*g_DM**2*m_n**2*mx**2/(3.1415927*mv**4*(m_n+mx)**2)
    #elif mediator_type == "AV":
    #  xsec = convert_to_cm2*3*0.32**2*g_q**2*g_DM**2*m_n**2*mx**2/(3.1415927*mv**4*(m_n+mx)**2)
#    print("r_limit_histo_copy->SetPoint("+str(i)+","+str(math.log10(mx))+","+str(math.log10(xsec))+","+str(r_limit)+");")
    if mx == 15: continue
    if mx == 25: continue
    if r_limit > 30.0: continue
    print("r_limit_histo_copy->SetPoint("+str(i)+","+str(math.log10(mx))+","+str(math.log10(xsec))+","+str(r_limit)+");")
    #print("r_limit_histo_copy->SetPoint("+str(i)+","+str(mx)+","+str(math.log10(xsec))+","+str(r_limit)+");")
    
    #print (" raman "+str(i)+" "+str(mv) + " " +str(mx)+ " " + str(xsec) + " " + str(r_limit) )
    already_printed.append((mv,mx))
    printed_limits[(mv,mx)] = r_limit
    i += 1


'''
for p in interpolations_on_line:
  if (int(p[0][0]),int(p[0][1])) in already_printed and (int(p[1][0]),int(p[1][1])) in already_printed:
    z_p1 = printed_limits[(int(p[0][0]),int(p[0][1]))]
    z_p2 = printed_limits[(int(p[1][0]),int(p[1][1]))]
    z_0_x = 99999.9
    z_0_y = 99999.9
    if p[1][0]-p[0][0] != 0.0:
      z_0_x = z_p1 + (z_p2-z_p1)*(p[2][0]-p[0][0])/(p[1][0]-p[0][0])
    if p[1][1]-p[0][1] != 0.0:
      z_0_y = z_p1 + (z_p2-z_p1)*(p[2][1]-p[0][1])/(p[1][1]-p[0][1])
    if z_0_x == 99999.9 and z_0_y == 99999.9:
      continue
    elif z_0_x == 99999.9:
      z_0_x = z_0_y
    elif z_0_y == 99999.9:
      z_0_y = z_0_x
    z_0 = (z_0_x + z_0_y)/2.0
    xsec = 0.0
    xsec = convert_to_cm2*9*g_q**2*g_DM**2*m_n**2*p[2][1]**2/(3.1415927*p[2][0]**4*(m_n+p[2][1])**2)
    #elif mediator_type == "AV":
    #  xsec = convert_to_cm2*3*0.32**2*g_q**2*g_DM**2*m_n**2*p[2][1]**2/(3.1415927*p[2][0]**4*(m_n+p[2][1])**2)
    print("r_limit_histo_copy->SetPoint("+str(i)+","+str(math.log10(p[2][1]))+","+str(math.log10(xsec))+","+str(z_0)+");")
    i += 1

'''
## do similar for observed limit also 


appendix = """
gStyle->SetLineStyleString(11,"15 20");

range_histo->GetYaxis()->SetRangeUser(-50.0,-33.0);
range_histo->GetXaxis()->SetRangeUser(log10(1.0),log10(900.0));
range_histo->SetTitle("");
range_histo->GetXaxis()->SetTickLength(0);
range_histo->GetXaxis()->SetLabelOffset(999);
range_histo->GetYaxis()->SetTickLength(0);
range_histo->GetYaxis()->SetLabelOffset(999);
range_histo->SetMarkerColor(kWhite);
range_histo->Draw("AP");

//lux_limit->SetLineColor(kPink+10);
//lux_limit->SetLineWidth(3);
//lux_limit->Draw("L SAME");
lux_limit_new->SetLineColor(kBlue+1);
lux_limit_new->SetLineWidth(3);
lux_limit_new->Draw("L SAME");
pandaX_limit->SetLineColor(kSpring-8);
pandaX_limit->SetLineWidth(3);
pandaX_limit->Draw("L SAME");
cdmslite_limit->SetLineColor(kPink+10);
cdmslite_limit->SetLineWidth(3);
cdmslite_limit->Draw("L SAME");
cresst_limit->SetLineColor(kOrange);
cresst_limit->SetLineWidth(3);
cresst_limit->Draw("L SAME");
//mp_limit->SetLineColor(kRed);
//mp_limit->SetLineWidth(3);
//mp_limit->Draw("L SAME");

r_limit_histo_expected50->SetNpx(100);
r_limit_histo_expected50->SetNpy(100);
r_limit_histo_expected50->SetMinimum(0.10);
TH2D* delaunay_expected50 = r_limit_histo_expected50->GetHistogram();
delaunay_expected50->SetMinimum(1);
delaunay_expected50->SetContour(1);
delaunay_expected50->SetLineWidth(3);
delaunay_expected50->SetLineStyle(11);
//gStyle->SetNumberContours(255);
delaunay_expected50->Draw("CONT3 SAME");

r_limit_histo_copy->SetNpx(100);
r_limit_histo_copy->SetNpy(100);
TH2D* delaunay_copy = r_limit_histo_copy->GetHistogram();
delaunay_copy->SetMinimum(1);
delaunay_copy->SetContour(1);
delaunay_copy->SetLineWidth(3);
delaunay_copy->Draw("CONT3 SAME");

delaunay_copy->SetLineColor(kBlack);
delaunay_expected50->SetLineColor(kBlack);

masspoints_histo->SetMarkerStyle(20);
masspoints_histo->SetMarkerSize(0.7);
masspoints_histo->SetMarkerColor(kBlue);
masspoints_histo->SetLineColor(kBlue);
//masspoints_histo->Draw("P SAME");

masspoints_remaining_histo->SetMarkerStyle(20);
masspoints_remaining_histo->SetMarkerSize(0.7);
masspoints_remaining_histo->SetMarkerColor(kBlack);
masspoints_remaining_histo->SetLineColor(kBlack);
//masspoints_remaining_histo->Draw("P SAME");

TLegend* leg = new TLegend(0.598711,0.740407,0.832665,0.862757,"");
leg->AddEntry(cresst_limit,"CRESST-II","L");
leg->AddEntry(cdmslite_limit,"CDMSLite 2015","L");
leg->AddEntry(pandaX_limit,"PandaX-II","L");
leg->AddEntry(lux_limit_new,"LUX 2016","L");
//leg->AddEntry(mp_limit,"Monophoton, #sqrt{s} = 8 TeV","L");
leg->SetFillColor(kWhite);
leg->SetFillStyle(0);
leg->SetTextSize(0.035);
leg->Draw();

TLegend* leg2 = new TLegend(0.358281,0.132561,0.742235,0.273279,"Vector, Dirac, #it{g}_{q} = 0.25, #it{g}_{DM} = 1");
leg2->AddEntry(delaunay_copy,"Observed 90% CL");
leg2->AddEntry(delaunay_expected50,"Median expected 90% CL");
leg2->SetFillColor(kWhite);
leg2->SetFillStyle(0);
leg2->SetTextSize(0.035);
leg2->Draw();

TLatex *texS = new TLatex(0.61023,0.907173,"35.9 fb^{-1} (13 TeV)");
texS->SetNDC();
texS->SetTextFont(42);
texS->SetTextSize(0.045);
texS->Draw();
TLatex *texS1 = new TLatex(0.19092,0.837173,"#bf{CMS}");
texS1->SetNDC();
texS1->SetTextFont(42);
texS1->SetTextSize(0.045);
texS1->Draw();

c->Update();

double xmin = c->GetUxmin();
double ymin = c->GetUymin();
double xmax = c->GetUxmax();
double ymax = c->GetUymax();

TGaxis *xaxis = new TGaxis(xmin,ymin,xmax,ymin,pow(10.0,xmin),pow(10.0,xmax),50510,"G");
xaxis->SetTitle("#it{m}_{DM} [GeV]");
xaxis->SetLabelFont(42);
xaxis->SetLabelSize(0.040);
xaxis->SetTitleFont(42);
xaxis->SetTitleSize(0.045);
xaxis->Draw("SAME");

TGaxis *xaxis_top = new TGaxis(xmin,ymax,xmax,ymax,pow(10.0,xmin),pow(10.0,xmax),50510,"G-");
xaxis_top->SetTitle("");
xaxis_top->SetLabelOffset(999);
xaxis_top->SetTitleOffset(999);
xaxis_top->Draw("SAME");

TGaxis *yaxis = new TGaxis(xmin,ymin,xmin,ymax,pow(10.0,ymin),pow(10.0,ymax),50510,"G");
yaxis->SetTitle("#it{#sigma}_{SI} (DM-nucleon) [cm^{2}]");
yaxis->SetLabelFont(42);
yaxis->SetLabelSize(0.040);
yaxis->SetTitleFont(42);
yaxis->SetTitleSize(0.045);
yaxis->SetTitleOffset(1.5);
yaxis->Draw("SAME");

TGaxis *yaxis_right = new TGaxis(xmax,ymin,xmax,ymax,pow(10.0,ymin),pow(10.0,ymax),50510,"G+");
yaxis_right->SetTitle("");
yaxis_right->SetLabelOffset(999);
yaxis_right->SetTitleOffset(999);
yaxis_right->Draw("SAME");

c->SaveAs("SI.png");
c->SaveAs("SI.pdf");
}"""

print(appendix)
