import ROOT as ROOT

luminosity_ = 41.5

def drawenergy1D(is2017, text_="Work in progress 2018", data=True):
    #pt = ROOT.TPaveText(0.0877181,0.9,0.9580537,0.96,"brNDC")
    pt = ROOT.TPaveText(0.0877181,0.95,0.9580537,0.96,"brNDC")
    pt.SetBorderSize(0)
    pt.SetTextAlign(12)
    pt.SetFillStyle(0)
    pt.SetTextFont(52)

    cmstextSize = 0.07
    preliminarytextfize = cmstextSize * 0.75
    lumitextsize = cmstextSize *0.85
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
    pt2 = ROOT.TPaveText(0.0877181,0.95,0.9080537,0.96,"brNDC")
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



def set_overflow(hist):
    bin_num = hist.GetXaxis().GetNbins()
    #print (bin_num)
    hist.SetBinContent(bin_num,hist.GetBinContent(bin_num+1)+hist.GetBinContent(bin_num)) #Add overflow bin content to last bin
    hist.SetBinContent(bin_num+1,0.)
    return hist

def setHistStyle(h_temp2,hist):
    dovarbin=False
    h_temp_=h_temp2
    return h_temp_

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
    c = ROOT.TCanvas("myCanvasName","The Canvas Title",650,600)
    c.SetBottomMargin(0.050)
    c.SetRightMargin(0.050)
    c.SetLeftMargin(0.050)
    c.SetTopMargin(0.050)
    return c

#def SetLegend(coordinate_=[.38,.7,.890,.87],ncol=2):
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
    legend.SetTextSize(0.046)
    return legend
