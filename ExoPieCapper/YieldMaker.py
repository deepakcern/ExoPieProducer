from ROOT import TFile, gROOT, kBlack,TH1F


Yield = open("Yield.txt","w")
Yield.write("bin"+"     "+"200-270 GeV"+"    "+"270-345 GeV" +"    "+"345-480 GeV"+"    "+"> 480 GeV"+'\n')
#f= TFile.Open('fitDiagnostics_merged_2017_data.root','READ')
f= TFile.Open('fitDiagnostics_resolved_2017_data.root','READ')

Hists = ["zjets","tt","wjets","singlet","diboson","smh",]
histLabel = ["Z+jets","tt","W+jets","Single t","Diboson","SM h"]
bin1Sum = 0 ; bin1SumErr = 0
bin2Sum = 0 ; bin2SumErr = 0
bin3Sum = 0 ; bin3SumErr = 0
bin4Sum = 0 ; bin4SumErr = 0

for  index,hist in enumerate(Hists):
    h_temp = f.Get("shapes_prefit/SR/"+hist)

    bin1 = h_temp.GetBinContent(1)*70
    bin1Sum+=bin1
    bin1Err = h_temp.GetBinError(1)*70
    bin1SumErr+=bin1Err

    bin2 = h_temp.GetBinContent(2)*75
    bin2Sum+=bin2
    bin2Err = h_temp.GetBinError(2)*75
    bin2SumErr+=bin2Err

    bin3 = h_temp.GetBinContent(3)*135
    bin3Sum+=bin3
    bin3Err = h_temp.GetBinError(3)*135
    bin3SumErr+=bin3Err

    bin4 = h_temp.GetBinContent(4)*520
    bin4Sum+=bin4
    bin4Err = h_temp.GetBinError(4)*520
    bin4SumErr+=bin4Err

    Yield.write(histLabel[index] + " & " + '%.2f'%bin1+" \pm "+'%.2f'%bin1Err + " & " + '%.2f'%bin2+" \pm "+'%.2f'%bin2Err+ " & " + '%.2f'%bin3+" \pm "+'%.2f'%bin3Err+ " & " + '%.2f'%bin4+" \pm "+'%.2f'%bin4Err+'\n')

Yield.write("SUM (SM)" + " & " + '%.2f'%bin1Sum+" \pm "+'%.2f'%bin1SumErr + " & " + '%.2f'%bin2Sum+" \pm "+'%.2f'%bin2SumErr+ " & " + '%.2f'%bin3Sum+" \pm "+'%.2f'%bin3SumErr+ " & " + '%.2f'%bin4Sum+" \pm "+'%.2f'%bin4SumErr+'\n')

h_signal = f.Get("shapes_prefit/SR/total_signal")

sig_bin1 = h_signal.GetBinContent(1)*70
sig_bin1Err = h_signal.GetBinError(1)*70

sig_bin2 = h_signal.GetBinContent(2)*75
sig_bin2Err = h_signal.GetBinError(2)*75

sig_bin3 = h_signal.GetBinContent(3)*135
sig_bin3Err = h_signal.GetBinError(3)*135

sig_bin4 = h_signal.GetBinContent(4)*520
sig_bin4Err = h_signal.GetBinError(4)*520

Yield.write("2HDM+a mA=1000,ma=150 GeV" + " & " + '%.2f'%sig_bin1+" \pm "+'%.2f'%sig_bin1Err + " & " + '%.2f'%sig_bin2+" \pm "+'%.2f'%sig_bin2Err+ " & " + '%.2f'%sig_bin3+" \pm "+'%.2f'%sig_bin3Err+ " & " + '%.2f'%sig_bin4+" \pm "+'%.2f'%bin4Err+'\n')
Yield.close()
