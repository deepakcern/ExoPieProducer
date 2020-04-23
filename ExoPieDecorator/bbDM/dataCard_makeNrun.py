import sys
import os

def createAndrunSetup(ma,mA):
    fout1b = open('datacards_bbDM/datacard_bbDM_2016_1b_Ma'+str(ma)+'_MChi1_MA'+str(mA)+'.txt', 'w')
    fout2b = open('datacards_bbDM/datacard_bbDM_2016_2b_Ma'+str(ma)+'_MChi1_MA'+str(mA)+'.txt', 'w')
    for iline in open('combine_maYY_MAXX_bbDM_sig1b.txt'):
        iline = iline.replace("_SIGMaPOINT", str(ma))
        iline = iline.replace("_SIGMAPOINT", str(mA))
        fout1b.write(iline)
    for iline in open('combine_maYY_MAXX_bbDM_sig2b.txt'):
        iline = iline.replace("_SIGMaPOINT", str(ma))
        iline = iline.replace("_SIGMAPOINT", str(mA))
        fout2b.write(iline)

ma_point_600 =[10,100,150,200,250,300,400,450]
ma_point_1200 =[10,50,100,150,200,250,300,400,450,500,750]

## add signal mass point histograms here if you want to extent the analysis

if len(sys.argv) > 1:
    if sys.argv[1]=='create':
        for sig in ma_point_600:
            createAndrunSetup(sig,600)
        for sig in ma_point_1200:
            createAndrunSetup(sig,1200)
    if sys.argv[1]=='run':
        listfiles = [f for f in os.listdir('datacards_bbDM')]
        for sig in listfiles:
            datacardname = 'datacards_bbDM/'+sig
            command_ = 'combine -M AsymptoticLimits --rAbsAcc 0 --rMax 30 -t -1 '+datacardname+' &> logfiles/'+sig.strip('datacard_bbDM_2016_')
            os.system(command_)
