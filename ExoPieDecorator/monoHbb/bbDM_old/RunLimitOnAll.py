import sys 
import os 

os.system('rm bin/limits_bbDM2016scalar.txt')
os.system('rm bin/limits_bbDM2016pseudo.txt')

def createAndrunSetup(signalStr):
    fout = open('datacard_bbDM_2016_'+signalStr+'.txt', 'w')
    
    for iline in open('comb_tmpl.txt'):
        iline = iline.replace("SIGNALPOINT", signalStr)
        fout.write(iline)


signal_ = ['NLO_pseudo_mphi_50_mchi_1', 'NLO_pseudo_mphi_100_mchi_1', 'NLO_pseudo_mphi_350_mchi_1', 'NLO_pseudo_mphi_400_mchi_1', 'NLO_pseudo_mphi_500_mchi_1', 'NLO_pseudo_mphi_1000_mchi_1', 'NLO_scalar_mphi_300_mchi_1', 'NLO_scalar_mphi_350_mchi_1','NLO_scalar_mphi_350_mchi_1','NLO_scalar_mphi_400_mchi_1','NLO_scalar_mphi_500_mchi_1', 'NLO_scalar_mphi_750_mchi_1']
## add signal mass point histograms here if you want to extent the analysis 

if len(sys.argv) > 1:
    if sys.argv[1]=='create':
        for sig in signal_:
            createAndrunSetup(sig)
    if sys.argv[1]=='run':
        for sig in signal_:
            datacardname = 'datacard_bbDM_2016_'+sig+'.txt'
            command_ = 'python scan.py '+datacardname
            os.system(command_)
