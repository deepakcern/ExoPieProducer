#!/usr/bin/python

import sys
import os
import sys
import datetime
import sys, optparse,argparse
import ROOT as rt
import array
import uproot
import numpy as np

usage = "python systLister.py -r path -b path "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-b", "--bPath",  dest="bPath",default=".")
parser.add_argument("-r", "--rPath",  dest="rPath",default=".")
args = parser.parse_args()


listFiles=''

dic={}

def HistValues(f_hist):
    return f_hist.numpy()[0]

def HistBins():
    return f_hist.numpy()[1]

for cat in ["B","R"]:
    print 'working on cat :: ',cat
    if cat=="B":path = args.bPath
    else:path = args.rPath
    for reg in ['SR','Zee','We','Tope','Zmumu','Wmu','Topmu']:
        for syst in ["btagweight","ewkweight","toppTweight","metTrigweight","puweight","jec","Res","En","lepweight"]:

            if reg=="SR":var="MET"
            else:var="Recoil"
            if reg=="SR" and 'lep' in syst:continue
            f_up = uproot.open(path+'/monoHROOT/'+'h_reg_'+reg+'_'+var+'_'+syst+'_up.root')
            f_cent = uproot.open(path+'/monoHROOT/'+'h_reg_'+reg+'_'+var+'.root')
            f_down=uproot.open(path+'/monoHROOT/'+'h_reg_'+reg+'_'+var+'_'+syst+'_down.root')
            array_up=HistValues(f_up["bkgSum"])
            array_cent=HistValues(f_cent["bkgSum"])
            array_down=HistValues(f_down["bkgSum"])
            # print array_up, array_cent,array_down
            updiff = np.abs(array_up-array_cent)
            downdiff=np.abs(array_cent-array_down)
            max=np.maximum(updiff,downdiff)
            systematics=np.divide(max,array_cent).tolist()
            dic[var+'_'+reg+'_'+cat+'_'+syst]=systematics
            # break


print "writing file"
f = open("syst_dic.py","w")
f.write( "dic={" )
for key , value in dic.items():
    print key,value
    f.write( '"'+str(key)+'":'+str(value)+','+'\n' )

f.write( "}" )
f.close()
print 'done'
