from ROOT import TFile

import os,traceback
import sys, optparse,argparse
import sample_xsec_2017_GenXSecAnalyser as sample_xsec
import glob

usage = "python DataframeToHist.py -F -inDir directoryName -D outputDir "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--inputfile",  dest="inputfile",default="myfiles.root")
parser.add_argument("-F", "--farmout", action="store_true",  dest="farmout")
parser.add_argument("-inDir", "--inputDir",  dest="inputDir",default=".")
parser.add_argument("-D", "--outputdir", dest="outputdir",default=".")

args = parser.parse_args()

if args.farmout==None:
    isfarmout = False
else:
    isfarmout = args.farmout

if args.inputDir and isfarmout:
    inDir=args.inputDir

outputdir = '.'
if args.outputdir:
    outputdir = str(args.outputdir)

infile  = args.inputfile


A=13.9
B=7.04
C=6.87
D=31.5
calclumi = 59.64
lumi = (C+D) * 1000

def scaleHists(infilename,outDir):
    xsec=1.0
    isDatafile = False
    rootFile   = infilename.split('/')[-1]
    if 'MET' in rootFile or 'SE' in rootFile:isDatafile=True
    print 'isDatafile',isDatafile

    f = TFile.Open(infilename, "READ")

    h_total_mcweight = f.Get("h_total_mcweight")
    totalEvents = h_total_mcweight.Integral()
    Keys = f.GetListOfKeys()

    if outDir == '.':outfilename = rootFile.replace(".root","_norm.root")

    else:outfilename = outDir+'/'+rootFile
    print 'outfilename',outfilename
    fout = TFile.Open(outfilename,"RECREATE")
    if not isDatafile:xsec = sample_xsec.getXsec(rootFile)
    print 'xsec',xsec,'lumi',lumi
    # xsec=1.0

    dirs = {}
    td = None

    for key in Keys:
        if key.GetClassName() == 'TDirectory':
        	td = key.ReadObj()
        	dirName = td.GetName()
        	print "found directory", dirName
        	dirs[dirName] = td

        elif key.GetClassName() == 'TH1F':
            hist = key.ReadObj()
            histName = hist.GetName()
            if 'h_total_mcweight' in histName or 'h_total' in histName:continue
            if not isDatafile:hist.Scale(lumi*xsec/totalEvents)

            fout.cd()
            hist.Write()
            fout.Write()

    	# print 'events',hist.Integral(),'total',totalEvents
        else:continue


if isfarmout:
    path=inDir
    files=glob.glob(path+'/'+'*.root')
    for inputFile in files:
        print 'running code for file:  ',inputFile
        scaleHists(inputFile,outputdir)

if not isfarmout:
    inputFile=infile
    print 'running code for file:  ',inputFile
    scaleHists(inputFile,outputdir)
