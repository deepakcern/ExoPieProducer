import os 
import sys

''' input file taken from Praveen at the moment from location
/afs/cern.ch/work/p/ptiwari/public/bbDM/WCr_Split/AllMETHistos.root                                                                                                                 '''
inputfile = "AllMETHistos.root"

''' For which analysis this is being done '''
analysisName = "_bbDM"

''' For which year data this model is being run, by default it is set to 2016, but can be changed here. '''
yearStr        = "_2016"

''' any other information needed to explain the details by rootfile name, by default it is version track, '''
postfix     = "_V0";

''' output file with the workspace including all information about the model including inputs'''
outputfile  = "combine" + analysisName + yearStr + postfix + ".root"

''' name of the workspace ''' 
wsname = "ws" + analysisName

''' fitting range '''
met_low = 200
met_hi = 2000


