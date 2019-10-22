import ROOT as r
import array, sys

from optparse import OptionParser

parser = OptionParser()
parser.add_option("","--skip",default=False,action='store_true',help="Ignore missing samples rather than failing")
(options,args) = parser.parse_args()


# Configurations Read in from Separate .py files
sys.path.append("configs")
configuration = args[0]
x = __import__(configuration) 

# book category should read list of samples and append them to as histograms 
# expect formats of 
# Region : Process : contributing MC samples
# 	   Process : contributing MC samples
# 	   ...
r.gROOT.SetBatch(1)
r.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ");
r.gSystem.AddIncludePath("-I$ROOFITSYS/include");
r.gSystem.Load("libRooFit.so")
r.gSystem.Load("libRooFitCore.so")
# All c++ functionalities
r.gROOT.ProcessLine('.L ./ModelBuilder_with_stat.cc+')
fout = r.TFile(x.out_file_name,'RECREATE')

# Loop and build components for categories
for cat_id,cat in enumerate(x.categories):
  
  fin  = r.TFile.Open(cat['in_file_name'])
  fout.cd(); 
  fdir = fout.mkdir("category_%s"%cat['name'])
  mb = r.ModelBuilder(cat_id,cat['name'])
  
  mb.fIn  = fin
  mb.fOut = fdir
  
  mb.cutstring  = cat['cutstring']
  mb.setvariable(cat['varstring'][0],cat['varstring'][1],cat['varstring'][2])
  mb.setweight(cat['weightname'])
  mb._pdfmodel=cat['pdfmodel']
  
  for avar in cat['additionalvars']: mb.addvariable(avar[0],avar[1],avar[2],avar[3])

  # create a template histogram from bins
  bins = cat["bins"]
  histo_base = r.TH1F("base_%d"%cat_id,"base"
        ,len(bins)-1
	,array.array('d',bins))

  mb.lTmp = histo_base.Clone()
  
  if "extra_cuts" in cat.keys():
   for ecut in cat["extra_cuts"]: mb.add_cut(ecut[0],ecut[1])

  # Run through regions and add MC/data processes for each 
  # Each region has 'signal' and 'backgrounds'
  samples = cat['samples'].keys()
  for sample in samples:
      entry = cat['samples'][sample]
      print sample
      if not options.skip:  mb.addSample(sample,entry[0],entry[1],entry[2],entry[3])  # name, region, process, is_mc, is_signal
      else : 
        try: 
	     sampentry = fin.Get(sample)
 	     sampentry.GetEntries()
	     mb.addSample(sample,entry[0],entry[1],entry[2],entry[3])  # name, region, process, is_mc, is_signal
        except : print " No Tree %s found, skipping "%sample

      try: 
        systematics = x.systematics
	for syst in systematics: 
	  extU = "_%sUp"%syst
	  extD = "_%sDown"%syst
          print sample+extU, fin.Get(sample+extU)
	  if fin.Get(sample+extU):
#	    print " Adding systematic variation %s for %s "%(syst,sample)
            mb.addSample(sample+extU,entry[0],entry[1]+extU,entry[2],entry[3],1)  # name, region, process, is_mc, is_signal
#	  else: print " No %s Up variation found for %s"%(syst,sample)
	  if fin.Get(sample+extD):
#	    print " Adding systematic variation %s for %s "%(syst,sample)
            mb.addSample(sample+extD,entry[0],entry[1]+extD,entry[2],entry[3],1)  # name, region, process, is_mc, is_signal
#	  else: print " No %s Down variation found for %s"%(syst,sample)
      except : 
          noSys = True
	  
  
  
  mb.save()
          
  print mb.cutstring

  # Add any 'cutstring' for future reference
  cstr = r.TNamed("cut_category_%s"%cat['name'],cat["cutstring"])
  fdir.cd(); cstr.Write()

# finally add the config used into the file 
config = r.TMacro("%s"%(args[0]))
config.ReadFile("configs/%s.py"%(args[0]))
fout.cd(); config.Write()

print "done!, Model saved in -> ", fout.GetName()

