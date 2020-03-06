# ExoPieProducer
Set of py modules to use ExoPieSlimmer output and produce results for AN/PAS/Papers

Directory structure: 

ExoPieProducer/ExoPieUtils: common utilities 

ExoPieProducer/ExoPieAnalyzer: analysis code to perform the selection and make tree for making plots

ExoPieProducer/ExoPieCapper: make all kind of plots needed for the AN/PAS/Slides/Paper 

ExoPieProducer/ExoPieDecorator: Limit model related stuff. 

## Setup framework 

```
cmsrel CMSSW_10_3_0
cd CMSSW_10_3_0/src
cmsenv
git clone https://github.com/deepakcern/ExoPieUtils.git
cd ExoPieUtils
git checkout monohbb_systematics
cd -
git clone https://github.com/deepakcern/ExoPieProducer.git
cd ExoPieProducer
git checkout monohbb
cd -
```

### Interactive Run
Note: change ```isCondor = False```  inside `monoHbbAnalyzer.py`
```
cd ExoPieProducer/ExoPieAnalyzer
python monoHbbAnalyzer.py -F -i pathOfInputTxtFile
```
txt files you can find here: ```https://github.com/deepakcern/CondorJobs/tree/master/Filelists```

### Submit Condor Jobs
Note: change ```isCondor = True```
```
cd ExoPieProducer/ExoPieAnalyzer
git clone git@github.com:deepakcern/CondorJobs.git
cd CondorJobs
. submitjobs_step2.sh
```
Note: Open `MultiSubmit_step2.py` and provide directory of Filelists name where all txt files of sample are saved. Make a directory and copy all txt files.
Fielists directory is already there you can use root files for now. Fielists will be updated soon for 2017/18

please change output path in the file ```runAnalysis_step2.sh```

#### For condor jobs on lxplus:
if you want to submit condor jobs on lxplus, please update three files ``` submit_multi_step2.sub``` , ```MultiSubmit_step2.py``` and ```runAnalysis_step2.sh```

add following lines in ```submit_multi_step2.sub```
```
Proxy_filename = x509up
Proxy_path = /afs/cern.ch/user/d/dekumar/private/$(Proxy_filename)
request_cpus = 4
+JobFlavour = "nextweek"
```

to get proxy file, open ```.bashrc``` file and add:
```
alias voms='voms-proxy-init --voms cms --valid 192:00 && cp -v /tmp/x509up_u104803 /afs/cern.ch/user/d/dekumar/private/x509up'
```
change username. 


open `MultiSubmit_step2.py` file add this string `$(Proxy_path) ` at last in line 19 for 5th arguments

```
export X509_USER_PROXY=$5
voms-proxy-info -all
voms-proxy-info -all -file $5
```
### Writting Histograms from Trees

```
cd ExoPieProducer/ExoPieAnalyzer
python DataFrameToHisto.py -F -inDir pathOfAnalyserRootFilesOutput -D OutputDirectory
```
combined data files into single file
```
cd OutputDirectory
hadd combined_data_SE.root SingleElectron-Run2017*.root
hadd combined_data_MET.root MET-Run2017*.root
```
### Making Control region plots

```
cd ExoPieProducer/ExoPieAnalyzer
wget https://raw.githubusercontent.com/deepakcern/ExoAnalysis/master/monoH/plottingTools/StackPlotter_2017_syst.py
wget https://raw.githubusercontent.com/deepakcern/ExoAnalysis/master/monoH/plottingTools/sample_xsec_2017.py
wget https://raw.githubusercontent.com/deepakcern/ExoAnalysis/master/monoH/plottingTools/samplelist_2017.txt

python StackPlotter_2017_syst.py -c B -d MET -m [muon region plots for boosted analysis]
python StackPlotter_2017_syst.py -c B -d MET -s [signal region]
python StackPlotter_2017_syst.py -c B -d SE -e [electron region]
```
change path of inputroot file inside ``` StackPlotter_2017_syst.py ``` file

