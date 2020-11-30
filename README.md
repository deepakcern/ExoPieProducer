# ExoPieProducer
Set of py modules to use ExoPieSlimmer output and produce results for AN/PAS/Papers

Directory structure: 

ExoPieProducer/ExoPieUtils: common utilities 

ExoPieProducer/ExoPieAnalyzer: analysis code to perform the selection and make tree for making plots

ExoPieProducer/ExoPieCapper: make all kind of plots needed for the AN/PAS/Slides/Paper 

ExoPieProducer/ExoPieDecorator: Limit model related stuff. 

## Setup framework 

```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_3_0
cd CMSSW_10_3_0/src
cmsenv
git clone https://github.com/ExoPie/ExoPieUtils.git
cd ExoPieUtils
git checkout test_systematics
cd -
git clone https://github.com/ExoPie/ExoPieProducer.git
cd ExoPieProducer
git checkout monohbb
cd -
```

### Interactive Run
Note: change ```isCondor = False```  inside `monoHbbAnalyzer.py`
```
cd ExoPieProducer/ExoPieAnalyzer
python monoHbbAnalyzer.py -y 2017 -F -i pathOfInputTxtFile
```
txt files you can find filelists (`Filelists_v07.04_noJER_updatedFilter_v2_2017`) here:

`https://github.com/deepakcern/ExoPieFileList/tree/master/ExoPieSlimmerFiles`

Important files to understand code:

1- `monoHbbAnalyzer.py` main file

2- `configs/eventSelector.py` event selection file

3- `configs/cutFlow.py` for cutflow

4-  scalefactors and other important functions are taken from ExopieUtils repo


### Submit Condor Jobs
Note: change ```isCondor = True```
```
cd ExoPieProducer/ExoPieAnalyzer
git clone https://github.com/deepakcern/CondorJobs.git
cd CondorJobs
mv AnalysisJobs/* .
. submitjobs.sh
```

Note: Open `MultiSubmit.py` and provide directory of Filelists name where all txt files of sample are saved. Make a directory and copy all txt files.
Fielists directory with txt files is already there.

please change output path in the file ```runAnalysis.sh``` before submitting condor jobs

To get proxy file, open ```.bashrc``` file and add:
```
alias voms='voms-proxy-init --voms cms --valid 192:00 && cp -v /tmp/x509up_u104803 /afs/cern.ch/user/d/dekumar/private/x509up'
```
change username. 

Once all jobs are done then merge output files to get single file per sample. You can use this shell script to merge the files.
`/eos/cms/store/group/phys_exotica/monoHiggs/monoHbb/2017_AnalyserOutput/mergeFiles_2017.sh` 

### Writting Histograms from Trees

#### For boosted category
```
cd ExoPieProducer/ExoPieAnalyzer
python DataFrameToHisto_B_syst.py -F -inDir pathOfAnalyserRootFilesOutput -D OutputDirectory
```
#### For resolved category
```
cd ExoPieProducer/ExoPieAnalyzer
python DataFrameToHisto_R_syst.py -F -inDir pathOfAnalyserRootFilesOutput -D OutputDirectory
```

#### For signale samples
```
cd ExoPieProducer/ExoPieAnalyzer
python DataFrameToHisto_signal.py -F -inDir pathOfAnalyserRootFilesOutput -D OutputDirectory
```


### Making Control region plots

```
cd ExoPieProducer/ExoPieCapper/

python StackPlotter_2017_syst.py -c B -d MET -s [muon region plots for boosted analysis, use R for resolved category]
python StackPlotter_2017_syst.py -c B -d MET -b [signal region]
python StackPlotter_2017_syst.py -c B -d MET -m [electron region]
python StackPlotter_2017_syst.py -c B -d MET -m [electron region]
```
change path of inputroot file inside ``` StackPlotter_2017_syst.py ``` file. Give path of DataFrameToHisto_B_syst.py or DataFrameToHisto_R_syst.py output root files

### Making of AllMETHistos.root file [input file of limit model]

```
cd ExoPieProducer/ExoPieDecorator
python monoH_combinedroot_v2.py 
```
give 3 paths inside the file:

1-path of stack plotter output root files for boosted category

2-path of stact plotter output root files for resolved ctegory

3-output of DataframeTohist_Signal.py for signal sample

Now run the command:
```
python monoH_combinedroot_v2.py
```
