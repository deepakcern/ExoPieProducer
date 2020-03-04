# ExoPieProducer
Set of py modules to use ExoPieSlimmer output and produce results for AN/PAS/Papers

Directory structure: 

ExoPieProducer/ExoPieUtils: common utilities 

ExoPieProducer/ExoPieAnalyzer: analysis code to perform the selection and make tree for making plots

ExoPieProducer/ExoPieCapper: make all kind of plots needed for the AN/PAS/Slides/Paper 

ExoPieProducer/ExoPieDecorator: Limit model related stuff. 

## Setup framework 

```
git clone git@github.com:deepakcern/ExoPieUtils.git
git checkout monohbb_systematics

git clone git@github.com:deepakcern/ExoPieProducer.git
git checkout monohbb
```

### Interactive Run
Note: keep ```isCondor = False``` 
```
cd ExoPieProducer/ExoPieAnalyzer
python monoHbbAnalyzer.py -F -i pathOfInputTxtFile
```

### Submit Condor Jobs
```
cd ExoPieProducer/ExoPieAnalyzer
git clone git@github.com:deepakcern/CondorJobs.git
cd CondorJobs
. Analyser_submitjobs.sh
```
Note: Open `submit_AnalyserJobs.py` and provide directory of Filelists name where all txt files of sample are saved. Make a directory and copy all txt files.

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

