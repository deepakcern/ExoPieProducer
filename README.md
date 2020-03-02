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
```
