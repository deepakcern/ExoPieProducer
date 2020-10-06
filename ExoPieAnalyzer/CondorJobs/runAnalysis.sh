#!/bin/sh
#### FRAMEWORK SANDBOX SETUP ####
# Load cmssw_setup function
export SCRAM_ARCH=slc6_amd64_gcc700
source ./cmssw_setup.sh

# Setup CMSSW Base
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# Setup framework from sandbox
scramv1 project CMSSW_10_3_0
cd CMSSW_10_3_0/src
eval `scramv1 runtime -sh`
cd ../../

python SkimTree.py -y 2017 -F -i "$1"
python bbMETAnalyzer.py -y 2017 -I -i "$2"
python DataFrameToHisto_bbDM_syst.py -i "$3"

if [ -e "$2" ]; then
  until xrdcp -f "$2" root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/bbMET/2017_SkimmedFiles/skim_v17_07-02/"$2"; do
    sleep 60
    echo "Retrying"
  done
fi

if [ -e "$3" ]; then
  until xrdcp -f "$3" root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/bbMET/2017_AnalysedFiles/analyser_v17_07_02_00/"$3"; do
    sleep 60
    echo "Retrying"
  done
fi

if [ -e "$4" ]; then
  until xrdcp -f "$4" root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/bbMET/2017_AnalysedHistos/df_Output_v17_07_02_00/"$4"; do
    sleep 60
    echo "Retrying"
  done
fi

exitcode=$?

if [ ! -e "$4" ]; then
  echo "Error: The python script failed, could not create the output file."

fi
exit $exitcode
