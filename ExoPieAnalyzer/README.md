## analyzer code here with systematics (using bbdm as base)

The analysis code is mainly written in the monoHbbAnalyzer.py file. This will create a .root file with trees for each CR/SR.


## To run interactive on a .txt file make following changes 

1. isCondor = False
2. runInteractive = False
3. year_file.write('era="2016"')  -- edit it accordingly 
4. keep the list of .root file in a .txt file (say ttskim.txt)
5. For 2017 and 18, change tight --> Tight in the zip statement 
   * st_nTau_discBased_looseEletightMuVeto --> st_nTau_discBased_looseEleTightMuVeto  and st_nTau_discBased_tightEletightMuVeto --> st_nTau_discBased_TightEleTightMuVeto
6. In config/varibale.py, change st_nTau_discBased_looseEletightMuVeto --> st_nTau_discBased_looseEleTightMuVeto  and st_nTau_discBased_tightEletightMuVeto --> st_nTau_discBased_TightEleTightMuVeto
7. run the code using: 
#### python bbMETAnalyzer.py -F -i ttskim.txt



### How to run
##for interactive run:

''' python bbMETAnalyzer.py -runOnTXT -F -inDir inputtextFileDir -D outputDir '''



The output tree can be used to make

1. the histogram and then make stack plots using the same framework used earlier.

2. use tree to make the stack plots using matplotlib (it is possible but slightly tricky to make divide canvas)

3. use tree to do the optimisation of one of the variable one might need. You have to remove cut on this variable and save in the rootfile.


For 1. i.e. make histogram use the script DataFrameToHisto_bbDM.py, which save histogram in a .root file.

### How to run

''' python DataFrameToHisto_bbDM.py -F -inDir Analyser_OutPut -D Output_Histograms '''

''' TreeSelector.py: it has the template to select the dataframe in one line. '''
