import os 
import  sys 



import argparse


class RunLimits:
    ''' class to perform all tasks related to the limits once datacards are prepared '''
    ''' this class exepcts that all the steps needed to prepare the datacards and prepration of its inputs are already performed '''
    
    ''' instantiation of the class is done here ''' 
    def __init__(self, datacardtemplatename):
        self.datacardtemplatename_ = datacardtemplatename
        print "class instantiation done"
        
    ''' get the full command to be run for a given datacards '''
    def getfullcommand(self, commandpre, datacard, command_, commandpost):
        return commandpre+datacard+command_+commandpost
        
    
    def makedatacards(self, templatecards, allparams):
        
        ma =str(allparams[0])
        mA =str(allparams[1])
        tb =(str(allparams[2])).replace(".","p")
        st =(str(allparams[3])).replace(".","p")
        mdm=str(allparams[4])
        
        datacardsname = self.datacardtemplatename_.replace("XXXMA", mA)
        datacardsname = datacardsname.replace("BBBMa",ma)
        datacardsname = datacardsname.replace("ZZZTB",tb)
        datacardsname = datacardsname.replace("YYYSP",st)
        datacardsname = datacardsname.replace("AAAMDM",mdm)
        
        os.system('rm '+datacardsname)
        fout = open(datacardsname,"a")
        for iline in open(templatecards):
            iline  = iline.replace("XXXMA", mA)
            ## add other params 
            iline = iline.replace("monoHbb2017_B","monoHbb2017_R")
            fout.write(iline)
        fout.close()
        return datacardsname
        

    def datacard_to_mparameters(self, name_):
        mparameters_ = ((name_.split("SR_ggF_")[1]).replace(".log","")).split("_")
        mparameters_ = [mp.replace("p",".") for mp in mparameters_]
        ## ma, mA, tb, st, mdm
        return ([mparameters_[9], mparameters_[7], mparameters_[3], mparameters_[1], mparameters_[5]])

    def LogToLimitList(self, logfile):
        expected25_="" 
        expected16_="" 
        expected50_="" 
        expected84_="" 
        expected975_=""
        observed_=""
        for ilongline in open(logfile):
            if "Observed Limit: r < " in ilongline:
                observed_ = ilongline.replace("Observed Limit: r < ","").rstrip()
            if "Expected  2.5%: r < " in ilongline:
                expected25_ = ilongline.replace("Expected  2.5%: r < ","").rstrip()
            if "Expected 16.0%: r < " in ilongline:
                expected16_ = ilongline.replace("Expected 16.0%: r < ","").rstrip()
            if "Expected 50.0%: r < " in ilongline:
                expected50_ = ilongline.replace("Expected 50.0%: r < ","").rstrip()
            if "Expected 84.0%: r < " in ilongline:
                expected84_ = ilongline.replace("Expected 84.0%: r < ","").rstrip()
            if "Expected 97.5%: r < " in ilongline:
                expected975_ = ilongline.replace("Expected 97.5%: r < ","").rstrip()
        
        allparameters  = self.datacard_to_mparameters(logfile)
        towrite =  str(allparameters[0])+" "+str(allparameters[1])+" "+expected25_+" "+expected16_+" "+ expected50_+" "+ expected84_+" "+ expected975_+" "+ observed_+"\n"
        
        print towrite
        outfile="bin/limits_monoH_R_2017.txt"
        #if args.merged: outfile = 'bin/limits_monoH_B_2017.txt'
        #if args.resolved: outfile = 'bin/limits_monoH_R_2017.txt'
        #if args.combined: outfile = 'bin/limits_monoH_Combo_2017.txt'
        
        fout = open(outfile,'a')
        fout.write(towrite)
        fout.close()


