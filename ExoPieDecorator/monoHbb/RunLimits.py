import os 
import  sys 



import argparse

from LimitHelper import RunLimits
#import params as parameters


usage = "run the script using python -i full/path/to/root/file "
parser = argparse.ArgumentParser(description=usage)

## string 
parser.add_argument("-i", "--inputdatacardpath",  dest="inputdatacardpath",default="monohbb2017_datacardslist.txt") ## this should be a .txt file which include the full path of the datacards


## booleans 
parser.add_argument("-B", "--runblind",  action="store_true", dest="runblind")
parser.add_argument("-A", "--runasimov",  action="store_true", dest="runasimov")
parser.add_argument("-L", "--runlimits",  action="store_true", dest="runlimits")
parser.add_argument("-D", "--rundiagonstics",  action="store_true", dest="rundiagonstics")

parser.add_argument("-c", "--createdatacards",  action="store_true", dest="createdatacards")

parser.add_argument("-m", "--merged",  action="store_true", dest="merged")
parser.add_argument("-r", "--resolved",  action="store_true", dest="resolved")
parser.add_argument("-C", "--combined",  action="store_true", dest="combined")


## integers 
parser.add_argument("-v", "--verbose",  dest="verbose", type=int, default=0)
parser.add_argument("-rmax", "--rmax",  dest="rmax", type=int, default=30)
parser.add_argument("-rmin", "--rmin",  dest="rmin", type=int, default=0.0000001)
parser.add_argument("-CL", "--CL",  dest="CL", type=int, default=0.95) ## can be used for SI interpretation but not using right now 

args = parser.parse_args()




''' all the defaults needed for rest of the class, and to execute the opetations are here ''' 
''' many of them are coming from the command line argumnet''' 

commandlist = []
commandpre=""

## add verbose everytime 
commandpost = " -v " + str(args.verbose) +  " --rMin "+str(args.rmin) + " --rMax "+str(args.rmax) + " "

''' collect all the options here ''' 
if args.runasimov: commandlist.append(" --noFitAsimov ")
if args.runblind: commandlist.append("  --run blind ")
command_=""
for ic in commandlist:  command_ = command_ + " "




''' which method to run '''
if args.runlimits and args.rundiagonstics:
    print "----- this is confusing, please keep only one of them true. only one can be run at one time. "
if args.runlimits: commandpre = "combine -M AsymptoticLimits "
if args.rundiagonstics: commandpre = "combine -M FitDiagnostics "

datacardtemplatename_ = 'datacards_monoHbb_2017/datacard_monoHbb2017_R_SR_ggF_sp_YYYSP_tb_ZZZTB_mXd_AAAMDM_mA_XXXMA_ma_BBBMa.txt'

####
####
####class RunLimits:
####    ''' class to perform all tasks related to the limits once datacards are prepared '''
####    ''' this class exepcts that all the steps needed to prepare the datacards and prepration of its inputs are already performed '''
####    
####    ''' instantiation of the class is done here ''' 
####    def __init__(self):
####        print "class instantiation done"
####        
####    ''' get the full command to be run for a given datacards '''
####    def getfullcommand(self, commandpre, datacard, command_, commandpost):
####        return commandpre+datacard+command_+commandpost
####        
####    def makedatacards(self, templatecards, MA):
####        
####        datacardsname = datacardtemplatename_.replace("XXXMA", str(MA))
####        os.system('rm '+datacardsname)
####        fout = open(datacardsname,"a")
####        for iline in open(templatecards):
####            iline  = iline.replace("XXXMA", str(MA))
####            #iline = iline.replace("monoHbb2017_B","monoHbb2017_R")
####            fout.write(iline)
####        fout.close()
####        return datacardsname
####        
####
####    def datacard_to_mparameters(self, name_):
####        mparameters_ = ((name_.split("SR_ggF_")[1]).replace(".log","")).split("_")
####        mparameters_ = [mp.replace("p",".") for mp in mparameters_]
####        ## ma, mA, tb, st, mdm
####        return ([mparameters_[9], mparameters_[7], mparameters_[3], mparameters_[1], mparameters_[5]])
####
####    def LogToLimitList(self, logfile):
####        expected25_="" 
####        expected16_="" 
####        expected50_="" 
####        expected84_="" 
####        expected975_=""
####        observed_=""
####        for ilongline in open(logfile):
####            if "Observed Limit: r < " in ilongline:
####                observed_ = ilongline.replace("Observed Limit: r < ","").rstrip()
####            if "Expected  2.5%: r < " in ilongline:
####                expected25_ = ilongline.replace("Expected  2.5%: r < ","").rstrip()
####            if "Expected 16.0%: r < " in ilongline:
####                expected16_ = ilongline.replace("Expected 16.0%: r < ","").rstrip()
####            if "Expected 50.0%: r < " in ilongline:
####                expected50_ = ilongline.replace("Expected 50.0%: r < ","").rstrip()
####            if "Expected 84.0%: r < " in ilongline:
####                expected84_ = ilongline.replace("Expected 84.0%: r < ","").rstrip()
####            if "Expected 97.5%: r < " in ilongline:
####                expected975_ = ilongline.replace("Expected 97.5%: r < ","").rstrip()
####        
####        allparameters  = self.datacard_to_mparameters(logfile)
####        towrite =  str(allparameters[0])+" "+str(allparameters[1])+" "+expected25_+" "+expected16_+" "+ expected50_+" "+ expected84_+" "+ expected975_+" "+ observed_+"\n"
####        
####        print towrite
####        outfile=""
####        if args.merged: outfile = 'bin/limits_monoH_B_2017.txt'
####        if args.resolved: outfile = 'bin/limits_monoH_R_2017.txt'
####        if args.combined: outfile = 'bin/limits_monoH_Combo_2017.txt'
####        
####        fout = open(outfile,'a')
####        fout.write(towrite)
####        fout.close()
####
####
## main code, 

def main():
    
    print "inside main"
    rl = RunLimits(datacardtemplatename_)
    
    ''' following is the syntax to get all the cards using template datacard ''' 
    
    fparam = open("params.txt","r")
    
    #MA=parameters.mA #[300, 400, 500, 600, 1000, 1200, 1600]
    if args.createdatacards:
        datacardtextfile = 'monohbb2017_datacardslist.txt'
        os.system('rm '+datacardtextfile)
        ftxt = open(datacardtextfile,'w')
        for iparam in fparam:
            ##allparams = iparam.split()
            ##ma =allparams[0]
            ##mA =allparams[1]
            ##tb =allparams[2]
            ##st =allparams[3]
            ##mdm=allparams[4]
            
            datacardname = rl.makedatacards('datacards_tmplate/combine_tmpl_sig2b_workspace.txt',iparam.split())
            print datacardname 
            ftxt.write(datacardname+' \n')
        ftxt.close()
        
        
    ''' following is the syntax to run all limits ''' 
    ''' datacards path in a text file is converted into a list''' 
    datacardnameslist=[]
    if args.runlimits:
        datacardnameslist = [iline.rstrip() for iline in open(args.inputdatacardpath)]
        
        for idc in datacardnameslist :
            print rl.getfullcommand(commandpre, idc, command_, commandpost)
            
            logfilename = "logs/"+idc.replace(".txt",".log")
            
            os.system(rl.getfullcommand(commandpre, idc, command_, commandpost) + " > "+logfilename)
            
            rl.LogToLimitList(logfilename)
            print "-----------------------------------------------------------------------------------------------------------------------"
        
if __name__ == "__main__":
    
    print "calling main"
    main()
    
    
    
    
