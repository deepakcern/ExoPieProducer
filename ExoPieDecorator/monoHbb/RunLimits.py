import os 
import  sys 



import argparse
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
parser.add_argument("-b", "--both",  action="store_true", dest="both")
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




class RunLimits:
    ''' class to perform all tasks related to the limits once datacards are prepared '''
    ''' this class exepcts that all the steps needed to prepare the datacards and prepration of its inputs are already performed '''
    
    ''' instantiation of the class is done here ''' 
    def __init__(self):
        print "class instantiation done"
        
    ''' get the full command to be run for a given datacards '''
    def getfullcommand(self, commandpre, datacard, command_, commandpost):
        return commandpre+datacard+command_+commandpost
        
    def makedatacards(self, templatecards, MA):
        
        datacardsname = datacardtemplatename_.replace("XXXMA", str(MA))
        os.system('rm '+datacardsname)
        fout = open(datacardsname,"a")
        for iline in open(templatecards):
            iline  = iline.replace("XXXMA", str(MA))
            #iline = iline.replace("monoHbb2017_B","monoHbb2017_R")
            fout.write(iline)
        fout.close()
        return datacardsname



## main code, 

def main():
    
    print "inside main"
    rl = RunLimits()
    
    ''' following is the syntax to get all the cards using template datacard ''' 
    
    
    MA=[300, 400, 500, 600, 1000, 1200, 1600]
    if args.createdatacards:
        datacardtextfile = 'monohbb2017_datacardslist.txt'
        os.system('rm '+datacardtextfile)
        ftxt = open(datacardtextfile,'w')
        for iMA in MA:
            datacardname = rl.makedatacards('datacards_tmplate/combine_tmpl_sig2b_workspace.txt',iMA)
            print datacardname 
            ftxt.write(datacardname+' \n')
        ftxt.close()
        
        
    ''' following is the syntax to run all limits ''' 
    ''' datacards path in a text file is converted into a list''' 
    datacardnameslist=[]
    if args.runlimits:
        datacardnameslist = [iline.rstrip() for iline in open(args.inputdatacardpath)]
        datacardtemplatename_ = 'datacards_monoHbb_2017/datacard_monoHbb2017_B_SR_ggF_sp_0p35_tb_1p0_mXd_10_mA_XXXMA_ma_150.txt'

        for idc in datacardnameslist :
            print rl.getfullcommand(commandpre, idc, command_, commandpost)
            os.system(rl.getfullcommand(commandpre, idc, command_, commandpost))
            print "-----------------------------------------------------------------------------------------------------------------------"
        
if __name__ == "__main__":
    
    print "calling main"
    main()
    
    
    
    
