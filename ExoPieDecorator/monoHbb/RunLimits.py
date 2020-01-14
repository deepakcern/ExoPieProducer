import os 
import  sys 



import argparse

from LimitHelper import RunLimits
#import params as parameters


usage = "run the script using python -i full/path/to/root/file "
parser = argparse.ArgumentParser(description=usage)

## string 
parser.add_argument("-i", "--inputdatacardpath",  dest="inputdatacardpath",default="monohbb2017_datacardslist.txt") ## this should be a .txt file which include the full path of the datacards

parser.add_argument("-model", "--model",  dest="model",default="params_2hdma.txt") 
parser.add_argument("-region", "--region",  dest="region",default="SR_default") ## default should lead to crash.  


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




## main code, 
def main():
    
    print "inside main"
    
    ## resolved/merged/combined analysis 
    analysis_tag = '_default' ## default is random, it will crash
    if args.merged: analysis_tag = '_B'
    if args.resolved: analysis_tag = '_R'
    if args.combined: analysis_tag = '_C'
    
    
    ## region tag: SR, ZEE, TOPE, WE, ZMUMU, TOPMU, WMU
    #if args.region == : region_tag = 'SR'
    
    regions = args.region.split(" ")
    print "region list: ", regions
    
    
    ## this can be different for each model. But for now lets keep it like this. 
    datacardtemplatename_ = 'datacards_monoHbb_2017/datacard_monoHbb2017'+analysis_tag+'_SR_ggF_sp_YYYSP_tb_ZZZTB_mXd_AAAMDM_mA_XXXMA_ma_BBBMa.txt'
    
    ## object of the RunLimits class
    rl = RunLimits(datacardtemplatename_)
    
    
    if args.createdatacards:
        fparam = open("params_"+args.model+".txt","r") 
        datacardtextfile = args.inputdatacardpath.replace(".txt", "_"+args.model+".txt")
        os.system('rm '+datacardtextfile)
        ftxt = open(datacardtextfile,'w')
        for iparam in fparam:
            ## collect datacard name from the function after preparing it. 
            ## rl.makedatacards should do in future: make all the SR and CR datacards and merge the datacards and return the merged datacard. 
            ## it should take the option of which CRs to use.
            ## to use boosted or resolved 
            datacards=[]
            for iregion in regions:
                datacardname = rl.makedatacards('datacards_tmplate/combine_tmpl_'+iregion+'_workspace.txt',iparam.split(), iregion)
                if iregion == "SR":
                    mergeddatacardmname = datacardname.replace("SR_ggF","Merged")
                    ftxt.write(mergeddatacardmname+' \n')
                datacards.append(" "+iregion+"="+datacardname)
            
            combostr = ""
            if len(datacards)>0:
                for idc in datacards:
                    combostr = combostr + idc
            print "datacards will me creating using categories: ",regions, " parameters: ", iparam.rstrip(), " datacard: ", mergeddatacardmname
            
            #ftxt.write("combineCards.py "+combostr+" > "+mergeddatacardmname+"\n")
            os.system("combineCards.py "+combostr+" > "+mergeddatacardmname+"\n")
            ## instead of datacard name write the combostr into .text file to make the combined card
        ftxt.close()
        


    
    ''' following is the syntax to run all limits ''' 
    ''' datacards path in a text file is converted into a list''' 
    datacardnameslist=[]
    if args.runlimits:
        datacardnameslist = [iline.rstrip() for iline in open(args.inputdatacardpath)]
        for idatacard in datacardnameslist :
            print rl.getfullcommand(commandpre, idatacard, command_, commandpost)
            
            logfilename = "logs/"+idatacard.replace(".txt",".log")
            
            os.system(rl.getfullcommand(commandpre, idatacard, command_, commandpost) + " > "+logfilename)
            
            rl.LogToLimitList(logfilename)
            print "-----------------------------------------------------------------------------------------------------------------------"
        
if __name__ == "__main__":
    
    print "calling main"
    main()
    
    
    
    
