from itertools import izip_longest
import os,sys,datetime
from glob import glob


def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def splitTxt(txtfile,dirName,maxFilePerTxt):
    n=maxFilePerTxt
    with open(txtfile) as f:
        newtxt=txtfile.split('/')[-1].replace('.txt','')
        for i, g in enumerate(grouper(n, f, fillvalue=''), 1):
            with open(dirName+'/'+newtxt+'_{0}.txt'.format(i), 'w') as fout:
                fout.writelines(g)

def submitjob(count,txtfile):
    #global count
    submittemp=open("submit_multi_temp.sub","w")
    submitfile=open("submit_multi.sub","r")
    for line in submitfile:
        if line.startswith('transfer_input_files'):
            submittemp.write(line.strip()+', '+txtfile+'\n')
        else:
            submittemp.write(line)
    submitfile.close()
    dummy='dummy'
    outfineName=txtfile.split('/')[-1].replace('.txt','.root')
    submittemp.write("arguments = "+txtfile.split('/')[-1]+" "+outfineName+" "+"Analysis_"+outfineName+"  "+"Output_Analysis_"+outfineName+'\nqueue')
    submittemp.close()
    print "\n===============================\nSubmitting jobs #"+str(count)+": "+ txtfile+"\n===============================\n"
    if not test: os.system("condor_submit submit_multi_temp.sub")
    
if __name__== "__main__":
    test=False
    count=0
    FilesToSubmit = "Filelists_test"
    FilesToResubmit = "Filelists_failed"

    if not test: os.system("chmod +x runAnalysis.sh")
    for outdirs in ['error','log','output']:
        os.system("mkdir -p "+outdirs)

    IsSkimJobs =False;submitJob=False;ResubmitJob=False
    if sys.argv[1]=="skim":
	IsSkimJobs=True
    if sys.argv[2]=="submit":
        submitJob=True
    elif sys.argv[2]=="resubmit":
        ResubmitJob=True
    else:
        print "Please provide correct arguments, check the code for usage+\n"
        sys.exit()

    maxfilesperjob=150
    if IsSkimJobs and not ResubmitJob: # This part is to submit skim jobs with multiple root file per job
        listfiles = [f for f in os.listdir(FilesToSubmit) if f.endswith('.txt')]
        dirName='tempFilelists_'+datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
	os.system("mkdir -p "+dirName)
	for txtfile in listfiles:
	    splitTxt(FilesToSubmit+'/'+txtfile,dirName,maxfilesperjob)
    elif ResubmitJob:
	#This part is to resubmit failed skim jobs, provide same input txt files which are failed,
	#copy those input txt files from tempFilelists* directory
	dirName=FilesToResubmit
    elif not IsSkimJobs: # this part is for analyser job: job with single root file. Each txt file is with single root file
	dirName=FilesToSubmit
    MytxtFiles = [f for f in os.listdir(dirName) if f.endswith('.txt')]

for ifile in MytxtFiles:
    submitjob(count,dirName+'/'+ifile)
    count+=1

print 'Total number of jobs: ',count
