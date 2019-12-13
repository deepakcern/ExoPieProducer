import sys
import os

# datacard=sys.argv[1]
#
# med=datacard.split("_")[-3]
#
# mchi=datacard.split("_")[-1].split(".")[0]
#
# model=datacard.split("_")[-5]
#
#
# logfile = "log_"+med+"_"+mchi+".del"
# command_  = "combine -M Asymptotic --rAbsAcc 0 --rMax 30 -t -1 "+ datacard + " >> " + logfile
#
# os.system(command_)
model = '2HDMa'
outfile = 'limits_bbDM2016_'+model+'.txt'
fout = open(outfile,'w')
logfileList= ['Ma50_MChi1_MA600.txt','Ma100_MChi1_MA600.txt', 'Ma150_MChi1_MA600.txt', 'Ma250_MChi1_MA600.txt', 'Ma300_MChi1_MA600.txt', 'Ma350_MChi1_MA600.txt', 'Ma400_MChi1_MA600.txt', ]
for ifile in logfileList:
    samp_file = open('logfiles/'+ifile, 'r')
    med = ifile.split('_')[0].strip('Ma')
    mchi = ifile.split('_')[1].strip('MChi')
    expected25_=""
    expected16_=""
    expected50_=""
    expected84_=""
    expected975_=""
    observed_=""
    for ilongline in samp_file:
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

    towrite =  med+" "+mchi+" "+expected25_+" "+expected16_+" "+ expected50_+" "+ expected84_+" "+ expected975_+" "+ observed_+"\n"

    print towrite
    fout.write(towrite)
fout.close()
