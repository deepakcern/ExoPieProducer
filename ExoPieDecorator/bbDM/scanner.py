import sys
import os

model = '2HDMa'

logfileList_600 = [f for f in os.listdir('logfiles150') if f.split('_')[2]==str('MA600.txt')]
logfileList_1200 = [f for f in os.listdir('logfiles150') if f.split('_')[2]==str('MA1200.txt')]
logfileList_600.sort(key = lambda x: int(x.split('_')[0].strip('Ma')))
logfileList_1200.sort(key = lambda x: int(x.split('_')[0].strip('Ma')))
outfile600 = 'limits_bbDM2016_'+model+'_MA600_150.txt'
fout600 = open(outfile600,'w')
for ifile in logfileList_600:
    samp_file = open('logfiles150/'+ifile, 'r')
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
    fout600.write(towrite)
fout600.close()
####################
outfile1200 = 'limits_bbDM2016_'+model+'_MA1200_150.txt'
fout1200 = open(outfile1200,'w')
for ifile in logfileList_1200:
    samp_file = open('logfiles150/'+ifile, 'r')
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
    fout1200.write(towrite)
fout1200.close()
