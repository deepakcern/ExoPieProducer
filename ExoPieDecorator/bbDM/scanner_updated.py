import sys
import os

model = '2HDMa'

#logfileList_600 = [f for f in os.listdir('logfiles') if f.split('_')[-1]==str('MA600.txt')]
logfileList_600_1b = [f for f in os.listdir('logfiles') if (f.split('_')[-1]==str('MA600.txt') and (f.split('_')[-4]==str('1b')))]
logfileList_600_2b = [f for f in os.listdir('logfiles') if (f.split('_')[-1]==str('MA600.txt') and (f.split('_')[-4]==str('2b')))]
#logfileList_1200 = [f for f in os.listdir('logfiles') if f.split('_')[-1]==str('MA1200.txt')]
logfileList_1200_1b = [f for f in os.listdir('logfiles') if (f.split('_')[-1]==str('MA1200.txt') and (f.split('_')[-4]==str('1b')))]
logfileList_1200_2b = [f for f in os.listdir('logfiles') if (f.split('_')[-1]==str('MA1200.txt') and (f.split('_')[-4]==str('2b')))]
logfileList_600_1b.sort(key = lambda x: int(x.split('_')[-3].strip('Ma')))
logfileList_600_2b.sort(key = lambda x: int(x.split('_')[-3].strip('Ma')))
logfileList_1200_1b.sort(key = lambda x: int(x.split('_')[-3].strip('Ma')))
logfileList_1200_2b.sort(key = lambda x: int(x.split('_')[-3].strip('Ma')))

#outfile600 = 'limits_bbDM2016_'+model+'_MA600_150.txt'
#fout600 = open(outfile600,'w')
cat_dict = {tuple(logfileList_600_1b):'MA600_1b',tuple(logfileList_600_2b):'MA600_2b',tuple(logfileList_1200_1b):'MA1200_1b',tuple(logfileList_1200_2b):'MA1200_2b'}
for logfiles in [logfileList_600_1b,logfileList_600_2b,logfileList_1200_1b,logfileList_1200_2b]:
    outfile = 'limits_bbDM2016_'+model+'_'+cat_dict[tuple(logfiles)]+'.txt'
    fout = open(outfile,'w')
    for ifile in logfiles:
        samp_file = open('logfiles/'+ifile, 'r')
        med = ifile.split('_')[-3].strip('Ma')
        mchi = ifile.split('_')[-2].strip('MChi')
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
