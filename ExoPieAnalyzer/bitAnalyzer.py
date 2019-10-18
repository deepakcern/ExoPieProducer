from ROOT import TH1F, TFile 

bits=[1,1,1,1,0,1,0]

x=1 ## always 1
y=0
z=0

for ibit in range(len(bits)):
    y = bits[ibit] << ibit
    z = z ^ y
    print bits[ibit], y, z
    
## now based on the value of z we can decide which cuts were passed. 
print "z  = ",z

print "binary output= ",bin(z)


def setBits(n):
    x = 1 ## always 1
    y = 0
    z = 0
    for i in range(n):
        y = x << i
        z = z ^ y
        
    return z

h = TH1F("h_cutflow", "h_cutflow", 10,0,10)        
for i in [1,2,3,4,5,6]:
    
    print i, setBits(i), setBits(i) & z 
    
    if setBits(i) == (setBits(i) & z) :
        print ("all selection upto ",i," passed")
        h.AddBinContent(i, 1)

f_out = TFile("test.root", "RECREATE")
f_out.cd()
h.Write()
f_out.Close()

        ## for those cuts, when this print statement is print on screen will be added in the cutflow histogram 
