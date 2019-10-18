bits=[1,1,1,1,0,1,0]

x=1 ## always 1
y=0
z=0

for ibit in range(len(bits)):
    y = bits[ibit] << ibit
    z = z ^ y
    print bits[ibit], y, z

    ## now based on the value of z we can decide which cuts were passed. 
    #print ibit
