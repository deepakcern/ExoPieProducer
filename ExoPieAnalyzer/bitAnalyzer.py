bits=[0,1,2,3,5]

x=1 ## always 1
y=0
z=0

for ibit in bits:
    y = x << ibit
    z = z ^ y
    print ibit, y, z

    ## now based on the value of z we can decide which cuts were passed. 
    #print ibit
