#main.py
#Generate HI profile line from different direction for subhalo in TNG simulation
#Haiyin Song, 2021/11
import sys
import os
import subhalo
import numpy as np



if len(sys.argv)!=3:
    raise ValueError("Expect argv = 3. Please provide snapNum and subhalo id")
try:
    snapNum = int(sys.argv[1])
    subhalo_id = int(sys.argv[2])
except:
    print('argv value error')

test = subhalo.subhalo(snapNum,subhalo_id)
test.drawPlot(np.pi,np.pi)
#test.LOSVelocity()
print(test.LOSVelocity)
print((test.LOSVelocity).shape)
#print((test.lower).shape)
