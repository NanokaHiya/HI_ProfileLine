#main.py
#Generate HI profile line from different direction for subhalo in TNG simulation
#Haiyin Song, 2021/11
import sys
import os
import subhalo
import numpy as np
from datetime import datetime


if len(sys.argv)!=3:
    raise ValueError("Expect argv = 3. Please provide snapNum and subhalo id")
try:
    snapNum = int(sys.argv[1])
    subhalo_id = int(sys.argv[2])
except:
    print('argv value error')

print('Loading subhalo data...')
time0 = datetime.now()
mysubhalo = subhalo.subhalo(snapNum,subhalo_id)
time1 = datetime.now()
print('Loading finished, time elapsed is '+str(time1-time0))

print('Plotting...(Total 32)')
thetaArr = np.linspace(0,np.pi,num=4,endpoint=False)
phiArr = np.linspace(0,2*np.pi,num=8,endpoint=False)


def get32Pics(flag):
    #if flag == False, all 32 line will be in a same plot
    runTime = datetime.now()
    counter = 0
    for i in thetaArr:
        for j in phiArr:
            mysubhalo.drawPlot(i,j,flag)
            counter += 1
            print(str(counter) + ' of 32 finished. Time elapsed is ' + str(datetime.now()-runTime))
            runTime = datetime.now()
    if flag == True:
        print('All 32 fishined. Total time: ' + str(datetime.now()-time0))
    else:
        print('Single finished. Total time: ' + str(datetime.now()-time0))
    return

if __name__ == "__main__":
    get32Pics(True)
    get32Pics(False)
    print('Job finished')



