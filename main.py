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

print('Plotting...(Total 26)')
thetaArr = np.linspace(0,np.pi,num=5,endpoint=True)
phiArr = np.linspace(0,2*np.pi,num=8,endpoint=False)


def get26Pics(flag):
    #if flag == False, all 32 line will be in a same plot
    runTime = datetime.now()
    counter = 0
    for i in thetaArr:
        for j in phiArr:
            counter += 1
            print(str(counter)+ ' of 26: theta = '+str("{:.2f}".format(i))+', phi = '+str("{:.2f}".format(j)) )
            mysubhalo.drawPlot(i,j,flag)
            #Time elapsed is ' + str(datetime.now()-runTime))
            runTime = datetime.now()
            if (i == 0)or(i == np.pi):
                break
    '''
    if flag == True:
        print('All 32 fishined. Total time: ' + str(datetime.now()-time0))
    else:
        print('Single finished. Total time: ' + str(datetime.now()-time0))
    '''
    return

if __name__ == "__main__":
    get26Pics(True)
    #get32Pics(False)
    #mysubhalo.drawPlot(0,0, center = False)
    print('Job finished')



