#!/usr/bin/env python

#main.py
#Generate HI profile line from different direction for subhalo in TNG simulation
#Haiyin Song, 2021/11
import sys
import os
import numpy as np
from datetime import datetime

import rotate_star
import subhalo


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


def get26Pics(checkSymmetry):
    #if checkSymmetry == False, all 32 line will be in a same plot
    runTime = datetime.now()
    counter = 0
    for i in thetaArr:
        for j in phiArr:
            counter += 1
            print(str(counter)+ ' of 26: theta = '+str("{:.2f}".format(i))+', phi = '+str("{:.2f}".format(j)) )
            mysubhalo.drawPlot(i,j,clean = checkSymmetry)
            #Time elapsed is ' + str(datetime.now()-runTime))
            runTime = datetime.now()
            if (i == 0)or(i == np.pi):
                break

def checkRatio(theta = np.pi/2):
    #Check M_simu/M_real ratio for a constant theta
    phi_12 = np.linspace(0, 2*np.pi, num = 12, endpoint = False)
    counter = 0
    for i in phi_12:
        counter += 1
        print(str(counter) + ' of 12: theta = ' + str(theta) + ' , phi = ' + str(i))
        mysubhalo.drawPlot(theta, i, ratioCheck = True)
    mysubhalo.drawRatioPlot()


def checkRatio2():
    #Check M_simu/M_real for 26Pics
    counter = 0
    for i in thetaArr:
        for j in phiArr:
            counter += 1
            print(str(counter)+ ' of 26: theta = '+str("{:.2f}".format(i))+', phi = '+str("{:.2f}".format(j)) )
            mysubhalo.drawPlot(i,j,center = True,ratioCheck = True)
            #Time elapsed is ' + str(datetime.now()-runTime))
            runTime = datetime.now()
            if (i == 0)or(i == np.pi):
                break
    mysubhalo.drawRatioPlot()
    

if __name__ == "__main__":
    #get26Pics(True)
    get26Pics(False)
    #checkRatio()
    #checkRatio2()
    print('Job finished. Time elapsed: '+ str(datetime.now()-time0))



