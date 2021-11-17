#util.py
#Methods to generate HI profile line for a given subhalo from given direction
#Haiyin Song, 2021/11
import h5py
import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt
import sys
import os

#No use yet
#class coord:
'''
    def __init__(self, a, b, c, flag='comoving'):
        #Init a point in a given coordinate system
        #These are used to check avaliablity.
        self.comoving = False
        self.recenter = False
        self.spherical = False
        if flag == 'comoving':
        #The original coordinates for particles in TNG
        #Default choice
            self.x = a
            self.y = b
            self.z = c
            self.comoving = True
            break
        else if flag == 'recenter':
        #Still cartisian; recentered at center of subhalo
            self.X = a
            self.Y = b
            self.Z = c
            self.recenter = True
            break
        else if flag == 'spherical':
        #Spherical coordinate centered at subhalo center
        #Pay attention to input of this init
            self.r = a
            self.theta = b
            self.phi = c
            self.spherical = True
            break
        else:
            raise ValueError('Problem assigining coordinate system')

    def recenter(self, center):
        #Perform a simple cartisian recenter action
        #comoving->recenter. Input center should be an array
        self.X = self.x - center[0]
        self.Y = self.y - center[1]
        self.Z = self.z - center[2]
        self.recenter = True

    def reToSpherical(self):
        #recenter->spherical
        #TODO
    '''


class subhalo:
    """Read given subhalo from a snapshot & corresponding group catalog.
       Generate HI profile line looking from (D,theta,phi) from the subhalo center.
    """
    basePath = '/public/furendeng/TNG-100/output'

    def __init__(self,snapNum,subhalo_id):
        self.snapNum = snapNum
        self.subhalo_id = subhalo_id
        self.outputName = str('snap'+str(self.snapNum)+'_subhalo'+str(self.subhalo_id).zfill(7))
        self.dir_name = 'Snap'+str(self.snapNum)+'/'
        #Create output path for this snapshot if not exist
        try:
            os.mkdir(self.dir_name)
        except:
            pass
        self.snapData = il.snapshot.loadSubhalo(basePath, self.snapNum, self.subhalo_id, 0, ['Velocities','CenterOfMass','NeutralHydrogenAbundance','GFM_Metals','Masses'])
        #Check if self contain Neutrual Hydrogen
        if self.snapData ['count'] == 0:
            #If no gas particle exist, 
            try:
                f = open(dir_name+str(outputName + '_NO_GAS'),'x')
                f.close()
            except FileExistsError:
                pass
            raise ValueError('Requested subhalo has no gas particle!')
        #Load information from first chunk of snapshot
        with h5py.File(il.snapshot.snapPath(self.basePath, self.snapNum),'r') as f:
            self.header = dict(f['Header'].attrs.items())
        self.H0 = self.header['HubbleParam']*100
        self.z = self.header['Redshift']
        self.a = self.header['Time']
        #Load information from group catagroy
        self.subhalo_g = il.groupcat.loadSingle(basePath, self.snapNum, subhaloID=self.subhalo_id)
        self.subhaloCM = self.subhalo_g['SubhaloCM']


    def luminosityDistance(self):
        #currently return simply distance from center of comoving coord
        #return luminosity distance in Mpc
        #TODO: define a luminosity distance by the size of this subhalo
        self.D = np.linalg.norm(self.subhaloCM*self.h*0.001,2) 
    
    def obsPoint(self, theta = 0, phi = 0):
        #with given angle (and luminosity distance) calculate the observation point
        #spherical->recentered->comoving
        self.theta = theta
        self.phi = phi
        self.obsPoint = np.add(self.subhaloCM,[self.D * np.cos(phi) * np.sin(theta),self.D * np.sin(phi) * np.sin(theta),self.D * np.cos(theta)])

    def LOSVelocity(self):
        #self.LOSVector is the new vector pointing from obsPoint to the subhaloCM (approximate)
        self.LOSVector = (self.subhaloCM-self.obsPoint)*3.086*np.power(10,16)*self.h
        self.LOSVelocity = np.transpose((((np.sqrt(self.a)*self.subhalo_g['Velocities']-self.subhaloVel)*self.LOSVector).sum(1))/(np.sqrt((self.LOSVector*self.LOSVector).sum(1))))

    def MassHI(self):
        self.MassHI = np.reshape(np.delete(self.snapData['GFM_Metals'],np.arange(1,10),axis=1),(len(list(self.snapData['Masses']))))* self.snapData['NeutralHydrogenAbundance']*h*np.power(10,10)*self.snapData['Masses']

    def fluxDen(self):
        self.fluxDen = MassHI*1000*(1+self.z)/(self.D*self.D*2.356*np.power(10,5)*delta_v)


    def bins(self):
        #calculate number of bins for histogram
        self.bins = int((np.nanmax(self.LOSVelocity) - np.nanmin(self.LOSVelocity))/delta_v)

    def drawPlot(self):
        #TODO: should save figure in subdir
        self.n, self.bins, self.patches = plt.hist(self.LOSVelocity, weights = self.fluxDen, bins=self.bins, histtype='step')
        plt.title('Snapshot '+str(self.snapNum)+' Subhalo '+str(self.subhalo_id))
        plt.xlabel('Velocity')
        plt.ylabel('FluxDen')
        plt.text(10, 10, '$D = '+str(self.D)+'\ntheta = '+str(self.theta)+'\nphi = '+str(self.phi) + '\nbins = ' + str(self.bins) + '$', fontsize=12)
        plt.savefig(self.dir_name + self.outputName + '_r'+str(self.D)+'_phi'+str(self.phi)+'_theta'+str(self.theta)+'.png' , format='png')



