#util.py
#Methods to generate HI profile line for a given subhalo from given direction
#Haiyin Song, 2021/11
import h5py
import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt
import sys
import os


#I don't think this is useful right now
'''
class coord:
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
    ############################################
    basePath = '/public/furendeng/TNG-100/output'
    ############################################
    outputName = str('snap'+str(snapNum)+'_subhalo'+str(subhalo_id).zfill(7))
    output = outputName + '.png'
    dir_name = 'Snap'+str(snapNum)+'/'
    #Create output path for this snapshot if not exist
    try:
        os.mkdir(dir_name)
    except:
        pass

    def __init__(self,snapNum,subhalo_id):
        self.snapNum = snapNum
        self.subhalo_id = subhalo_id
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
    
    def luminosityDistance(self):
        #TODO: define a luminosity distance by the size of this subhalo

    
    def obsPoint(self,theta,phi):
        #TODO: 



    def LOSVelocity(self):
        #TODO: I'm not sure if I should keep this function but I'll let it be here for a while


    def drawPlot(self):
        #TODO: should save figure in subdir
