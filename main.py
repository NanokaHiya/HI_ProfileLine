#Generate a simple HI line profile from TNG-100
#Haiyin Song, 2021/10
import h5py
import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt

#Constants
delta_v = 1.4 #As compatible resolution of Arecibo, unit=km/s

def LOSVelocity(a,b):
    #TODO: check why returning always positive velocity!
    #TODO: the units of this velocity is INCORRECT
    #(velocities, centerofmass)
    #Takes velocities onto centerofmass direction to get LOS velocity
    #Should return array shape (N,)
    upper = (a*b).sum(1)
    lower = np.sqrt((b*b).sum(1))
    return np.transpose(upper/lower)

def MassHI(HiAbundance,GFM_Metals,Masses):
    #TODO: change number 19710 with length number of particles
    #return the mass of HI in unit of 10^10M_dot/h
    #Caculated by gas mass*HydrogenAbundance(from GFM_Metals)*NeutralHydrogenAbundance
    #Return array shape (N,)
    remove = np.arange(1,10)
    M_x = np.reshape(np.delete(GFM_Metals,remove,axis=1),(19710))
    return M_x*HiAbundance*Masses

#Load gas particles of single subhalo from TNG-100 snapshot 99 with redshift 0
basePath = '/public/furendeng/TNG-100/output'
fields = ['Velocities','CenterOfMass','NeutralHydrogenAbundance','GFM_Metals','Masses']
subhalo_id = 42
snapNum = 99
z=0
partType = 0

subhalo42 = il.snapshot.loadSubhalo(basePath, snapNum, subhalo_id, partType, fields)
rad_velocity = LOSVelocity(subhalo42['Velocities'],subhalo42['CenterOfMass'])
HIMASS = MassHI(subhalo42['NeutralHydrogenAbundance'],subhalo42['GFM_Metals'],subhalo42['Masses'])
v_max = np.nanmax(rad_velocity)
v_min = np.nanmin(rad_velocity)


#First attemption of histogram
#Note: There is no meaning of the actual numbers here!
n, bins, patches = plt.hist(rad_velocity, weights = HIMASS, bins=1000)
plt.savefig('test1000.png',format='png')

