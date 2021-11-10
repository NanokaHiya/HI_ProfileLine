#Generate a simple HI line profile from TNG-100
#Haiyin Song, 2021/10
import h5py
import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt
import sys
import os
import astropy
from astropy.cosmology import LambdaCDM

######################################################################################################
#Loading information & Checking
if len(sys.argv)!=3:
    raise ValueError("Expect argv = 3. Please provide snapNum and subhalo id")
try:
    snapNum = int(sys.argv[1])
    subhalo_id = int(sys.argv[2])
except:
    print('argv value error')
basePath = '/public/furendeng/TNG-100/output'
outputName = str('snap'+str(snapNum)+'_subhalo'+str(subhalo_id).zfill(7))
output = outputName + '.png'
dir_name = 'Snap'+str(snapNum)+'/'
os.mkdir(dir_name)
######################################################################################################

#Load constants
fields = ['Velocities','CenterOfMass','NeutralHydrogenAbundance','GFM_Metals','Masses']
partType = 0 #Gas particle
delta_v = 1.4 #As compatible resolution of Arecibo, unit = km/s

#Load subhalo info from snapshot
subhalo = il.snapshot.loadSubhalo(basePath, snapNum, subhalo_id, partType, fields)

#Check if this subhalo has no gas
if subhalo['count']==0:
    try:
        f = open(dir_name+str(outputName + '_NO_GAS'),'x')
        f.close()
    except FileExistsError:
        pass
    raise ValueError('Requested subhalo has no gas particle!')

def snapInfo():
    #Load snapshot Header from first chunk
    #Copied from il.snapshot.loadSubset()
    #Get cosmological constants
    with h5py.File(il.snapshot.snapPath(basePath, snapNum),'r') as f:
        header = dict(f['Header'].attrs.items())
    return header

#Set up cosmology
snapInfo = snapInfo()
#print(snapInfo)
h = snapInfo['HubbleParam']
H0 = 100*h
#print(h0)
z = snapInfo['Redshift']
a = snapInfo['Time']
#Previously I use AstroPy to calculate luminosity distance, however I believe redshift here
#is a time indicator instead of the distance of light traveled. Gas particles are always
#at their "current" state.
"""
Omega0 = snapInfo['Omega0']
OmegaLambda = snapInfo['OmegaLambda']
OmegaBaryon = snapInfo['OmegaBaryon']
#cosmo = astropy.cosmology.FLRW(H0 = H0,Om0 = Omega0,Ode0 = OmegaLambda,Ob0 = OmegaBaryon)
#Get luminosity distance
cosmo = LambdaCDM(H0 = H0,Om0 = Omega0,Ode0 = OmegaLambda)
D = cosmo.luminosity_distance(z)
"""

def groupcatInfo():
    #Read subhalo position from groupcat
    #Return: 1. luminosity distance in Mpc and subhalo velocity in km/s
    #        2. Peculiar velocity of group in km/s
    subhalo_g = il.groupcat.loadSingle(basePath, snapNum, subhaloID=subhalo_id)
    #print(subhalo_g)
    #subhalo center position [Mpc]
    subhaloCM=subhalo_g['SubhaloCM']*h*0.001
    return np.linalg.norm(subhaloCM,2),subhalo_g['SubhaloVel']
    
D,subhaloVel = groupcatInfo()
#print('subhaloVel=',subhaloVel)

######################################################################################################
#Calculation utilities
def LOSVelocity():
    #Takes velocities onto centerofmass direction to get LOS velocity
    #Return array shape (N,) with unit km/s
    CM_km = subhalo['CenterOfMass']*3.086*np.power(10,16)*h
    upper = ((np.sqrt(a)*subhalo['Velocities']-subhaloVel)*CM_km).sum(1)
    lower = np.sqrt((CM_km*CM_km).sum(1))
    return np.transpose(upper/lower)

def MassHI():
    #Caculated by gas mass*HydrogenAbundance(from GFM_Metals)*NeutralHydrogenAbundance
    #Return array shape (N,) with unit M_dot
    M_x = np.reshape(np.delete(subhalo['GFM_Metals'],np.arange(1,10),axis=1),(len(list(subhalo['Masses']))))
    return M_x*subhalo['NeutralHydrogenAbundance']*h*np.power(10,10)*subhalo['Masses']

def fluxDensity(MassHI,D):
    #convert MassHI into flux density (Catinella et al. 2010)
    #Set N(v)=0
    #Unit is mJy
    return MassHI*1000*(1+z)/(D*D*2.356*np.power(10,5)*delta_v)

######################################################################################################

MassHI = MassHI()
LOSVelocity = LOSVelocity()
print('LOSVelocity=',LOSVelocity)
bins = int((np.nanmax(LOSVelocity) - np.nanmin(LOSVelocity))/delta_v)
print('bins=',bins)
fluxDensity = fluxDensity(MassHI,D)

######################################################################################################
#Save histogram as PNG file
#TODO: tags
n, bins, patches = plt.hist(LOSVelocity, weights = fluxDensity, bins=bins, histtype='step')
plt.title('Snapshot '+str(snapNum)+' Subhalo '+str(subhalo_id))
plt.xlabel('Velocity')
plt.ylabel('FluxDen')
plt.savefig(output,format='png')
