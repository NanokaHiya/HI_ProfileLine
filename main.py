#Generate a simple HI line profile from TNG-100
#Haiyin Song, 2021/10
import h5py
import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt
import sys

#TODO: take redshift into account!
######################################################################################################
#Loading information
if len(sys.argv)!=3:
    raise ValueError('Expect argv = 3. Please provide snapNum and subhalo id')
try:
    snapNum = int(sys.argv[1])
    subhalo_id = int(sys.argv[2])
except:
    print('argv value error')
basePath = '/public/furendeng/TNG-100/output'
output = str('snap'+str(snapNum)+'_subhalo'+str(subhalo_id)+'.png')
######################################################################################################

#Load constants
fields = ['Velocities','CenterOfMass','NeutralHydrogenAbundance','GFM_Metals','Masses']
partType = 0 #Gas particle
delta_v = 1.4 #As compatible resolution of Arecibo, unit = km/s
a=1
def snapInfo():
    #Load snapshot Header from first chunk
    #Copied from il.snapshot.loadSubset()
    #Get cosmological constants
    with h5py.File(il.snapshot.snapPath(basePath, snapNum),'r') as f:
        header = dict(f['Header'].attrs.items())
    return header

snapInfo = snapInfo()
h = snapInfo['HubbleParam']
z = snapInfo['Redshift']
#print(sHeader)

def groupcatInfo():
    #Read subhalo position from groupcat
    #Return: 1. luminosity distance in Mpc and subhalo velocity in km/s
	#		 2. Peculiar velocity of group in km/s
    subhalo_g = il.groupcat.loadSingle(basePath, snapNum, subhaloID=subhalo_id)
    #subhalo center position [Mpc]
    subhaloCM=subhalo_g['SubhaloCM']*h*0.001
    return np.linalg.norm(subhaloCM,2),subhalo_g['SubhaloVel']
    
D, subhaloVel = groupcatInfo()
print('subhaloVel=',subhaloVel)

######################################################################################################

def LOSVelocity():
    #TODO: check why returning always positive velocity!
    #Takes velocities onto centerofmass direction to get LOS velocity
    #Return array shape (N,) with unit km/s
    CM_km = subhalo['CenterOfMass']*3.086*np.power(10,16)*h
    #CM_km = []
    #CM_km.append(
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
    return MassHI/(D*D*2.356*np.power(10,5)*delta_v)

######################################################################################################
subhalo = il.snapshot.loadSubhalo(basePath, snapNum, subhalo_id, partType, fields)
MassHI = MassHI()
LOSVelocity = LOSVelocity()
print('LOSVelocity=',LOSVelocity)
fluxDensity = fluxDensity(MassHI,D)

'''
Depreciated
rad_velocity = LOSVelocity(subhalo42['Velocities'],subhalo42['CenterOfMass'])
HIMASS = MassHI(subhalo42['NeutralHydrogenAbundance'],subhalo42['GFM_Metals'],subhalo42['Masses'])
fluxDen = fluxDensity(HIMASS,subhalo42['CenterOfMass'])
v_max = np.nanmax(rad_velocity)
v_min = np.nanmin(rad_velocity)
bins=int((v_max-v_min)/delta_v)
#print(bins)
'''

######################################################################################################
#Save histogram as PNG file
#TODO: tags
n, bins, patches = plt.hist(LOSVelocity, weights = fluxDensity, bins=300)
plt.savefig(output,format='png')
