#Generate a simple HI line profile from TNG-100
#Haiyin Song, 2021/10
import h5py
import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt

#TODO: take redshift into account!
######################################################################################################
#Loading information
subhalo_id = 42
snapNum = 99
basePath = '/public/furendeng/TNG-100/output'
######################################################################################################

#Load constants
fields = ['Velocities','CenterOfMass','NeutralHydrogenAbundance','GFM_Metals','Masses']
partType = 0 #Gas particle
delta_v = 1.4 #As compatible resolution of Arecibo, unit = km/s

def snapHeader(basePath,snapNum):
    #Load Header from first chunk
    #Copied from il.snapshot.loadSubset()
    #Get cosmological constants
    with h5py.File(il.snapshot.snapPath(basePath, snapNum),'r') as f:
        header = dict(f['Header'].attrs.items())
    return header

sHeader = snapHeader(basePath,snapNum)
h = sHeader['HubbleParam']
#print(sHeader)
z = sHeader['Redshift']

def luminosityDis():
    #Read subhalo position from groupcat
    #Return luminosity distance in Mpc
    subhalo_g = il.groupcat.loadSingle(basePath, snapNum, subhaloID=subhalo_id)
    #subhalo center position [Mpc]
    subhaloCM=subhalo_g['SubhaloCM']*h*0.001
    return np.linalg.norm(subhaloCM,2)
    
D = luminosityDis()

######################################################################################################

def LOSVelocity():
    #TODO: check why returning always positive velocity!
    #Takes velocities onto centerofmass direction to get LOS velocity
    #Return array shape (N,) with unit km/s
    CM_km = subhalo['CenterOfMass']*3.086*np.power(10,16)*h
    upper = (subhalo['Velocities']*CM_km).sum(1)
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
plt.savefig('fluxDensity.png',format='png')
