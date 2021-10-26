#Generate a simple HI line profile from TNG-100
#Haiyin Song, 2021/10
import h5py
import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt


#Loading information
basePath = '/public/furendeng/TNG-100/output'
fields = ['Velocities','CenterOfMass','NeutralHydrogenAbundance','GFM_Metals','Masses']
subhalo_id = 42
snapNum = 99
z=0
partType = 0 #Gas particle

#Load constants
def snapHeader(basePath,snapNum):
    #Load Header from first chunk
    #Redundant to some extent. Copied from il.snapshot.loadSubset()
    #Target: get cosmological constants
    with h5py.File(il.snapshot.snapPath(basePath, snapNum),'r') as f:
        header = dict(f['Header'].attrs.items())
    return header

header42 = snapHeader(basePath,snapNum)
h = header42['HubbleParam'] #Hubble constant, value subject to change, unit = H0/(100km/s/Mpc)
delta_v = 1.4 #As compatible resolution of Arecibo, unit = km/s


def LOSVelocity(velocities,centerofmass):
    #TODO: check why returning always positive velocity!
    #Takes velocities onto centerofmass direction to get LOS velocity
    #Return array shape (N,) with unit km/s
    centerofmass_km = centerofmass*3.086*np.power(10,16)*h
    upper = (velocities*centerofmass_km).sum(1)
    lower = np.sqrt((centerofmass_km*centerofmass_km).sum(1))
    return np.transpose(upper/lower)

def MassHI(HiAbundance,GFM_Metals,Masses):
    #TODO: change number 19710 with length number of particles
    #Caculated by gas mass*HydrogenAbundance(from GFM_Metals)*NeutralHydrogenAbundance
    #Return array shape (N,) with unit M_dot
    Masses_Mdot = Masses*h*np.power(10,10)
    remove = np.arange(1,10)
    M_x = np.reshape(np.delete(GFM_Metals,remove,axis=1),(19710))
    return M_x*HiAbundance*Masses_Mdot

def fluxDensity(MassHI,centerofmass):
    #convert MassHI into flux density (Catinella et al. 2010)
    #Set N(v)=0
    #use distance to center of mass as luminosity distance (z=0 case)
    #should Return array shape(N,) with unit Jy
    centerofmass_Mpc = centerofmass*h*0.001
    D = np.sqrt((centerofmass_Mpc*centerofmass_Mpc).sum(1))
    return MassHI/(D*D*2.356*np.power(10,5)*delta_v)



subhalo42 = il.snapshot.loadSubhalo(basePath, snapNum, subhalo_id, partType, fields)
rad_velocity = LOSVelocity(subhalo42['Velocities'],subhalo42['CenterOfMass'])
HIMASS = MassHI(subhalo42['NeutralHydrogenAbundance'],subhalo42['GFM_Metals'],subhalo42['Masses'])
fluxDen = fluxDensity(HIMASS,subhalo42['CenterOfMass'])
v_max = np.nanmax(rad_velocity)
v_min = np.nanmin(rad_velocity)
bins=int((v_max-v_min)/delta_v)
#print(bins)


#Save histogram as PNG file
#Note: There is no meaning of the actual numbers here!
n, bins, patches = plt.hist(fluxDen, weights = HIMASS, bins=300)
plt.savefig('conversionTest.png',format='png')

