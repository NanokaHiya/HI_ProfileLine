#subhalo.py
#Methods to generate HI profile line for a given subhalo from given direction
#Haiyin Song, 2021/11
import h5py
import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt
import os
import rotate_star
#from astropy import units

class subhalo:
    """Read given subhalo from a snapshot & corresponding group catalog.
       Generate HI profile line looking from (D,theta,phi) from the subhalo center.
    """
    #Change basePath for different simulation
    basePath = '/public/furendeng/TNG-100/output'
    delta_v = 1.4 #unit: km/s

    def __init__(self,snapNum,subhalo_id):
        self.snapNum = snapNum
        self.subhalo_id = subhalo_id
        self.dir_name = 'Snap'+str(self.snapNum)+'/' + 'Subhalo' + str(self.subhalo_id)
        self.gasData = il.snapshot.loadSubhalo(self.basePath, self.snapNum, self.subhalo_id, 0, ['Velocities','CenterOfMass','NeutralHydrogenAbundance','GFM_Metals','Masses'])
        #Check if self contain Neutrual Hydrogen
        #Create output path for this snapshot & subhalo if not exist
        if self.gasData ['count'] == 0:
            #If no gas particle exist, 
            try:
                os.makedirs(self.dir_name + '_NO_GAS')
            except:
                pass
            raise ValueError('Requested subhalo has no gas particle!')
        else:
            try:
                os.makedirs(self.dir_name)
            except FileExistsError:
                pass
        #Load information from first chunk of snapshot
        with h5py.File(il.snapshot.snapPath(self.basePath, self.snapNum),'r') as f:
            self.header = dict(f['Header'].attrs.items())
        #self.h = self.header['HubbleParam']*100
        self.h = self.header['HubbleParam']
        self.z = self.header['Redshift']
        self.a = self.header['Time']
        #Load information from group catagroy
        self.subhalo_g = il.groupcat.loadSingle(self.basePath, self.snapNum, subhaloID=self.subhalo_id)
        self.subhaloCM = self.subhalo_g['SubhaloCM']     #unit: ckpc/h
        self.subhaloVel = self.subhalo_g['SubhaloVel']   #unit: km/s
        #Observation distance, which is distance to the center of subhalo,  in kpc/h unit
        #This should NOT be used as Luminosity distance
        #self.kpcD = np.linalg.norm(self.subhaloCM,2)
        self.ratioDict = {}
        print('Subhalo CM: '+str(self.subhaloCM))

    def __luminosityDistance(self):
        #currently return simply distance from center of comoving coord
        #return luminosity distance in Mpc
        #TODO: define a luminosity distance by the size of this subhalo
        self.arrayD = np.linalg.norm(self.obsPoint - self.gasData['CenterOfMass'], axis = 1)*self.h/1000
        #print('Average luminosity distance ratio: {:.2f}'.format(np.average(self.arrayD)*1000/(self.h*self.kpcD)))

    def axis(self):
        #Use rotate_star to get subhalo shape information
        self.starData = il.snapshot.loadSubhalo(self.basePath, self.snapNum, self.subhalo_id, 4, ['Coordinates', 'Velocities','Masses'])
        self.xstar = self.starData['Coordinates']*self.h
        self.vstar = self.starData['Velocities']*np.sqrt(self.a)
        self.mstar = self.starData['Masses']*self.h
        self.ba, self. ca, self.angle, self.Tiv = rotate_star.get_shape(self.xstar, self.mstar)
        self.xstar_axis = np.dot(self.Tiv, self.xstar.T).T
        self.vstar_axis = np.dot(self.Tiv, self.vstar.T).T


    def __obsPoint(self, theta, phi, kpcD):
        #with given angle (and luminosity distance) calculate the observation point
        #spherical->recentered->comoving
        #unit is ckpc/h
        self.theta = theta
        self.phi = phi
        self.kpcD = kpcD
        self.MpcD = self.kpcD*self.h/1000
        self.obsPoint = np.add(self.subhaloCM, [self.kpcD* np.cos(self.phi)* np.sin(self.theta), self.kpcD* np.sin(self.phi)* np.sin(self.theta), self.kpcD* np.cos(self.theta)])
        print('obsPoint: ' + str(self.obsPoint))

    def __LOSVelocity(self, center = True):
        #self.LOSVector is the new vector pointing from obsPoint to particles
        #When center == True, total subhalo velocity is subtracted from particle velocity
        #All units here km/s
        self.LOSVector = (self.gasData['CenterOfMass']-self.obsPoint)*3.086*np.power(10,16)*self.h*1000
        self.absLOSVector = np.linalg.norm(self.LOSVector, axis = 1)
        if center == True:
            self.LOSVelocity = (((np.sqrt(self.a)*self.gasData['Velocities']-self.subhaloVel)*self.LOSVector).sum(1))/self.absLOSVector
        else:    
            self.LOSVelocity = (((np.sqrt(self.a)*self.gasData['Velocities'])*self.LOSVector).sum(1))/self.absLOSVector

    def __MassHI(self):
        #unit: M_sun
        self.MassHI = np.reshape(np.delete(self.gasData['GFM_Metals'],np.arange(1,10),axis=1),(len(list(self.gasData['Masses']))))* self.gasData['NeutralHydrogenAbundance']*self.h*np.power(10,10)*self.gasData['Masses']
        self.actualHIMass = np.sum(self.MassHI)

    def __fluxDensity(self):
        #unit: Jy
        self.__MassHI()
        self.__luminosityDistance()
        self.fluxDen = self.MassHI*(1+self.z)/(self.arrayD*self.arrayD*2.356*np.power(10,5)*self.delta_v)


    def __bins(self):
        #calculate number of bins for histogram
        self.nBins = int((np.nanmax(self.LOSVelocity) - np.nanmin(self.LOSVelocity))/self.delta_v)


    def drawPlot(self, theta = 0, phi = 0, kpcD = 1000000, clean = True, center = True, ratioCheck = False):
        self.__obsPoint(theta,phi,kpcD)
        self.__LOSVelocity(center)
        self.__fluxDensity()
        self.__bins()

        print('Doing theta = {:2f}, phi = {:2f}'.format(self.theta,self.phi))
        self.n, self.bins, self.patches = plt.hist(self.LOSVelocity, weights = 1000*self.fluxDen, bins=self.nBins, histtype='step')
        #Title and file name
        if clean == False:
            plt.title( 'Snapshot '+str(self.snapNum)+' Subhalo '+str(self.subhalo_id)+' Symmetry Check')
            self.outputName = 'symmetry_check'
        else:
            plt.title('Snapshot '+str(self.snapNum)+' Subhalo '+str(self.subhalo_id)+r', Mock HI Profile line, $ \theta = ' + str("{:.2f}".format(self.theta))+', \phi = '+str("{:.2f}".format(self.phi))+'$')
            self.outputName = 'theta_'+str("{:.2f}".format(self.theta))+'_phi_'+str("{:.2f}".format(self.phi))
        if center == True:
            self.outputName += '_center.png'
        else:
            self.outputName += '.png'

        #Plotting details
        plt.xlabel('Velocity / $km*s^{-1}$')
        plt.ylabel('FluxDen / $mJy$')
        
        #Calculate simulated HI mass from plt
        if clean == True:
            self.actual_delta_v = self.bins[1]-self.bins[0]
            self.simu_real_mass_ratio = np.sum(self.n/1000)*self.actual_delta_v*2.356*np.power(10,5)*(self.MpcD**2)/self.actualHIMass
            print('M_simu/M_real = {:5f}'.format(self.simu_real_mass_ratio))
            plt.annotate(r'$\frac{M_{simu}}{M_{real}} = $'+str(self.simu_real_mass_ratio), xy=(0.05, 0.85), xycoords='axes fraction')
            if ratioCheck == True:
                #Warning: ratio check should only be used with constant theta!
                self.ratioDict[str(self.theta)+', '+str(self.phi)] = self.simu_real_mass_ratio-1
        
        #fig, axs = plt.subplots(2, 1)
        #self.n, self.bins, self.patches = axs[0].hist(self.LOSVelocity, weights = self.fluxDen, bins=self.bins, histtype='step')
        #fig.suptitle('Snapshot '+str(self.snapNum)+' Subhalo '+str(self.subhalo_id))
        #axs[0].set_title(r'Mock HI Profile line, $ \theta = ' + str("{:.2f}".format(self.theta))+', \phi = '+str("{:.2f}".format(self.phi))+'$')
        #axs[0].set_xlabel('Velocity / $km*s^{-1}$')
        #axs[0].set_ylabel('FluxDen / mJy')
        #axs[1].barh(0.2, self.actualHIMass)
        #axs[1].get_yaxis().set_visible(False)
        #plt.subplots_adjust(hspace=0.2)
        #plt.text(10, 10, '$D = '+str(self.arrayD)+'\ntheta = '+str(self.theta)+'\nphi = '+str(self.phi) + '\nbins = ' + str(self.bins) + '$', fontsize=12)
        
        plt.savefig(self.dir_name +'/'+ self.outputName, format='png')
        if clean == True:
            plt.clf()

    def drawRatioPlot(self):
        #Warning: This function is not well implented for performance and should be redone
        try:
            plt.clf()
        except:
            pass
        #self.ratioNames = ['0',r'$\frac{\pi}{6}$',r'$\frac{\pi}{3}$',r'$\frac{\pi}{2}$',r'$\frac{2\pi}{3}$',r'$\frac{5\pi}{6}$',r'$\pi$', r'$\frac{7\pi}{6}$',r'$\frac{4\pi}{3}$',r'$\frac{3\pi}{2}$',r'$\frac{5\pi}{3}$',r'$\frac{11\pi}{6}$']
        self.ratioNames = list(self.ratioDict.keys())
        self.ratioData = list(self.ratioDict.values())
        fig, ax = plt.subplots()
        ax.bar(self.ratioNames,self.ratioData)
        ax.set_title('Snapshot '+str(self.snapNum)+' Subhalo '+str(self.subhalo_id)+' Mass Ratio Check')
        #ax.set_title(r'Percentage difference between $M_{simu}$ and $M_{real}$ at $\theta = $'+str("{:.2f}".format(self.ratioCheckTheta)))
        plt.savefig(self.dir_name + '/ratio2.png', format = 'png')
        plt.clf()



