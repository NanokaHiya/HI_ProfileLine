# HI_ProfileLine
Plot a HI profile line from TNG-100 simulation

# Logs
## 11/18/2021
One-key generate 32 + 1 plot for a single subhalo.

## 11/17/2021
Main implentation of class subhalo completed.

## 11/16/2021
Reconstructing.
Looks like currently there is no need for coordinate class.

Create class structure for subhalo reading/ploting.
Utility for coordinate system change
## 11/3/2021
Reason for last error: no gas particle exist in some subhalos. Add check if gas exists.
## 11/1/2021
Known error for snapShot 99, subhalos with higher ordinal returns KeyError: 'GFM_Metals'. Need to check the shape of GFM_Metals for 43
Unexpected shape for subhalo: 
1. Probabaly too large?
2. Symmetry problem
7. Seems like not much HI present here 
42. Symmetry problem
43. Seems like not much HI present here 
