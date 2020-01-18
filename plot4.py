#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


plt.figure(figsize=(15,10))    
plt.title('TITLE')
ax=plt.gca()

ax.set_xlim(500,5000)
ax.set_ylim(-30,4)

plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
plt.ylabel('Energy Density log($\lambda F_\lambda$) ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)



plt.xticks(range(0,8000, 250))


age = []

step = 40   # how many ages to skip (20)
gap = 1    # how many wavelengths to skip (30)



### NO DUST ###
    
file_number = "31"
fileName = "/Users/lester/BEAGLE/BEAGLE-general/results/Lester_" + file_number + "/my_first_test.fits"
catalogue = fits.open(fileName)
N_Objects = len(catalogue['FULL SED'].data)

for obj in range(0, N_Objects, step):
    
    age_s = int(catalogue['GALAXY PROPERTIES'].data[obj][7])
    
    if age_s >= 1000000:
        wl = (catalogue['FULL SED WL'].data[0][0]) ### angstroms
        spec = catalogue['FULL SED'].data[obj,:]     
        plt.plot((wl[::gap]), np.log10((spec)[::gap]), color='b')
        age.append('$10^{%s}$' % (str(int(np.log10(age_s)))))
        
        
### Calzetti, tau = 1 ###

file_number = "2"
fileName = "/Users/lester/BEAGLE/BEAGLE-general/results/Lester2_" + file_number + "/my_first_test.fits"
catalogue = fits.open(fileName)
N_Objects = len(catalogue['FULL SED'].data)

for obj in range(0, N_Objects, step):
    
    age_s = int(catalogue['GALAXY PROPERTIES'].data[obj][7])
    
    if age_s >= 1000000: 
        wl = (catalogue['FULL SED WL'].data[0][0]) ### angstroms
        spec = catalogue['FULL SED'].data[obj,:]
        plt.plot((wl[::gap]), np.log10((spec)[::gap]), color='r')   
        age.append('$10^{%s}$' % (str(int(np.log10(age_s)))))
        
    
### Calzetti, tau = 10 ###

file_number = "3"
fileName = "/Users/lester/BEAGLE/BEAGLE-general/results/Lester2_" + file_number + "/my_first_test.fits"
catalogue = fits.open(fileName)
N_Objects = len(catalogue['FULL SED'].data)

for obj in range(0, N_Objects, step):
    
    age_s = int(catalogue['GALAXY PROPERTIES'].data[obj][7])
    
    if age_s >= 1000000: 
        wl = (catalogue['FULL SED WL'].data[0][0]) ### angstroms
        spec = catalogue['FULL SED'].data[obj,:]
        plt.plot((wl[::gap]), np.log10((spec)[::gap]), color='g')   
        age.append('$10^{%s}$' % (str(int(np.log10(age_s)))))


### Calzetti, tau = 20 ###

file_number = "7"
fileName = "/Users/lester/BEAGLE/BEAGLE-general/results/Lester2_" + file_number + "/my_first_test.fits"
catalogue = fits.open(fileName)
N_Objects = len(catalogue['FULL SED'].data)

for obj in range(0, N_Objects, step):
    
    age_s = int(catalogue['GALAXY PROPERTIES'].data[obj][7])
    
    if age_s >= 1000000: 
        wl = (catalogue['FULL SED WL'].data[0][0]) ### angstroms
        spec = catalogue['FULL SED'].data[obj,:]
        plt.plot((wl[::gap]), np.log10((spec)[::gap]), color='k')   
        age.append('$10^{%s}$' % (str(int(np.log10(age_s)))))
    

plt.legend(age, prop={'size': 14}) 
plt.show()
    










