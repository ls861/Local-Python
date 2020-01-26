#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import time
import numpy as np
from astropy.io import fits

file_number = "31"
fileName = "/Users/lester/BEAGLE/BEAGLE-general/results/Lester_" + file_number + "/my_first_test.fits"
catalogue = fits.open(fileName)

print catalogue.info()
N_Objects = len(catalogue['FULL SED'].data)
print(catalogue['FULL SED'].data[9][:])

print('test')

#print(catalogue['GALAXY PROPERTIES'].data[7][obj])

#    age = int(catalogue['GALAXY PROPERTIES'].data[obj][7])
#    plt.title(age)

plt.figure(figsize=(15,10))    
plt.title('test')
ax=plt.gca()



### pre log10
ax.set_xlim(500,5000)
ax.set_ylim(-4,4)

#ax.set_xlim(0,10000)
#ax.set_ylim(0,7)

plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
#plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)

#plt.xlabel('Wavelength log($\lambda$) ($\AA$)', fontsize=14)
plt.ylabel('Energy Density log($\lambda F_\lambda$) ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)

age = []

step = 20   # how many ages to skip (20)
gap = 1    # how many wavelengths to skip (30)

for obj in range(0, N_Objects, step):
    
    age_s = int(catalogue['GALAXY PROPERTIES'].data[obj][7])
    
    if age_s >= 1000000:
    
    ### FULL SED ###    
        wl = (catalogue['FULL SED WL'].data[0][0]) ### angstroms
        spec = catalogue['FULL SED'].data[obj,:]
        
        plt.plot((wl[::gap]), np.log10((spec)[::gap]))
        
        age.append('10^'+str(int(np.log10(age_s))))
        
    #    plt.plot(wl, np.log10(spec))
    #    plt.plot(wl, spec_cont)
    
plt.xticks(range(0,5000, 250))
plt.legend(age) 
plt.show()
    
print(N_Objects)
print(age)
