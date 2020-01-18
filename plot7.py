#!/usr/bin/env python
# -*- coding: utf-8 -*-

### DUST COMPARISON PRE 10^7 ###

import matplotlib.pyplot as plt
from astropy.io import fits


plt.figure(figsize=(13,13))    
plt.title('BEAGLE models')
ax=plt.gca()

ax.set_xlim(10**2,10**5)
ax.set_ylim(10**1,10**7)

plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)

plt.xscale('log')
plt.yscale('log')

#plt.xticks(range(0,6000, 250))

age = []

step = 20   # how many ages to skip (20)
gap = 1    # how many wavelengths to skip (30)


file_number     = ['2_9', '2_10', '2_11']
c               = ['b', 'r', 'g']
lab             = ['No Dust', 
                   r'Calzetti, $\hat{\tau}_V$=1',
                   r'Charlot&Fall, $\hat{\tau}_V$=1, $\mu$=0.3']


for i in range(len(file_number)):
    
    fileName = "/Users/lester/BEAGLE/BEAGLE-general/results/Lester" + file_number[i] + "/my_first_test.fits"
    catalogue = fits.open(fileName)
    N_Objects = len(catalogue['FULL SED'].data)
    
    for obj in range(0, N_Objects, step):
        
        age_s = catalogue['GALAXY PROPERTIES'].data[obj][7]
        
        if age_s >= 1000000 and age_s <= 100000000:
            
            age_s = '%.1e' % catalogue['GALAXY PROPERTIES'].data[obj][7]
            wl = (catalogue['FULL SED WL'].data[0][0]) ### angstroms
            spec = wl*(catalogue['FULL SED'].data[obj,:])  ### multiplied by wl
            plt.plot((wl[::gap]), spec[::gap], color=c[i], label=(age_s + ' years, ' + lab[i]))
            

plt.legend(prop={'size': 14})
plt.show()










