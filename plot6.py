#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


plt.figure(figsize=(13,5))    
plt.title('BEAGLE models')
ax=plt.gca()

ax.set_xlim(10**2,10**5)
ax.set_ylim(10**-5,10**9)

plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)

plt.xscale('log')
plt.yscale('log')

#plt.xticks(range(0,6000, 250))


age = []

step = 40   # how many ages to skip (20)
gap = 1    # how many wavelengths to skip (30)


file_number     = ['_31', '2_2', '2_3', '2_neb_3', '2_8']
c               = ['b', 'k', 'g', 'r', 'c']
lab             = ['No Dust', 
                   'Calzetti, tau =1', 
                   'Calzetti, tau =10', 
                   'Calzetti, tau =1, neb=0.3', 
                   'No Dust, Z/Zs = 3']



for i in range(len(file_number)):
    
    fileName = "/Users/lester/BEAGLE/BEAGLE-general Oct/results/Lester" + file_number[i] + "/my_first_test.fits"
    catalogue = fits.open(fileName)
    N_Objects = len(catalogue['FULL SED'].data)
    
    for obj in range(0, N_Objects, step):
        
        age_s = int(catalogue['GALAXY PROPERTIES'].data[obj][7])
        
        if age_s >= 1000000:
            wl = (catalogue['FULL SED WL'].data[0][0]) ### angstroms
            spec = wl*(catalogue['FULL SED'].data[obj,:])  ### multiplied by wl
            plt.plot((wl[::gap]), spec[::gap], color=c[i], label=('$10^{%s} yrs, $' % (str(int(np.log10(age_s))))) + lab[i])



plt.legend(prop={'size': 4})
plt.show()










