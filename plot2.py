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



for obj in range(N_Objects):

    
    
### FULL SED ###    
    wl = (catalogue['FULL SED WL'].data[0][0]) / 10000 ### microns
    spec = catalogue['FULL SED'].data[obj,:]
    spec_cont = catalogue['CONTINUUM SED'].data[obj,:]
    
    
    plt.figure(figsize=(15,5))    
#    plt.plot(wl, np.log10(spec))
    plt.plot(wl, spec)
#    plt.plot(wl, spec_cont)
    
    
    ax=plt.gca()
    ax.set_xlim(0,3)
#    ax.set_ylim(-5,5)
    
    plt.xlabel('Microns', fontsize=14)
    plt.ylabel('Energy Density', fontsize=14)

    age = int(catalogue['GALAXY PROPERTIES'].data[obj][7])
    plt.title(age)
        
    plt.show()
    

    
    
    
    time.sleep(0.01)
    
print(N_Objects)

