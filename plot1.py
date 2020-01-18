#!/usr/bin/env python
# -*- coding: utf-8 -*-



import matplotlib.pyplot as plt

from astropy.io import fits

file_number = ["40", "42", "44", "45", "46", "48", "47"]

### 40  sfh_type = constant
### 42  sfh_type = constant
### 44  sfh_type = exponential tau = 7
### 45  sfh_type = exponential tau = 8
### 46  sfh_type = exponential tau = 9
### 48  sfh_type = exponential tau = 9.5
### 47  sfh_type = exponential tau = 10


for i in file_number:
    

    fileName = "/Users/lester/BEAGLE/BEAGLE-general/results/Old/BEAGLE_Docker_lester_" + i + "/my_first_test.fits"
    print(fileName)
    catalogue = fits.open(fileName)
    ### print catalogue.info()
    
    
    wl = (catalogue['FULL SED WL'].data[0][0]) / 10000 ### microns
    obj = 0
    spec = catalogue['FULL SED'].data[obj,:]
    
    
    
    plt.plot(wl, spec)
    
    
    ax=plt.gca()
    ax.set_xlim(0,1)
    #ax.set_xlim(0.01,0.5)
    #ax.set_ylim(0,1.5)
    
    
    plt.xlabel('Microns', fontsize=14)
    plt.ylabel('Energy Density', fontsize=14)
    
    
    plt.show()
