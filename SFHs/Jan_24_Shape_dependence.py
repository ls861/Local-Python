#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 00:10:37 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Jan_20_1_018/mock_catalogue.fits'
data_fits = fits.open(fileName)


wl  = data_fits['FULL SED WL'].data['wl'][0]                                    # A
f   = data_fits['FULL SED'].data[:26]                                           # erg s-1 cm-2 A-1
tau = data_fits['STAR FORMATION BINS'].data['bin_tau'][:26]                     # yrs

data_fits.close()

### ### ### ### ###

A = 1
t = 1.12 * 1E9 # time of observation / age of universe at z=5
msa = 5E8  # 10**8.70
t0 = t - msa

### ### ### ### ###

ind = (wl > 500) & (wl < 20000)

shape = np.zeros(len(f))

for i in range(len(f)):
    
    i0 = f[i][ind]/np.mean(f[i][ind])
    
    if i == 0:
        ip = f[i][ind]/np.mean(f[i][ind])
        shape[i] = sum(ip/i0)
        
    else:
        im = f[i-1][ind]/np.mean(f[i-1][ind])
        shape[i] = sum(i0/im)
        

            
shape = shape / max(shape)

### ### ### ### ###



shape2 = np.zeros((len(f), len(f[0])))
shape3 = np.zeros(len(f))

for i in range(len(f)):
    if i == 0:
        k = +1
    else:
        k = -1
        
    for j in range(len(f[0])):
        shape2[i][j] = max(f[i+k][j], f[i][j]) / min(f[i+k][j], f[i][j])
    
    shape3[i] = 5000 * sum(shape2[i]) / abs(tau[i]-tau[i+k])
    
### ### ### ### ###
    
    
    
tau_arr = 10**np.arange(7.5, 10, 0.1)
mass_arr = A*(-tau_arr * ( ((np.exp(-(t-t0)/tau_arr))*(-t0+tau_arr+t)) - tau_arr) )


plt.figure(figsize=(10, 5))


plt.plot(tau_arr, mass_arr / max(mass_arr), label='OFFSET')
#plt.plot(tau, shape)
plt.plot(tau, shape3, label='SHAPE')
plt.plot((msa, msa), (0, 1.1), color='k', linestyle=':')   

plt.xlim(0-t0, 1E10-t0)
plt.ylim(0, 1.1)
plt.xlabel('TAU')
plt.ylabel('MASS')
plt.legend()
plt.show()


### ### ### ### ###


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

