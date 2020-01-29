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


### ### ### ### ### 28th Jan, shape using normalisation


'''
normalisation for one SED onto the other: (1/N) sum (Ai/Bi)
shape = (norm * Bi -Ai)^2
'''

shape4 = np.zeros(len(f))

for i in range(len(f)):
    
    if i == 0:
        k = +1
    else:
        k = -1

    norm = (1./len(f[0])) * sum(f[i+k]/f[i])
    d_tau = abs(tau[i]-tau[i+k])
    
    shape4[i] = np.mean((norm*f[i] - f[i+k]) ** 2)
    shape4[i] /= d_tau
    

plt.figure(figsize=(10, 5))

plt.plot(tau_arr, mass_arr / max(mass_arr), label='OFFSET')
#plt.plot(tau, shape)
plt.plot(tau[1:], shape4[1:] / 4E-48, label='SHAPE')
plt.plot((msa, msa), (0, 1.1), color='k', linestyle=':')   

plt.xlim(0-t0, 1E10-t0)
plt.ylim(0, 1.1)
plt.xlabel('TAU')
plt.ylabel('MASS')
plt.legend()
plt.show()
        



### ### ### ### ### 28th Jan, shape using normalisation against a single SED, and messing with indicies


'''
normalisation for one SED onto the other: (1/N) sum (Ai/Bi)
shape = (norm * Bi -Ai)^2
'''


fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Jan_20_1_018/mock_catalogue.fits'
data_fits = fits.open(fileName)

wl  = data_fits['FULL SED WL'].data['wl'][0]                                    # A
f   = data_fits['FULL SED'].data[:26]                                           # erg s-1 cm-2 A-1
tau = data_fits['STAR FORMATION BINS'].data['bin_tau'][:26]                     # yrs

data_fits.close()

c = int(1+len(f)/2) - 1 # index of roughly central SED ) 13 for N=26
ind = (wl > 500) & (wl < 20000)
wl = wl[ind]
f1 = np.zeros((len(f), len(wl)))
shape5 = np.zeros(len(f))
shape6 = np.zeros(len(f))
shape7 = np.zeros(len(f))
shape8 = np.zeros(len(f))
norm = np.zeros(len(f))

for i in range(len(f)):
    f1[i] = f[i][ind]

for i in range(len(norm)):
    norm[i] = np.median(f1[c] / f1[i])

for i in range(len(f1)):
    
    if i == 0:
        k = +1
    else:
        k = -1

    s1 = max( (norm[i]*f1[i] - norm[i+k]*f1[i+k]) ** 2)
    s2 = min( (norm[i]*f1[i] - norm[i+k]*f1[i+k]) ** 2)
    s3 = np.mean( (norm[i]*f1[i] - norm[i+k]*f1[i+k]) ** 2)
    s4 = np.median( (norm[i]*f1[i] - norm[i+k]*f1[i+k]) ** 2)
    
    d_tau = abs(tau[i]-tau[i+k])
    
    shape5[i] = s3 / d_tau
    shape6[i] = s4 / d_tau
    shape7[i] = s3
    shape8[i] = s4  



shape5 /= 10**-1
shape6 /= 10**-4
shape7 /= 10**-39.7
shape8 /= 10**-42.7
    
    
plt.figure(figsize=(10, 10))

plt.plot(tau_arr, mass_arr / max(mass_arr), label='OFFSET')
#plt.plot(tau, shape)
plt.plot(tau, shape5 / 4E-48, label='SHAPE')
plt.plot(tau, shape6 / 4E-48, label='SHAPE')
plt.plot(tau, shape7, label='SHAPE')
plt.plot(tau, shape8, label='SHAPE')
plt.plot((msa, msa), (0, 1.1), color='k', linestyle=':')   

plt.xlim(0-t0, 1E10-t0)
plt.ylim(0, 1.1)
plt.xlabel('TAU')
plt.ylabel('MASS')
plt.legend()
plt.show()
    




























        
        
        
        
        

