#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 14:42:28 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Jan_20_1_015/mock_catalogue.fits'
data_fits = fits.open(fileName)

#info = data_fits.info()
#header = data_fits['GALAXY PROPERTIES'].header

#wl_spec = data_fits[6].data[0][0]
#redshift = data_fits[1].data


z   = data_fits['GALAXY PROPERTIES'].data['redshift']
wl  = data_fits['FULL SED WL'].data['wl'][0]                                    # A
f   = data_fits['FULL SED'].data                                                # erg s-1 cm-2 A-1
tau = data_fits['STAR FORMATION BINS'].data['bin_tau']                          # yrs
msa = data_fits['STAR FORMATION'].data['max_stellar_age']                       # yrs

tau_size = 14
msa_size = 7

ind = (wl > 500) & (wl < 20000)
xlim = (0, 20000)
ylim = (-1, 1)

### ### constant t, vary tau ### ### 11

fig11, axs11 = plt.subplots(2, 4, figsize=(20,10))
fig12, axs12 = plt.subplots(2, 4, figsize=(20,10))

a = 0 
b = tau_size

c = 0
d = 0

for i in range(msa_size):
    
    r = np.arange(a, b)[7] # residual    
    
    for j in np.arange(a, b):
        if tau[j] <= msa[a]:
            axs11[c,d].plot(wl[ind], np.log10(f[j][ind]/f[r][ind]), label=r'$\tau$ = %.1g' % (tau[j]) )
        else:
            axs12[c,d].plot(wl[ind], np.log10(f[j][ind]/f[r][ind]), label=r'$\tau$ = %.1g' % (tau[j]) )

    axs11[c,d].set_xlim(xlim)
    axs11[c,d].set_ylim(-1, 0.1)
    axs11[c,d].set_title('msa = %.1g' % (msa[a]))
    axs11[c,d].legend()  
    
    axs12[c,d].set_xlim(xlim)
    axs12[c,d].set_ylim(-0.1, 0.15)
    axs12[c,d].set_title('msa = %.1g' % (msa[a]))
    axs12[c,d].legend()  

    a += tau_size
    b += tau_size

    if d == 3:
        c += 1

    d = (d+1)%4

axs11[1,3].axis('off') 
fig11.show()  

axs12[1,3].axis('off') 
fig12.show() 



### ### constant tau, vary t ### ### SHAPE OF PLOT to 4x4
    
fig21, axs21 = plt.subplots(4, 4, figsize=(20,20))
fig22, axs22 = plt.subplots(4, 4, figsize=(20,20))

a = 0
b = tau_size

c = 0
d = 0

for i in range(tau_size):
    
    r = np.arange(a, len(z), b)[3] # residual
    
    for j in np.arange(a, len(z), b):
        if msa[j] <= tau[i]:
            axs21[c,d].plot(wl[ind], np.log10(f[j][ind]/f[r][ind]), label=r'msa = %.1g' % (msa[j]) )
        else:
            axs22[c,d].plot(wl[ind], np.log10(f[j][ind]/f[r][ind]), label=r'msa = %.1g' % (msa[j]) )

    axs21[c,d].set_xlim(xlim)
    axs21[c,d].set_ylim(ylim)
    axs21[c,d].set_title(r'$\tau$ = %.1g' % (tau[i]))
    axs21[c,d].legend()
    
    axs22[c,d].set_xlim(xlim)
    axs22[c,d].set_ylim(ylim)
    axs22[c,d].set_title(r'$\tau$ = %.1g' % (tau[i]))
    axs22[c,d].legend() 
    
    a += 1
    
    if d == 3:
        c += 1

    d = (d+1)%4
    
axs21[3,2].axis('off')
axs21[3,3].axis('off') 
fig21.show()   

axs22[3,2].axis('off')
axs22[3,3].axis('off') 
fig22.show()   

### ### constant tau, vary t ### ### SHAPE OF PLOT to 2x7
'''
fig2, axs2 = plt.subplots(7, 2, figsize=(10,30))

a = 0
b = tau_size

c = 0
d = 0

for i in range(tau_size):
    
    for j in np.arange(a, len(z), b):
        
        axs2[c,d].set_xlim(0, 20000)
        axs2[c,d].set_ylim(-24, -16)
        axs2[c,d].plot(wl, np.log10(f[j]), label=r't = %.1g' % (msa[j]) )
        axs2[c,d].set_title(r'$\tau$ = %.1g' % (tau[i]))

    a += 1
    
    c += d
    d = (d+1)%2
   
axs2[6,1].legend()    
fig2.show()   
'''



data_fits.close()































