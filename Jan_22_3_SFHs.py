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

wl_spec = data_fits[6].data[0][0]
redshift = data_fits[1].data



z   = data_fits['GALAXY PROPERTIES'].data['redshift']
wl  = data_fits['FULL SED WL'].data['wl'][0]                                    # A
f   = data_fits['FULL SED'].data                                                # erg s-1 cm-2 A-1
tau = data_fits['STAR FORMATION BINS'].data['bin_tau']                          # yrs
msa = data_fits['STAR FORMATION'].data['max_stellar_age']                       # yrs

tau_size = 14
msa_size = 7

ind = (wl > 500) & (wl < 20000)


time = 1E9*np.linspace(0, 3, len(msa))
time_now = 1.12 * 1E9



### ### constant t, vary tau, plotting SFH ### ###

fig1, axs1 = plt.subplots(2, 4, figsize=(20,10))

a = 0 
b = tau_size

c = 0
d = 0

max_y = 0

for i in range(msa_size):
    
    for j in np.arange(a, b):
        
        t0 = time_now-msa[a]
        tau0 = tau[j]
        
        norm = max((time-t0) * np.exp(-(time-t0)/tau0))
        
        axs1[c,d].plot(time, ( (time-t0) * np.exp(-(time-t0)/tau0) )/ norm, label=r'$\tau$ = %.1g' % (tau[j]))
        
        axs1[c,d].set_title('t = %.1g' % (msa[a]))
        
        if max((time-t0) * np.exp(-(time-t0)/tau0)) > max_y:
            max_y = max((time-t0) * np.exp(-(time-t0)/tau0))

        axs1[c,d].set_xlim(0, 3E9)        
        axs1[c,d].set_ylim(0, 1.1)

    axs1[c,d].plot((time_now, time_now), (0, 1), color='k', linestyle=':')            
        
    a += tau_size
    b += tau_size

    if d == 3:
        c += 1

    d = (d+1)%4

axs1[1,3].axis('off')
axs1[0,3].legend()    
fig1.show()  





### ### constant tau, vary t, plotting SFH ### ### SHAPE OF PLOT to 4x4
    
fig2, axs2 = plt.subplots(4, 4, figsize=(20,20))

a = 0
b = tau_size

c = 0
d = 0

max_y = 0

for i in range(tau_size):
    
    for j in np.arange(a, len(z), b):
        
        t0 = time_now-msa[j]
        tau0 = tau[i]
        axs2[c,d].plot(time, (time-t0) * np.exp(-(time-t0)/tau0), label=r'msa = %.1g' % (msa[j]))

        axs2[c,d].set_title(r'$\tau$ = %.1g' % (tau[i]))
        
        if max((time-t0) * np.exp(-(time-t0)/tau0)) > max_y:
            max_y = max((time-t0) * np.exp(-(time-t0)/tau0))

        axs2[c,d].set_xlim(0, 3E9)        
        axs2[c,d].set_ylim(0, 1.1*max_y)

    axs2[c,d].plot((time_now, time_now), (0, 2*max_y), color='k', linestyle=':')

    a += 1
    
    if d == 3:
        c += 1

    d = (d+1)%4

    
axs2[3,2].axis('off')
axs2[3,3].axis('off')
axs2[0,3].legend()    
fig2.show()   





### ### constant tau, vary t ### ###
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































