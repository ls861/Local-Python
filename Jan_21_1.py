#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 14:42:28 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Jan_20_1_012/mock_catalogue.fits'
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

print(tau) # 14 
print(msa) # 7

tau_size = 14
msa_size = 7


### ### constant t, vary tau ### ###
a = 0 
b = tau_size

for i in range(msa_size):
    

    #plt.figure(figsize=(10, 10))
    plt.title('t = %.1g' % (msa[a]))
    plt.xlim(0, 20000)
    plt.ylim(-24, -16)
    
    for j in np.arange(a, b):
        plt.plot(wl, np.log10(f[j]), label=r'$\tau$ = %.1g' % (tau[j]) )

    plt.legend()
    plt.show()

    a += tau_size
    b += tau_size



### ### constant tau, vary t ### ### SHAPE OF PLOT to 4x4
    
fig2, axs2 = plt.subplots(4, 4, figsize=(20,20))

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
    
    if d == 3:
        c += 1

    d = (d+1)%4
    
axs2[3,2].axis('off')
axs2[3,3].axis('off')
axs2[2,3].legend()    
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































