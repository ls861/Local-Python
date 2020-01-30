#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 16:06:01 2020

@author: lester
"""

### ### DELAYED SFH OVERVIEW ### ###


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


'''
fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Jan_20_1_015/mock_catalogue.fits'
data_fits = fits.open(fileName)

#info = data_fits.info()
#header = data_fits['GALAXY PROPERTIES'].header

z   = data_fits['GALAXY PROPERTIES'].data['redshift']
wl  = data_fits['FULL SED WL'].data['wl'][0]                                    # A
f   = data_fits['FULL SED'].data                                                # erg s-1 cm-2 A-1
tau = data_fits['STAR FORMATION BINS'].data['bin_tau']                          # yrs
msa = data_fits['STAR FORMATION'].data['max_stellar_age']                       # yrs
data_fits.close()

tau_size = 14
msa_size = 7

ind = (wl > 500) & (wl < 20000)
xlim = (0, 20000)
ylim = (-21, -17)

time = 1E9*np.linspace(0, 3, len(msa))
time_now = 1.12 * 1E9
'''

### ### SFHs ### ###
### ### constant t, vary tau, plotting SFH ### ###

fig1, axs1 = plt.subplots(2, 4, figsize=(10,5))

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
        if tau[j] > msa[a]:
            axs1[c,d].plot(time, ( (time-t0) * np.exp(-(time-t0)/tau0) )/ norm, label=r'$\tau$ = %.1g' % (tau[j]))
        
        axs1[c,d].set_title('msa = %.1g' % (msa[a]))
        
        if max((time-t0) * np.exp(-(time-t0)/tau0)) > max_y:
            max_y = max((time-t0) * np.exp(-(time-t0)/tau0))

        axs1[c,d].set_xlim(0, 3E9)        
        axs1[c,d].set_ylim(0, 1.1)
    axs1[c,d].legend() 
    axs1[c,d].plot((time_now, time_now), (0, 1), color='k', linestyle=':')            
        
    a += tau_size
    b += tau_size

    if d == 3:
        c += 1

    d = (d+1)%4

axs1[1,3].axis('off')
   
fig1.show()  


### ### constant tau, vary t, plotting SFH ### ### SHAPE OF PLOT to 4x4
''' 
fig2, axs2 = plt.subplots(4, 4, figsize=(10,10))

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

'''
### ### SEDs ### ###
### ### constant t, vary tau ### ###
'''
fig1, axs1 = plt.subplots(2, 4, figsize=(10,5))

a = 0 
b = tau_size

c = 0
d = 0

for i in range(msa_size):
    
    for j in np.arange(a, b):
        
        r = np.arange(a, b)[7] # residual
        
        axs1[c,d].set_xlim(xlim)
        axs1[c,d].set_ylim(ylim)
        axs1[c,d].plot(wl[ind], np.log10(f[j][ind]), label=r'$\tau$ = %.1g' % (tau[j]) )
        axs1[c,d].set_title('msa = %.1g' % (msa[a]))
    
    a += tau_size
    b += tau_size

    if d == 3:
        c += 1

    d = (d+1)%4

axs1[1,3].axis('off')
axs1[0,3].legend()    
fig1.show()  

### ### constant tau, vary t ### ### SHAPE OF PLOT to 4x4
    
fig2, axs2 = plt.subplots(4, 4, figsize=(10,10))

a = 0
b = tau_size

c = 0
d = 0

for i in range(tau_size):
    
    for j in np.arange(a, len(z), b):
        r = np.arange(a, len(z), b)[3] # residual
        
        axs2[c,d].set_xlim(xlim)
        axs2[c,d].set_ylim(ylim)
        axs2[c,d].plot(wl[ind], np.log10(f[j][ind]), label=r'msa = %.1g' % (msa[j]) )
        axs2[c,d].set_title(r'$\tau$ = %.1g' % (tau[i]))

    a += 1
    
    if d == 3:
        c += 1

    d = (d+1)%4
    
axs2[3,2].axis('off')
axs2[3,3].axis('off')
axs2[0,3].legend()    
fig2.show()   

'''
### ### RESIDUALS ### ###
'''
ylim = (-1, 1)
'''
### ### constant t, vary tau ### ### 11
'''
fig11, axs11 = plt.subplots(2, 4, figsize=(10,5))
fig12, axs12 = plt.subplots(2, 4, figsize=(10,5))

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
    
fig21, axs21 = plt.subplots(4, 4, figsize=(10,10))
fig22, axs22 = plt.subplots(4, 4, figsize=(10,10))

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
'''

### ### Shape Dependence ### ###
'''
fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Jan_20_1_018/mock_catalogue.fits'
data_fits = fits.open(fileName)

wl  = data_fits['FULL SED WL'].data['wl'][0]                                    # A
f   = data_fits['FULL SED'].data[:26]                                           # erg s-1 cm-2 A-1
tau = data_fits['STAR FORMATION BINS'].data['bin_tau'][:26]                     # yrs

data_fits.close()

A = 1
t = 1.12 * 1E9 # time of observation / age of universe at z=5
msa = 5E8  # 10**8.70
t0 = t - msa

tau_arr = 10**np.arange(7.5, 10, 0.1)
mass_arr = A*(-tau_arr * ( ((np.exp(-(t-t0)/tau_arr))*(-t0+tau_arr+t)) - tau_arr) )

c = int(1+len(f)/2) - 1 # index of roughly central SED ) 13 for N=26
ind = (wl > 500) & (wl < 20000)
wl1 = wl[ind]
f1 = np.zeros((len(f), len(wl1)))
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
    
    
plt.figure(figsize=(5, 5))

plt.plot(tau_arr, mass_arr / max(mass_arr), label='OFFSET')
plt.plot(tau, shape5 / 4E-48, label='mean')
plt.plot(tau, shape6 / 4E-48, label='median')
plt.plot(tau, shape7, label='mean')
plt.plot(tau, shape8, label='median')
plt.plot((msa, msa), (0, 1.1), color='k', linestyle=':')   

plt.xlim(0-t0, 1E10-t0)
plt.ylim(0, 1.1)
plt.xlabel('TAU')
plt.ylabel('MASS')
plt.legend()
plt.show()
'''



























































