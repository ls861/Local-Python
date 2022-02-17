#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 29 21:01:05 2022

@author: lester
"""




import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from astropy.io import fits
import random

# might want to use 33? or 34?
filename = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_34_clusters_z1p25-6p0.fits'



cat = fits.open(filename)

M = cat[1].data['mass_BEAGLE_stellar']
z = cat[1].data['redshift_BEAGLE']

# plt.hist(M)
# plt.show()

# plt.hist(z)
# plt.show()



M1 = M[(z>1.25)&(z<2.0)]
M2 = M[(z>2.0)&(z<3.0)]
M3 = M[(z>3.0)&(z<4.0)]
M4 = M[(z>4.0)&(z<5.0)]
M5 = M[(z>5.0)&(z<6.0)]

plt.hist(M1, histtype='step', density=False)
plt.hist(M2, histtype='step', density=False)
plt.hist(M3, histtype='step', density=False)
plt.hist(M4, histtype='step', density=False)
plt.hist(M5, histtype='step', density=False)
plt.show()


Ms = [M1, M2, M3, M4, M5]

for i in range(len(Ms)):
    print('')
    for j in range(len(Ms)):
        ks = ks_2samp(Ms[i], Ms[j])
        print('M{} M{} stat={:.6f} p={:.6f}'.format(i+1,j+1,ks[0],ks[1]))



#%%
m1_1 = random.choices(M1, k=1000)
m1_2 = random.choices(M1, k=1000000)
m1_3 = random.choices(M1, k=1000000)

k1 = ks_2samp(m1_2, m1_2)
k2 = ks_2samp(m1_2, m1_3)
k3 = ks_2samp(m1_1, m1_3)
k4 = ks_2samp(m1_1, m1_2)

print(k1)
print(k2)
print(k3)
print(k4)

plt.hist(m1_1, histtype='step', density=True, bins=50)
plt.hist(m1_2, histtype='step', density=True, bins=50)
plt.hist(m1_3, histtype='step', density=True, bins=50)
plt.show()


#%%
# WORKING WITH EMMA



M1sub = random.choices(M1, k=len(M5))
k = ks_2samp(M1, M1sub)

plt.hist(M1, density=True)
plt.hist(M1sub, density=True)
plt.hist(M5, density=True)
plt.show()

print(k)



print(M3)















