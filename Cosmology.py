#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 10:24:35 2019

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def Hubble_P(H0, z, O, w):
    '''
    Take H0 (km s-1 Mpc-1) and z to calculate the Hubble constant for a given redshift. 
    O is an array of density parameters such as: O = [Or, Ob, Odm, Ode]
    w is an array of indicies such as: w = [1/3, 0, 0, -1]
    '''
    
    H0 = H0 * (10**-19) / 3.086             # s-1
    
    E = 0
    
    for i in range(len(O)):
        E += (O[i] * ((1+z)**(3*(1+w[i])))) + ((1-sum(O)) * ((1+z)**2))

    return H0 * (E**0.5)

def D_p(H0, z, O, w):
    '''
    '''
    c = 299792458.
    n = 10000
    total = 0.
    dz = z/n
    
    for i in range(n):
        total += (c * dz) / Hubble_P(H0, (i+0.5)*dz, O, w)
        
    return total / ((10**6) * (3.086 * (10**16)))   # Mpc

Or = 0.
Ob = 0.
Odm = 0.3
Ode = 0.7

H0 = 70 # km s-1 Mpc-1
z = np.arange(0, 10, 0.03, dtype=float)
z = np.array([0, 2, 4, 6, 8, 10], dtype=float)
O = [Or, Ob, Odm, Ode]
w = [1./3., 0., 0., -1.]

H = np.zeros(len(z))
proper_dist = np.zeros(len(z))


start = datetime.now()

for i in range(len(z)):
    H[i] = Hubble_P(H0, z[i], O, w)
    proper_dist[i] = D_p(H0, z[i], O, w)
    
finish = datetime.now()

plt.plot(z, proper_dist)
plt.show()


plt.figure(figsize=(10, 10))

# plotting luminosity distance

plt.plot(z, proper_dist*(1+z) + 2000, label='manual calculation', linestyle=':', zorder=3)

# comparison with http://www.astro.ucla.edu/~wright/CosmoCalc.html


lum = np.array([0, 15537.2, 35842.1, 57707.7, 80454.2, 103791.7], dtype=float)
plt.plot(z,lum, label='online calculation', zorder=1)


# comparison with astropy   https://docs.astropy.org/en/stable/cosmology/

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


ast = np.array(cosmo.luminosity_distance(z))
plt.plot(z, ast - 2000, label='astropy calculation', linestyle='--', zorder=2)


plt.legend()
plt.show()





print(sum(O))


