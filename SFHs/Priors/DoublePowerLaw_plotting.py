#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 13:20:17 2020

@author: lester
"""




import numpy as np
import matplotlib.pyplot as plt
import cosmolopy.distance as cd
import cosmolopy.constants as cc
from scipy.special import erf
from scipy.integrate import quad
import matplotlib.colors as mcolors

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)

size        = np.int(1E1)
z           = np.random.uniform(low=4.5, high=5.5, size=size)                   # redshift

age_z_4p5   = cd.age(4.5, **cosmo)/cc.yr_s                                      # yr, (log, 9.107)
age_z_15    = cd.age(15, **cosmo)/cc.yr_s                                       # yr, age of universe at z=15
age_galaxy  = (10**np.random.uniform(low=6, high=10, size=size))                # yr, age of galaxy

t_arr       = cd.age(z, **cosmo)/cc.yr_s                                        # yr, age of Universe
t0_arr      = t_arr-age_galaxy                                                  # yr, start of star formation
tau_arr     = np.random.uniform(low=0.1, high=t_arr, size=size)                 # yr, width of function

m_arr       = np.random.uniform(low=5, high=12, size=size)                      # log, mass of galaxy
massArr     = np.arange(5, 13, 1)

alpha_arr   = (10**np.random.uniform(low=-1, high=3, size=size))
beta_arr    = (10**np.random.uniform(low=-1, high=3, size=size))



tau     = 8E9
t       = (1E9)*np.arange(0.1, 15, 0.1)
alpha   = 10
beta    = 10
m       = 10**10 

integrand = lambda T: 1 / (((T/tau)**alpha)+((T/tau)**-beta))
integral  = quad(integrand, 0, 2E9)
A = m / integral[0] 

sfr    = A / (((t/tau)**alpha)+((t/tau)**-beta))
plt.plot(t, sfr)










    
    
















