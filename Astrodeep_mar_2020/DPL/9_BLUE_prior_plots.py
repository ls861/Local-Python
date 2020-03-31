#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 13:46:49 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import truncnorm
import cosmolopy.distance as cd
import cosmolopy.constants as cc
import math


np.seterr(invalid='warn', divide='warn', over='warn', under='warn')

L=1.0
Harr=1.0

#Larr = np.linspace(0.2, 1.0, num=17)
#Harr = np.linspace(1.0, 1.8, num=17)

#for H in Harr:
nObj = 10000

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
cosmo = cd.set_omega_k_0(cosmo)
zObs=2.0
ageUniv = cd.age(zObs, **cosmo)/cc.yr_s

mass = 10**np.random.uniform(low=7, high=12, size=nObj)                         # solar masses
alpha = 10**np.random.uniform(low=-1, high=3, size=nObj)
beta = 10**np.random.uniform(low=-1, high=3, size=nObj)


#    tau = np.random.uniform(low=0.1, high=ageUniv, size=nObj)                       # yrs
tau = 10**np.random.uniform(low=7, high=10, size=nObj)                          # yrs
#    tau = 10**np.random.uniform(low=9, high=9.7, size=nObj)                         # yrs
#    tau = np.random.uniform(low=L*ageUniv, high=H*ageUniv, size=nObj)               # yrs

#new beta and tau - weighted distributions
#beta = 10**(truncnorm.rvs(-2, 0, scale=2, loc=3, size=nObj))
#tau = truncnorm.rvs(-2, 0, scale=ageUniv/2, loc=ageUniv, size=nObj)
    
A = np.zeros(nObj)
sfr = np.zeros(nObj)
integral = np.zeros(nObj)

for i in range(nObj):
    
#        er=0
    integrand = lambda T: 1 / (((T/tau[i])**alpha[i])+((T/tau[i])**-beta[i]))
    
    integral[i]  = quad(integrand, 0, ageUniv)[0]

    

    A[i] = mass[i] / integral[i]



    
    sfr[i] = A[i] / (((ageUniv/tau[i])**alpha[i])+((ageUniv/tau[i])**-beta[i]))
#    print(integral[i], A[i], sfr[i])
    
    

#        print(er)
    if math.isnan(sfr[i]):
        print(alpha[i], beta[i], tau[i], mass[i], A[i], integral[i], sfr[i])

#    idx = sfr != 0
idx = sfr > 1e-4

# =============================================================================
# make PLOT sfr vs mass
# =============================================================================

#    plt.figure(figsize=(7, 5))
##    plt.title('tau from {} ageUniv to {} ageUniv'.format(L, H))
#    plt.hist2d(np.log10(mass[idx]), np.log10(sfr[idx]), bins=[50,2500], cmap='Blues', normed=True)
#    c=plt.colorbar()
#    c.set_label("weighting in prior")
#    plt.xlabel("log mass")
#    plt.ylabel("log sfr")
#    
#    plt.xlim(7,11)
#    plt.ylim(-4,3)
#    
#    ageUniv_2p0 = cd.age(2.0, **cosmo)/cc.yr_s                                      # yr
#    ageUniv_5p5 = cd.age(5.5, **cosmo)/cc.yr_s                                      # yr
#    ageUniv_4p5 = cd.age(4.5, **cosmo)/cc.yr_s                                      # yr
#    
#    plt.plot(np.log10(mass), np.log10(mass/(1E6)), label='t=1E6')
#    plt.plot(np.log10(mass), np.log10(mass/ageUniv_2p0), label='z=2.0')
#    #plt.plot(np.log10(mass), np.log10(mass/ageUniv_4p5), label='z=4.5')
#    #plt.plot(np.log10(mass), np.log10(mass/ageUniv_5p5), label='z=5.5')
#    
#    slope2 = 0.84 - 0.026*ageUniv*1e-9 # 0.7560501633562806
#    intercept2 = -(6.51 - 0.11*ageUniv*1e-9) # -6.1548276141996485
#    xlin = np.array([7, 12])
#    ylin2 = slope2*xlin + intercept2
#    plt.plot(xlin, ylin2, color='g', label='Speagle, z=2') # plot speagle MS @ z=2
#    
#    plt.legend()
#    plt.show()






















































