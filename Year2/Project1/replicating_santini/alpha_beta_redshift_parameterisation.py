#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 17:07:08 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from astropy.cosmology import FlatLambdaCDM

norm = 9.7
#norm = 0

# =============================================================================
# Santini+17 'True' values - delayed SFH, SFR from UV slope
# =============================================================================
# logSFR = alpha log(M / M_9.7) + beta
z_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))

# converting normalisation
alpha_san = B_san - 9.7*A_san
alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
beta_san = A_san
beta_err_san = A_err_san

alpha_san_n = alpha_san + (norm*beta_san) # santini normalised
alpha_err_san_n = (alpha_err_san**2 - (norm*beta_err_san)**2) ** 0.5

# =============================================================================
# Santini+17 Original values - obtained by eye - delayed SFH, SFR from UV slope
# =============================================================================
# logSFR = alpha log(M / M_9.7) + beta
z_san0 = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san0 = np.array((1.05, 1.1, 0.9, 0.75, 0.55))
A_err_san0 = np.array((0.03, 0.03, 0.04, 0.05, 0.18))
B_san0 = np.array((1.0, 1.15, 1.25, 1.2, 1.5))
B_err_san0 = np.array((0.05, 0.03, 0.03, 0.06, 0.12))

# converting normalisation
alpha_san0 = B_san0 - 9.7*A_san0
alpha_err_san0 = (B_err_san0**2 + (9.7*A_err_san0)**2) ** 0.5
beta_san0 = A_san0
beta_err_san0 = A_err_san0

alpha_san0_n = alpha_san0 + (norm*beta_san0) # santini normalised
alpha_err_san0_n = (alpha_err_san0**2 - (norm*beta_err_san0)**2) ** 0.5

# =============================================================================
# scenario 7_7, NO BOOTSTRAPPING, just with and without clipping
# =============================================================================

alpha_scenario_7_7_wo = np.array([-6.38,-5.37,-3.77,-5.67,-7.01])
beta_scenario_7_7_wo = np.array([0.744,0.638,0.482,0.762,0.942])
alpha_scenario_7_7_wo_n = alpha_scenario_7_7_wo + (norm*beta_scenario_7_7_wo) # santini normalised

alpha_scenario_7_7 = np.array([-5.14,-5.66,-5.09,-6.32,-6.76])
beta_scenario_7_7 = np.array([0.585,0.667,0.64,0.819,0.881])
alpha_scenario_7_7_n = alpha_scenario_7_7 + (norm*beta_scenario_7_7) # santini normalised

#alpha_err_scenario_7_7 = np.array((alpha_scenario_7_7-np.array([,,,,]),np.array([,,,,])-alpha_scenario_7_7))
#beta_err_scenario_7_7 = np.array((beta_scenario_7_7-np.array([,,,,]),np.array([,,,,])-beta_scenario_7_7))
#alpha_err_scenario_7_7_n = abs(((((alpha_err_scenario_7_7[1]+alpha_err_scenario_7_7[0])/2))**2 - (norm*(((beta_err_scenario_7_7[1]+beta_err_scenario_7_7[0])/2)))**2)) ** 0.5

# =============================================================================
# STUFF NEEDED
# =============================================================================
z_lester = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
x = np.linspace(1.3, 6.0, 1000)
inc = 0.03 # increment in offset
color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
         
var_s = [beta_san, beta_san0, beta_scenario_7_7_wo, beta_scenario_7_7, alpha_san, alpha_san0, alpha_scenario_7_7_wo, alpha_scenario_7_7, alpha_san_n, alpha_san0_n, alpha_scenario_7_7_wo_n, alpha_scenario_7_7_n]

# =============================================================================
# functional form fitting
# =============================================================================

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html

def func1(x, a, b):
    return a + b*x

def func2(x, a, b, c):
    return a + b*x + c*x**2

def func3(x, a, b, c, d):
    return a + b*x + c*x**2 + d*x**3   

def func4(x, a, b, c, d, e):
    return a + b*x + c*x**2 + d*x**3 + e*x**4

def func5(x, a, b, c):
    return a + b*(1/x) + c*(1/x)**2

def func6(x, a, b, c, d):
    return a + b*(1/x)

#def func7(x, a, b, c):
#    return a + b*(1/x) + c*(1/x)**2

#def func8(x, a, b, c):
#    return a * np.exp(-b * x) + c

#functions = [func1, func2, func3, func4, func5, func7]
functions = [func1, func2]
functions = [func3, func4]
functions = [func5]
functions = [func6]

# =============================================================================
# PLOT slope
# =============================================================================
plt.figure(figsize=(15, 6))
plt.title('slope')

plt.scatter(z_san, beta_san, label='Santini+17')
plt.errorbar(z_san, beta_san, yerr=beta_err_san, ls='none')

plt.scatter(z_san0, beta_san0, label='Santini+17 Raw')
plt.errorbar(z_san0, beta_san0, yerr=beta_err_san0, ls='none')

plt.scatter(z_lester+ 0*inc, beta_scenario_7_7_wo, label='Scenario 7 7 Without Clipping')
plt.scatter(z_lester+ 0*inc, beta_scenario_7_7, label='Scenario 7 7 With Clipping')

for f in functions:
    for i in range(4):
        popt, pcov = curve_fit(f, z_lester, var_s[i])
        plt.plot(x, f(x, *popt), color=color[i])
    
plt.ylim(0.0, 1.5)
plt.legend()
plt.show()

# =============================================================================
# PLOT intercept
# =============================================================================
plt.figure(figsize=(15, 6))
plt.title('intercept')

plt.scatter(z_san, alpha_san, label='Santini+17')
plt.errorbar(z_san, alpha_san, yerr=alpha_err_san, ls='none')

plt.scatter(z_san0, alpha_san0, label='Santini+17 Raw')
plt.errorbar(z_san0, alpha_san0, yerr=alpha_err_san0, ls='none')

plt.scatter(z_lester+ 0*inc, alpha_scenario_7_7_wo, label='Scenario 7 7 Without Clipping')
plt.scatter(z_lester+ 0*inc, alpha_scenario_7_7, label='Scenario 7 7 With Clipping')

for f in functions:
    for i in range(4,8):
        popt, pcov = curve_fit(f, z_lester, var_s[i])
        plt.plot(x, f(x, *popt), color=color[i-4])
        
plt.ylim(-11, -2)
plt.legend()
plt.show()

# =============================================================================
# PLOT normalised as in Santini - slope unaffected
# =============================================================================
plt.figure(figsize=(15, 6))
plt.title('santini normalised intercept')

plt.scatter(z_san, alpha_san_n, label='Santini+17')
plt.errorbar(z_san, alpha_san_n, yerr=alpha_err_san_n, ls='none')

plt.scatter(z_san0, alpha_san0_n, label='Santini+17 Raw')
plt.errorbar(z_san0, alpha_san0_n, yerr=alpha_err_san0_n, ls='none')

plt.scatter(z_lester+ 0*inc, alpha_scenario_7_7_wo_n, label='Scenario 7 7 Without Clipping')
plt.scatter(z_lester+ 0*inc, alpha_scenario_7_7_n, label='Scenario 7 7 With Clipping')

for f in functions:
    for i in range(8,12):
        popt, pcov = curve_fit(f, z_lester, var_s[i])
        plt.plot(x, f(x, *popt), color=color[i-8])
    
plt.ylim(0, 3)
plt.legend()
plt.show()

# =============================================================================
# Also confirming that 1/z is similar to time:
# =============================================================================
#%%
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

t = cosmo.age(x).value

plt.plot(x, 1/x, label='1 / z')
plt.plot(x, t/6, label='t / 6Gyr')
plt.xlabel('z')
plt.legend()
plt.show()





























