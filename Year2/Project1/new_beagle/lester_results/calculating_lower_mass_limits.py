#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 16:39:28 2021

@author: lester
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/mass_completeness_limits_0p95.fits'
mass_completeness_limits_0p95 = fits.open(fileName)
#print(data_fits.info())
#print(data_fits[1].header)

fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/mass_completeness_limits_0p95_ap_0p8.fits'
mass_completeness_limits_0p95_ap_0p8 = fits.open(fileName)
#print(data_fits[1].header)

fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/mass_completeness_limits_0p95_new.fits'
mass_completeness_limits_0p95_new = fits.open(fileName)
#print(data_fits.info())
#print(data_fits[1].header)

fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/mass_completeness_limits_0p90_new.fits'
mass_completeness_limits_0p90_new = fits.open(fileName)
#print(data_fits.info())
#print(data_fits[1].header)

mcl_mass = mass_completeness_limits_0p95[1].data['mass']
mcl_z = mass_completeness_limits_0p95[1].data['redshift']

mcl_ap_0p8_mass = mass_completeness_limits_0p95_ap_0p8[1].data['mass']
mcl_ap_0p8_z = mass_completeness_limits_0p95_ap_0p8[1].data['redshift']

mcl_mass_new = mass_completeness_limits_0p95_new[1].data['mass']
mcl_z_new = mass_completeness_limits_0p95_new[1].data['redshift']

mcl_mass_new_90 = mass_completeness_limits_0p90_new[1].data['mass']
mcl_z_new_90 = mass_completeness_limits_0p90_new[1].data['redshift']

santini_mass = [8.3, 8.5, 8.8, 8.8, 8.8]
santini_z = [1.65, 2.5, 3.5, 4.5, 5.5]

fit = np.polyfit(mcl_z, mcl_mass, 1)
fit_ap_0p8 = np.polyfit(mcl_ap_0p8_z, mcl_ap_0p8_mass, 1)
fit_new = np.polyfit(mcl_z_new, mcl_mass_new, 1)
fit_new_90 = np.polyfit(mcl_z_new_90, mcl_mass_new_90, 1)


# =============================================================================
# determining the limits to use
# =============================================================================
'''
I think you could safely fit the first part with a straight line
 from z~2-4.5 or so, put a minimum at 8.1 and plateau at ~9.1 
 (wherever your fitted line crosses 9.1).
'''

fit_new_2_to_4p5 = np.polyfit(mcl_z_new[(mcl_z_new>2.0)&(mcl_z_new<4.5)], mcl_mass_new[(mcl_z_new>2.0)&(mcl_z_new<4.5)], 1)

mass_low = 8.1
mass_high = 9.1

z_low = (mass_low - fit_new_2_to_4p5[1]) / fit_new_2_to_4p5[0]
z_high = (mass_high - fit_new_2_to_4p5[1]) / fit_new_2_to_4p5[0]


fit_new_90_2_to_4p5 = np.polyfit(mcl_z_new_90[(mcl_z_new_90>2.0)&(mcl_z_new_90<4.5)], mcl_mass_new_90[(mcl_z_new_90>2.0)&(mcl_z_new_90<4.5)], 1)

mass_low_90 = 8.0
mass_high_90 = 9.0

z_low_90 = (mass_low_90 - fit_new_90_2_to_4p5[1]) / fit_new_90_2_to_4p5[0]
z_high_90 = (mass_high_90 - fit_new_90_2_to_4p5[1]) / fit_new_90_2_to_4p5[0]

# =============================================================================
# plot
# =============================================================================

x = np.linspace(0, 9, 2)
plt.xlim(x)
plt.scatter(mcl_z, mcl_mass, label='MCL')
plt.plot(x, fit[1]+x*fit[0])

plt.scatter(mcl_ap_0p8_z, mcl_ap_0p8_mass, label='MCL ap 0.8')
plt.plot(x, fit_ap_0p8[1]+x*fit_ap_0p8[0])

plt.scatter(santini_z, santini_mass, label='Santini', color='r')

x = np.linspace(z_low, z_high, 2)
plt.scatter(mcl_z_new, mcl_mass_new, label='MCL new', color='tab:green')
plt.plot(x, fit_new_2_to_4p5[1]+x*fit_new_2_to_4p5[0], color='tab:green')
plt.plot((0, z_low), (mass_low, mass_low), color='tab:green')
plt.plot((z_high, 9), (mass_high, mass_high), color='tab:green')

x = np.linspace(z_low_90, z_high_90, 2)
plt.scatter(mcl_z_new_90, mcl_mass_new_90, label='MCL new 90', color='k')
plt.plot(x, fit_new_90_2_to_4p5[1]+x*fit_new_90_2_to_4p5[0], color='k')
plt.plot((0, z_low_90), (mass_low_90, mass_low_90), color='k')
plt.plot((z_high_90, 9), (mass_high_90, mass_high_90), color='k')

x = np.linspace(z_low, z_high, 2)
plt.scatter(mcl_z_new, mcl_mass_new-0.1, label='MCL new - 0.1', color='tab:pink')
plt.plot(x, fit_new_2_to_4p5[1]+x*fit_new_2_to_4p5[0] - 0.1, color='tab:pink')
plt.plot((0, z_low), (mass_low - 0.1, mass_low - 0.1), color='tab:pink')
plt.plot((z_high, 9), (mass_high - 0.1, mass_high - 0.1), color='tab:pink')



plt.xlabel('Redshift')
plt.ylabel('log(Mass)')
plt.legend()
plt.show()

print(fit)
print(fit_ap_0p8)
print(fit_new_2_to_4p5)
print(z_low, z_high)
print(fit_new_90_2_to_4p5)
print(z_low_90, z_high_90)
















