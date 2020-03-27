#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:17:13 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from astropy.io import fits

param1 = 'DPL'
revision1 = '010'
title1 = param1 + ' ' + revision1

param2 = 'DPL'
revision2 = '011'
title2 = param2 + ' ' + revision2

size = 13
fsize = 8

# =============================================================================
# INPUT - get "real" parameters
# =============================================================================

fileName = '/Users/lester/Documents/PhD/param_DPL/mock_013/mock_catalogue_DPL_001.fits'
data_fits = fits.open(fileName)

mtot_r = np.log10(data_fits['GALAXY PROPERTIES'].data['m_tot'])
sfr_r = data_fits['STAR FORMATION'].data['SFR']
tau_r = np.log10(data_fits['STAR FORMATION BINS'].data['bin_tau'])

data_fits.close()

# =============================================================================
# OUTPUT - get BEAGLE parameters 1
# =============================================================================

fileName = '/Users/lester/Documents/PhD/param_{}/fit_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(param1, revision1)
data_fits = fits.open(fileName)

id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1
mtot_b1 = data_fits['POSTERIOR PDF'].data['mass_mean']
mtot_68_b1 = data_fits['POSTERIOR PDF'].data['mass_68.00']
sfr_b1 = data_fits['STAR FORMATION'].data['SFR_mean']
sfr_68_b1 = data_fits['STAR FORMATION'].data['SFR_68.00']

data_fits.close()

# =============================================================================
# OUTPUT - get BEAGLE parameters 2
# =============================================================================

fileName = '/Users/lester/Documents/PhD/param_{}/fit_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(param2, revision2)
data_fits = fits.open(fileName)

id_b2 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1
mtot_b2 = data_fits['POSTERIOR PDF'].data['mass_mean']
mtot_68_b2 = data_fits['POSTERIOR PDF'].data['mass_68.00']
sfr_b2 = data_fits['STAR FORMATION'].data['SFR_mean']
sfr_68_b2 = data_fits['STAR FORMATION'].data['SFR_68.00']

data_fits.close()

# =============================================================================
# calculate input gradient for rising or falling
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_MS_parameters_004.fits'
data_fits = fits.open(fileName)
#print(data_fits[1].header)
alpha = data_fits[1].data['dpl_alpha']
beta = data_fits[1].data['dpl_beta']
tau = 10**(data_fits[1].data['tau'])
A = data_fits[1].data['A']
sfr = data_fits[1].data['sfr']
data_fits.close()

ageUniv2 = 3228839870.9122815
xlin = np.linspace(1, 1e10, 100000)
grad = np.empty(len(A))

# nice trick to find index in xlin which has value closest to ageUniv2
idx = (np.abs(xlin - ageUniv2)).argmin()

for i in range(len(A)):
    sfr_calc = A[i] / (((xlin/tau[i])**alpha[i])+((xlin/tau[i])**-beta[i]))
    grad[i] = np.gradient(sfr_calc, xlin)[idx]

idx_rising = grad >= 0
idx_falling = grad < 0

# =============================================================================
# PLOT - input mass vs output mass 1
# =============================================================================

print(len(mtot_r))
print(len(mtot_b1))

idx_r1 = idx_rising[id_b1]
idx_f1 = idx_falling[id_b1]

plt.figure(figsize=(fsize, fsize))
plt.title('Input Mass (DPL) vs Output Mass ({})'.format(title1), size=size)
plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)
plt.plot((7.5, 11), (7.5, 11))
plt.scatter(mtot_r[id_b1][idx_r1], mtot_b1[idx_r1], s=10, zorder=1, color='r', label='rising')
plt.scatter(mtot_r[id_b1][idx_f1], mtot_b1[idx_f1], s=10, zorder=1, color='g', label='falling')
plt.errorbar(mtot_r[id_b1], mtot_b1, yerr=[mtot_b1 - mtot_68_b1[:, 0], mtot_68_b1[:, 1] - mtot_b1], linestyle="None", elinewidth=0.5, color='k')

plt.xlim(7.5, 11)
plt.ylim(7.5, 11)
plt.legend()
plt.show()


# =============================================================================
# PLOT - input mass vs output mass 2
# =============================================================================

print(len(mtot_r))
print(len(mtot_b2))

idx_r2 = idx_rising[id_b2]
idx_f2 = idx_falling[id_b2]

plt.figure(figsize=(fsize, fsize))
plt.title('Input Mass (DPL) vs Output Mass ({})'.format(title2), size=size)
plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)
plt.plot((7.5, 11), (7.5, 11))
plt.scatter(mtot_r[id_b2][idx_r2], mtot_b2[idx_r2], s=10, zorder=1, color='r', label='rising')
plt.scatter(mtot_r[id_b2][idx_f2], mtot_b2[idx_f2], s=10, zorder=1, color='g', label='falling')
plt.errorbar(mtot_r[id_b2], mtot_b2, yerr=[mtot_b2 - mtot_68_b2[:, 0], mtot_68_b2[:, 1] - mtot_b2], linestyle="None", elinewidth=0.5, color='k')

plt.xlim(7.5, 11)
plt.ylim(7.5, 11)
plt.legend()
plt.show()

# =============================================================================
# PLOT - input sfr vs output sfr 1
# =============================================================================

plt.figure(figsize=(fsize, fsize))
plt.title('Input SFR (DE) vs Output SFR ({})'.format(title1), size=size)
plt.xlabel(r'$\text{Input - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.ylabel(r'$\text{Output - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.plot((-1, 3.5), (-1, 3.5))
plt.scatter(np.log10(sfr_r[id_b1])[idx_r1], np.log10(sfr_b1)[idx_r1], s=10, zorder=1, color='r', label='rising')
plt.scatter(np.log10(sfr_r[id_b1])[idx_f1], np.log10(sfr_b1)[idx_f1], s=10, zorder=1, color='g', label='falling')
plt.errorbar(np.log10(sfr_r[id_b1]), np.log10(sfr_b1), yerr=[np.log10(sfr_b1 / sfr_68_b1[:, 0]), np.log10(sfr_68_b1[:, 1] / sfr_b1)], linestyle="None", elinewidth=1, color='k', zorder=0)

plt.xlim(-1, 3.5)
plt.ylim(-1, 3.5)
plt.legend()
plt.show()


# =============================================================================
# PLOT - input sfr vs output sfr 2
# =============================================================================

plt.figure(figsize=(fsize, fsize))
plt.title('Input SFR (DE) vs Output SFR ({})'.format(title2), size=size)
plt.xlabel(r'$\text{Input - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.ylabel(r'$\text{Output - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.plot((-1, 3.5), (-1, 3.5))
plt.scatter(np.log10(sfr_r[id_b2])[idx_r2], np.log10(sfr_b2)[idx_r2], s=10, zorder=1, color='r', label='rising')
plt.scatter(np.log10(sfr_r[id_b2])[idx_f2], np.log10(sfr_b2)[idx_f2], s=10, zorder=1, color='g', label='falling')
plt.errorbar(np.log10(sfr_r[id_b2]), np.log10(sfr_b2), yerr=[np.log10(sfr_b2 / sfr_68_b2[:, 0]), np.log10(sfr_68_b2[:, 1] / sfr_b2)], linestyle="None", elinewidth=1, color='k', zorder=0)

plt.xlim(-1, 3.5)
plt.ylim(-1, 3.5)
plt.legend()
plt.show()

# =============================================================================
# matching indicies for 1 & 2
# =============================================================================

test = np.intersect1d(id_b1, id_b2)

b_boo1 = []
for i in range(len(id_b1)):
    if id_b1[i] in test:
        b_boo1.append(True)
    else:
        b_boo1.append(False)

b_boo2 = []
for i in range(len(id_b2)):
    if id_b2[i] in test:
        b_boo2.append(True)
    else:
        b_boo2.append(False)
        
        
print(len(b_boo1), len(b_boo2), len(test))

# =============================================================================
# PLOT - input mass vs output mass for 1 & 2
# =============================================================================

plt.figure(figsize=(1.6*fsize, fsize))
plt.title('Mass ({}) vs Mass ({})'.format(title1, title2), size=size)
plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)
plt.plot((7.5, 11), (7.5, 11))
#plt.scatter(mtot_r[id_b1][b_boo1], mtot_b1[b_boo1], s=10, c=tau_r[id_b1][b_boo1], zorder=10)
plt.scatter(mtot_r[id_b1][b_boo1], mtot_b1[b_boo1], s=10, zorder=10)
#plt.errorbar(mtot_r[id_b1], mtot_b1, yerr=[mtot_b1 - mtot_68_b1[:, 0], mtot_68_b1[:, 1] - mtot_b1], linestyle="None", elinewidth=0.5, color='k')

#plt.colorbar()

for i in range(len(test)):
    plt.plot((mtot_r[id_b1][b_boo1][i], mtot_r[id_b2][b_boo2][i]), (mtot_b1[b_boo1][i], mtot_b2[b_boo2][i]))

plt.xlim(7.5, 11)
plt.ylim(7.5, 11)
plt.show()

