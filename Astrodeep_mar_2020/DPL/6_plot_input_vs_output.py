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
revision1 = '002'
title1 = param1 + ' ' + revision1

param2 = '006'
revision2 = '006'
title2 = param2 + ' ' + revision2

size = 15
fsize = 10

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
# OUTPUT1 - get BEAGLE parameters
# =============================================================================

fileName = '/Users/lester/Documents/PhD/param_{}/fit_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(param1, revision1)
data_fits = fits.open(fileName)

id_b = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1
mtot_b = data_fits['POSTERIOR PDF'].data['mass_mean']
mtot_68_b = data_fits['POSTERIOR PDF'].data['mass_68.00']
sfr_b = data_fits['STAR FORMATION'].data['SFR_mean']
sfr_68_b = data_fits['STAR FORMATION'].data['SFR_68.00']

data_fits.close()

print(len(mtot_r))
print(len(mtot_b))

# =============================================================================
# PLOT - input mass vs output mass
# =============================================================================

plt.figure(figsize=(fsize, fsize))
plt.title('Input Mass (DE) vs Output Mass ({})'.format(title1), size=size)
plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)
plt.plot((7.5, 11), (7.5, 11))
plt.scatter(mtot_r[id_b], mtot_b, s=10, c=tau_r[id_b], zorder=10)
plt.errorbar(mtot_r[id_b], mtot_b, yerr=[mtot_b - mtot_68_b[:, 0], mtot_68_b[:, 1] - mtot_b], linestyle="None", elinewidth=0.5, color='k')

plt.colorbar()

plt.xlim(7.5, 11)
plt.ylim(7.5, 11)
plt.show()





'''
# =============================================================================
# OUTPUT - get BEAGLE parameters FOR NO REL ERROR
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/param_{}/astrodeep_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(param2, revision2)
data_fits = fits.open(fileName)

id_bn = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1
mtot_bn = data_fits['POSTERIOR PDF'].data['mass_mean']
mtot_68_bn = data_fits['POSTERIOR PDF'].data['mass_68.00']
sfr_bn = data_fits['STAR FORMATION'].data['SFR_mean']
sfr_68_bn = data_fits['STAR FORMATION'].data['SFR_68.00']

data_fits.close()

# =============================================================================
# PLOT - input mass vs output mass FOR NO REL ERROR
# =============================================================================

plt.figure(figsize=(fsize, fsize))
plt.title('Input Mass (DE) vs Output Mass ({})'.format(title2), size=size)
plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)
plt.plot((7.5, 11), (7.5, 11))
plt.scatter(mtot_r[id_bn], mtot_bn, s=10, c=tau_r[id_bn], zorder=10)
plt.errorbar(mtot_r[id_bn], mtot_bn, yerr=[mtot_bn - mtot_68_bn[:, 0], mtot_68_bn[:, 1] - mtot_bn], linestyle="None", elinewidth=0.5, color='k')

plt.colorbar()

plt.xlim(7.5, 11)
plt.ylim(7.5, 11)
plt.show()

# =============================================================================
# PLOT - input mass vs output for with and without REL ERROR
# =============================================================================

test = np.intersect1d(id_b, id_bn)

b_boo = []
for i in range(len(id_b)):
    if id_b[i] in test:
        b_boo.append(True)
    else:
        b_boo.append(False)

bn_boo = []
for i in range(len(id_bn)):
    if id_bn[i] in test:
        bn_boo.append(True)
    else:
        bn_boo.append(False)
        
        
print(len(b_boo), len(bn_boo), len(test))
    
# =============================================================================
# PLOT - input mass vs output mass
# =============================================================================

plt.figure(figsize=(1.6*fsize, fsize))
plt.title('Mass ({}) vs Mass ({})'.format(title1, title2), size=size)
plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)
plt.plot((7.5, 11), (7.5, 11))
#plt.scatter(mtot_r[id_b][b_boo], mtot_b[b_boo], s=10, c=tau_r[id_b][b_boo], zorder=10)
plt.scatter(mtot_r[id_b][b_boo], mtot_b[b_boo], s=10, zorder=10)
#plt.errorbar(mtot_r[id_b], mtot_b, yerr=[mtot_b - mtot_68_b[:, 0], mtot_68_b[:, 1] - mtot_b], linestyle="None", elinewidth=0.5, color='k')

#plt.colorbar()

for i in range(len(test)):
    plt.plot((mtot_r[id_b][b_boo][i], mtot_r[id_bn][bn_boo][i]), (mtot_b[b_boo][i], mtot_bn[bn_boo][i]))

plt.xlim(7.5, 11)
plt.ylim(7.5, 11)
plt.show()







# =============================================================================
# PLOT - input sfr vs output sfr
# =============================================================================

plt.figure(figsize=(fsize, fsize))
plt.title('Input SFR (DE) vs Output SFR ({})'.format(title1), size=size)
plt.xlabel(r'$\text{Input - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.ylabel(r'$\text{Output - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.plot((-1, 3.5), (-1, 3.5))
plt.scatter(np.log10(sfr_r[id_b]), np.log10(sfr_b), s=10, zorder=1)
plt.errorbar(np.log10(sfr_r[id_b]), np.log10(sfr_b), yerr=[np.log10(sfr_b / sfr_68_b[:, 0]), np.log10(sfr_68_b[:, 1] / sfr_b)], linestyle="None", elinewidth=0.5, color='k', zorder=0)

plt.xlim(-1, 3.5)
plt.ylim(-1, 3.5)
plt.show()

# =============================================================================
# PLOT - input sfr vs output sfr FOR NO REL ERROR
# =============================================================================

plt.figure(figsize=(fsize, fsize))
plt.title('Input SFR (DE) vs Output SFR ({})'.format(title2), size=size)
plt.xlabel(r'$\text{Input - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.ylabel(r'$\text{Output - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.plot((-1, 3.5), (-1, 3.5))
plt.scatter(np.log10(sfr_r[id_bn]), np.log10(sfr_bn), s=10, zorder=1)
plt.errorbar(np.log10(sfr_r[id_bn]), np.log10(sfr_bn), yerr=[np.log10(sfr_bn / sfr_68_bn[:, 0]), np.log10(sfr_68_bn[:, 1] / sfr_bn)], linestyle="None", elinewidth=0.5, color='k', zorder=10)

plt.xlim(-1, 3.5)
plt.ylim(-1, 3.5)
plt.show()





'''





