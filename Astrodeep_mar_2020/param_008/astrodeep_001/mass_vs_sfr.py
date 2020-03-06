#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:26:47 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

size = 15
fsize = 15

# =============================================================================
# INPUT - get "real" parameters
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/param_008/astrodeep_001/mock_catalogue_005_010.fits'
data_fits = fits.open(fileName)

mtot_r = np.log10(data_fits['GALAXY PROPERTIES'].data['m_tot'])
sfr_r = data_fits['STAR FORMATION'].data['SFR']

data_fits.close()

# =============================================================================
# OUTPUT - get BEAGLE parameters
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/param_008/astrodeep_001/pyp-beagle/data/BEAGLE_summary_catalogue.fits'
data_fits = fits.open(fileName)

id_b = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1
mtot_b = data_fits['POSTERIOR PDF'].data['mass_mean']
mtot_68_b = data_fits['POSTERIOR PDF'].data['mass_68.00']
sfr_b = data_fits['STAR FORMATION'].data['SFR_mean']
sfr_68_b = data_fits['STAR FORMATION'].data['SFR_68.00']

data_fits.close()

# =============================================================================
# PLOT - input mass vs output mass
# =============================================================================


plt.figure(figsize=(fsize, fsize/2))
plt.title('Input Mass (DE) vs Output Mass (DPL)', size=size)
plt.xlabel(r'$\text{log}(M_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(M_{tot}/M_{\odot})$', size=size)
plt.plot((7.5, 11), (7.5, 11))
plt.scatter(mtot_r[id_b], mtot_b, s=10)
plt.errorbar(mtot_r[id_b], mtot_b, yerr=[mtot_b - mtot_68_b[:, 0], mtot_68_b[:, 1] - mtot_b], linestyle="None", elinewidth=0.5, color='k')

plt.xlim(7.5, 11)
plt.ylim(7.5, 11)
plt.show()

# =============================================================================
# PLOT - input sfr vs output sfr
# =============================================================================


plt.figure(figsize=(fsize, fsize/2))
plt.title('Input SFR (DE) vs Output SFR (DPL)', size=size)
plt.xlabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.plot((0, 4000), (0, 4000))
plt.scatter(sfr_r[id_b], sfr_b, s=10)
plt.errorbar(sfr_r[id_b], sfr_b, yerr=[sfr_b - sfr_68_b[:, 0], sfr_68_b[:, 1] - sfr_b], linestyle="None", elinewidth=0.5, color='k')

plt.xlim(0, 2000)
plt.ylim(0, 4000)
plt.show()



# =============================================================================
# PLOT - main sequence as a 2d histogram
# =============================================================================

# input values
plt.figure(figsize=(1.2*fsize, fsize))
plt.hist2d(mtot_r, np.log10(sfr_r), range=[[7.5, 11], [-1, 3.5]], bins=100)
plt.colorbar()
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()

massh = np.empty(0)
sfrh = np.empty(0)

for i in range(len(id_b)):
#for i in range(10):

    beagleData = fits.open('/Users/lester/Documents/PhD/param_008/astrodeep_001/{}_BEAGLE.fits'.format(id_b[i]+1))
    
    #needs float64 to provide precision needed for the random.choice weights
    temp_probs = np.float64(beagleData['POSTERIOR PDF'].data['probability'])
    temp_probs = temp_probs/np.sum(temp_probs)

    #here's the key line - take weighted samples from the multinest output!
    idx = np.random.choice(len(temp_probs), size=10000, p=temp_probs)
    massh = np.append(massh, np.log10(beagleData['GALAXY PROPERTIES'].data['M_tot'][idx]))
    sfrh = np.append(sfrh, np.log10(beagleData['STAR FORMATION'].data['sfr'][idx]))


plt.figure(figsize=(1.2*fsize, fsize))
plt.hist2d(massh, sfrh, range=[[7.5, 11], [-1, 3.5]], bins=100)
plt.colorbar()
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()


























































