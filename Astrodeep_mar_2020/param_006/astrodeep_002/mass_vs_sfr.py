#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:26:47 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from astropy.io import fits

title = 'DE'
param = '006'
revision = '002'

#title = 'DE NO MinRelError'
#param = '006'
#revision = '003'

#title = 'LE'
#param = '007'
#revision = '001'

#title = 'DPL'
#param = '008'
#revision = '001'

size = 15
fsize = 10

# =============================================================================
# INPUT - get "real" parameters
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/param_{}/astrodeep_{}/mock_catalogue_005_010.fits'.format(param, revision)
data_fits = fits.open(fileName)

mtot_r = np.log10(data_fits['GALAXY PROPERTIES'].data['m_tot'])
sfr_r = data_fits['STAR FORMATION'].data['SFR']
msa_r = np.log10(data_fits['STAR FORMATION'].data['max_stellar_age'])
tau_r = np.log10(data_fits['STAR FORMATION BINS'].data['bin_tau'])

data_fits.close()

# =============================================================================
# OUTPUT - get BEAGLE parameters
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/param_{}/astrodeep_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(param, revision)
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

plt.figure(figsize=(fsize, fsize))
plt.title('Input Mass (DE) vs Output Mass ({})'.format(title), size=size)
plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)
plt.plot((7.5, 11), (7.5, 11))
plt.scatter(mtot_r[id_b], mtot_b, s=10, c=msa_r[id_b]-tau_r[id_b], zorder=10)
plt.errorbar(mtot_r[id_b], mtot_b, yerr=[mtot_b - mtot_68_b[:, 0], mtot_68_b[:, 1] - mtot_b], linestyle="None", elinewidth=0.5, color='k')

plt.colorbar()

plt.xlim(7.5, 11)
plt.ylim(7.5, 11)
plt.show()

# =============================================================================
# PLOT - input sfr vs output sfr
# =============================================================================

plt.figure(figsize=(fsize, fsize))
plt.title('Input SFR (DE) vs Output SFR ({})'.format(title), size=size)
plt.xlabel(r'$\text{Input - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.ylabel(r'$\text{Output - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.plot((-1, 3.5), (-1, 3.5))
plt.scatter(np.log10(sfr_r[id_b]), np.log10(sfr_b), s=10)
plt.errorbar(np.log10(sfr_r[id_b]), np.log10(sfr_b), yerr=[np.log10(sfr_b / sfr_68_b[:, 0]), np.log10(sfr_68_b[:, 1] / sfr_b)], linestyle="None", elinewidth=0.5, color='k')

plt.xlim(-1, 3.5)
plt.ylim(-1, 3.5)
plt.show()

# =============================================================================
# PLOT - main sequence as a 2d histogram
# =============================================================================

# input values
plt.figure(figsize=(1.2*fsize, fsize))
plt.title('INPUT (DE) SFR vs Mass', size=size)
plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.hist2d(mtot_r, np.log10(sfr_r), range=[[7.5, 11], [-1, 3.5]], bins=100)
plt.colorbar()
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()

massh = np.empty(0)
sfrh = np.empty(0)

#for i in range(len(id_b)):
for i in range(10):

    beagleData = fits.open('/Users/lester/Documents/PhD/param_{}/astrodeep_{}/{}_BEAGLE.fits'.format(param, revision, id_b[i]+1))
    
    #needs float64 to provide precision needed for the random.choice weights
    temp_probs = np.float64(beagleData['POSTERIOR PDF'].data['probability'])
    temp_probs = temp_probs/np.sum(temp_probs)

    #here's the key line - take weighted samples from the multinest output!
    idx = np.random.choice(len(temp_probs), size=1000, p=temp_probs)
    massh = np.append(massh, np.log10(beagleData['GALAXY PROPERTIES'].data['M_tot'][idx]))
    sfrh = np.append(sfrh, np.log10(beagleData['STAR FORMATION'].data['sfr'][idx]))


plt.figure(figsize=(1.2*fsize, fsize))
plt.title('FITTED ({}) SFR vs Mass'.format(title), size=size)
plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.hist2d(massh, sfrh, range=[[7.5, 11], [-1, 3.5]], bins=100, norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
plt.colorbar()

cmap = cm.get_cmap('viridis')
rgba = cmap(0)
ax = plt.axes()
ax.set_facecolor(rgba)

plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()


















































