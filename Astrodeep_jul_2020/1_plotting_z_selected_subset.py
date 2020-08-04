#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 16:50:52 2020

@author: lester
"""

import numpy as np
from numpy import errstate,isneginf
import matplotlib.pyplot as plt
from astropy.io import fits

# =============================================================================
# Plotting redshift selected subset for new BEAGLE test (logU changes)
# =============================================================================

field = '0A2744C'

# =============================================================================
# GET DATA from entire field, previously fitted
# =============================================================================

directory = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jun_2020/from_cluster/'

# BEAGLE OUTPUT SUMMARY
fileName = directory+'{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(field[0])
data_fits = fits.open(fileName)
#        print(data_fits.info())
#        print(data_fits[1].header)
id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)    

with errstate(divide='ignore'):
    sfr_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_median'])
    sfr_68_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_68.00'])

sfr_b1[isneginf(sfr_b1)]=-40
sfr_68_b1[isneginf(sfr_68_b1)]=-40
  
mass_b1 = data_fits['POSTERIOR PDF'].data['mass_median']    
mass_68_b1 = data_fits['POSTERIOR PDF'].data['mass_68.00']

data_fits.close()        

# BEAGLE INPUT FLUXES - need to compare IDs
fileName = directory+'{}/astrodeep_{}_{}_subset_RF1_001.fits'.format(field[0], field[1:-1], field[-1].lower()) 
data_fits = fits.open(fileName)
#        print(data_fits.info())
#        print(data_fits[1].header)
id_input = np.asarray(data_fits[1].data['ID'][id_b1-1], dtype=int)
field_original = np.asarray(data_fits[1].data['field'][id_b1-1], dtype=int)
id_original = np.asarray(data_fits[1].data['ID_original'][id_b1-1], dtype=int)
zbest = data_fits[1].data['ZBEST'][id_b1-1]
data_fits.close()

# =============================================================================
# match IDs
# =============================================================================

ff = "astrodeep_A2744_c_subset_RF1_002"
D_RF1 = np.load(ff+'.npy')
print(D_RF1['ID'])

idx = np.isin(id_original, D_RF1['ID'])


# =============================================================================
# PLOT MAIN SEQUENCE
# =============================================================================

x = mass_b1[idx]
y = sfr_b1[idx]

plt.scatter(x, y)








































