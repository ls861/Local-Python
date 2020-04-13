#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 12:56:05 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# =============================================================================
# CALCULATE FLUXES FROM MAGNITUDE FITS FILE
# =============================================================================

##column name from apparent mag table in BEAGLE mock (same as left column in config + _APP)
#filters = ['F090W_NRC_and_OTE_APP', 'F115W_NRC_and_OTE_APP', 'F150W_NRC_and_OTE_APP', 'F200W_NRC_and_OTE_APP', 'F277W_NRC_and_OTE_APP', 'F356W_NRC_and_OTE_APP', 'F444W_NRC_and_OTE_APP']
#
##column name from JADES config file + NOT ACTUALLY USED HERE
#filter_label = ['F090W', 'F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F444W']

#column name from apparent mag table in BEAGLE mock (same as left column in config + _APP)
filters = ['HST_ACS_WFC_F435W_APP', 'HST_ACS_WFC_F606W_APP', 'HST_ACS_WFC_F814W_APP', 'HST_WFC3_IR_F105W_APP', 'HST_WFC3_IR_F125W_APP', 'HST_WFC3_IR_F140W_APP', 'HST_WFC3_IR_F160W_APP', 'Paranal_HAWKI_Ks_APP', 'Spitzer_IRAC_I1_APP', 'Spitzer_IRAC_I2_APP']

#column name from ASTRODEEP config file + NOT ACTUALLY USED HERE
filter_label = ['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2']

fileName = '/Users/lester/Documents/PhD/param_100/mock_100_LE/mock_catalogue_100_LE.fits'
data_fits = fits.open(fileName)

m = []
for i in range(len(filters)):
    m.append(np.array(data_fits['APPARENT MAGNITUDES'].data[filters[i]]))
    
data_fits.close()

f = []
for i in range(len(filters)):
    f.append(10**( (23.9 - m[i]) / 2.5 )) #mJy (in frequency not wavelength)

# =============================================================================
# FIND ERRORS
# =============================================================================
    
#sig_mag5 = np.array([29.4, 29.6, 29.7, 29.8, 29.4, 29.4, 29.1])    
#sig_mag5 = np.full(len(filters), 30.0) # original for 012_001
sig_mag5 = np.array([28.95, 29.35, 28.84, 28.45, 28.34, 28.34, 28.16, 26.45, 26.52, 26.25])    # CANDELS for 012_010, NOTE, 140 wasn't in table, hence duplicated 125.

sig_flux5 = 10**( (23.9 - sig_mag5) / 2.5 )     # microJy
sig_flux = sig_flux5 / 5.

print(sig_flux)
    
# =============================================================================
# Create ASCII file    
# =============================================================================

#header_string = '#ID F090W F115W F150W F200W F277W F356W F444W errF090W errF115W errF150W errF200W errF277W errF356W errF444W\n'

header_string = '#ID b_B435 b_V606 b_I814 b_Y105 b_J125 b_JH140 b_H160 b_Ks b_CH1 b_CH2 b_errB435 b_errV606 b_errI814 b_errY105 b_errJ125 b_errJH140 b_errH160 b_errKs b_errCH1 b_errCH2\n'

# =============================================================================
#file1 = open("DPL_mock_fluxes.txt","w+") # original 012_001
#file1 = open("100_LE_mock_fluxes_CANDELS.txt","w+") # CANDELS 012_010
#file1.write(header_string)
#
#for i in range(len(m[0])):
#    row = str(i+1)
#    
#    for j in range(len(filters)):
#        row = row + ' ' + str(f[j][i] + np.random.normal(loc=0.0, scale=sig_flux[j]))
#    
#    for j in range(len(sig_flux)):
#        row = row + ' ' + str(sig_flux[j])
#
#    row = row + '\n'    
#    file1.write(row)
#
#file1.close()
# =============================================================================














