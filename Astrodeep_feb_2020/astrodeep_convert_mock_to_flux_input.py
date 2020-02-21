#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 12:56:05 2020

@author: lester
"""

import numpy as np

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/fit/mock_catalogue_005_009.fits'
data_fits = fits.open(fileName)
id_B = np.asarray(data_fits[1].data['ID'], dtype=int)
z_B = data_fits[1].data['redshift_median']
z_Berr = data_fits[1].data['redshift_68.00']
#z_B = data_fits[1].data['redshift_mean']
data_fits.close()














# =============================================================================
# filters = ['HST_ACS_WFC_F435W', 'HST_ACS_WFC_F606W', 'HST_ACS_WFC_F814W', 'HST_WFC3_IR_F105W', 'HST_WFC3_IR_F125W', 'HST_WFC3_IR_F140W', 'HST_WFC3_IR_F160W', 'Paranal_HAWKI_Ks', 'Spitzer_IRAC_I1', 'Spitzer_IRAC_I2']
# =============================================================================

header = ['b_B435', 'b_errB435', 'b_V606', 'b_errV606', 'b_I814', 'b_errI814', 'b_Y105', 'b_errY105', 'b_J125', 'b_errJ125', 'b_JH140', 'b_errJH140', 'b_H160', 'b_errH160', 'b_Ks', 'b_errKs', 'b_CH1', 'b_errCH1', 'b_CH2', 'b_errCH2', 'ZBEST', 'field', 'ID']

header_string = '#ID b_B435 b_errB435 b_V606 b_errV606 b_I814 b_errI814 b_Y105 b_errY105 b_J125 b_errJ125 b_JH140 b_errJH140 b_H160 b_errH160 b_Ks b_errKs b_CH1 b_errCH1 b_CH2 b_errCH2 ZBEST field ID_original\n'



#f= open("astrodeep_ascii.txt","w+")  this was the original ~300, which doesn't include the ID to determine which galaxy
#f= open("astrodeep_ascii2.txt","w+")   this was a subset of the above of just first 3 to fit locally

f= open("astrodeep_ascii_002.txt","w+")  # this includes data to determine source galaxy, not used anywhere, can overwrite for now

f.write(header_string)

for i in range(len(D)):
    row = str(i+1)
    
    for j in range(len(header)):
        row = row + ' ' + str(D[header[j]][i])
        
    row = row + '\n'    
    f.write(row)

f.close()

















