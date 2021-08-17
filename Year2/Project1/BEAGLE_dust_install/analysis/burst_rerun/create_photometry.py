#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 22:12:05 2021

@author: lester
"""

import numpy as np
from astropy.io import fits

# =============================================================================
# get redshift investigation BEAGLE IDs for field 2
# =============================================================================
fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation/subset_to_investigate.fits'
d = fits.open(fileName)
#print(d.info())
#print(d['Joined'].columns)
subset_IDs = d['Joined'].data['id_BEAGLE']
d.close()

# =============================================================================
# get redshift investigation BEAGLE IDs for field 2 for ON RELATION SUBSET
# =============================================================================
fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation/subset_to_investigate_on_relation.fits'
d_on_relation = fits.open(fileName)
#print(d.info())
#print(d['Joined'].columns)
subset_IDs_on_relation = d_on_relation['Joined'].data['id_BEAGLE']
d_on_relation.close()

# =============================================================================
# open original photometry
# =============================================================================
fileName = './astrodeep_field_2_subset_RF1_001.fits'
original_photometry = fits.open(fileName)
#print(original_photometry.info())
#print(original_photometry[1].columns)
#print(original_photometry[1].data['b_B435'])
original_photometry = original_photometry[1].data

# =============================================================================
# output into an ASCII file which can then be converted to FITS
# =============================================================================
header = ['b_B435', 'b_errB435', 'b_V606', 'b_errV606', 'b_I814', 'b_errI814', 'b_Y105', 'b_errY105', 'b_J125', 'b_errJ125', 'b_JH140', 'b_errJH140', 'b_H160', 'b_errH160', 'b_Ks', 'b_errKs', 'b_CH1', 'b_errCH1', 'b_CH2', 'b_errCH2', 'ZBEST', 'field', 'ID_original']

header_string = '#ID b_B435 b_errB435 b_V606 b_errV606 b_I814 b_errI814 b_Y105 b_errY105 b_J125 b_errJ125 b_JH140 b_errJH140 b_H160 b_errH160 b_Ks b_errKs b_CH1 b_errCH1 b_CH2 b_errCH2 ZBEST field ID_original\n'

#f = open('./BEAGLE_inputs/astrodeep_field_2_subset_RF1_001_exponential_redshift_prior.txt','w+')
#g = open('./BEAGLE_inputs/astrodeep_field_2_subset_RF1_001_tauv_eff_0_to_6.txt','w+')
#h = open('./BEAGLE_inputs/astrodeep_field_2_subset_RF1_001_excluding_IRAC_filters.txt','w+')

#f.write(header_string)
#g.write(header_string)
#h.write(header_string)

for i in range(len(original_photometry['ID'])):

    if np.isin(original_photometry['ID'][i], subset_IDs):
        row_f = str(i+1)
        row_g = str(i+1)
        row_h = str(i+1)

        for j in range(len(header)):
            row_f = row_f + ' ' + str(original_photometry[header[j]][i])
            row_g = row_g + ' ' + str(original_photometry[header[j]][i])

            if header[j] in ['b_CH1', 'b_errCH1', 'b_CH2', 'b_errCH2']:
                row_h = row_h + ' -67.0'
            else:
                row_h = row_h + ' ' + str(original_photometry[header[j]][i])


        row_f = row_f + '\n'
        row_g = row_g + '\n'
        row_h = row_h + '\n'

#        f.write(row_f)
#        g.write(row_g)
#        h.write(row_h)

#f.close()
#g.close()
#h.close()


# =============================================================================
# ON RELATION
# =============================================================================
#o = open('./BEAGLE_inputs/astrodeep_field_2_subset_RF1_001_exponential_redshift_prior_on_relation.txt','w+')
#o.write(header_string)

for i in range(len(original_photometry['ID'])):

    if np.isin(original_photometry['ID'][i], subset_IDs_on_relation):
        row_o = str(i+1)

        for j in range(len(header)):
            row_o = row_o + ' ' + str(original_photometry[header[j]][i])

        row_o = row_o + '\n'
#        o.write(row_o)

#o.close()
