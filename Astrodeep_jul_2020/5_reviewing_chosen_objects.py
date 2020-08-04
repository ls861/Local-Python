#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:39:00 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P', '4M0717C', '5M0717P', '6M1149C', '7M1149P']
#fields = ['0A2744C']
runs = ['001']

directory = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/'

# ASTRODEEP CATALOG - needed for magnification etc
catalog = np.load(directory+'astrodeep_rawfile_1234_ABCZ.npy')
print(catalog.dtype.names)
        

plt.hist(catalog['H160'], bins=50, range=[15,35])
plt.show()

GMMH160_arr = []

massLim             = '8p5'
dimension           = '2d' # 2d or 3d - 3d includes redshift in the GMM fit

fsize = 5
size = 8

for field in fields:

    for run in runs:

        # =============================================================================
        # GET DATA
        # =============================================================================
        
        # BEAGLE OUTPUT SUMMARY
        fileName = '{}fields/{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(directory, field[0])
        data_fits = fits.open(fileName)
#        print(data_fits.info())
#        print(data_fits[1].header)
        id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)        
        data_fits.close()

        # BEAGLE INPUT FLUXES - need to compare IDs
        fileName = '{}data/astrodeep_{}_{}_subset_RF1_001.fits'.format(directory, field[1:-1], field[-1])
        data_fits = fits.open(fileName)
        id_input = np.asarray(data_fits[1].data['ID'][id_b1-1], dtype=int)
        field_original = np.asarray(data_fits[1].data['field'][id_b1-1], dtype=int)
        id_original = np.asarray(data_fits[1].data['ID_original'][id_b1-1], dtype=int)
        zbest = data_fits[1].data['ZBEST'][id_b1-1]
        data_fits.close()

        # GMM INPUTS (sorted by ID)
        fileName = '{}linmix_inputs/linmix_inputs_GMM_{}_{}/{}_id.npy'.format(directory, dimension, field[0], massLim)
        GMMid       = np.asarray(np.load(fileName), dtype=int)             # BEAGLE ID

        print(len(GMMid))


        '''
        catalog['ID'] is AD ID for every object, ALL IMAGES
        id_original is AD ID for fitted objects
        
        id_input is BEAGLE ID for fitted objects
        id_b1 is BEAGLE ID for fitted objects
        GMMid is BEAGLE ID for GMM objects, inc mass cut
        
        field_original is field for fitted objects
        catalog['field'] is field for every object
        
        I want id_original which corresponds to GMMid
        
        id_original and id_input are from same file
        
        Need:
            1. location of GMMid in id_input
            2. take same subset of id_original
            3. get H band fluxes based on these from catalog
        
        '''
        
#        print('LENGTHS', len(catalog['ID']), len(id_original), len(id_input), len(id_b1), len(GMMid), field_original, catalog['field'])
        
        test1 = np.isin(id_input, GMMid)
        test2 = id_original[test1]

        f = (catalog['field'] == int(field[0]))
        
        test3 = np.isin(catalog['ID'][f], test2)
             
#        plt.hist(catalog['H160'][f][test3], bins=50, range=[15,35])
#        plt.show()

        for i in range(len(catalog['H160'][f][test3])):
            GMMH160_arr.append(catalog['H160'][f][test3][i])
        
plt.hist(GMMH160_arr, bins=50, range=[15,35])
plt.show()

        
        
plt.hist(catalog['H160'], bins=50, range=[15,35]) 
plt.hist(GMMH160_arr, bins=50, range=[15,35])
plt.show() 
        
        
'''
      
        
        # =============================================================================
        # filter by redshift
        # =============================================================================

        idx_GMM = np.isin(id_b1, GMMid)
        GMMredshift = redshift_b1[idx_GMM]

        idx_z = (abs(GMMredshift - redshift_center) < redshift_width)

        GMMid = GMMid[idx_z]
        pi_err = pi_err[idx_z]
        GMMx = GMMx[idx_z]
        GMMy = GMMy[idx_z]
        GMMxsig = GMMxsig[idx_z]
        GMMysig = GMMysig[idx_z]
        GMMxycov = GMMxycov[idx_z]
        GMMredshift = GMMredshift[idx_z]

        GMMmass = mass_b1[idx_GMM]
        GMMmass = GMMmass[idx_z]

        GMMsfr = sfr_b1[idx_GMM]
        GMMsfr = GMMsfr[idx_z]
        
        if dimension == '3d':
            GMMz            = GMMz[idx_z]
            GMMzsig         = GMMzsig[idx_z]
            GMMxzcov        = GMMxzcov[idx_z]
            GMMyzcov        = GMMyzcov[idx_z]


        # =============================================================================
        # filter by chi2 (id_b1 == id_chi2 == id_input)
        # =============================================================================

        idx_GMM = np.isin(id_chi2, GMMid)
        GMMchi2 = chi2[idx_GMM]

        idx_chi2 = (GMMchi2 < chi2_max)

        GMMid = GMMid[idx_chi2]
        pi_err = pi_err[idx_chi2]
        GMMx = GMMx[idx_chi2]
        GMMy = GMMy[idx_chi2]
        GMMxsig = GMMxsig[idx_chi2]
        GMMysig = GMMysig[idx_chi2]
        GMMxycov = GMMxycov[idx_chi2]
        GMMredshift = GMMredshift[idx_chi2]
        GMMchi2 = GMMchi2[idx_chi2]

        # adding mass
        GMMmass = GMMmass[idx_chi2]
        GMMsfr = GMMsfr[idx_chi2]
        
        if dimension == '3d':
            GMMz = GMMz[idx_chi2]
            GMMzsig = GMMzsig[idx_chi2]
            GMMxzcov = GMMxzcov[idx_chi2]
            GMMyzcov = GMMyzcov[idx_chi2]

        # =============================================================================
        # include magnification
        # =============================================================================

        idx_GMM = np.isin(id_input, GMMid)
        GMMid_original = id_original[idx_GMM]
        GMMfield_original = field_original[idx_GMM]

        GMMcatalog = catalog[catalog['field']==GMMfield_original[0]]

        idx_mag = np.isin(GMMcatalog['ID'], GMMid_original)
        GMMmag = np.log10(GMMcatalog['MAGNIF'][idx_mag])

        GMMmass = GMMmass - GMMmag
        GMMsfr = GMMsfr - GMMmag

        GMMmag = np.array([GMMmag]*3).transpose() # this is genius

        GMMx = GMMx - GMMmag
        GMMy = GMMy - GMMmag

        # =============================================================================
        # Combining all fields into single linmix input
        # =============================================================================

        GMMid_arr.append(GMMid)
        pi_err_arr.append(pi_err)
        GMMx_arr.append(GMMx)
        GMMy_arr.append(GMMy)
        GMMxsig_arr.append(GMMxsig)
        GMMysig_arr.append(GMMysig)
        GMMxycov_arr.append(GMMxycov)

        GMMredshift_arr.append(GMMredshift)
        GMMchi2_arr.append(GMMchi2)
        GMMmass_arr.append(GMMmass)
        GMMsfr_arr.append(GMMsfr)
        
        if dimension == '3d':
            GMMz_arr.append(GMMz)
            GMMzsig_arr.append(GMMzsig)
            GMMxzcov_arr.append(GMMxzcov)
            GMMyzcov_arr.append(GMMyzcov)     
            

        # for visual inspection
        print(id_input[np.isin(id_input, GMMid[GMMsfr > 2.1])])
        print(id_original[np.isin(id_input, GMMid[GMMsfr > 2.1])])

GMMid       = np.concatenate(GMMid_arr)
pi_err      = np.concatenate(pi_err_arr)
GMMx        = np.concatenate(GMMx_arr)
GMMy        = np.concatenate(GMMy_arr)
GMMxsig     = np.concatenate(GMMxsig_arr)
GMMysig     = np.concatenate(GMMysig_arr)
GMMxycov    = np.concatenate(GMMxycov_arr)

GMMchi2     = np.concatenate(GMMchi2_arr)
GMMmass     = np.concatenate(GMMmass_arr)
GMMsfr      = np.concatenate(GMMsfr_arr)
GMMredshift = np.concatenate(GMMredshift_arr)

if dimension == '3d':
    GMMz = np.concatenate(GMMz_arr)
    GMMzsig = np.concatenate(GMMzsig_arr)
    GMMxzcov = np.concatenate(GMMxzcov_arr)
    GMMyzcov = np.concatenate(GMMyzcov_arr)

# =============================================================================
# RANDOM CUTS
# =============================================================================

plt.scatter(GMMmass, GMMsfr)
plt.show()

# removing SFR < -2
idx1 = (GMMsfr < -2.0)

# removing post magnification mass > 10.5
idx2 = (GMMmass > 10.5)

# removing anomalously high sfr at high mass
idx3 = (GMMsfr > 2.1)

print(len(GMMx[idx1]), len(GMMx[idx2]), len(GMMx[idx3]))

random_cuts = np.invert(np.logical_or.reduce((idx1, idx2, idx3)))

GMMid       = GMMid[random_cuts]
pi_err      = pi_err[random_cuts]
GMMx        = GMMx[random_cuts]
GMMy        = GMMy[random_cuts]
GMMxsig     = GMMxsig[random_cuts]
GMMysig     = GMMysig[random_cuts]
GMMxycov    = GMMxycov[random_cuts]

GMMredshift = GMMredshift[random_cuts]
GMMchi2     = GMMchi2[random_cuts]
GMMmass     = GMMmass[random_cuts]
GMMsfr      = GMMsfr[random_cuts]

if dimension == '3d':
    GMMz = GMMz[random_cuts]
    GMMzsig = GMMzsig[random_cuts]
    GMMxzcov = GMMxzcov[random_cuts]
    GMMyzcov = GMMyzcov[random_cuts]

print(len(GMMx))


plt.scatter(GMMmass, GMMsfr)
plt.show()

plt.hist(GMMredshift)
plt.show()

'''
