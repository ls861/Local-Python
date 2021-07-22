#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:39:00 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from numpy import errstate,isneginf

'''
Takes the GMM fits files created from the BEAGLE output
Selects them based on a mass cut (already calculated from BEAGLE posterior medians)
Filters into a redshift bin
Filters by a maximum chi squared value (calculated elsewhere)
Amends GMMx and GMMy according to magnification in astrodeep catalog
Combines all fields into single linmix inputs and saves as .npy files
'''

#fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P']
fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P', '4M0717C', '5M0717P', '6M1149C', '7M1149P']
#fields = ['0A2744C']
runs = ['001']

massType = '_mStar'
#massType = '_mTot'

sfrType = '_delayed'
#sfrType = ''

#massLim = '8p5'
#redshift_center     = 2.0
#redshift_width      = 0.5 # either side, so 0.5 is 1.5 to 2.5
#chi2_max            = 9.5

massLim = '8p4'
redshift_center     = 1.65 # 2.5, 3.5, 4.5
redshift_width      = 0.35 # either side
chi2_max            = 9.5 # 2.5, 9.5
dimension           = '2d' # 2d or 3d - 3d includes redshift in the GMM fit
comment             = '201'


#massLim = '7p0'
#redshift_center     = 7.0
#redshift_width      = 10.0 # either side, so 0.5 is 1.5 to 2.5
#chi2_max            = 9999

#Hmm, it does look like the high mass objects might be affecting the sig0. This is where we would visually inspect the 4(?) Objects above the relation at high mass to see whether we believe them. Anything below sfr -2 could probably be removed from the fit. We would overplot other people's measured relations and check alpha and beta look decent (you should compare to santini). I'm worried that the objects above 10.5 at low see are where the relation starts to bend. For now I would not fit the 3 objs at high mass, low sfr, as well as all objs below -2 which we can categorically say are not in the relation and I would probably remove the 4 above the relation but that one's iffier and would benefit from visually inspecting the objects. Unfortunately we can't show the results for the scatter as they are...

zLim = str(redshift_center).replace(".","p")
wLim = str(redshift_width).replace(".","p")
chi2Lim = str(chi2_max).replace(".","p")

GMMid_arr           = []
pi_err_arr          = []
GMMx_arr            = []
GMMy_arr            = []
GMMxsig_arr         = []
GMMysig_arr         = []
GMMxycov_arr        = []

GMMchi2_arr         = []
GMMmass_arr         = []
GMMsfr_arr          = []
GMMredshift_arr     = []

if dimension == '3d':
    GMMz_arr            = []
    GMMzsig_arr         = []
    GMMxzcov_arr        = []
    GMMyzcov_arr        = []

fsize = 5
size = 8

for field in fields:

    for run in runs:

        # =============================================================================
        # GET DATA
        # =============================================================================

        directory = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/'

        # BEAGLE OUTPUT SUMMARY
        fileName = directory+'fields/{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(field[0])
        data_fits = fits.open(fileName)
#        print(data_fits.info())
#        print(data_fits[1].header)
        id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)

        with errstate(divide='ignore'):
            sfr_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_median'])
            sfr_68_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_68.00'])

        sfr_b1[isneginf(sfr_b1)]=-50
        sfr_68_b1[isneginf(sfr_68_b1)]=-50

        ssfr_b1 = data_fits['STAR FORMATION'].data['sSFR_median']
        ssfr_68_b1 = data_fits['STAR FORMATION'].data['sSFR_68.00']
        
        redshift_b1 = data_fits['POSTERIOR PDF'].data['redshift_median']
        mass_b1 = data_fits['POSTERIOR PDF'].data['mass_median']
        msa_b1 = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_median']
        tauV_eff_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_median']
        metallicity_b1 = data_fits['POSTERIOR PDF'].data['metallicity_median']
        tau_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_median']
        
        redshift_68_b1 = data_fits['POSTERIOR PDF'].data['redshift_68.00']
        mass_68_b1 = data_fits['POSTERIOR PDF'].data['mass_68.00']
        msa_68_b1 = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_68.00']
        tauV_eff_68_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_68.00']
        metallicity_68_b1 = data_fits['POSTERIOR PDF'].data['metallicity_68.00']
        tau_68_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_68.00']
        
        nebular_xi_b1 = np.full(len(id_b1), 0.3)
        nebular_xi_68_b1 = np.full((len(id_b1),2), 0.3)
        
        data_fits.close()

        # BEAGLE OUTPUT CHI SQUARED
        fileName = directory+'chi2/{}_chi2.fits'.format(field)
        data_fits = fits.open(fileName)
#        print(data_fits.info())
#        print(data_fits[1].header)
        id_chi2 = np.asarray(data_fits[1].data['id'], dtype=int)
        chi2 = data_fits[1].data['chi2']
        data_fits.close()

        # BEAGLE INPUT FLUXES - need to compare IDs
        fileName = directory+'data/astrodeep_{}_{}_subset_RF1_001.fits'.format(field[1:-1], field[-1])
        data_fits = fits.open(fileName)
#        print(data_fits.info())
#        print(data_fits[1].header)
        id_input = np.asarray(data_fits[1].data['ID'][id_b1-1], dtype=int)
        field_original = np.asarray(data_fits[1].data['field'][id_b1-1], dtype=int)
        id_original = np.asarray(data_fits[1].data['ID_original'][id_b1-1], dtype=int)
        zbest = data_fits[1].data['ZBEST'][id_b1-1]
        data_fits.close()

        # ASTRODEEP CATALOG - needed for magnification etc
        catalog = np.load(directory+'astrodeep_rawfile_1234_ABCZ.npy')
#        print(catalog.dtype.names)

        # GMM INPUTS (sorted by ID)
        sbf = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/linmix_inputs/linmix_inputs_GMM_{}_{}{}{}/'.format(dimension, field[0], massType, sfrType)

        GMMid       = np.asarray(np.load(sbf+massLim+'_'+'id.npy'), dtype=int)             # BEAGLE ID

        print(sbf+massLim+'_'+'id.npy')

        idx = np.argsort(GMMid)
        GMMid       = GMMid[idx]

        pi_err      = np.load(sbf+massLim+'_'+'pi_err.npy')[idx]            # 3x probability of each posterior gaussian
        GMMx        = np.load(sbf+massLim+'_'+'GMMx.npy')[idx]              # 3x posterior means per mass
        GMMy        = np.load(sbf+massLim+'_'+'GMMy.npy')[idx]              # 3x posterior means per sfr
        GMMxsig     = np.load(sbf+massLim+'_'+'GMMxsig.npy')[idx]           # 3x posterior sigmas per mass
        GMMysig     = np.load(sbf+massLim+'_'+'GMMysig.npy')[idx]           # 3x posterior sigmas per sfr
        GMMxycov    = np.load(sbf+massLim+'_'+'GMMxycov.npy')[idx]          # 3x posterior covar per mass-sfr pair

        nK          = np.load(sbf+massLim+'_'+'nK.npy')                     # 3 #gaussians modelling xi
        nGauss      = np.load(sbf+massLim+'_'+'nGauss.npy')                 # 3 #gaussians modelling BEAGLE posterior
        nChains     = np.load(sbf+massLim+'_'+'nChains.npy')                # 2
        minIter     = np.load(sbf+massLim+'_'+'minIter.npy')                # 3000
        maxIter     = np.load(sbf+massLim+'_'+'maxIter.npy')                # 3000

        if dimension == '3d':
            GMMz            = np.load(sbf+massLim+'_'+'GMMz.npy')[idx]
            GMMzsig         = np.load(sbf+massLim+'_'+'GMMzsig.npy')[idx]
            GMMxzcov        = np.load(sbf+massLim+'_'+'GMMxzcov.npy')[idx]
            GMMyzcov        = np.load(sbf+massLim+'_'+'GMMyzcov.npy')[idx]
        
        print(len(GMMx))


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

        print(len(GMMx))
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
            
        print(len(GMMx))

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

print(len(GMMx))
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
from os import mkdir
from os import path
sbf = './linmix_AD_combined/linmix_npy_files_z{}_w{}_chi{}_dim{}{}{}_{}/'.format(zLim, wLim, chi2Lim, dimension, massType, sfrType, comment)

if not path.exists(sbf):
    mkdir(sbf)

sbf = '{}{}_'.format(sbf, massLim)


# =============================================================================
# adding field information
# =============================================================================
k = 0
field_indicator = 0
field_test = []
for i in range(len(GMMid)):
    if GMMid[i] < k:
        field_indicator += 1
    field_test.append(field_indicator)
    k = GMMid[i]
field_test = np.array(field_test)
# =============================================================================
# 
# =============================================================================

np.save(sbf+'id', GMMid)
np.save(sbf+'field', field_test)
np.save(sbf+'pi_err', pi_err)
np.save(sbf+'GMMx', GMMx)
np.save(sbf+'GMMy', GMMy)
np.save(sbf+'GMMxsig', GMMxsig)
np.save(sbf+'GMMysig', GMMysig)
np.save(sbf+'GMMxycov', GMMxycov)
np.save(sbf+'nK', nK)
np.save(sbf+'nGauss', nGauss)
np.save(sbf+'nChains', nChains)
np.save(sbf+'minIter', minIter)
np.save(sbf+'maxIter', maxIter)

np.save(sbf+'GMMredshift', GMMredshift)
np.save(sbf+'GMMchi2', GMMchi2)
np.save(sbf+'GMMmass', GMMmass)
np.save(sbf+'GMMsfr', GMMsfr)

if dimension == '3d':
    np.save(sbf+'GMMz', GMMz)
    np.save(sbf+'GMMzsig', GMMzsig)
    np.save(sbf+'GMMxzcov', GMMxzcov)
    np.save(sbf+'GMMyzcov', GMMyzcov)



#%%
x = np.array([7.0,11.0])
m = 1.0
c = -7.01
y = m*x + c

plt.figure(figsize=(10,10))
plt.scatter(GMMmass, GMMsfr)
plt.plot(x, y)
plt.show()

plt.hist(GMMredshift)
plt.show()

plt.hist(GMMmass)
plt.show()




idx_line = GMMsfr > (m*GMMmass + c)
plt.scatter(GMMmass[idx_line], GMMsfr[idx_line])
plt.show()

# =============================================================================
# investigating the line
# =============================================================================

k = 0
field_indicator = 0
field_test = []

for i in range(len(GMMid)):
    if GMMid[i] < k:
        field_indicator += 1
    field_test.append(field_indicator)
    k = GMMid[i]
      
field_test = np.array(field_test)

print(GMMid)
print(field_test)
print(GMMid[idx_line])
print(field_test[idx_line])



# =============================================================================
# 
# =============================================================================
idx_line = GMMsfr < (m*GMMmass + c)
plt.scatter(GMMmass[idx_line], GMMsfr[idx_line])
plt.show()



















