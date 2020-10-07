#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:29:30 2020

@author: lester
"""


'''
THIS IS ALL JUST FOR FIELD 1 AT THE MOMENT

'''


import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import pickle
from numpy import errstate,isneginf
import sys


local = True

if local == True:
    fields = ['1A2744P']
else:
    fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P', '4M0717C', '5M0717P', '6M1149C', '7M1149P']
    
if local == True:
    output = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/scripts/sfr_instant_and_m_star/instant_sfr_medians_and_sample_outputs/'
else:
    output = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/kelly/instant_sfr_medians_and_sample_outputs/'

mass_sfr_options = ['_mTot', '_mTot_delayed', '_mStar', '_mStar_delayed']


samples = 100

xlow = 6.0
xhigh = 11.0
ylow = -3
yhigh = 3


# =============================================================================
# MS medians (BEAGLE SFR, instant SFR, Santini SFR)
# =============================================================================


for field in fields:
    
    
    if local == True:
        folder = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/fields/1/'
    else:
        folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/fit_001/'.format(field)

    if local == True:
        santini_folder = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/scripts/sfr_santini/results/{}_'.format(field[0])
    else:
        print('NOT CODED FOR THIS')
        sys.exit() 
        
        
    
    # BEAGLE OUTPUT SUMMARY
    fileName = folder+'pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(field[0])
    data_fits = fits.open(fileName)
    # print(data_fits.info())
    # print(data_fits[1].header)
    ID = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    
    with errstate(divide='ignore'):
        sfr = np.log10(data_fits['STAR FORMATION'].data['SFR_median'])
        sfr_68 = np.log10(data_fits['STAR FORMATION'].data['SFR_68.00'])

    sfr[isneginf(sfr)]=-30
    sfr_68[isneginf(sfr_68)]=-30

    mTot = np.log10(np.array(data_fits['GALAXY PROPERTIES'].data['M_tot_median'], np.float64))
    mTot_68 = np.log10(np.array(data_fits['GALAXY PROPERTIES'].data['M_tot_68.00'], np.float64))
    
    mStar = np.log10(np.array(data_fits['GALAXY PROPERTIES'].data['M_star_median'], np.float64))
    mStar_68 = np.log10(np.array(data_fits['GALAXY PROPERTIES'].data['M_star_68.00'], np.float64))    

    redshift = data_fits['POSTERIOR PDF'].data['redshift_median']
    
    data_fits.close()     
        
    ID_instant = np.array(pickle.load(open(output+'{}_instant_sfr_medians.p'.format(field[0]),"r"))['ID'], int)

    idx_sort = np.argsort(ID_instant)

    ID_instant = ID_instant[idx_sort]
    
    sfr_instant = pickle.load(open(output+'{}_instant_sfr_medians.p'.format(field[0]),"r"))['log_instant_sfr_median'][idx_sort]
    sfr_instant_68 = pickle.load(open(output+'{}_instant_sfr_medians.p'.format(field[0]),"r"))['log_instant_sfr_median_68_00'][idx_sort]

    ID_santini = np.load(santini_folder+'sfr_ID_santini.npy')
    sfr_santini = np.load(santini_folder+'sfr_santini.npy')
    sfr_santini_68 = np.load(santini_folder+'sfr_err_santini.npy')
    idx_santini = np.isin(ID_instant, ID_santini)

    plt.figure(figsize=(14, 14))
    plt.xlim(xlow, xhigh)
    plt.ylim(ylow, yhigh)
    plt.scatter(mStar, sfr, marker='x')
    plt.scatter(mStar[idx_santini], sfr[idx_santini], marker='x')
    plt.show()
    
    plt.figure(figsize=(14, 14))
#    plt.xlim(xlow, xhigh)
#    plt.ylim(ylow, yhigh)
    plt.scatter(mStar[idx_santini], sfr[idx_santini], marker='x', label='averaged')
    plt.scatter(mStar[idx_santini], sfr_instant[idx_santini], marker='x', label='instant')
    plt.scatter(mStar[idx_santini], sfr_santini, marker='x', label='santini')
    plt.legend()
    plt.show()
    

        
#    fileList = os.listdir(folder)
    


#    BEAGLE_samples = pickle.load(open(output+'{}_{}_BEAGLE_samples.p'.format(field[0], ID),"r"))
#    GMM_samples = pickle.load(open(output+'{}_{}_GMM_samples{}.p'.format(field[0], ID, mass_sfr_option),"r"))
    
    
    '''
    
    pickle.dump({'ID':np.array(ID_arr), 'log_instant_sfr_mean':np.array(mean_arr), 'log_instant_sfr_median.npy':np.array(median_arr), 'log_instant_sfr_median_68_00.npy':np.array(interval_arr)}, open(output+'{}_instant_sfr_medians.p'.format(field[0]),'w'))
    
    pickle.dump({'mTot':log_mTot, 'mStar':log_mStar, 'sfr':log_sfr, 'sfr_instant':log_sfr_instant}, open(output+'{}_{}_BEAGLE_samples.p'.format(field[0], ID),'w'))
    
    pickle.dump({'mass':hp_x, 'sfr':hp_y}, open(output+'{}_{}_GMM_samples{}.p'.format(field[0], ID, mass_sfr_option),'w'))
    
    '''
    


# =============================================================================
# get data for AD vs BEAGLE vs (Santini SFR)
# =============================================================================

if local == True:
    AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'
else:
    print('NOT CODED FOR THIS')
    sys.exit() 

AD = np.load(AD_location)

AD1 = AD[AD['field']==1]
print(AD1.dtype.names)

# BEAGLE INPUT FLUXES - need to compare IDs
fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/data/astrodeep_A2744_p_subset_RF1_001.fits'
data_fits = fits.open(fileName)
# print(data_fits.info())
# print(data_fits[1].header)
id_input = np.asarray(data_fits[1].data['ID'][ID-1], dtype=int)
field_original = np.asarray(data_fits[1].data['field'][ID-1], dtype=int)
id_original = np.asarray(data_fits[1].data['ID_original'][ID-1], dtype=int)
zbest = data_fits[1].data['ZBEST'][ID-1]
data_fits.close()

idx_AD1 = np.isin(AD1['ID'], id_original) # finds indicies where BEAGLE fitted values are

with errstate(divide='ignore'):
    mass_AD = np.log10(AD1['MSTAR']*1e9)[idx_AD1]
    mass_AD_neb = np.log10(AD1['MASTAR_NEB']*1e9)[idx_AD1]
    sfr_AD = np.log10(AD1['SFR'])[idx_AD1]
    sfr_AD_neb = np.log10(AD1['SFR_NEB'])[idx_AD1]

mass_AD[isneginf(mass_AD)]=-30
mass_AD_neb[isneginf(mass_AD_neb)]=-30
sfr_AD[isneginf(sfr_AD)]=-30
sfr_AD_neb[isneginf(sfr_AD_neb)]=-30

redshift_AD = AD1['ZBEST'][idx_AD1]
mag_AD = np.log10(AD1['MAGNIF'][idx_AD1])

# =============================================================================
# plot AD vs BEAGLE vs (Santini SFR)
# =============================================================================

idx_redshift = np.isin(ID_santini, (ID[idx_santini][abs(redshift[idx_santini]-1.65)<0.35]))
# all sfrs, masses and redshifts are matched to BEAGLE fitted subset, except santini sfr

'''
sfr
sfr_instant
sfr_AD
sfr_AD_neb
sfr_santini
'''
'''
mTot
mStar
mass_AD
mass_AD_neb
'''
'''
redshift
redshift_AD
zbest # same as AD
'''

xs = [sfr_instant, 
      sfr_instant, 
      sfr_instant, 
      sfr_instant[idx_santini],
      sfr_instant[idx_santini][idx_redshift],
      sfr_instant[idx_santini][idx_redshift],
      sfr_instant[idx_santini][idx_redshift],
      sfr_instant[idx_santini][idx_redshift],
      mass_AD,
      mass_AD[idx_santini][idx_redshift],
      (mass_AD-mag_AD)[idx_santini][idx_redshift],
      redshift_AD,
      redshift_AD[idx_santini][idx_redshift],
      sfr_AD_neb[idx_santini], 
      sfr_AD_neb[idx_santini][idx_redshift],
      redshift_AD,
      mass_AD_neb,
      sfr_AD_neb]

ys = [sfr,
      sfr_AD,
      sfr_AD_neb,
      sfr_santini,
      sfr[idx_santini][idx_redshift],
      sfr_AD[idx_santini][idx_redshift],
      sfr_AD_neb[idx_santini][idx_redshift],
      sfr_santini[idx_redshift],
      mStar,
      mStar[idx_santini][idx_redshift],
      (mStar-mag_AD)[idx_santini][idx_redshift],
      redshift,
      redshift[idx_santini][idx_redshift],
      sfr_santini, 
      sfr_santini[idx_redshift],
      redshift,
      mStar,
      sfr_instant]

print(mStar)
print(mStar-mag_AD)
print(mag_AD)
xlow = [-4, -4, -4, -4, -4, -4, -4, -4, 6, 6, 6, 0, 0, -4, -4, 0, 6, -4]
xhigh = [4, 4, 4, 4, 4, 4, 4, 4, 10, 10, 10, 10, 10, 4, 4, 10, 10, 4]
ylow = [-4, -4, -4, -4, -4, -4, -4, -4, 6, 6, 6, 0, 0, -4, -4, 0, 6, -4]
yhigh = [4, 4, 4, 4, 4, 4, 4, 4, 10, 10, 10, 10, 10, 4, 4, 10, 10, 4]

xlabel = ['sfr instant',
          'sfr instant',
          'sfr instant',
          'sfr instant',
          'sfr instant',
          'sfr instant',
          'sfr instant',
          'sfr instant',
          'massAD',
          'massAD',
          'massAD',
          'redshiftAD',
          'redshiftAD',
          'sfr AD neb',
          'sfr AD neb',
          'redshift AD',
          'mass AD neb',
          'sfr AD neb']

ylabel = ['sfr',
          'sfr AD',
          'sfr AD neb',
          'sfr santini',
          'sfr',
          'sfr AD',
          'sfr AD neb',
          'sfr santini',
          'mStar',
          'mStar',
          'mStar',
          'redshift',
          'redshift',
          'sfr santini',
          'sfr santini',
          'redshift',
          'mStar',
          'sfr instant']

title = ['sfr',
         'sfr',
         'sfr',
         'sfr',
         'sfr 1.3$<$z$<$2.0',
         'sfr 1.3$<$z$<$2.0',
         'sfr 1.3$<$z$<$2.0',
         'sfr 1.3$<$z$<$2.0',
         'mStar',
         'mStar 1.3$<$z$<$2.0',
         'mStar w MAG 1.3$<$z$<$2.0',
         'redshift',
         'redshift 1.3$<$z$<$2.0', 
         'sfr',
         'sfr 1.3$<$z$<$2.0',
         'redshift',
         'mStar',
         'sfr']


for i in range(len(xs)):

    plt.hist2d(xs[i], ys[i], bins=30, range=[[xlow[i], xhigh[i]], [ylow[i], yhigh[i]]])
    plt.xlabel(xlabel[i])
    plt.ylabel(ylabel[i])
    plt.title(title[i])
    plt.plot((xlow[i], xhigh[i]), (ylow[i], yhigh[i]), color='k')
    plt.show()
    
    plt.scatter(xs[i], ys[i])
    plt.xlim(xlow[i], xhigh[i])
    plt.ylim(ylow[i], yhigh[i])
    plt.xlabel(xlabel[i])
    plt.ylabel(ylabel[i])
    plt.plot((xlow[i], xhigh[i]), (ylow[i], yhigh[i]), color='k')
    plt.title(title[i])
    plt.show()







# =============================================================================
# BEAGLE MS heatplot
# =============================================================================














# =============================================================================
# GMM MS heatplot
# =============================================================================








# =============================================================================
# object sample, BEAGLE and GMM
# =============================================================================
















