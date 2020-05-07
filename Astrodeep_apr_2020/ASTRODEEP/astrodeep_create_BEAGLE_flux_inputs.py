#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 23:24:03 2020

@author: lester
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('figure', figsize=(4,4))

'''

A2744_c    0
A2744_p    1 
M0416_c    2
M0416_p    3
M0717_c    4
M0717_p    5
M1149_c    6
M1149_p    7

'''

# dtype=[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]
    
#D = np.load('astrodeep_rawfile.npy')

#('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8'), ('c_ID', '<f8'), ('c_X', '<f8'), ('c_Y', '<f8'), ('c_XMIN', '<f8'), ('c_YMIN', '<f8'), ('c_XMAX', '<f8'), ('c_YMAX', '<f8'), ('c_CLASS_STAR', '<f8'), ('c_SEXFLAG', '<f8'), ('c_RESFLAG', '<f8'), ('c_VISFLAG', '<f8'), ('c_TPHOTFLAG_Ks', '<f8'), ('c_COVMAX_Ks', '<f8'), ('c_TPHOFLAG_IRAC1', '<f8'), ('c_COVMAX_IRAC1', '<f8'), ('c_TPHOTFLAG_IRAC2', '<f8'), ('c_COVMAX_IRAC2', '<f8')])

# =============================================================================
# defining some useful lists
# =============================================================================

fluxes = ['b_B435','b_V606','b_I814','b_Y105','b_J125','b_JH140','b_H160','b_Ks','b_CH1','b_CH2']
flux_errs = ['b_errB435','b_errV606','b_errI814','b_errY105','b_errJ125','b_errJH140','b_errH160','b_errKs','b_errCH1','b_errCH2']
flux_covs = ['b_Ks','b_CH1','b_CH2']
flux_cov_errs = ['b_errKs','b_errCH1','b_errCH2']
covs = ['c_COVMAX_Ks', 'c_COVMAX_IRAC1', 'c_COVMAX_IRAC2']

# =============================================================================
# read in the file + get some basic info
# =============================================================================

D = np.load('astrodeep_rawfile_ABCZ.npy')

print('length of input catalog', len(D))
print('relflag == 1', len(D[D['RELFLAG']==1]))
print('relflag == 0', len(D[D['RELFLAG']==0]))

# =============================================================================
# begin to select a subset
# =============================================================================

# select RELFLAG ==1
D_RF1 = D[D['RELFLAG']==1]

# select a field as a subset if needed 
D_RF1 = D_RF1[D_RF1 ['field']==0]

print('length of subset (field & relflag)', len(D_RF1))

for i, flux in enumerate(fluxes):
    flux_err = flux_errs[i]

    # if flux == 0 make flux_err == -66
    D_RF1[flux_err][D_RF1[flux]==0] = -66

    # if flux == 0 make flux == -66
    D_RF1[flux][D_RF1[flux]==0] = -66
    

for i, flux_cov in enumerate(flux_covs):
    cov = covs[i]
    flux_cov_err = flux_cov_errs[i]
    
    # test to see whether Ks, CH1 or CH2 had any flux==0, they did not in first 4 fields
    print('COV TEST', flux_cov)
    print('# -66s, flux==0 rejections', len(D_RF1[D_RF1[flux_cov]==-66]))    

    # if covmax >= 1 make flux_err == -67
    D_RF1[flux_cov_err][D_RF1[cov]>=1] = -67

    # if covmax >= 1 make flux == -67
    D_RF1[flux_cov][D_RF1[cov]>=1] = -67

    print('# -67s, covmax rejections', len(D_RF1[D_RF1[flux_cov]==-67]))    

# =============================================================================
# plotting some histograms
# =============================================================================

#for i, flux in enumerate(fluxes):
#    flux_err = flux_errs[i]
#
#    plt.title(flux)
#    plt.hist(D_RF1[flux][D_RF1[flux]>-65], bins=100, range=(-12, 12))
#    plt.show()
#    plt.title(flux_err)
#    plt.hist(D_RF1[flux_err][D_RF1[flux_err]>-65], bins=100)
#    plt.show()


# =============================================================================
# min and max fluxes and errors remaining in subset
# =============================================================================
    
for i, flux in enumerate(fluxes):
    flux_err = flux_errs[i]

#    print(min(D_RF1[flux_err][D_RF1[flux_err]>-65]))

for i, flux in enumerate(fluxes):
    flux_err = flux_errs[i]

#    print(min(D_RF1[flux][D_RF1[flux]>-65]))
    
for i, flux in enumerate(fluxes):
    flux_err = flux_errs[i]

#    print(max(D_RF1[flux_err][D_RF1[flux_err]>-65]))

for i, flux in enumerate(fluxes):
    flux_err = flux_errs[i]

#    print(max(D_RF1[flux][D_RF1[flux]>-65]))

# =============================================================================
# seeing how many objects have -67 set for ks, ch1 and ch2
# =============================================================================
    
print(D_RF1['ID'][(D_RF1[flux_covs[0]]==-67)&(D_RF1[flux_covs[1]]==-67)&(D_RF1[flux_covs[2]]==-67)])
print(D_RF1['field'][(D_RF1[flux_covs[0]]==-67)&(D_RF1[flux_covs[1]]==-67)&(D_RF1[flux_covs[2]]==-67)])



# =============================================================================
# output into an ASCII file which can then be converted to FITS
# =============================================================================

header = ['b_B435', 'b_errB435', 'b_V606', 'b_errV606', 'b_I814', 'b_errI814', 'b_Y105', 'b_errY105', 'b_J125', 'b_errJ125', 'b_JH140', 'b_errJH140', 'b_H160', 'b_errH160', 'b_Ks', 'b_errKs', 'b_CH1', 'b_errCH1', 'b_CH2', 'b_errCH2', 'ZBEST', 'field', 'ID']

header_string = '#ID b_B435 b_errB435 b_V606 b_errV606 b_I814 b_errI814 b_Y105 b_errY105 b_J125 b_errJ125 b_JH140 b_errJH140 b_H160 b_errH160 b_Ks b_errKs b_CH1 b_errCH1 b_CH2 b_errCH2 ZBEST field ID_original\n'


#f= open("astrodeep_A2744_c_subset_RF1_001.txt","w+")
#f.write(header_string)
#
#for i in range(len(D_RF1)):
#    row = str(i+1) # ID from 1 to len(D_RF1)
#    
#    for j in range(len(header)):
#        row = row + ' ' + str(D_RF1[header[j]][i])
#        
#    row = row + '\n'    
#    f.write(row)
#
#f.close()

#np.save("astrodeep_A2744_c_subset_RF1_001", D_RF1)






























mpl.rcdefaults()