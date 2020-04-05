#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:57:57 2019

@author: lester
"""

#import sys
import numpy as np
import matplotlib.pyplot as plt

#np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(threshold=25, edgeitems=10)


### Z ###

catalogs_z = np.empty(8)

col_1 = ('ID', 'ZBEST', 'ZBEST_SIQR', 'MAGNIF', 'ZSPECFLAG', 'Chi2', 'MSTAR', 'MSTAR_MIN', 'MSTAR_MAX', 'SFR', 'SFR_MIN', 'SFR_MAX', 'Chi2_NEB', 'MSTAR_NEB', 'MSTAR_MIN_NEB', 'MSTAR_MAX_NEB', 'SFR_NEB', 'SFR_MIN_NEB', 'SFR_MAX_NEB', 'RELFLAG')

col_2 = ('ID', 'ZBEST', 'ZBEST_SIQR', 'MAGNIF', 'ZSPECFLAG', 'chi2', 'MSTAR', 'MASTAR_MIN', 'MSTAR_MAX', 'SFR', 'SFR_MIN', 'SFR_MAX', 'chi2_NEB', 'MASTAR_NEB', 'MASS_MIN_NEB', 'MASS_MAX_NEB', 'SFR_NEB', 'SFR_MIN_NEB', 'SFR_MAX_NEB', 'RELFLAG')

col_z = ('ID', 'ZBEST', 'ZBEST_SIQR', 'MAGNIF', 'ZSPECFLAG', 'chi2', 'MSTAR', 'MASTAR_MIN', 'MSTAR_MAX', 'SFR', 'SFR_MIN', 'SFR_MAX', 'chi2_NEB', 'MASTAR_NEB', 'MASS_MIN_NEB', 'MASS_MAX_NEB', 'SFR_NEB', 'SFR_MIN_NEB', 'SFR_MAX_NEB', 'RELFLAG')

A2744_c_z = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'A2744cl_ZPHOT.cat', names=True, dtype=float, usecols=col_1)
A2744_p_z = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'A2744PAR_ZPHOT.cat', names=True, dtype=float, usecols=col_1)
M0416_c_z = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'M0416cl_ZPHOT.cat', names=True, dtype=float, usecols=col_1)
M0416_p_z = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'M0416PAR_ZPHOT.cat', names=True, dtype=float, usecols=col_1)
M0717_c_z = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS0717cl_ZPHOT_complete.cat', names=True, dtype=float, usecols=col_2)
M0717_p_z = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS0717par_ZPHOT_complete.cat', names=True, dtype=float, usecols=col_2)
M1149_c_z = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS1149cl_ZPHOT_complete.cat', names=True, dtype=float, usecols=col_2)
M1149_p_z = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS1149par_ZPHOT_complete.cat', names=True, dtype=float, usecols=col_2)
    
catalogs_z        = [A2744_c_z,
                     A2744_p_z,
                     M0416_c_z,
                     M0416_p_z,
                     M0717_c_z,
                     M0717_p_z,
                     M1149_c_z,
                     M1149_p_z]

### A ###

catalogs_a = np.empty(8)

col_3 = ('ID', 'RA', 'DEC', 'MAG_B435', 'MAG_V606', 'MAG_I814', 'MAG_Y105', 'MAG_J125', 'MAG_JH140', 'MAG_H160', 'MAG_Ks', 'MAG_IRAC1', 'MAG_IRAC2', 'MAGERR_B435', 'MAGERR_V606', 'MAGERR_I814', 'MAGERR_Y105', 'MAGERR_J125', 'MAGERR_JH140', 'MAGERR_H160', 'MAGERR_Ks', 'MAGERR_IRAC1', 'MAGERR_IRAC2')

col_4 = ('ID', 'RA', 'DEC', 'B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2', 'errB435', 'errV606', 'errI814', 'errY105', 'errJ125', 'errJH140', 'errH160', 'errKs', 'errCH1', 'errCH2')

col_a = ('a_ID', 'a_RA', 'a_DEC', 'a_B435', 'a_V606', 'a_I814', 'a_Y105', 'a_J125', 'a_JH140', 'a_H160', 'a_Ks', 'a_CH1', 'a_CH2', 'a_errB435', 'a_errV606', 'a_errI814', 'a_errY105', 'a_errJ125', 'a_errJH140', 'a_errH160', 'a_errKs', 'a_errCH1', 'a_errCH2')

A2744_c_a = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'A2744cl_A.cat', names=True, dtype=float, usecols=col_3)
A2744_p_a = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'A2744PAR_A.cat', names=True, dtype=float, usecols=col_3)
M0416_c_a = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'M0416cl_A.cat', names=True, dtype=float, usecols=col_3)
M0416_p_a = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'M0416PAR_A.cat', names=True, dtype=float, usecols=col_3)
M0717_c_a = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS0717cl_A.cat', names=True, dtype=float, usecols=col_4)
M0717_p_a = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS0717par_A.cat', names=True, dtype=float, usecols=col_4)
M1149_c_a = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS1149cl_A.cat', names=True, dtype=float, usecols=col_4)
M1149_p_a = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS1149par_A.cat', names=True, dtype=float, usecols=col_4)
    
catalogs_a        = [A2744_c_a,
                     A2744_p_a,
                     M0416_c_a,
                     M0416_p_a,
                     M0717_c_a,
                     M0717_p_a,
                     M1149_c_a,
                     M1149_p_a]

### B ###

catalogs_b = np.empty(8)

col_5 = ('ID', 'FLUX_B435', 'FLUX_V606', 'FLUX_I814', 'FLUX_Y105', 'FLUX_J125', 'FLUX_JH140', 'FLUX_H160', 'FLUX_Ks', 'FLUX_IRAC1', 'FLUX_IRAC2', 'FLUXERR_B435', 'FLUXERR_V606', 'FLUXERR_I814', 'FLUXERR_Y105', 'FLUXERR_J125', 'FLUXERR_JH140', 'FLUXERR_H160', 'FLUXERR_Ks', 'FLUXERR_IRAC1', 'FLUXERR_IRAC2')

col_6 = ('ID', 'B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2', 'errB435', 'errV606', 'errI814', 'errY105', 'errJ125', 'errJH140', 'errH160', 'errKs', 'errCH1', 'errCH2')

col_b = ('b_ID', 'b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2', 'b_errB435', 'b_errV606', 'b_errI814', 'b_errY105', 'b_errJ125', 'b_errJH140', 'b_errH160', 'b_errKs', 'b_errCH1', 'b_errCH2')

A2744_c_b = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'A2744cl_B.cat', names=True, dtype=float, usecols=col_5)
A2744_p_b = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'A2744PAR_B.cat', names=True, dtype=float, usecols=col_5)
M0416_c_b = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'M0416cl_B.cat', names=True, dtype=float, usecols=col_5)
M0416_p_b = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'M0416PAR_B.cat', names=True, dtype=float, usecols=col_5)
M0717_c_b = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS0717cl_B.cat', names=True, dtype=float, usecols=col_6)
M0717_p_b = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS0717par_B.cat', names=True, dtype=float, usecols=col_6)
M1149_c_b = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS1149cl_B.cat', names=True, dtype=float, usecols=col_6)
M1149_p_b = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS1149par_B.cat', names=True, dtype=float, usecols=col_6)
    
catalogs_b        = [A2744_c_b,
                     A2744_p_b,
                     M0416_c_b,
                     M0416_p_b,
                     M0717_c_b,
                     M0717_p_b,
                     M1149_c_b,
                     M1149_p_b]

###  ###  ###

for i in range(len(catalogs_z)):  
    
    # remove all IDs < 100000 and give arrays equal dimensions
    catalogs_z[i] = catalogs_z[i][catalogs_z[i]['ID'] < 100000]
    catalogs_a[i] = catalogs_a[i][catalogs_a[i]['ID'] < 100000]
    catalogs_b[i] = catalogs_b[i][catalogs_b[i]['ID'] < 100000]


### ### ###
# need to append C or P onto the z catalog of each field (x8)
# and also merge catalogs z, a and b...

for i in range(len(catalogs_z)):
    
    # converts from structured array to basic headerless numpy array
    catalogs_z[i] = catalogs_z[i].view(np.float64).reshape(catalogs_z[i].shape + (-1,))
    catalogs_a[i] = catalogs_a[i].view(np.float64).reshape(catalogs_a[i].shape + (-1,))
    catalogs_b[i] = catalogs_b[i].view(np.float64).reshape(catalogs_b[i].shape + (-1,))

    # creates an array of identical numbers (0 - 7) depending on the field
    field_flag = np.full((len(catalogs_z[i]), 1), i)
    
    # smashes the above together
    catalogs_z[i] = np.append(field_flag, catalogs_z[i], axis=1)
    catalogs_z[i] = np.append(catalogs_z[i], catalogs_a[i], axis=1)
    catalogs_z[i] = np.append(catalogs_z[i], catalogs_b[i], axis=1)
    
    # converts numpy array into suitable input to recreate structured array
    catalogs_z[i] = np.core.records.fromarrays(catalogs_z[i].transpose())
    
    # creates a super useful structured array
    catalogs_z[i] = np.array(catalogs_z[i], dtype=[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')])
    
    
cat       = np.hstack(catalogs_z)
cat       = cat[cat['ZBEST'] > 0] # removes values with -1

#np.save("cat", cat)

















    
    
    
    