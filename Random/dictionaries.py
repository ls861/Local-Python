#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:08:05 2019

@author: lester
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:57:57 2019

@author: lester
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

#np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(threshold=10)


                  
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

for i in range(len(catalogs_z)):
    catalogs_z[i].dtype.names = col_z
    catalogs_z[i] = catalogs_z[i][np.where(catalogs_z[i]['ID'] < 100000)]
    

### A ###

catalogs_a = np.empty(8)

col_3 = ('ID', 'RA', 'DEC', 'MAG_B435', 'MAG_V606', 'MAG_I814', 'MAG_Y105', 'MAG_J125', 'MAG_JH140', 'MAG_H160', 'MAG_Ks', 'MAG_IRAC1', 'MAG_IRAC2', 'MAGERR_B435', 'MAGERR_V606', 'MAGERR_I814', 'MAGERR_Y105', 'MAGERR_J125', 'MAGERR_JH140', 'MAGERR_H160', 'MAGERR_Ks', 'MAGERR_IRAC1', 'MAGERR_IRAC2')

col_4 = ('ID', 'RA', 'DEC', 'B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2', 'errB435', 'errV606', 'errI814', 'errY105', 'errJ125', 'errJH140', 'errH160', 'errKs', 'errCH1', 'errCH2')

col_a = ('ID', 'RA', 'DEC', 'B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2', 'errB435', 'errV606', 'errI814', 'errY105', 'errJ125', 'errJH140', 'errH160', 'errKs', 'errCH1', 'errCH2')

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

for i in range(len(catalogs_a)):
    catalogs_a[i].dtype.names = col_a
    catalogs_a[i] = catalogs_a[i][np.where(catalogs_a[i]['ID'] < 100000)]
    

### B ###

catalogs_b = np.empty(8)

col_5 = ('ID', 'FLUX_B435', 'FLUX_V606', 'FLUX_I814', 'FLUX_Y105', 'FLUX_J125', 'FLUX_JH140', 'FLUX_H160', 'FLUX_Ks', 'FLUX_IRAC1', 'FLUX_IRAC2', 'FLUXERR_B435', 'FLUXERR_V606', 'FLUXERR_I814', 'FLUXERR_Y105', 'FLUXERR_J125', 'FLUXERR_JH140', 'FLUXERR_H160', 'FLUXERR_Ks', 'FLUXERR_IRAC1', 'FLUXERR_IRAC2')

col_6 = ('ID', 'B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2', 'errB435', 'errV606', 'errI814', 'errY105', 'errJ125', 'errJH140', 'errH160', 'errKs', 'errCH1', 'errCH2')

col_b = ('ID', 'B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2', 'errB435', 'errV606', 'errI814', 'errY105', 'errJ125', 'errJH140', 'errH160', 'errKs', 'errCH1', 'errCH2')

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

for i in range(len(catalogs_b)):
    catalogs_b[i].dtype.names = col_b
    catalogs_b[i] = catalogs_b[i][np.where(catalogs_b[i]['ID'] < 100000)]

###  ###  ###


total_z         = np.hstack(catalogs_z)
clusters_z      = np.hstack((A2744_c_z, M0416_c_z, M0717_c_z, M1149_c_z))
parallels_z     = np.hstack((A2744_p_z, M0416_p_z, M0717_p_z, M1149_p_z))

total_a         = np.hstack(catalogs_a)
clusters_a      = np.hstack((A2744_c_a, M0416_c_a, M0717_c_a, M1149_c_a))
parallels_a     = np.hstack((A2744_p_a, M0416_p_a, M0717_p_a, M1149_p_a))

total_b         = np.hstack(catalogs_b)
clusters_b      = np.hstack((A2744_c_b, M0416_c_b, M0717_c_b, M1149_c_b))
parallels_b     = np.hstack((A2744_p_b, M0416_p_b, M0717_p_b, M1149_p_b))


print(np.array_equal(total_z['ID'], total_a['ID']))
print(np.array_equal(total_z['ID'], total_b['ID']))


###  ###  ###













