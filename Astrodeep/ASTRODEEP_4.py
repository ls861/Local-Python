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



###  ###  ###

for i in range(len(catalogs_z)):
    
    # make column names consistent
    catalogs_z[i].dtype.names = col_z
    catalogs_a[i].dtype.names = col_a
    catalogs_b[i].dtype.names = col_b    

    # remove all IDs < 100000
    catalogs_z[i] = catalogs_z[i][catalogs_z[i]['ID'] < 100000]
    catalogs_a[i] = catalogs_a[i][catalogs_a[i]['ID'] < 100000]
    catalogs_b[i] = catalogs_b[i][catalogs_b[i]['ID'] < 100000]

    # remove all lines corresponding to -ve redshift
    indices = np.where(catalogs_z[i]['ZBEST'] > 0)
    
    catalogs_z[i] = catalogs_z[i][indices]
    catalogs_a[i] = catalogs_a[i][indices]
    catalogs_b[i] = catalogs_b[i][indices]


total_z         = np.hstack(catalogs_z)
clusters_z      = np.hstack(catalogs_z[0::2])
parallels_z     = np.hstack(catalogs_z[1::2])

total_a         = np.hstack(catalogs_a)
clusters_a      = np.hstack(catalogs_a[0::2])
parallels_a     = np.hstack(catalogs_a[1::2])

total_b         = np.hstack(catalogs_b)
clusters_b      = np.hstack(catalogs_b[0::2])
parallels_b     = np.hstack(catalogs_b[1::2])

_z = [total_z, clusters_z, parallels_z]
_a = [total_a, clusters_a, parallels_a]
_b = [total_b, clusters_b, parallels_b]
title = ['Total', 'Clusters', 'Parallels']

print(np.array_equal(total_z['ID'], total_a['ID']))
print(np.array_equal(total_z['ID'], total_b['ID']))


###  ###  ###



### PLOTTING number per redshift ###

for i in range(len(_z)):
    plt.hist(_z[i]['ZBEST'], bins=np.linspace(0, 10.00000001, 11), histtype='step', label=title[i]) # this puts 10.0 max value in bin 9-10.

plt.legend()
plt.show()


### MAG PER REDSHIFT Plotting spitzer 3.6 (13) per redshift bin (combined catalog)

for i in range(len(_a)):
    for j in range(10):
        # , bins=np.linspace(15, 33, 40)  
        plt.hist(_a[i]['CH1'][ (_z[i]['ZBEST'] > j) & (_z[i]['ZBEST'] <= j+1)  & (_a[i]['CH1'] > 0) ], bins=50, label='%i $<$ z $\leq$ %i' % (j, j+1)) # this strips all -ve mags, MUST amend
    plt.title(title[i] + ' ' + 'CH1')
    plt.xlim(15, 31)
    plt.legend()
    plt.show()
    
### MAG PER REDSHIFT Plotting H160 per redshift bin (combined catalog)

for i in range(len(_a)):
    for j in range(10):
        # , bins=np.linspace(15, 33, 40)  
        plt.hist(_a[i]['H160'][ (_z[i]['ZBEST'] > j) & (_z[i]['ZBEST'] <= j+1)  & (_a[i]['H160'] > 10)  & (_a[i]['H160'] < 50) ], bins=50, label='%i $<$ z $\leq$ %i' % (j, j+1)) # this strips all mags outside of 10 and 50, MUST amend
    plt.title(title[i] + ' ' + 'H160')
    plt.xlim(15, 31)
    plt.legend()
    plt.show()


### Finding filter which includes rest frame 1500 UV, yet excludes 1215.
    # will determine 'in' filter by fwhm opposed to min and max. Can recode if necessary.
    
### FILTERS ###

filter_label = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])

filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])

filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])

filter_hwhm = filter_fwhm / 2


### ### ###

#total_z = total_z[total_z['ZBEST'] > 1]

uv = (1+total_z['ZBEST']) * 1500    # UV rest frame
la = (1+total_z['ZBEST']) * 1216    # lyman alpha rest frame

filter_firstguess = np.zeros(len(uv), dtype=int)  # if UV in only 1 filter, this is that filter, if UV in > 1 filter, this is the nearest to central filter, and if UV in 0 filters, this is next reddest filter. 


for i in range(len(uv)):
    
    in_filter = np.zeros(len(filter_fwhm_centre), dtype=float) 
    
    
    # determine which filters UV is IN
    for j in range(len(filter_fwhm_centre)):
        
        
        # UV in filter?
        if np.abs(uv[i] - filter_fwhm_centre[j]) < filter_hwhm[j]: 
            in_filter[j] = 1
        else:
            in_filter[j] = 0
            
    # UV in only 1 filter, select that filter 
    if sum(in_filter) == 1:
        
        filter_firstguess[i] = np.nonzero(in_filter)[0]
        
        
    # UV in > 1 filter, select nearest filter
    elif sum(in_filter) > 1:
        
        # select filter based on UV inside and nearest to centre of filter
        for j in range(len(filter_fwhm_centre)):
            
            in_filter[j] *= np.abs(uv[i] - filter_fwhm_centre[j])


        # this is finding the nearest central wavelength to UV[i]
        for j in np.nonzero(in_filter)[0]:
            
            if in_filter[j] <= in_filter[np.nonzero(in_filter)[0][0]]:
                
                filter_firstguess[i] = j
            

    # UV in NO filters, select nearest filter    
    elif sum(in_filter) == 0:

        # select filter nearest to centre of filter
        for j in range(len(filter_fwhm_centre)):
            
            # distance to each filter
            in_filter[j] = np.abs(uv[i] - filter_fwhm_centre[j])


        # this is finding the nearest central wavelength to UV[i]
        for j in range(len(in_filter)):
            
            if in_filter[j] <= in_filter[0]:
                
                filter_firstguess[i] = j
        
    #### comparing with LA 1216, if LA in firstguess, then try redder filter and retest. 
    
    while filter_firstguess[i] <= len(filter_fwhm_centre):
        
        #print('while', i)
        #print(filter_firstguess[i])
        
        
        if filter_firstguess[i] == len(filter_fwhm_centre):
            print("there aren't enough red filters you muppet")
            break
            
        elif filter_fwhm_centre[filter_firstguess[i]] - filter_hwhm[filter_firstguess[i]]  > la[i]:
            # first guess is suitable
            break
        
        else:
            # must try redder filter
            filter_firstguess[i] += 1
        


# plot quantity of filters
plt.hist(filter_firstguess, bins=[0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
plt.show()


# plot distribution of differences between 1500 and 'deredshifted' central filter wavelengths
plt.hist(filter_fwhm_centre[filter_firstguess] / (1+total_z['ZBEST']), bins=50) # amazed that this worked as expected, gives array of chosen central wavelengths converted to rest frame

plt.show()



### Calculating Absolute mags from chosen filters ###
print(filter_firstguess)


filter_chosen = filter_label[filter_firstguess]


print(filter_chosen)


#plt.hist(total_b[filter_chosen])
#plt.show()


print(type(total_b))
print(len(filter_chosen))
print(len(total_b))
# again, not obvious how to avoid a loop

uv_flux = np.zeros(len(filter_chosen))

for i in range(len(filter_chosen)):
    uv_flux[i] = total_b[filter_chosen[i]][i]
    


logbins = np.geomspace(0.0000001, 25000, 100)
#logbins = np.concatenate((-logbins[::-1], np.array([0]), logbins))

# reason behind the decision to cut anything under 0.0000001
print(len(uv_flux[uv_flux > 0.0000001]))
print(len(uv_flux[(uv_flux < 0.0000001) & (uv_flux > 0)]))
print(len(uv_flux[uv_flux == 0]))
print(len(uv_flux[uv_flux < 0]))


plt.xscale('log')
plt.hist(uv_flux[uv_flux > 0.0000001], bins=logbins)
plt.show()

m_app = -2.5 * np.log10(  (uv_flux[uv_flux > 0.0000001] * (10**-6))  /  3631) # calculated apparent mags from flux

# same as uv_flux but for uv_mag to check my calculations were correct
uv_mag = np.zeros(len(filter_chosen))

for i in range(len(filter_chosen)):
    uv_mag[i] = total_a[filter_chosen[i]][i]

uv_mag = uv_mag[uv_flux > 0.0000001] # from catalog

print(m_app)
print(uv_mag)
# seems legit, they match where uv_mag is +ve, still not sure what the -ve is about



### Converting fluxes to absolute mags using redshift ###

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z = total_z['ZBEST'][uv_flux > 0.0000001]
D_l = cosmo.luminosity_distance(z) 

m_abs = m_app - ( 5 * np.log10(D_l.value * (10**5)))

plt.hist(m_abs, bins=30)
























