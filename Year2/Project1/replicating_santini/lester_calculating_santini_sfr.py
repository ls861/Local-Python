#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 22:03:39 2020

@author: lester
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

### TO CONSIDER
# CHECK WAVELENGTH LIMITS

### DONE
# COMPARE WITH FLUX CALCULATIONS (as remember mags were strange...), fluxes are in microJy - using fluxes, checked that usually gives same result as ABmag
# REMEMBER TO TRY MAGNIFICATION BEFORE AND AFTER - CHECKED ITS FINE
# CHECK LOG BASE 10 for BETA CALC - SURELY based on definition of mag
# CHECK I'VE APPLIED EXTINCTION CORRECTLY (just subtracted from AB mag at 1600) - that's also what emma did



AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/astrodeep/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
#print(AD.dtype.names)

# SUBSET FOR TESTING
#AD = AD[:20]


#find filters spanning 1280-2600A
filter_label_flux = np.array(['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2'])
filter_label_ABmag = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])
filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])
filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])   
filter_low = filter_fwhm_centre - filter_fwhm/2
filter_high = filter_fwhm_centre + filter_fwhm/2


redshift_med = AD['ZBEST']

sfrArr = []
betaArr = []

for i in range(len(redshift_med)):
    
#    redshift_med = [1.3]
    lambda_1280 = 1280*(redshift_med[i]+1.0)
    lambda_2600 = 2600*(redshift_med[i]+1.0) # 3100 ensures 2 filters available at z==1.3 using edges, 2600 is 2 filters at about z=1.28 when using filter centers
#    print(lambda_1280, lambda_2600)
    
    appMagArr = []
#    appMagArr2 = [] # to test is ABmag in catalogue was same as my calculation from flux
    appMagErr = []
    lambdaArr = []
    
    
    for f in range(len(filter_label_flux)):
        
#        if lambda_1280  < filter_low[f] and lambda_2600 > filter_high[f]:
        if lambda_1280  < filter_fwhm_centre[f] and lambda_2600 > filter_fwhm_centre[f]:
#            print(filter_label_flux[f])
            temp_fnu = AD[filter_label_flux[f]][i] # uJy, F_nu
            
            if temp_fnu > 0:
                
                temp_fnu_err = AD[filter_label_flux[f].replace('_', '_err')][i]
                appMagArr.append(-2.5 * np.log10(temp_fnu*1e-6) + 8.90)
#                appMagArr2.append(AD[filter_label_ABmag[f]][i])
                lambdaArr.append(filter_fwhm_centre[f])

    if len(appMagArr) >= 2:
        
        # CASTELLANO 2012 - CALCULATE BETA
        fit = np.polyfit(np.log10(lambdaArr), appMagArr, 1)    

#        for m in range(len(filter_label_ABmag)):
#            plt.scatter(np.log10(filter_fwhm_centre[m]), AD[filter_label_ABmag[m]][i], marker='x', color='k')       
#        plt.scatter(np.log10(lambdaArr), appMagArr, marker='x', color='r')
#        x = np.array((np.log10(lambdaArr)[0],np.log10(lambdaArr)[-1]))
#        plt.plot(x,fit[0]*x+fit[1])
#        plt.show()
        
        beta = (fit[0]/-2.5)-2.0 # Mi = −2.5(β + 2.0) log(λi) + c
#        print(beta)
        # MEURER 1999 - CORRECT 1600A for dust extinction
        A_1600 = 4.43 + (1.99*beta)
        appMag_1600 = (-2.5*(beta+2.0)*np.log10(1600*(1+redshift_med[i]))) + fit[1] - A_1600

        F_v = 10**((appMag_1600 - 8.90)/-2.5) # Jy
#        print(appMag_1600)
#        print(F_v)
        F_v = F_v * 1e-26 # W Hz-1 m-2
        
        D_l = cosmo.luminosity_distance(redshift_med[i]).value # Mpc
#        print(D_l)
        D_l = D_l * 1e6 * 3.086e16 # m
        
        L_v = F_v * 4 * np.pi * (D_l**2) / (1+redshift_med[i]) # W Hz-1
        
        vL_v = (299792458 / (1600*1e-10)) * L_v # W
        vL_v = vL_v * 1e7 # erg s-1

        # KENNICUTT 2012
        sfr = np.log10(vL_v) - 43.35 # log (solar masses per year)
#        print(L_v)
#        print(sfr)

        sfrArr.append(sfr)
        betaArr.append(beta)
    else:
        sfrArr.append(-99.0)
        betaArr.append(-99.0)
  

sbf = './calculate_santini_sfr_results/'
#np.save(sbf+'sfr_santini_lester_2600_median.npy', sfrArr)
#np.save(sbf+'sfr_santini_beta_lester_2600_median.npy', betaArr)
#np.save(sbf+'sfr_santini_temp.npy', sfrArr)
#np.save(sbf+'sfr_santini_beta_temp.npy', betaArr)

#%%


#sfr_SAN = np.load(sbf+'sfr_santini.npy')
#sfr_SAN_lester = np.load(sbf+'sfr_santini_lester.npy')
#sfr_SAN_lester_2600 = np.load(sbf+'sfr_santini_beta_lester_2600_median.npy')
#
#
#plt.scatter(sfr_SAN, sfr_SAN_lester)
#plt.xlim(-4, 4)
#plt.ylim(-4, 4)
#plt.show()

#print((sfr_SAN-sfr_SAN_lester))




'''
appMag = 10

# trying to get F_v

LESTER = 10**((appMag - 8.90)/-2.5) * 1e-26 # W Hz-1 m-2

EMMA = 10**((appMag+48.6)/(-2.5)) # erg s−1 cm−2 Hz−1

#print(LESTER, EMMA, LESTER/EMMA, EMMA/LESTER)
'''

#plt.hist(betaArr, bins=100, range=(-20, 20))
#plt.show()











