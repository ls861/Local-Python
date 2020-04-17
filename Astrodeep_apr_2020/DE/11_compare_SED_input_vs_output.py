#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 16:22:23 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

param = 'DE'
revision = '100'

samples = 10

IDs=[22, 61, 63]

fsize=10
c = 299792458 # m s^-2


for ID in IDs:
    
    title = '{}-{} ID{} {} samples'.format(param, revision.replace('_', '-'), str(ID), str(samples))
    
    # =============================================================================
    # ASTRODEEP filter info
    # =============================================================================
    
    # column name from apparent mag table in BEAGLE mock (same as left column in config + _APP)
    filters = ['HST_ACS_WFC_F435W_APP', 'HST_ACS_WFC_F606W_APP', 'HST_ACS_WFC_F814W_APP', 'HST_WFC3_IR_F105W_APP', 'HST_WFC3_IR_F125W_APP', 'HST_WFC3_IR_F140W_APP', 'HST_WFC3_IR_F160W_APP', 'Paranal_HAWKI_Ks_APP', 'Spitzer_IRAC_I1_APP', 'Spitzer_IRAC_I2_APP']
    
    # column name from ASTRODEEP config file
    filter_label = ['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2'] 
    
    # for taking errors of the perturbed input fluxes
    filter_err = ['b_errB435', 'b_errV606', 'b_errI814', 'b_errY105', 'b_errJ125', 'b_errJH140', 'b_errH160', 'b_errKs', 'b_errCH1', 'b_errCH2'] 
    
    filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])
    filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])
    
    # =============================================================================
    # from INPUT fits file get SEDs and apparent ABMAG
    # =============================================================================
    
    fileName = '/Users/lester/Documents/PhD/param_100/mock_100_{}/mock_catalogue_100_{}.fits'.format(param, param)
    data_fits = fits.open(fileName)
    #print(data_fits.info())
    
    z_mock = data_fits['GALAXY PROPERTIES'].data['redshift'][ID-1]
    wl_spec_mock = data_fits['FULL SED WL'].data['wl'][0]*(1+z_mock)
    f_spec_mock = data_fits['FULL SED'].data[ID-1]/(1+z_mock)
    lfl_spec_mock = f_spec_mock*wl_spec_mock
    
    # PHOTOMETRY 
    appmag_phot_mock = np.zeros(len(filters))
    for i in range(len(filters)):
        appmag_phot_mock[i] = data_fits['APPARENT MAGNITUDES'].data[filters[i]][ID-1]
        
    data_fits.close()
    
    f_phot_mock_mJy = (10**( (23.9 - appmag_phot_mock) / 2.5 )) # not used [mJy]
    lfl_test = (3631 * 10**( -appmag_phot_mock / 2.5 ) )   /    ( 3.34e4 * filter_fwhm_centre) # not used [erg cm-2 s-1]
    lfl_phot_mock = (c / filter_fwhm_centre) * (10 ** (-(appmag_phot_mock+23.6)/2.5)) # [erg cm-2 s-1]
    
    
    # =============================================================================
    # get INPUT perturbed fluxes and errors
    # =============================================================================
    
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/{}/100_{}_mock_fluxes_CANDELS.fits'.format(param, param)
    
    data_fits = fits.open(fileName)
    # print(data_fits.info())
    # print(data_fits[1].header)
    
    # PHOTOMETRY 
    ptbf_phot_mock = np.zeros(len(filter_label))
    for i in range(len(filters)):
        ptbf_phot_mock[i] = data_fits[1].data[filter_label[i]][ID-1] # uJy
        
    ptbferr_phot_mock = np.zeros(len(filter_label))
    for i in range(len(filters)):
        ptbferr_phot_mock[i] = data_fits[1].data[filter_err[i]][ID-1] # uJy    
        
    # adding min rel error to errors: 
    min_rel_error = np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.05, 0.1, 0.1])
    ptbferr_phot_mock = ( (ptbferr_phot_mock**2) + ((min_rel_error*ptbf_phot_mock)**2) ) ** 0.5
    
    ptblfl_phot_mock = (ptbf_phot_mock/filter_fwhm_centre)*(1e-10 / 3.34) # erg cm-2 s-1
    ptblflerr_phot_mock = (ptbferr_phot_mock/filter_fwhm_centre)*(1e-10 / 3.34) # erg cm-2 s-1
    
    data_fits.close()
    
    # =============================================================================
    # Get OUTPUT SEDs
    # =============================================================================
    
    data_fits = fits.open('/Users/lester/Documents/PhD/param_100/fit_{}_{}/{}_BEAGLE.fits.gz'.format(revision, param, ID))
    # print(data_fits.info())
    # print(data_fits['POSTERIOR PDF'].header)
    
    chi2_fit_total = data_fits['POSTERIOR PDF'].data['chi_square']
    
    #needs float64 to provide precision needed for the random.choice weights
    temp_probs = np.float64(data_fits['POSTERIOR PDF'].data['probability'])
    temp_probs = temp_probs/np.sum(temp_probs)
      
    z_fit_arr = []
    wl_spec_fit_arr = []
    lfl_spec_fit_arr = []
    
    lfl_phot_fit_arr = []
    
    chi2_fit_arr = []
    
    for i in range(samples):
        
        #here's the key line - take weighted samples from the multinest output!
        idx = np.random.choice(len(temp_probs), size=1, p=temp_probs)
    
        z_fit_arr.append(data_fits['GALAXY PROPERTIES'].data['redshift'][idx])
        wl_spec_fit_arr.append(data_fits['FULL SED WL'].data['wl'][0]*(1+z_fit_arr[i]))
        lfl_spec_fit_arr.append(data_fits['FULL SED'].data[idx][0]/(1+z_fit_arr[i]) * wl_spec_fit_arr[i])
        
    
        # PHOTOMETRY 
        appmag_phot_fit = np.zeros(len(filters))
        for i in range(len(filters)):
            appmag_phot_fit[i] = data_fits['APPARENT MAGNITUDES'].data[filters[i]][idx]
            
            
        data_fits.close()
        
        lfl_phot_fit = (c / filter_fwhm_centre) * (10 ** (-(appmag_phot_fit+23.6)/2.5)) # [erg cm-2 s-1]
        lfl_phot_fit_arr.append(lfl_phot_fit)
        
        # CHI SQUARED
        chi2_fit_arr.append(data_fits['POSTERIOR PDF'].data['chi_square'][idx][0])
    
    # =============================================================================
    # averaging the samples to assist in plotting
    # =============================================================================
    
    # spec
    lfl_spec_fit_min = lfl_spec_fit_arr[0]
    lfl_spec_fit_max = lfl_spec_fit_arr[0]
    
    # phot
    lfl_phot_fit_min = lfl_phot_fit_arr[0]
    lfl_phot_fit_max = lfl_phot_fit_arr[0]
    
    
    for i in range(samples-1):
        lfl_spec_fit_min = np.minimum(lfl_spec_fit_min, lfl_spec_fit_arr[i+1])
        lfl_spec_fit_max = np.maximum(lfl_spec_fit_max, lfl_spec_fit_arr[i+1])
    
        lfl_phot_fit_min = np.minimum(lfl_phot_fit_min, lfl_phot_fit_arr[i+1])
        lfl_phot_fit_max = np.maximum(lfl_phot_fit_max, lfl_phot_fit_arr[i+1])
    
    
    # =============================================================================
    # PLOT with fitted SEDs shaded
    # =============================================================================
    
    plt.figure(figsize=(2*fsize, 1*fsize))
    plt.title(title + ', chi2 = {0:.3g}'.format(np.average(chi2_fit_arr)), fontsize=14)
    plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
    plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)

    
    plt.plot(wl_spec_mock, lfl_spec_mock, linewidth=0.5, zorder=2, label='Real SED')
    
    plt.fill_between(wl_spec_mock, lfl_spec_fit_min, lfl_spec_fit_max, alpha=0.4, color='k', linewidth=0, zorder=1, label='Output SED samples') 
    
    plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='c', zorder=3, label='Real Phot')
    
    plt.errorbar(filter_fwhm_centre, lfl_phot_fit_min, yerr=[np.zeros(len(lfl_phot_fit_min)), lfl_phot_fit_max - lfl_phot_fit_min], linestyle="None", linewidth=10, color='r', zorder=4, label='Output Phot samples')
    
    plt.scatter(filter_fwhm_centre, ptblfl_phot_mock, color='k', marker='x', zorder=5, label='Perturbed Input Fluxes')
    plt.errorbar(filter_fwhm_centre, ptblfl_phot_mock, yerr= ptblflerr_phot_mock, linestyle="None", color='k', zorder=5)

    idx_lim1 = (np.abs(wl_spec_mock - 45000)).argmin()
    y_lim1 = lfl_spec_mock[idx_lim1]
    
    idx_lim2 = (np.abs(wl_spec_mock - 12500)).argmin()
    y_lim2 = lfl_spec_mock[idx_lim2]
    
    plt.xlim(0, 50000)
#    plt.xlim(10000, 15000)
#    plt.ylim(0.9*y_lim1, 1.5*y_lim2)
    plt.ylim(0.8*min(ptblfl_phot_mock), 1.5*max(ptblfl_phot_mock))
    
    plt.yscale('log')
    plt.legend()
    plt.show()
    

    # =============================================================================
    # chi_squared random comparisons
    # =============================================================================
    
#    chi2_mock = np.sum(((ptblfl_phot_mock - lfl_phot_mock)**2) / (ptblflerr_phot_mock**2) )
#    
#    #print(chi2_fit_total)
#    plt.figure(figsize=(0.2*fsize, 0.2*fsize))
#    plt.hist(chi2_fit_total, bins=20)
#    plt.show()
#    
#    #print(chi2_fit_arr)
#    plt.figure(figsize=(0.2*fsize, 0.2*fsize))
#    plt.hist(chi2_fit_arr, bins=20)
#    plt.show()
#    
#    # just choose a number from the samples
#    i=0
#    
#    test = np.sum(((lfl_phot_fit_arr[i] - ptblfl_phot_mock)**2) / (ptblflerr_phot_mock**2) )
#    print('lester')
#    print(test)
#    print(chi2_fit_arr[i])
    







# =============================================================================
# PLOT mock
# =============================================================================

#plt.figure(figsize=(2*fsize, fsize))
#plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
#plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)
#plt.xlim(0, 50000)
#plt.ylim(1e-15, 1e-14)
#plt.yscale('log')
#plt.plot(wl_spec_mock, f_spec_mock*wl_spec_mock, linewidth=0.5, zorder=0)
##plt.scatter(filter_fwhm_centre, lfl_phot_mock, marker='x', color='r', zorder=2)
#plt.scatter(filter_fwhm_centre, lfl_phot_fit, marker='x', color='k', zorder=2)
#plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='r', zorder=1)
##plt.legend()
#plt.show()



# =============================================================================
# PLOT a few SEDs
# =============================================================================

#plt.figure(figsize=(6*fsize, 3*fsize))
#plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
#plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)
#plt.xlim(0, 50000)
#plt.ylim(1e-15, 1e-14)
#plt.yscale('log')
#plt.plot(wl_spec_mock, f_spec_mock*wl_spec_mock, linewidth=0.5, zorder=0)
#for i in range(samples):
#    plt.plot(wl_spec_fit_arr[i], f_spec_fit_arr[i]*wl_spec_fit_arr[i], linewidth=0.5, zorder=1)
#    
##plt.scatter(filter_fwhm_centre, lfl_phot_mock, marker='x', color='r', zorder=3)
#plt.scatter(filter_fwhm_centre, lfl_phot_fit, marker='x', color='k', zorder=3)
#plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='r', zorder=2)
##plt.legend()
#plt.show()

# =============================================================================
# PLOT with fitted SEDs shaded
# =============================================================================

#plt.figure(figsize=(6*fsize, 2*fsize))
#plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
#plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)
#plt.xlim(0, 50000)
#plt.ylim(1e-15, 1e-14)
#
#plt.xlim(15000, 20000)
#plt.ylim(2.5e-15, 3.5e-15)
#
#plt.plot(wl_spec_mock, f_spec_mock*wl_spec_mock, linewidth=0.5, zorder=1)
#
#
##plt.fill_between(wl_spec_mock, f_spec_fit_min*wl_spec_mock, f_spec_fit_max*wl_spec_mock) 
#plt.plot(wl_spec_mock, f_spec_fit_min*wl_spec_mock, alpha=0.3, color='k', zorder=0) 
#plt.plot(wl_spec_mock, f_spec_fit_max*wl_spec_mock, alpha=0.3, color='k', zorder=0) 
#
##plt.scatter(filter_fwhm_centre, lfl_phot_mock, marker='x', color='r', zorder=3)
##plt.scatter(filter_fwhm_centre, lfl_phot_fit, marker='x', color='k', zorder=3)
##plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='r', zorder=2)
##plt.legend()
#
#for i in range(samples):
#    plt.plot(wl_spec_fit_arr[i], f_spec_fit_arr[i]*wl_spec_fit_arr[i], linewidth=0.5, zorder=1)
#    
#plt.yscale('log')
#plt.show()


# =============================================================================
# PLOT with fitted SEDs shaded
# =============================================================================

#plt.figure(figsize=(6*fsize, 2*fsize))
#plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
#plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)
#plt.xlim(0, 50000)
#plt.ylim(1e-15, 1e-14)
#
##plt.xlim(10000, 20000)
##plt.ylim(2.5e-15, 3.5e-15)
#
#plt.plot(wl_spec_mock, f_spec_mock*wl_spec_mock, linewidth=0.5, zorder=1)
#
#
##plt.fill_between(wl_spec_mock, f_spec_fit_min*wl_spec_mock, f_spec_fit_max*wl_spec_mock) 
#plt.plot(wl_spec_mock, f_spec_fit_min*wl_spec_mock, alpha=0.3, color='k', zorder=0) 
#plt.plot(wl_spec_mock, f_spec_fit_max*wl_spec_mock, alpha=0.3, color='k', zorder=0) 
#
##plt.scatter(filter_fwhm_centre, lfl_phot_mock, marker='x', color='r', zorder=3)
##plt.scatter(filter_fwhm_centre, lfl_phot_fit, marker='x', color='k', zorder=3)
#plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='c', zorder=2)
#plt.errorbar(filter_fwhm_centre, lfl_phot_fit_min, yerr=[np.zeros(len(lfl_phot_fit_min)), lfl_phot_fit_max - lfl_phot_fit_min], linestyle="None", linewidth=10, color='r', zorder=5)
#
#    
#
##plt.legend()
#
#plt.yscale('log')
#plt.show()









