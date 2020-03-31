#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 13:00:38 2020

@author: lester
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

param1 = 'DPL'
revisions = ['004', '005', '006', '007', '008', '009', '010', '011']
    
for revision1 in revisions:
    
    title1 = param1 + ' ' + revision1
    
    # =============================================================================
    # calculate input gradient for rising or falling + get INPUT params
    # =============================================================================
    
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_MS_parameters_004.fits'
    data_fits = fits.open(fileName)
    #print(data_fits[1].header)
    
    mass = data_fits[1].data['mass']
    dpl_alpha = data_fits[1].data['dpl_alpha']
    dpl_beta = data_fits[1].data['dpl_beta']
    tauV_eff = data_fits[1].data['tauV_eff']
    metallicity = data_fits[1].data['metallicity']
    nebular_logU = data_fits[1].data['nebular_logU']
    tau = 10**(data_fits[1].data['tau'])
    nebular_xi = data_fits[1].data['nebular_xi']
    
    A = data_fits[1].data['A']
    sfr = data_fits[1].data['sfr']
    
    data_fits.close()
    
    # =============================================================================
    # OUTPUT - get BEAGLE parameters
    # =============================================================================
    
    fileName = '/Users/lester/Documents/PhD/param_{}/fit_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(param1, revision1)
    data_fits = fits.open(fileName)
    
    id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1
    mtot_b1 = data_fits['POSTERIOR PDF'].data['mass_mean']
    mtot_68_b1 = data_fits['POSTERIOR PDF'].data['mass_68.00']
    sfr_b1 = data_fits['STAR FORMATION'].data['SFR_mean']
    sfr_68_b1 = data_fits['STAR FORMATION'].data['SFR_68.00']
    
    mass_b1 = data_fits['POSTERIOR PDF'].data['mass_mean']
    dpl_alpha_b1 = data_fits['POSTERIOR PDF'].data['dpl_alpha_mean']
    dpl_beta_b1 = data_fits['POSTERIOR PDF'].data['dpl_beta_mean']
    tauV_eff_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_mean']
    metallicity_b1 = data_fits['POSTERIOR PDF'].data['metallicity_mean']
    nebular_logU_b1 = data_fits['POSTERIOR PDF'].data['nebular_logu_mean']
    tau_b1 = data_fits['POSTERIOR PDF'].data['tau_mean']
    nebular_xi_b1 = data_fits['POSTERIOR PDF'].data['nebular_xi_mean']
    
    
    data_fits.close()
    
    # =============================================================================
    # comparing individual parameters
    # =============================================================================
    
    params_names = ['mass', 'alpha', 'beta', 'tauVeff', 'metallicity', 'nebularlogU', 'tau', 'nebularxi']
    params = [mass, np.log10(dpl_alpha), np.log10(dpl_beta), tauV_eff, metallicity, nebular_logU, np.log10(tau), nebular_xi]
    params_b1 = [mass_b1, np.log10(dpl_alpha_b1), np.log10(dpl_beta_b1), tauV_eff_b1, metallicity_b1, nebular_logU_b1, tau_b1, nebular_xi_b1]
    
    fig, axs = plt.subplots(2, len(params)/2, figsize=(15, 8))
    fig.suptitle(title1)
    for j in [0, 1]:
        for i in range(len(params)/2):
            axs[j,i].set_title(params_names[i+4*j])
            axs[j,i].scatter(params[i+4*j][id_b1], params_b1[i+4*j])
            min_ax = min(min(params[i+4*j][id_b1]), min(params_b1[i+4*j]))
#            min_ax = 9.2
            max_ax = max(max(params[i+4*j][id_b1]), max(params_b1[i+4*j]))
            axs[j,i].set_xlim(min_ax, max_ax)
            axs[j,i].set_ylim(min_ax, max_ax)
    plt.show()
    



















