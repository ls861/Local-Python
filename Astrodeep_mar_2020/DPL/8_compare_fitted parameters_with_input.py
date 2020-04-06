#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 13:00:38 2020

@author: lester
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import cosmolopy.distance as cd
import cosmolopy.constants as cc
from scipy.integrate import quad

param1 = 'DPL'
revisions = ['004', '005', '006', '007', '008', '009', '010', '010_002', '011', '012_001']
revisions = ['012_005']
    
for revision1 in revisions:
    
    title1 = param1 + ' ' + revision1.replace('_', ' ')
    
    # =============================================================================
    # get INPUT params (len 100)
    # =============================================================================
    
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_MS_parameters_005.fits'
    data_fits = fits.open(fileName)
#    print(data_fits[1].header)

    id_i = data_fits[1].data['id']
    
    mass = data_fits[1].data['mass']
    dpl_alpha = data_fits[1].data['dpl_alpha']
    dpl_beta = data_fits[1].data['dpl_beta']
    tauV_eff = data_fits[1].data['tauV_eff']
    metallicity = data_fits[1].data['metallicity']
    nebular_logU = data_fits[1].data['nebular_logU']
    tau = 10**(data_fits[1].data['tau'])
    nebular_xi = data_fits[1].data['nebular_xi']
    
    A = data_fits[1].data['A']
    sfr = np.log10(data_fits[1].data['sfr'])
    
    data_fits.close()
    
    # =============================================================================
    # calculate input gradient for rising or falling (len 100)
    # =============================================================================
    
    cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
    cosmo = cd.set_omega_k_0(cosmo)
    ageUniv2 = cd.age(2.0, **cosmo)/cc.yr_s
    ageUniv999 = cd.age(999.0, **cosmo)/cc.yr_s
    msa = ageUniv2 - ageUniv999

    xlin = np.linspace(1, 1e10, 100000)
    grad_in = np.empty(len(A))
    
    # nice trick to find index in xlin which has value closest to msa
    idx = (np.abs(xlin - msa)).argmin()
    
    for i in range(len(A)):
        sfr_calc = A[i] / (((xlin/tau[i])**dpl_alpha[i])+((xlin/tau[i])**-dpl_beta[i]))
        grad_in[i] = np.gradient(sfr_calc, xlin)[idx]
    
    idx_rising_in = grad_in >= 0
    idx_falling_in = grad_in < 0
    
    # =============================================================================
    # OUTPUT - get BEAGLE parameters (<100)
    # =============================================================================
    
    fileName = '/Users/lester/Documents/PhD/param_{}/fit_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(param1, revision1)
    data_fits = fits.open(fileName)
    
    id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1

    sfr_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_mean'])
    sfr_68_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_68.00'])
    
    ssfr_b1 = data_fits['STAR FORMATION'].data['sSFR_mean']
    ssfr_68_b1 = data_fits['STAR FORMATION'].data['sSFR_68.00']
    
    mass_b1 = data_fits['POSTERIOR PDF'].data['mass_mean']
    dpl_alpha_b1 = data_fits['POSTERIOR PDF'].data['dpl_alpha_mean']
    dpl_beta_b1 = data_fits['POSTERIOR PDF'].data['dpl_beta_mean']
    tauV_eff_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_mean']
    metallicity_b1 = data_fits['POSTERIOR PDF'].data['metallicity_mean']
    nebular_logU_b1 = data_fits['POSTERIOR PDF'].data['nebular_logu_mean']
    tau_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_mean']
    nebular_xi_b1 = data_fits['POSTERIOR PDF'].data['nebular_xi_mean']
    
    mass_68_b1 = data_fits['POSTERIOR PDF'].data['mass_68.00']
    dpl_alpha_68_b1 = data_fits['POSTERIOR PDF'].data['dpl_alpha_68.00']
    dpl_beta_68_b1 = data_fits['POSTERIOR PDF'].data['dpl_beta_68.00']
    tauV_eff_68_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_68.00']
    metallicity_68_b1 = data_fits['POSTERIOR PDF'].data['metallicity_68.00']
    nebular_logU_68_b1 = data_fits['POSTERIOR PDF'].data['nebular_logu_68.00']
    tau_68_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_68.00']
    nebular_xi_68_b1 = data_fits['POSTERIOR PDF'].data['nebular_xi_68.00']
    
    data_fits.close()
    
    # =============================================================================
    # calculate output gradient for rising or falling (<100)
    # =============================================================================
    
    grad_out = np.empty(len(mass_b1))
    A_b1 = np.empty(len(mass_b1))
    integral = np.empty(len(mass_b1))
    
    # nice trick to find index in xlin which has value closest to msa
    idx = (np.abs(xlin - msa)).argmin()
    
    for i in range(len(mass_b1)):
        integrand = lambda T: 1 / (((T/tau_b1[i])**dpl_alpha_b1[i])+((T/tau_b1[i])**-dpl_beta_b1[i]))    
        integral[i]  = quad(integrand, 0, msa)[0]
        A_b1[i] = mass_b1[i] / integral[i]
        # A_b1[i] = 1    # noticed most of the above (integration etc) was pointless as A cancels out in gradient
        sfr_at_msa = A_b1[i] / (((msa/tau_b1[i])**dpl_alpha_b1[i])+((msa/tau_b1[i])**-dpl_beta_b1[i]))
        sfr_at_msa_plus1 = A_b1[i] / ((((1.01*msa)/tau_b1[i])**dpl_alpha_b1[i])+(((1.01*msa)/tau_b1[i])**-dpl_beta_b1[i]))
        grad_out[i] = sfr_at_msa_plus1 - sfr_at_msa
    
    idx_rising_out = grad_out >= 0
    idx_falling_out = grad_out < 0
    
    idx_r1_out = idx_rising_out
    idx_f1_out = idx_falling_out
    
    idx_r1_in = idx_rising_in[id_b1]
    idx_f1_in = idx_falling_in[id_b1]
    
    idx_rr = np.logical_and(idx_r1_in, idx_r1_out)
    idx_ff = np.logical_and(idx_f1_in, idx_f1_out)
    idx_rf = np.logical_and(idx_r1_in, idx_f1_out)
    idx_fr = np.logical_and(idx_f1_in, idx_r1_out)
    
    sum_rr = np.sum(np.logical_and(idx_r1_in, idx_r1_out))
    sum_ff = np.sum(np.logical_and(idx_f1_in, idx_f1_out))
    sum_rf = np.sum(np.logical_and(idx_r1_in, idx_f1_out))
    sum_fr = np.sum(np.logical_and(idx_f1_in, idx_r1_out))
    
    # =======================================================================f======
    # PLOT comparing individual parameters
    # =============================================================================
    
    params_names = ['mass', 'alpha', 'beta', 'tauVeff', 'metallicity', 'nebularlogU', 'tau', 'nebularxi']
    params = [mass, np.log10(dpl_alpha), np.log10(dpl_beta), tauV_eff, metallicity, nebular_logU, np.log10(tau), nebular_xi]
    params_b1 = [mass_b1, np.log10(dpl_alpha_b1), np.log10(dpl_beta_b1), tauV_eff_b1, metallicity_b1, nebular_logU_b1, np.log10(tau_b1), nebular_xi_b1]
    params_68_b1 = [mass_68_b1, np.log10(dpl_alpha_68_b1), np.log10(dpl_beta_68_b1), tauV_eff_68_b1, metallicity_68_b1, nebular_logU_68_b1, np.log10(tau_68_b1), nebular_xi_68_b1]
    
    fig, axs = plt.subplots(2, len(params)/2, figsize=(15, 8))
    fig.suptitle(title1)
    for j in [0, 1]:
        for i in range(len(params)/2):
            axs[j,i].set_title(params_names[i+4*j])
            
            axs[j,i].scatter(params[i+4*j][id_b1][idx_rr], params_b1[i+4*j][idx_rr], s=10, zorder=2, color='r', label='RR {}'.format(sum_rr))
            axs[j,i].scatter(params[i+4*j][id_b1][idx_ff], params_b1[i+4*j][idx_ff], s=10, zorder=2, color='b', label='FF {}'.format(sum_ff))      
            
            axs[j,i].scatter(params[i+4*j][id_b1][idx_rf], params_b1[i+4*j][idx_rf], s=10, zorder=2, color='c', label='RF {}'.format(sum_rf), marker='o')
            axs[j,i].scatter(params[i+4*j][id_b1][idx_fr], params_b1[i+4*j][idx_fr], s=10, zorder=2, color='m', label='FR {}'.format(sum_fr), marker='o')   
            
            axs[j,i].errorbar(params[i+4*j][id_b1], params_b1[i+4*j], yerr=[params_b1[i+4*j] - params_68_b1[i+4*j][:, 0], params_68_b1[i+4*j][:, 1] - params_b1[i+4*j]], linestyle="None", elinewidth=0.5, color='k', zorder=1)
            
            min_ax = min(min(params[i+4*j][id_b1]), min(params_b1[i+4*j]))
#            min_ax = 9.3
            max_ax = max(max(params[i+4*j][id_b1]), max(params_b1[i+4*j]))
            
            axs[j,i].plot((min_ax, max_ax), (min_ax, max_ax), color='k', zorder=0)
            
            axs[j,i].set_xlim(min_ax, max_ax)
            axs[j,i].set_ylim(min_ax, max_ax)
            axs[j,i].legend()
    plt.show()

    # =============================================================================
    # PLOT - input mass vs output mass 1
    # =============================================================================
    
    fsize = 6
    size = 13
    
    plt.figure(figsize=(fsize, fsize))
    plt.title('Input Mass (DPL) vs Output Mass ({})'.format(title1), size=size)
    plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
    plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)

    plt.scatter(mass[id_b1][idx_rr], mass_b1[idx_rr], s=10, zorder=2, color='r', label='RR {}'.format(sum_rr))
    plt.scatter(mass[id_b1][idx_ff], mass_b1[idx_ff], s=10, zorder=2, color='b', label='FF {}'.format(sum_ff))      
    
    plt.scatter(mass[id_b1][idx_rf], mass_b1[idx_rf], s=10, zorder=2, color='c', label='RF {}'.format(sum_rf), marker='o')
    plt.scatter(mass[id_b1][idx_fr], mass_b1[idx_fr], s=10, zorder=2, color='m', label='FR {}'.format(sum_fr), marker='o')   
    
    plt.errorbar(mass[id_b1], mass_b1, yerr=[mass_b1 - mass_68_b1[:, 0], mass_68_b1[:, 1] - mass_b1], linestyle="None", elinewidth=0.5, color='k', zorder=1)
    
    min_ax = 7.5
    max_ax = 11.0
    
    plt.plot((min_ax, max_ax), (min_ax, max_ax), color='k', zorder=0)
    plt.xlim(min_ax, max_ax)
    plt.ylim(min_ax, max_ax)
    plt.legend()
    plt.show()

    # =============================================================================
    # PLOT - input sfr vs output sfr 1
    # =============================================================================

    plt.figure(figsize=(fsize, fsize))
    plt.title('Input SFR (DE) vs Output SFR ({})'.format(title1), size=size)
    plt.xlabel(r'$\text{Input - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
    plt.ylabel(r'$\text{Output - log}(\Psi / M_{\odot} yr^{-1})$', size=size)

    plt.scatter(sfr[id_b1][idx_rr], sfr_b1[idx_rr], s=10, zorder=2, color='r', label='RR {}'.format(sum_rr))
    plt.scatter(sfr[id_b1][idx_ff], sfr_b1[idx_ff], s=10, zorder=2, color='b', label='FF {}'.format(sum_ff))      
    
    plt.scatter(sfr[id_b1][idx_rf], sfr_b1[idx_rf], s=10, zorder=2, color='c', label='RF {}'.format(sum_rf), marker='o')
    plt.scatter(sfr[id_b1][idx_fr], sfr_b1[idx_fr], s=10, zorder=2, color='m', label='FR {}'.format(sum_fr), marker='o')   
    
    plt.errorbar(sfr[id_b1], sfr_b1, yerr=[sfr_b1 - sfr_68_b1[:, 0], sfr_68_b1[:, 1] - sfr_b1], linestyle="None", elinewidth=0.5, color='k', zorder=1)
    
    min_ax = -1.0
    max_ax = 3.5
    
    plt.plot((min_ax, max_ax), (min_ax, max_ax), color='k', zorder=0)
    plt.xlim(min_ax, max_ax)
    plt.ylim(min_ax, max_ax)
    plt.legend()
    plt.show()


    # =============================================================================
    # PLOT - input sfr vs output ssfr 1
    # =============================================================================

    ssfr = sfr - mass # log space    
   
    plt.figure(figsize=(fsize, fsize))
    plt.title('Input SSFR (DE) vs Output SSFR ({})'.format(title1), size=size)
    plt.xlabel(r'$\text{Input - log}(\text{SSFR} / yr^{-1})$', size=size)
    plt.ylabel(r'$\text{Output - log}(SSFR / yr^{-1})$', size=size)
    
    plt.scatter(ssfr[id_b1][idx_rr], ssfr_b1[idx_rr], s=10, zorder=2, color='r', label='RR {}'.format(sum_rr))
    plt.scatter(ssfr[id_b1][idx_ff], ssfr_b1[idx_ff], s=10, zorder=2, color='b', label='FF {}'.format(sum_ff))      
    
    plt.scatter(ssfr[id_b1][idx_rf], ssfr_b1[idx_rf], s=10, zorder=2, color='c', label='RF {}'.format(sum_rf), marker='o')
    plt.scatter(ssfr[id_b1][idx_fr], ssfr_b1[idx_fr], s=10, zorder=2, color='m', label='FR {}'.format(sum_fr), marker='o')   
    
    plt.errorbar(ssfr[id_b1], ssfr_b1, yerr=[ssfr_b1 - ssfr_68_b1[:, 0], ssfr_68_b1[:, 1] - ssfr_b1], linestyle="None", elinewidth=0.5, color='k', zorder=1)

    min_ax = -10.0
    max_ax = -7.0
    
    plt.plot((min_ax, max_ax), (min_ax, max_ax), color='k', zorder=0)
    plt.xlim(min_ax, max_ax)
    plt.ylim(min_ax, max_ax)
    plt.legend()
    plt.show()    
