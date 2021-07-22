#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 17:09:07 2021

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle


def get_data(subfolder):
    AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/data/astrodeep_rawfile_1234_ABCZ.npy'
    AD = np.load(AD_location)
    mag_GMM = np.array([np.log10(AD['MAGNIF'])]*3).transpose() # this is genius
    #print(AD.dtype.names)
    sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/' # real Santini values
    astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/{}/{}_pickle.p'.format(subfolder, subfolder),'r'))    
    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/sfr_santini_temp3.npy' # my 1500, central filter, quoted wavelength method
    sfr_SAN_beta_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/sfr_santini_beta_temp3.npy'    
    
    with np.errstate(divide='ignore', invalid='ignore'):
        data = {    'field_AD':             AD['field'],
                    'id_AD':                AD['ID'], 
                    'mag_AD':               np.log10(AD['MAGNIF']), 
                    'redshift_AD':          AD['ZBEST'], 
                    'mass_AD':              np.log10(AD['MSTAR']*1e9) - np.log10(AD['MAGNIF']), 
                    'mass_AD_neb':          np.log10(AD['MASTAR_NEB']*1e9) - np.log10(AD['MAGNIF']), 
                    'sfr_AD':               np.log10(AD['SFR']) - np.log10(AD['MAGNIF']), 
                    'sfr_AD_neb':           np.log10(AD['SFR_NEB']) - np.log10(AD['MAGNIF']), 
                    'relflag_AD':           AD['RELFLAG'],
    
                    'sfr_SAN':              np.load(sfr_SAN_location) - np.log10(AD['MAGNIF']), 
                    'sfr_SAN_beta':         np.load(sfr_SAN_beta_location), 
    
                    'id_BEAGLE':            astrodeep_pickle['id_BEAGLE'],                     
                    'mass_BEAGLE_tot':      np.log10(astrodeep_pickle['mass_BEAGLE_tot']) - np.log10(AD['MAGNIF']), 
                    'mass_BEAGLE_stellar':  np.log10(astrodeep_pickle['mass_BEAGLE_stellar']) - np.log10(AD['MAGNIF']), 
                    'sfr_BEAGLE_instant':   astrodeep_pickle['sfr_BEAGLE_instant'] - np.log10(AD['MAGNIF']), 
                    'redshift_BEAGLE':      astrodeep_pickle['redshift_BEAGLE'], 
                    'tau_BEAGLE':           astrodeep_pickle['tau_BEAGLE'], 
                    'tauv_BEAGLE':          astrodeep_pickle['tauv_BEAGLE'], 
                    'msa_BEAGLE':           astrodeep_pickle['msa_BEAGLE'], 
                    'metallicity_BEAGLE':   astrodeep_pickle['metallicity_BEAGLE'],
                    'min_chi2_BEAGLE':      astrodeep_pickle['min_chi2_BEAGLE'],
                    
                    'id_GMM_2d':            astrodeep_pickle['id_GMM_2d'],
                    'x_GMM_2d':             astrodeep_pickle['x_GMM_2d'] - mag_GMM,
                    'y_GMM_2d':             astrodeep_pickle['y_GMM_2d'] - mag_GMM,
                    'xsig_GMM_2d':          astrodeep_pickle['xsig_GMM_2d'],
                    'ysig_GMM_2d':          astrodeep_pickle['ysig_GMM_2d'],
                    'xycov_GMM_2d':         astrodeep_pickle['xycov_GMM_2d'],
                    'amp_GMM_2d':           astrodeep_pickle['amp_GMM_2d'],
    
                    'id_GMM_3d':            astrodeep_pickle['id_GMM_3d'],
                    'x_GMM_3d':             astrodeep_pickle['x_GMM_3d'] - mag_GMM,
                    'y_GMM_3d':             astrodeep_pickle['y_GMM_3d'] - mag_GMM,
                    'z_GMM_3d':             astrodeep_pickle['z_GMM_3d'],
                    'xsig_GMM_3d':          astrodeep_pickle['xsig_GMM_3d'],
                    'ysig_GMM_3d':          astrodeep_pickle['ysig_GMM_3d'],
                    'zsig_GMM_3d':          astrodeep_pickle['zsig_GMM_3d'],
                    'xycov_GMM_3d':         astrodeep_pickle['xycov_GMM_3d'],
                    'xzcov_GMM_3d':         astrodeep_pickle['xzcov_GMM_3d'],
                    'yzcov_GMM_3d':         astrodeep_pickle['yzcov_GMM_3d'],
                    'amp_GMM_3d':           astrodeep_pickle['amp_GMM_3d'],
                    
                    'id_SANTINI':           np.load(sbf+'id_SANTINI.npy', allow_pickle=True).astype(float),
                    'mass_SANTINI':         np.load(sbf+'mass_SANTINI.npy', allow_pickle=True).astype(float),
                    'sfr_SANTINI':          np.load(sbf+'sfr_SANTINI.npy', allow_pickle=True).astype(float),
                    'redshift_SANTINI':     np.load(sbf+'redshift_SANTINI.npy', allow_pickle=True).astype(float),
                    'mag_SANTINI':          np.log10(np.load(sbf+'mag_SANTINI.npy', allow_pickle=True).astype(float)) # -103 -> nan
    
                    }
    return data


def get_plots(run1, param1, extra1, run2, param2, extra2, low, high):

    data1 = get_data(run1)
    data2 = get_data(run2)
    
    plt.scatter(data1[param1+extra1], data2[param2+extra2], alpha=0.01)
    plt.xlabel((run1+' '+param1+extra1).replace('_',' '))
    plt.ylabel((run2+' '+param2+extra2).replace('_',' '))
    plt.xlim(low, high)
    plt.ylim(low, high)
    plt.plot((low, high),(low, high),color='k')
    plt.show()
    
    return 


get_plots('delayed_uniform_logtau', 'redshift', '_BEAGLE', 'danger_delayed_uniform_logtau', 'redshift', '_BEAGLE', 0, 10)
get_plots('delayed_uniform_logtau', 'redshift', '_BEAGLE', 'danger_constant', 'redshift', '_BEAGLE', 0, 10)
get_plots('danger_delayed_uniform_logtau', 'redshift', '_BEAGLE', 'danger_constant', 'redshift', '_BEAGLE', 0, 10)

get_plots('delayed_uniform_logtau', 'redshift', '_AD', 'delayed_uniform_logtau', 'redshift', '_BEAGLE', 0, 10)
get_plots('danger_delayed_uniform_logtau', 'redshift', '_AD', 'danger_delayed_uniform_logtau', 'redshift', '_BEAGLE', 0, 10)
get_plots('danger_constant', 'redshift', '_AD', 'danger_constant', 'redshift', '_BEAGLE', 0, 10)


get_plots('delayed_uniform_logtau', 'mass', '_BEAGLE_stellar', 'danger_delayed_uniform_logtau', 'mass', '_BEAGLE_stellar', 5, 13)
get_plots('delayed_uniform_logtau', 'mass', '_BEAGLE_stellar', 'danger_constant', 'mass', '_BEAGLE_stellar', 5, 13)
get_plots('danger_delayed_uniform_logtau', 'mass', '_BEAGLE_stellar', 'danger_constant', 'mass', '_BEAGLE_stellar', 5, 13)

get_plots('delayed_uniform_logtau', 'mass', '_AD', 'delayed_uniform_logtau', 'mass', '_BEAGLE_stellar', 5, 13)
get_plots('danger_delayed_uniform_logtau', 'mass', '_AD', 'danger_delayed_uniform_logtau', 'mass', '_BEAGLE_stellar', 5, 13)
get_plots('danger_constant', 'mass', '_AD', 'danger_constant', 'mass', '_BEAGLE_stellar', 5, 13)


get_plots('delayed_uniform_logtau', 'sfr', '_BEAGLE_instant', 'danger_delayed_uniform_logtau', 'sfr', '_BEAGLE_instant', -4, 4)
get_plots('delayed_uniform_logtau', 'sfr', '_BEAGLE_instant', 'danger_constant', 'sfr', '_BEAGLE_instant', -4, 4)
get_plots('danger_delayed_uniform_logtau', 'sfr', '_BEAGLE_instant', 'danger_constant', 'sfr', '_BEAGLE_instant', -4, 4)

get_plots('delayed_uniform_logtau', 'sfr', '_AD', 'delayed_uniform_logtau', 'sfr', '_BEAGLE_instant', -4, 4)
get_plots('danger_delayed_uniform_logtau', 'sfr', '_AD', 'danger_delayed_uniform_logtau', 'sfr', '_BEAGLE_instant', -4, 4)
get_plots('danger_constant', 'sfr', '_AD', 'danger_constant', 'sfr', '_BEAGLE_instant', -4, 4)


get_plots('danger_constant', 'redshift', '_AD', 'danger_constant', 'redshift', '_AD', 0, 10)























