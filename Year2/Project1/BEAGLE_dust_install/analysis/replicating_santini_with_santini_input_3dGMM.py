#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import copy

# =============================================================================
# NOTES
# =============================================================================

'''
NOTE sfr_SAN has -99.0 when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)

NOTE BEAGLE PARAMS have -101 AD OBJECT WAS NOT A BEAGLE INPUT, -102 AD OBJECT WAS NOT FITTED BY BEAGLE, 
sfr_BEAGLE_instant can also have -30 from during creation of instant sfr

THESE BECOME NAN WHEN TAKING LOG: BEAGLE had -101.0 when the object was not a BEAGLE input and -102.0 if the object was a BEAGLE input, but was NOT fitted

Objects not included by SANTINI have -103.0 for all params
Objects not fitted by GMM 2d or 3d are also set to -103.0 just for GMM params

NOTE the -30s, -101s and -102s aren't strict as magnification was added to them!

# log(0) -> -inf (mass_AD, mass_AD_neb, sfr_AD, sfr_AD_neb)
# lof(-ve) -> nan (mass_BEAGLE_tot and )

'''

# =============================================================================
# SCENARIOS
# =============================================================================

scenarioA = 28

# =============================================================================
# get data
# =============================================================================
AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/data/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
mag_GMM = np.array([np.log10(AD['MAGNIF'])]*3).transpose() # this is genius
print(AD.dtype.names)

#    sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/npy_files_matching_AD_and_BEAGLE_8_fields/'
sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/' # real Santini values

#    astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/saved_astrodeep_pickle/astrodeep_pickle.p','r'))
#    astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/{}/{}_pickle.p'.format(subfolder, subfolder),'r'))

# for new dust beagle only, ignore new chi2 value:
astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/astrodeep_pickle.p','r'))

#ORIGINAL
#astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/recalc_chi2/astrodeep_pickle.p','r'))

#    MY SANTINI VALUES
#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy' # emma technique I think

#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_temp3.npy' # my 1500, central filter, quoted wavelength method
sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/sfr_santini_temp3.npy' # my 1500, central filter, quoted wavelength method

#    sfr_SAN_beta_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_beta_temp3.npy'
sfr_SAN_beta_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/sfr_santini_beta_temp3.npy'    

# =============================================================================
# make combined file
# =============================================================================
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
                'RA_AD':                AD['RA'],
                'DEC_AD':               AD['DEC'],
                
                'b_CH1_AD':             AD['b_CH1'],
                'b_errCH1_AD':          AD['b_errCH1'],
                'b_CH2_AD':             AD['b_CH2'],
                'b_errCH2_AD':          AD['b_errCH2'],

                'sfr_SAN':              np.load(sfr_SAN_location) - np.log10(AD['MAGNIF']), 
                'sfr_SAN_beta':         np.load(sfr_SAN_beta_location), 

                'id_BEAGLE':            astrodeep_pickle['id_BEAGLE'],                     
                'mass_BEAGLE_tot':      np.log10(astrodeep_pickle['mass_BEAGLE_tot']) - np.log10(AD['MAGNIF']), 
                'mass_BEAGLE_stellar':  np.log10(astrodeep_pickle['mass_BEAGLE_stellar']) - np.log10(AD['MAGNIF']), 
                'sfr_BEAGLE_instant':   astrodeep_pickle['sfr_BEAGLE_instant'] - np.log10(AD['MAGNIF']), 
                'redshift_BEAGLE':      astrodeep_pickle['redshift_BEAGLE'], 
                'redshift_BEAGLE_mean':      astrodeep_pickle['redshift_BEAGLE_mean'], 
                'tau_BEAGLE':           astrodeep_pickle['tau_BEAGLE'], 
                'tauv_BEAGLE':          astrodeep_pickle['tauv_BEAGLE'], 
                'msa_BEAGLE':           astrodeep_pickle['msa_BEAGLE'], 
                'metallicity_BEAGLE':   astrodeep_pickle['metallicity_BEAGLE'],
                'min_chi2_BEAGLE':      astrodeep_pickle['min_chi2_BEAGLE'],
                'new_min_chi2_BEAGLE':  astrodeep_pickle['new_min_chi2_BEAGLE'],
                'Ks_BEAGLE_input':      astrodeep_pickle['Ks'],
                'CH1_BEAGLE_input':     astrodeep_pickle['CH1'],
                'CH2_BEAGLE_input':     astrodeep_pickle['CH2'],

                'ch1_beagle_mag_median':     astrodeep_pickle['ch1_beagle_mag_median'],
                'ch1_beagle_mag_lower':     astrodeep_pickle['ch1_beagle_mag_lower'],
                'ch1_beagle_mag_upper':     astrodeep_pickle['ch1_beagle_mag_upper'],
                'ch2_beagle_mag_median':     astrodeep_pickle['ch2_beagle_mag_median'],
                'ch2_beagle_mag_lower':     astrodeep_pickle['ch2_beagle_mag_lower'],
                'ch2_beagle_mag_upper':     astrodeep_pickle['ch2_beagle_mag_upper'],

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

#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/data.p','w')) # all redshifts   


# =============================================================================
# SCENARIO A
# =============================================================================


if scenarioA == 23: 

    idx1 = (AD['field']%2.0==0.0) # clusters
    idx2 = (AD['RELFLAG']==1.0) # relflag
    idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut

    idx3_IRAC = np.logical_and(data['CH1_BEAGLE_input']<-60.0, data['CH2_BEAGLE_input']<-60.0)
    idx3_IRAC = np.logical_and(idx3_IRAC, data['redshift_BEAGLE']>4.0)
    idx3_IRAC = ~idx3_IRAC
    
    MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
    for i in range(len(MCLmassLow)):
        if data['redshift_BEAGLE'][i] <2.1789654:
            MCLmassLow[i] = 8.0
        elif data['redshift_BEAGLE'][i] > 4.195:
            MCLmassLow[i] = 9.0
        else:
            MCLmassLow[i] = 6.91926521 + 0.49598529*data['redshift_BEAGLE'][i]

    idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > MCLmassLow) 
#        print(MCLmassLow, MCLmassLow-0.1)
#        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > -10.0) 
    
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
    idx5_1 = (abs(data['redshift_BEAGLE']-((0.5+6.5)/2.0)) < (((0.5+6.5)/2.0) - 0.5)) # 0.5<redshift_BEAGLE<6.5 (only affects up to visual inspection of z=3.5)
    idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
    idx5_3 = np.logical_and((data['redshift_BEAGLE'] < 3.5), (data['redshift_AD'] < 3.5)) 
    
    idx5_zlt3p5 = np.logical_and(idx5_1,idx5_2)
    idx5_zlt3p5 = np.logical_and(idx5_zlt3p5,idx5_3)
    
    vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_4/investigation_3_and_4_selection.csv', delimiter=",", skip_header=1)

    idx5_zgt3p5 = np.full(len(data['id_AD']), False)
    for i in range(len(data['id_AD'])):
        idx5_zgt3p5_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)])
        if idx5_zgt3p5_temp:
            idx5_zgt3p5[i] = True        
            
    idx5_zgt3p5 = np.logical_and(idx5_zgt3p5, data['redshift_BEAGLE'] < 6.5)

    idx5_z = np.logical_or(idx5_zlt3p5,idx5_zgt3p5)
    
    TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
    idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 

    idx7 = (data['new_min_chi2_BEAGLE']>0) & (data['new_min_chi2_BEAGLE']<13.28) # chi-squared
    idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
    
    idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103

    
    
    
    idx = np.logical_and(idx1,idx2)

    idx = np.logical_and(idx,idx3)

    idx = np.logical_and(idx,idx4)

    idx = np.logical_and(idx,idx6)

    idx = np.logical_and(idx,idx5_z)

    idx = np.logical_and(idx,idx3_IRAC)

    idx = np.logical_and(idx,idx7)

    idx = np.logical_and(idx,idx8)
       
    idx = np.logical_and(idx,idx9)

    print(len(idx), sum(idx))
    
if scenarioA == 24:
    
    idx1 = (data['field_AD']%2.0==1.0) # == 0 clusters, == 1 parallels
    idx2 = (data['relflag_AD']==1.0) # relflag
    idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
    
    idx3_IRAC = np.logical_and(data['CH1_BEAGLE_input']<-60.0, data['CH2_BEAGLE_input']<-60.0)
    idx3_IRAC = ~idx3_IRAC
    
    idx3_KIRAC = np.logical_and(data['CH1_BEAGLE_input']<-60.0, data['CH2_BEAGLE_input']<-60.0)
    idx3_KIRAC = np.logical_and(idx3_KIRAC, data['Ks_BEAGLE_input']<-60.0)
    idx3_KIRAC = ~idx3_KIRAC
    
    MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
    for i in range(len(MCLmassLow)):
        if data['redshift_BEAGLE'][i] <2.1789654:
            MCLmassLow[i] = 8.0
        elif data['redshift_BEAGLE'][i] > 4.195:
            MCLmassLow[i] = 9.0
        else:
            MCLmassLow[i] = 6.91926521 + 0.49598529*data['redshift_BEAGLE'][i]
    idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow - 0.2) )
    
    idx5_z1 = (data['redshift_BEAGLE'] > 4) & (data['redshift_BEAGLE'] < 5)
    idx5_z2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
    vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis1/vis1_selection.csv', delimiter=",", skip_header=1)
    idx5_z3 = np.full(len(data['id_AD']), False)
    for i in range(len(data['id_AD'])):
        idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)])
        if idx5_z3_temp:
            idx5_z3[i] = True
    idx5_z4 = (data['redshift_BEAGLE'] > 3.5) | (data['redshift_AD'] > 3.5)        

    TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies
    idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 
    
    idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<13.28) # chi-squared
    idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
    idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
    
    # clusters, relflag and H<27.5
    print(sum(idx1))
    idx = np.logical_and(idx1,idx2) #clusters+relflag
    print(sum(idx))
    idx = np.logical_and(idx,idx3) #H<27.5
    print(sum(idx))
    
    # redshift bin 
#    idx = np.logical_and(idx,idx5_z1) #4<z<5
    print(sum(idx))
#    idx = np.logical_and(idx,idx5_z2) #beagle within 1 from AD
    print(sum(idx))
#    idx = np.logical_and(idx,idx5_z3) #vis1 visual inspection
    print(sum(idx))
    idx = np.logical_and(idx,idx5_z4) #beagle z > 3.5  or AD z > 3.5
    print(sum(idx))
    
    # upper and lower mass, chi2, arbitrary sfr & 3d GMM
    idx = np.logical_and(idx,idx6) #higher
    print(sum(idx))
    idx = np.logical_and(idx,idx4) #lower
    print(sum(idx))
    idx = np.logical_and(idx,idx7) #chi2
    print(sum(idx))
    idx = np.logical_and(idx,idx8) #arbitrary sfr
    print(sum(idx))
    idx = np.logical_and(idx,idx9) #3d GMM
    print(sum(idx))
    
    # IRAC
#    idx_IRAC_old = np.logical_and(idx,idx3_IRAC)
#    idx_KIRAC_old = np.logical_and(idx,idx3_KIRAC)

    idx = np.logical_and(idx,idx3_KIRAC)   
#    idx = np.logical_and(idx,idx3_IRAC)
    print(sum(idx))
#    print(sum(idx_KIRAC_old))
#    print(sum(idx_IRAC_old))

    

    #i_IRAC = ~((data['CH1_BEAGLE_input']<-60.0)&(data['CH2_BEAGLE_input']<-60.0))
    #i_KIRAC = ~((data['Ks_BEAGLE_input']<-60.0)&(data['CH1_BEAGLE_input']<-60.0)&(data['CH2_BEAGLE_input']<-60.0))
    #i_upper_mass = (data['mass_SANTINI'] < (9.244 + (0.753*4.0) - (0.090*(4.0**2)))) # all galaxies (for z>4)
    #i_lower_mass = (data['mass_SANTINI'] + data['mag_SANTINI'] > MCLmassLow_SANTINI) 
    #i_chi2_beagle_old = (data['new_min_chi2_BEAGLE']>0) & (data['new_min_chi2_BEAGLE']<13.28) # chi-squared
    #i_z_beagle_new = (data['redshift_median'] > 4) & (data['redshift_median'] < 5)


if scenarioA == 25:
    
    idx1 = (data['field_AD']%2.0==1.0) # == 0 clusters, == 1 parallels

    vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis3/vis3_selection.csv', delimiter=",", skip_header=1)
    idx5_z3 = np.full(len(data['id_AD']), False)
    for i in range(len(data['id_AD'])):
        idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,4]==1)])
        if idx5_z3_temp:
            idx5_z3[i] = True

    print(sum(idx1))
    idx = np.logical_and(idx1,idx5_z3) 
    print(sum(idx))

  
if scenarioA == 26:
    
    idx1 = (data['field_AD']%1.0==0.0) # == 0 clusters, == 1 parallels
    idx2 = (data['relflag_AD']==1.0) # relflag
    idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
    
    idx3_KIRAC = np.logical_and(data['CH1_BEAGLE_input']<-60.0, data['CH2_BEAGLE_input']<-60.0)
    idx3_KIRAC = np.logical_and(idx3_KIRAC, data['Ks_BEAGLE_input']<-60.0)
    idx3_KIRAC = ~idx3_KIRAC
    
    MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
    for i in range(len(MCLmassLow)):
        if data['redshift_BEAGLE'][i] <2.1789654:
            MCLmassLow[i] = 8.0
        elif data['redshift_BEAGLE'][i] > 4.195:
            MCLmassLow[i] = 9.0
        else:
            MCLmassLow[i] = 6.91926521 + 0.49598529*data['redshift_BEAGLE'][i]
    idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow) )
    
    idx5_z1 = (data['redshift_BEAGLE'] > 1.25) & (data['redshift_BEAGLE'] < 6.0)
    idx5_z2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
    vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis3/vis3_selection.csv', delimiter=",", skip_header=1)
    idx5_z3 = np.full(len(data['id_AD']), False)
    for i in range(len(data['id_AD'])):
        idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)])
        if idx5_z3_temp:
            idx5_z3[i] = True
    idx5_z5 = (idx5_z3) | ((data['redshift_BEAGLE'] < 3.5) & (data['redshift_AD'] < 3.5) & (idx5_z2))   

    TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies
    idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 
    
    idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<13.28) # chi-squared
    idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
    idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
    
    # clusters, relflag and H<27.5
    print(sum(idx1))
    idx = np.logical_and(idx1,idx2) #clusters+relflag
    print(sum(idx))
    idx = np.logical_and(idx,idx3) #H<27.5
    print(sum(idx))
    
    # redshift bin 
    idx = np.logical_and(idx,idx5_z1) #redshift bin
    print(sum(idx))
    idx = np.logical_and(idx,idx5_z5) #beagle within 1 from AD OR visual inspection
    print(sum(idx))
    
    # upper and lower mass, chi2, arbitrary sfr & 3d GMM
    idx = np.logical_and(idx,idx6) #higher
    print(sum(idx))
    idx = np.logical_and(idx,idx4) #lower
    print(sum(idx))
    idx = np.logical_and(idx,idx7) #chi2
    print(sum(idx))
    idx = np.logical_and(idx,idx8) #arbitrary sfr
    print(sum(idx))
    idx = np.logical_and(idx,idx9) #3d GMM
    print(sum(idx))
    
    idx = np.logical_and(idx,idx3_KIRAC)   
    print(sum(idx))

if scenarioA == 27:
    
    fields = ['clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels']
    z_bins = ['z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0', 'z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0']
    idx_clusters_parallels = [2.0,2.0,2.0,2.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0,1.0]
    z_lower = [1.25, 1.25, 2.0, 3.0, 4.0, 5.0, 1.25, 1.25, 2.0, 3.0, 4.0, 5.0]
    z_upper = [6.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    
    for s in range(len(fields)):
        print(s, len(data['field_AD']))
        idx1 = (data['field_AD']%idx_clusters_parallels[s]==0.0) # 2 == 0 clusters, 2 == 1 parallels, 1 == 0 both
        print(len(idx1))
        idx2 = (data['relflag_AD']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        idx3_KIRAC = np.logical_and(data['CH1_BEAGLE_input']<-60.0, data['CH2_BEAGLE_input']<-60.0)
        idx3_KIRAC = np.logical_and(idx3_KIRAC, data['Ks_BEAGLE_input']<-60.0)
        idx3_KIRAC = ~idx3_KIRAC
        
        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] <2.1789654:
                MCLmassLow[i] = 8.0
            elif data['redshift_BEAGLE'][i] > 4.195:
                MCLmassLow[i] = 9.0
            else:
                MCLmassLow[i] = 6.91926521 + 0.49598529*data['redshift_BEAGLE'][i]
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow) )

        idx5_z1 = (data['redshift_BEAGLE'] > z_lower[s]) & (data['redshift_BEAGLE'] < z_upper[s])
        idx5_z2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
        vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis3/vis3_selection.csv', delimiter=",", skip_header=1)
        idx5_z3 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,5]==1)])
            if idx5_z3_temp:
                idx5_z3[i] = True
        idx5_z5 = (idx5_z3) | ((data['redshift_BEAGLE'] < 3.5) & (data['redshift_AD'] < 3.5) & (idx5_z2))   
    
        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 
        
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<13.28) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        # clusters, relflag and H<27.5
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2) #clusters+relflag
        print(sum(idx))
        idx = np.logical_and(idx,idx3) #H<27.5
        print(sum(idx))
        
        # redshift bin 
        idx = np.logical_and(idx,idx5_z1) #redshift bin
        print(sum(idx))
        idx = np.logical_and(idx,idx5_z5) #beagle within 1 from AD OR visual inspection
        print(sum(idx))
        
        # upper and lower mass, chi2, arbitrary sfr & 3d GMM
        idx = np.logical_and(idx,idx6) #higher
        print(sum(idx))
        idx = np.logical_and(idx,idx4) #lower
        print(sum(idx))
        idx = np.logical_and(idx,idx7) #chi2
        print(sum(idx))
        idx = np.logical_and(idx,idx8) #arbitrary sfr
        print(sum(idx))
        idx = np.logical_and(idx,idx9) #3d GMM
        print(sum(idx))
        
        idx = np.logical_and(idx,idx3_KIRAC)   
        print(sum(idx))
    
        data_new = copy.deepcopy(data)
        for key in data_new.keys():
            data_new[key] = data_new[key][idx] 
        
        pickle.dump(data_new, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, fields[s], z_bins[s], ),'w')) 
        

# VISUAL INSPECTION 4
if scenarioA == 28:
    
    fields = ['clusters+parallels']
    z_bins = ['z1p25-6p0']
    idx_clusters_parallels = [1.0]
    z_lower = [1.25]
    z_upper = [6.0]
    
    for s in range(len(fields)):
        print(s, len(data['field_AD']))
        idx1 = (data['field_AD']%idx_clusters_parallels[s]==0.0) # 2 == 0 clusters, 2 == 1 parallels, 1 == 0 both
        print(len(idx1))
        idx2 = (data['relflag_AD']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        # =============================================================================
        # new IRAC test (do 68 beagle fitted credible intervals overlap with +-1sigma measured input data)
        # =============================================================================
        ch1_beagle_flux_median = (10**6) * (10**((data['ch1_beagle_mag_median'] - 8.9)/(-2.5)))
        ch1_beagle_flux_lower = (10**6) * (10**((data['ch1_beagle_mag_lower'] - 8.9)/(-2.5)))
        ch1_beagle_flux_upper = (10**6) * (10**((data['ch1_beagle_mag_upper'] - 8.9)/(-2.5)))
        
        ch2_beagle_flux_median = (10**6) * (10**((data['ch2_beagle_mag_median'] - 8.9)/(-2.5)))
        ch2_beagle_flux_lower = (10**6) * (10**((data['ch2_beagle_mag_lower'] - 8.9)/(-2.5)))
        ch2_beagle_flux_upper = (10**6) * (10**((data['ch2_beagle_mag_upper'] - 8.9)/(-2.5)))
        
        b_CH1_AD = data['b_CH1_AD']
        b_errCH1_AD = data['b_errCH1_AD']
        b_CH2_AD = data['b_CH2_AD']
        b_errCH2_AD = data['b_errCH2_AD']
        
        ### recap:
        ### reject z>4 objects with CH1<-60 & CH2<-60 if BOTH input+fit don't overlap
        ### reject z>4 objects with CH1<-60 & CH2<-60 & K<-60
        idx_new_IRAC = ((((b_CH1_AD-b_errCH1_AD > ch1_beagle_flux_upper) | (b_CH1_AD+b_errCH1_AD < ch1_beagle_flux_lower)) | \
                       ((b_CH2_AD-b_errCH2_AD > ch2_beagle_flux_upper) | (b_CH2_AD+b_errCH2_AD < ch2_beagle_flux_lower))) & \
                       ((data['CH1_BEAGLE_input'] < -60) & (data['CH2_BEAGLE_input'] < -60))) | \
                       ((data['CH1_BEAGLE_input'] < -60) & (data['CH2_BEAGLE_input'] < -60) & (data['Ks_BEAGLE_input'] < -60))
        idx_new_IRAC = ~(idx_new_IRAC & (data['redshift_BEAGLE']>4.0))
        print(sum(idx_new_IRAC), len(idx_new_IRAC))
        # =============================================================================

        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] < 2.1789654:
                MCLmassLow[i] = 8.0
            elif data['redshift_BEAGLE'][i] > 4.195:
                MCLmassLow[i] = 9.0
            else:
                MCLmassLow[i] = 6.91926521 + 0.49598529*data['redshift_BEAGLE'][i]
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow) )

#        idx5_z1 = (data['redshift_BEAGLE'] > z_lower[s]) & (data['redshift_BEAGLE'] < z_upper[s])
#        idx5_z2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
#        vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis3/vis3_selection.csv', delimiter=",", skip_header=1)
#        idx5_z3 = np.full(len(data['id_AD']), False)
#        for i in range(len(data['id_AD'])):
#            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,5]==1)])
#            if idx5_z3_temp:
#                idx5_z3[i] = True
        idx5_z5 = ((data['redshift_BEAGLE'] > 3.5) | (data['redshift_AD'] > 3.5))
    
        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 
        
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<13.28) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        # clusters, relflag and H<27.5
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2) #clusters+parallels+relflag
        print(sum(idx))
        idx = np.logical_and(idx,idx3) #H<27.5
        print(sum(idx))
        
        # redshift bin 
#        idx = np.logical_and(idx,idx5_z1) #redshift bin
#        print(sum(idx))
        idx = np.logical_and(idx,idx5_z5) 
        print(sum(idx))
        
        # upper and lower mass, chi2, arbitrary sfr & 3d GMM
        idx = np.logical_and(idx,idx6) #higher
        print(sum(idx))
        idx = np.logical_and(idx,idx4) #lower
        print(sum(idx))
        idx = np.logical_and(idx,idx7) #chi2
        print(sum(idx))
        idx = np.logical_and(idx,idx8) #arbitrary sfr
        print(sum(idx))
        idx = np.logical_and(idx,idx9) #3d GMM
        print(sum(idx))
        
        idx = np.logical_and(idx,idx_new_IRAC)   
        print(sum(idx))
    
        data_new = copy.deepcopy(data)
        for key in data_new.keys():
            data_new[key] = data_new[key][idx] 

#        pickle.dump(data_new, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_vis4.p'.format(scenarioA),'w')) 
        
        # SORTING OUT THE 35 NEW OBJECTS WHICH ALSO NEED VISUAL INSPECTION
        
        idx_vis4 = ((data_new['CH1_BEAGLE_input'] < -60) & (data_new['CH2_BEAGLE_input'] < -60) & (data_new['Ks_BEAGLE_input'] < -60))
        print(sum(idx_vis4))
        
        data_x35 = copy.deepcopy(data_new)
        for key in data_x35.keys():
            data_x35[key] = data_x35[key][idx_vis4] 

#        pickle.dump(data_x35, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_vis4_x35.p'.format(scenarioA),'w'))



            
if scenarioA == 29:
    
    fields = ['clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels']
    z_bins = ['z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0', 'z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0']
    idx_clusters_parallels = [2.0,2.0,2.0,2.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0,1.0]
    z_lower = [1.25, 1.25, 2.0, 3.0, 4.0, 5.0, 1.25, 1.25, 2.0, 3.0, 4.0, 5.0]
    z_upper = [6.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0, 2.0, 3.0, 4.0, 5.0, 6.0]

    fields = ['clusters']
    z_bins = ['z3p0-6p0']
    idx_clusters_parallels = [2.0]
    z_lower = [3.0]
    z_upper = [6.0]
    
    for s in range(len(fields)):
        print(s, len(data['field_AD']))
        idx1 = (data['field_AD']%idx_clusters_parallels[s]==0.0) # 2 == 0 clusters, 2 == 1 parallels, 1 == 0 both
        print(len(idx1))
        idx2 = (data['relflag_AD']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] <2.1789654:
                MCLmassLow[i] = 8.0
            elif data['redshift_BEAGLE'][i] > 4.195:
                MCLmassLow[i] = 9.0
            else:
                MCLmassLow[i] = 6.91926521 + 0.49598529*data['redshift_BEAGLE'][i]
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow) )

        idx5_z1 = (data['redshift_BEAGLE'] > z_lower[s]) & (data['redshift_BEAGLE'] < z_upper[s])
        idx5_z2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) & (data['redshift_BEAGLE'] < 3.5) & (data['redshift_AD'] < 3.5)
        vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)
        idx5_z3 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)])

            if idx5_z3_temp:
                idx5_z3[i] = True
        idx5_z = (idx5_z1) & (idx5_z2 | idx5_z3)
    
        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 
        
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<13.28) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        # clusters, relflag and H<27.5
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2) #clusters+relflag
        print(sum(idx))
        idx = np.logical_and(idx,idx3) #H<27.5
        print(sum(idx))
        
        # redshift bin 
        idx = np.logical_and(idx,idx5_z) #redshift bin
        print(sum(idx))
        
        # upper and lower mass, chi2, arbitrary sfr & 3d GMM
        idx = np.logical_and(idx,idx6) #higher
        print(sum(idx))
        idx = np.logical_and(idx,idx4) #lower
        print(sum(idx))
        idx = np.logical_and(idx,idx7) #chi2
        print(sum(idx))
        idx = np.logical_and(idx,idx8) #arbitrary sfr
        print(sum(idx))
        idx = np.logical_and(idx,idx9) #3d GMM
        print(sum(idx))
        

        data_new = copy.deepcopy(data)
        for key in data_new.keys():
            data_new[key] = data_new[key][idx] 
        
        
#        pickle.dump(data_new, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, fields[s], z_bins[s], ),'w')) 
        

# low slope issue, feeding in arbitrary ids and making fits file with details to make SEDs
if scenarioA == 30:
    
    low_slope_file = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis5/low_slope_issue.fits'

    low_slope_data = fits.open(low_slope_file)
    #print(low_slope_data.info())
    #print(low_slope_data[1].columns)

    data_uid = []
    low_slope_data_uid = []
    
    for i in range(len(low_slope_data[1].data['field_AD'])):
        low_slope_data_uid.append('{}_{}'.format(low_slope_data[1].data['field_AD'][i], low_slope_data[1].data['id_AD'][i]))
        
    for i in range(len(data['field_AD'])):
        data_uid.append('{}_{}'.format(int(data['field_AD'][i]), int(data['id_AD'][i])))

    idx = np.isin(data_uid, low_slope_data_uid)

    data_new = copy.deepcopy(data)
    for key in data_new.keys():
        data_new[key] = data_new[key][idx] 
    
    
#    pickle.dump(data_new, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/vis5_low_slope_issue.p','w')) 
    
    

#    print(data['field_AD'], data['id_AD'], data['id_BEAGLE'], data['redshift_BEAGLE'], data['redshift_AD'])

#    plt.scatter(data['mag_AD'], data['mag_SANTINI'], alpha=0.2)
#    plt.plot((-0.3,1.5),(-0.3,1.5), color='k')
#    plt.title('Scenario 23')
#    plt.xlabel('AD mag')
#    plt.ylabel('SANTINI mag')
#    plt.show()


# =============================================================================
# Create INPUT for KELLY
# =============================================================================

#     SCENARIO 
#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/{}/scenario_{}_data_z{}.p'.format(subfolder, scenarioA, str(zLow).replace('.','p')),'w'))
#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/{}/scenario_{}_data_z0p5.p'.format(subfolder, scenarioA),'w')) # all redshifts
   
#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/{}/scenario_{}_subset_zgt5.p'.format(subfolder, scenarioA),'w')) # all redshifts
#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_14.p'.format(scenarioA),'w')) 

#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_vis1.p'.format(scenarioA),'w')) 

#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_vis2.p'.format(scenarioA),'w')) 

#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_vis3.p'.format(scenarioA),'w')) 

#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_vis3_check_clusters.p'.format(scenarioA),'w')) 

#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_vis3_check_parallels.p'.format(scenarioA),'w')) 










