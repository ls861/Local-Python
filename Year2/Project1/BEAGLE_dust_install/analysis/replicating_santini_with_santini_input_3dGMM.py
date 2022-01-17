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
from astropy.table import Table
import os

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

scenarioA = 35

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


from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import matplotlib as mpl

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)



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
        ### reject BEAGLE z>4 objects with CH1<-60 & CH2<-60 if BOTH input+fit don't overlap
        ### reject BEAGLE z>4 objects with CH1<-60 & CH2<-60 & K<-60
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

#    fields = ['clusters']
#    z_bins = ['z3p0-6p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [3.0]
#    z_upper = [6.0]

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
        # vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)
        vis = np.genfromtxt('/Users/lester/Documents/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)


        idx5_z3 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)]) # prior to z3 above MS visual inspection
#            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,4]==1)])

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

        print('FINAL: {}'.format(len(data_new[key])))
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


#same as 29, but with 16 objects removed from z3 bin by further visual inspection of objects above MS (4/8/21)
if scenarioA == 31:

    fields = ['clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels']
    z_bins = ['z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0', 'z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0']
    idx_clusters_parallels = [2.0,2.0,2.0,2.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0,1.0]
    z_lower = [1.25, 1.25, 2.0, 3.0, 4.0, 5.0, 1.25, 1.25, 2.0, 3.0, 4.0, 5.0]
    z_upper = [6.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0, 2.0, 3.0, 4.0, 5.0, 6.0]

#    fields = ['clusters']
#    z_bins = ['z3p0-6p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [3.0]
#    z_upper = [6.0]

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
        # vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)
        vis = np.genfromtxt('/Users/lester/Documents/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)


        idx5_z3 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
#            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)]) # prior to z3 above MS visual inspection
            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,4]==1)])

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

        print('FINAL: {}'.format(len(data_new[key])))
        pickle.dump(data_new, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, fields[s], z_bins[s], ),'w'))




#same as 31, but trying to include all redshifts to investigate redshift distribution (ie from z=0, not 1.25)
if scenarioA == 32:

    fields = ['clusters','clusters+parallels']
    z_bins = ['z0p0-10p0','z0p0-10p0']
    idx_clusters_parallels = [2.0,1.0]
    z_lower = [0.0,0.0]
    z_upper = [10.0,10.0]

#    z_lower = [1.25,1.25]
#    z_upper = [6.0,6.0]

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
        # vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)
        vis = np.genfromtxt('/Users/lester/Documents/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)


        idx5_z3 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
#            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)]) # prior to z3 above MS visual inspection
            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,4]==1)])

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

        print('FINAL: {}'.format(len(data_new[key])))
#        pickle.dump(data_new, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, fields[s], z_bins[s]),'w'))

        outputTable = Table(data_new)
        outputTable.write('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.fits'.format(scenarioA, fields[s], z_bins[s]), overwrite=True)





#same as 31, but with large sigma sfr objects removed per z bin, and then removed small sfr sig sfr GMM peaks
if scenarioA == 33:

    fields = ['clusters', 'clusters', 'clusters', 'clusters', 'clusters']
    z_bins = ['z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0']
    idx_clusters_parallels = [2.0,2.0,2.0,2.0,2.0]
    z_lower = [1.25, 2.0, 3.0, 4.0, 5.0]
    z_upper = [2.0, 3.0, 4.0, 5.0, 6.0]

#    fields = ['clusters']
#    z_bins = ['z1p25-6p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [1.25]
#    z_upper = [6.0]

#    fields = ['clusters']
#    z_bins = ['z1p25-2p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [1.25]
#    z_upper = [2.0]
#
#    fields = ['clusters']
#    z_bins = ['z2p0-3p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [2.0]
#    z_upper = [3.0]
#
#    fields = ['clusters']
#    z_bins = ['z5p0-6p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [5.0]
#    z_upper = [6.0]


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
        # vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)
        vis = np.genfromtxt('/Users/lester/Documents/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)


        idx5_z3 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
#            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)]) # prior to z3 above MS visual inspection
            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,4]==1)])

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



        ### REMOVE OBJECTS WITH LARGE SFR STD DEV

#        z_med_hp = (z_lower[s]+z_upper[s])/2.0
#        z_med_hp_gap = (z_lower[s]+z_upper[s])/2.0 - z_lower[s]

        z_med_hp = np.empty(len(data_new['id_AD']))
        z_med_hp_gap = np.empty(len(data_new['id_AD']))

        for i in range(len(data_new['id_AD'])):
            if (data_new['redshift_BEAGLE'][i] > 1.25) & (data_new['redshift_BEAGLE'][i] < 2.0):
                z_med_hp[i] = (1.25 + 2.0)/2.0
                z_med_hp_gap[i] = (1.25 + 2.0)/2.0 - 1.25
            elif (data_new['redshift_BEAGLE'][i] > 2.0) & (data_new['redshift_BEAGLE'][i] < 3.0):
                z_med_hp[i] = (2.0 + 3.0)/2.0
                z_med_hp_gap[i] = (2.0 + 3.0)/2.0 - 2.0
            elif (data_new['redshift_BEAGLE'][i] > 3.0) & (data_new['redshift_BEAGLE'][i] < 4.0):
                z_med_hp[i] = (3.0 + 4.0)/2.0
                z_med_hp_gap[i] = (3.0 + 4.0)/2.0 - 3.0
            elif (data_new['redshift_BEAGLE'][i] > 4.0) & (data_new['redshift_BEAGLE'][i] < 5.0):
                z_med_hp[i] = (4.0 + 5.0)/2.0
                z_med_hp_gap[i] = (4.0 + 5.0)/2.0 - 4.0
            elif (data_new['redshift_BEAGLE'][i] > 5.0) & (data_new['redshift_BEAGLE'][i] < 6.0):
                z_med_hp[i] = (5.0 + 6.0)/2.0
                z_med_hp_gap[i] = (5.0 + 6.0)/2.0 - 5.0


        x_hp = np.array([])
        y_hp = np.array([])
        z_hp = np.array([])

        n_hp = 3000 # number of samples to take from GMM in total

        fig, ax = plt.subplots(1, 1, figsize=(15, 15))

        sfr_std = []

        count_blue = 1

        for i in range(len(data_new['id_AD'])):
            if i <100000:

                x_temp = np.array([])
                y_temp = np.array([])
                z_temp = np.array([])

                for G in range(3):

                    mean = np.array([data_new['x_GMM_3d'][i,G],data_new['y_GMM_3d'][i,G],data_new['z_GMM_3d'][i,G]])
                    cov = np.array([[np.power(data_new['xsig_GMM_3d'][i,G],2), data_new['xycov_GMM_3d'][i,G], data_new['xzcov_GMM_3d'][i,G]],[data_new['xycov_GMM_3d'][i,G], np.power(data_new['ysig_GMM_3d'][i,G],2), data_new['yzcov_GMM_3d'][i,G]],[data_new['xzcov_GMM_3d'][i,G], data_new['yzcov_GMM_3d'][i,G], np.power(data_new['zsig_GMM_3d'][i,G],2)]])

                    xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*data_new['amp_GMM_3d'][i,G]))


                    z_idx = (abs(xyz[:,2] - z_med_hp[i]) < z_med_hp_gap[i])

                    if True: # only keep samples within z bin
                        x_temp = np.concatenate((x_temp,xyz[:,0][z_idx]))
                        y_temp = np.concatenate((y_temp,xyz[:,1][z_idx]))
                        z_temp = np.concatenate((z_temp,xyz[:,2][z_idx]))

                        x_hp = np.concatenate((x_hp,xyz[:,0][z_idx]))
                        y_hp = np.concatenate((y_hp,xyz[:,1][z_idx]))
                        z_hp = np.concatenate((z_hp,xyz[:,2][z_idx]))

                    else:
                        x_hp = np.concatenate((x_hp,xyz[:,0]))
                        y_hp = np.concatenate((y_hp,xyz[:,1]))
                        z_hp = np.concatenate((z_hp,xyz[:,2]))

                        x_temp = np.concatenate((x_temp,xyz[:,0]))
                        y_temp = np.concatenate((y_temp,xyz[:,1]))
                        z_temp = np.concatenate((z_temp,xyz[:,2]))

                # ax.scatter(x_temp, y_temp, s=0.5)
                # confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='firebrick', linewidth=1)

                sfr_std.append(np.std(y_temp))

                sfr_limit = 2.5
                mass_limit = 9.5

                if np.std(y_temp) > 2:
                    confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='blue', linewidth=1)
                    print('MMMMMMM', count_blue)
                    count_blue += 1

                elif (data_new['x_GMM_3d'][i,0] > mass_limit and data_new['y_GMM_3d'][i,0] > sfr_limit) or (data_new['x_GMM_3d'][i,1] > mass_limit and data_new['y_GMM_3d'][i,1] > sfr_limit) or (data_new['x_GMM_3d'][i,2] > mass_limit and data_new['y_GMM_3d'][i,2] > sfr_limit):

                    confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='green', linewidth=1)
                    print(data_new['field_AD'][i], data_new['id_AD'][i], data_new['id_BEAGLE'][i])
                    print(data_new['x_GMM_3d'][i], data_new['y_GMM_3d'][i])
                else:
                    confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='firebrick', linewidth=1)




        ax.set_xlim(6, 11)
        # ax.set_ylim(-3, 4)
        ax.set_ylim(-20, 4)

        # ax.set_xlim(8, 11)
        # ax.set_ylim(0, 4)

        # ax.set_xlim(6.5, 10.5)
        # ax.set_ylim(-5.5, 2.5)

        plt.show()


        plt.hist(sfr_std, bins=30)
        plt.show()

        sfr_std = np.array(sfr_std)
        print(len(sfr_std))
        print(len(sfr_std[sfr_std>2]))


        plt.scatter(data_new['mass_BEAGLE_stellar'][sfr_std<2], data_new['sfr_BEAGLE_instant'][sfr_std<2], label=len(data_new['mass_BEAGLE_stellar'][sfr_std<2]))
        plt.scatter(data_new['mass_BEAGLE_stellar'][sfr_std>2], data_new['sfr_BEAGLE_instant'][sfr_std>2], label=len(data_new['mass_BEAGLE_stellar'][sfr_std>2]))
        plt.title('scenario{} {} {}'.format(scenarioA, fields[s], z_bins[s]))
        plt.xlabel('stellar mass')
        plt.ylabel('instant sfr')
        plt.legend()
        plt.show()


        data_new_sfr = copy.deepcopy(data_new)
        for key in data_new_sfr.keys():
            data_new_sfr[key] = data_new_sfr[key][sfr_std<2]



        #%%
        # =============================================================================
        # remove 'bad' gmm peaks
        # =============================================================================

        plt.scatter(data_new_sfr['y_GMM_3d'][:,0], data_new_sfr['ysig_GMM_3d'][:,0])
        plt.scatter(data_new_sfr['y_GMM_3d'][:,1], data_new_sfr['ysig_GMM_3d'][:,1])
        plt.scatter(data_new_sfr['y_GMM_3d'][:,2], data_new_sfr['ysig_GMM_3d'][:,2])
        x_tmp = np.array([-40.0, 5.0])
        plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
        plt.show()


        idx_sort = np.argsort(data_new_sfr['amp_GMM_3d'].flatten())

        #MASS
        plt.title('mass, each GMM, sigma vs mean')
        plt.scatter(data_new_sfr['x_GMM_3d'].flatten()[idx_sort], data_new_sfr['xsig_GMM_3d'].flatten()[idx_sort], c=data_new_sfr['amp_GMM_3d'].flatten()[idx_sort])
        # x_tmp = np.array([-40.0, 5.0])
        # plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
        plt.colorbar()
        plt.show()

        cmap = mpl.cm.get_cmap("tab20", 20)
        #SFR
        plt.title('SFR, each GMM, sigma vs mean')
        plt.scatter(data_new_sfr['y_GMM_3d'].flatten()[idx_sort], data_new_sfr['ysig_GMM_3d'].flatten()[idx_sort], c=data_new_sfr['amp_GMM_3d'].flatten()[idx_sort], cmap=cmap)
        x_tmp = np.array([-40.0, 5.0])
        plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
        # plt.xlim(-5,5)
        # plt.ylim(0,6)
        plt.colorbar()
        plt.show()

        #Z
        plt.title('redshift, each GMM, sigma vs mean')
        plt.scatter(data_new_sfr['z_GMM_3d'].flatten()[idx_sort], data_new_sfr['zsig_GMM_3d'].flatten()[idx_sort], c=data_new_sfr['amp_GMM_3d'].flatten()[idx_sort])
        # x_tmp = np.array([-40.0, 5.0])
        # plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
        # plt.xlim(0.5,3)
        # plt.ylim(0,0.5)
        plt.colorbar()
        plt.show()



        #%% for objects where I remove a peak, show full GMM weighting + mean SFR + sig SFR

        print('bad gmm peaks')
        for i in range(len(data_new_sfr['amp_GMM_3d'])):
            if data_new_sfr['ysig_GMM_3d'][i,0] < (-9.0/26.0)*data_new_sfr['y_GMM_3d'][i,0] - (9.0/13.0) or data_new_sfr['ysig_GMM_3d'][i,1] < (-9.0/26.0)*data_new_sfr['y_GMM_3d'][i,1] - (9.0/13.0) or data_new_sfr['ysig_GMM_3d'][i,2] < (-9.0/26.0)*data_new_sfr['y_GMM_3d'][i,2] - (9.0/13.0):
                print(data_new_sfr['amp_GMM_3d'][i])
                print(data_new_sfr['y_GMM_3d'][i])
                print(data_new_sfr['ysig_GMM_3d'][i])
                print('')


        #%% actually remove bad peaks
        data_new_sfr['amp_GMM_3d'] = np.where(data_new_sfr['ysig_GMM_3d'] < (-9.0/26.0)*data_new_sfr['y_GMM_3d'] - (9.0/13.0), 0.0, data_new_sfr['amp_GMM_3d']) # set bad gmms to p==0
        data_new_sfr['amp_GMM_3d'] = data_new_sfr['amp_GMM_3d'] / np.array([np.sum(data_new_sfr['amp_GMM_3d'], axis=1)]*3).T # renormalise - GENIUS

        print('FINAL: {}'.format(len(data_new_sfr[key])))

#        pickle.dump(data_new_sfr, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, fields[s], z_bins[s]),'w'))
#        outputTable = Table(data_new_sfr)
#        outputTable.write('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.fits'.format(scenarioA, fields[s], z_bins[s]), overwrite=True)
#
#

#%%
# applying pick a peak to scenario 33 full redshift for clusters
if scenarioA == 34:

    field = 'clusters'
    z_bin = 'z1p25-6p0'
    z_lower = 1.25
    z_upper = 6.0

#    field = 'clusters'
#    z_bin = 'z1p25-6p0' # NEED TO MAKE     z_bin = 'z1p25-5p0' LATER ON
#    z_lower = 1.25
#    z_upper = 5.0

#    field = 'clusters'
#    z_bin = 'z1p25-6p0'
#    z_lower = 1.25
#    z_upper = 3.0

#    field = 'clusters'
#    z_bin = 'z1p25-6p0'
#    z_lower = 2.0
#    z_upper = 4.0
#
#    field = 'clusters'
#    z_bin = 'z1p25-6p0'
#    z_lower = 3.0
#    z_upper = 5.0
#
#    field = 'clusters'
#    z_bin = 'z1p25-6p0'
#    z_lower = 1.25
#    z_upper = 4.0
##
#    field = 'clusters'
#    z_bin = 'z1p25-6p0'
#    z_lower = 2.0
#    z_upper = 5.0








    with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format('33', field, z_bin), 'rb') as f:
        data = pickle.load(f)

    filename = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/pick_a_peak/from_emma/scenario_29_clusters_z1p25-6p0_post_pap.fits'
    pap = fits.open(filename)
#    print(data_pap.info())
#    print(data_pap[1].data)
    pap = pap[1].data
#    print(data_pap['id_AD'])

    uid_pap_arr = []
    for i in range(len(pap['id_AD'])):
        uid_pap_arr.append(pap['field_AD'][i].astype(int).astype(str) + '_' + pap['id_AD'][i].astype(int).astype(str))
    uid_pap_arr = np.array(uid_pap_arr)


    idx_pap = np.full(len(data['id_AD']), True)

    for i in range(len(data['id_AD'])):

        uid = str(int(data['field_AD'][i]))+'_'+str(int(data['id_AD'][i]))

        idx_isin = np.isin(uid_pap_arr, uid)



        if sum(pap['amp_GMM_3d_post_pap'][idx_isin][0]) == 0.0: # 0 0 0 amplitude in pick a peak, means remove object
            print('                                  PAP REMOVAL {}'.format(uid))
            idx_pap[i] = False

        if np.prod(data['amp_GMM_3d'][i]) == 0: # before looking at pick a peak, if 1 or more amplitudes = 0, then I've removed a peak
            print('                 PAP CHECK {}'.format(uid))
            print(data['amp_GMM_3d'][i])

        if np.prod(pap['amp_GMM_3d_post_pap'][idx_isin][0]) == 0: # now we take the pick a peak amplitudes, and then check whether they clash with previously removed ones
#            print('')
#            print('PAP UPDATE {}'.format(uid))
#            print(data['amp_GMM_3d'][i])
            data['amp_GMM_3d'][i] = pap['amp_GMM_3d_post_pap'][idx_isin][0]
#            print(data['amp_GMM_3d'][i])

            boo1 = data['ysig_GMM_3d'][i,0] < (-9.0/26.0)*data['y_GMM_3d'][i,0] - (9.0/13.0)
            boo2 = data['ysig_GMM_3d'][i,1] < (-9.0/26.0)*data['y_GMM_3d'][i,1] - (9.0/13.0)
            boo3 = data['ysig_GMM_3d'][i,2] < (-9.0/26.0)*data['y_GMM_3d'][i,2] - (9.0/13.0)
#            print(boo1, boo2, boo3)
#            print(data['redshift_BEAGLE'][i])
#            print(data['x_GMM_3d'][i])
#            print(data['y_GMM_3d'][i])
#            print(data['ysig_GMM_3d'][i])
#            print(data['z_GMM_3d'][i])
#            print(data['id_BEAGLE'][i])
#            print('')



        data_pap = copy.deepcopy(data)
        for key in data_pap.keys():
            data_pap[key] = data_pap[key][idx_pap]





    # =============================================================================
    # remove 'bad' gmm peaks AGAIN
    # =============================================================================

    plt.scatter(data_pap['y_GMM_3d'][:,0], data_pap['ysig_GMM_3d'][:,0])
    plt.scatter(data_pap['y_GMM_3d'][:,1], data_pap['ysig_GMM_3d'][:,1])
    plt.scatter(data_pap['y_GMM_3d'][:,2], data_pap['ysig_GMM_3d'][:,2])
    x_tmp = np.array([-40.0, 5.0])
    plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
    plt.show()


    idx_sort = np.argsort(data_pap['amp_GMM_3d'].flatten())

    #MASS
    plt.title('mass, each GMM, sigma vs mean')
    plt.scatter(data_pap['x_GMM_3d'].flatten()[idx_sort], data_pap['xsig_GMM_3d'].flatten()[idx_sort], c=data_pap['amp_GMM_3d'].flatten()[idx_sort])
    # x_tmp = np.array([-40.0, 5.0])
    # plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
    plt.colorbar()
    plt.show()

    cmap = mpl.cm.get_cmap("tab20", 20)
    #SFR
    plt.title('SFR, each GMM, sigma vs mean')
    plt.scatter(data_pap['y_GMM_3d'].flatten()[idx_sort], data_pap['ysig_GMM_3d'].flatten()[idx_sort], c=data_pap['amp_GMM_3d'].flatten()[idx_sort], cmap=cmap)
    x_tmp = np.array([-40.0, 5.0])
    plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
    # plt.xlim(-5,5)
    # plt.ylim(0,6)
    plt.colorbar()
    plt.show()

    #Z
    plt.title('redshift, each GMM, sigma vs mean')
    plt.scatter(data_pap['z_GMM_3d'].flatten()[idx_sort], data_pap['zsig_GMM_3d'].flatten()[idx_sort], c=data_pap['amp_GMM_3d'].flatten()[idx_sort])
    # x_tmp = np.array([-40.0, 5.0])
    # plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
    # plt.xlim(0.5,3)
    # plt.ylim(0,0.5)
    plt.colorbar()
    plt.show()



    #%% for objects where I remove a peak, show full GMM weighting + mean SFR + sig SFR

    print('bad gmm peaks')
    for i in range(len(data_pap['amp_GMM_3d'])):
        if (data_pap['ysig_GMM_3d'][i,0] < (-9.0/26.0)*data_pap['y_GMM_3d'][i,0] - (9.0/13.0) and  data_pap['amp_GMM_3d'][i,0] > 0.0)  \
        or (data_pap['ysig_GMM_3d'][i,1] < (-9.0/26.0)*data_pap['y_GMM_3d'][i,1] - (9.0/13.0) and  data_pap['amp_GMM_3d'][i,1] > 0.0)  \
        or (data_pap['ysig_GMM_3d'][i,2] < (-9.0/26.0)*data_pap['y_GMM_3d'][i,2] - (9.0/13.0) and  data_pap['amp_GMM_3d'][i,2] > 0.0)  :
            print(data_pap['amp_GMM_3d'][i])
            print(data_pap['y_GMM_3d'][i])
            print(data_pap['ysig_GMM_3d'][i])
            print(data_pap['field_AD'][i], data_pap['id_AD'][i])
            print('')


    #%% actually remove bad peaks
    data_pap['amp_GMM_3d'] = np.where(data_pap['ysig_GMM_3d'] < (-9.0/26.0)*data_pap['y_GMM_3d'] - (9.0/13.0), 0.0, data_pap['amp_GMM_3d']) # set bad gmms to p==0
    data_pap['amp_GMM_3d'] = data_pap['amp_GMM_3d'] / np.array([np.sum(data_pap['amp_GMM_3d'], axis=1)]*3).T # renormalise - GENIUS


#    for i in range(len(data_pap['amp_GMM_3d'])):
#        if np.prod(data_pap['amp_GMM_3d'][i]) == 0:
#            print(data_pap['field_AD'][i], data_pap['id_AD'][i], data_pap['amp_GMM_3d'][i])
#
#


    # =============================================================================
    # REMOVE OBJECTS WITH MEDIANS WHICH ARE NOW OUTSIDE Z RANGE
    # =============================================================================
    #%%
    z_new_medians = []
    for i in range(len(data_pap['id_BEAGLE'])):

        n_hp = 30000 # number of samples to take from GMM in total
        z_temp = np.array([])

        for G in range(3):

            mean = data_pap['z_GMM_3d'][i,G]
            std = data_pap['zsig_GMM_3d'][i,G]

            z = np.random.normal(mean, std, size=int(n_hp*data_pap['amp_GMM_3d'][i,G]))
            z_temp = np.concatenate((z_temp, z))

        z_new_medians.append(np.median(z_temp))

        print('')
        print('here', len(z_temp))
        print(data_pap['id_BEAGLE'][i], data_pap['redshift_BEAGLE'][i], np.median(z_temp))
        #%%

    z_new_medians = np.array(z_new_medians)

    print(z_new_medians[z_new_medians<z_lower])
    print(z_new_medians[z_new_medians>z_upper])

    print(data_pap['id_AD'][z_new_medians<z_lower])

    plt.hist(z_new_medians, bins=30)
    plt.show()


    print('FINAL: {}'.format(len(data_pap[key])))

    idx_pap_2 = (z_new_medians<z_lower) | (z_new_medians>z_upper)

    data_pap_2 = copy.deepcopy(data_pap)
    for key in data_pap_2.keys():
        data_pap_2[key] = data_pap_2[key][~idx_pap_2]


    print('FINAL: {}'.format(len(data_pap_2[key])))


#    z_bin = 'z1p25-5p0'
#    z_bin = 'z1p25-3p0'
#    z_bin = 'z2p0-4p0'
#    z_bin = 'z3p0-5p0'
#    z_bin = 'z1p25-4p0'
#    z_bin = 'z2p0-5p0'


#    pickle.dump(data_pap_2, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, field, z_bin),'w'))
#    outputTable = Table(data_pap_2)
#    outputTable.write('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.fits'.format(scenarioA, field, z_bin), overwrite=True)

#



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









#same as 33, but with iterative sigma clipping before sfr object removal
if scenarioA == 35:

#    fields = ['clusters', 'clusters', 'clusters', 'clusters', 'clusters']
#    z_bins = ['z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0']
#    idx_clusters_parallels = [2.0,2.0,2.0,2.0,2.0]
#    z_lower = [1.25, 2.0, 3.0, 4.0, 5.0]
#    z_upper = [2.0, 3.0, 4.0, 5.0, 6.0]

#    fields = ['clusters']
#    z_bins = ['z1p25-6p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [1.25]
#    z_upper = [6.0]

#    fields = ['clusters']
#    z_bins = ['z1p25-2p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [1.25]
#    z_upper = [2.0]
#
#    fields = ['clusters']
#    z_bins = ['z2p0-3p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [2.0]
#    z_upper = [3.0]
#
#    fields = ['clusters']
#    z_bins = ['z5p0-6p0']
#    idx_clusters_parallels = [2.0]
#    z_lower = [5.0]
#    z_upper = [6.0]


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
        # vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)
        vis = np.genfromtxt('/Users/lester/Documents/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis4/vis4_selection.csv', delimiter=",", skip_header=1)


        idx5_z3 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
#            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)]) # prior to z3 above MS visual inspection
            idx5_z3_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,4]==1)])

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


        data_temp = copy.deepcopy(data)
        for key in data_temp.keys():
            data_temp[key] = data_temp[key][idx]



        # =============================================================================
        # SIGMA CLIPPING
        # =============================================================================

        from scipy.optimize import curve_fit
        
        def straight_line(x, A, B): # this is your 'straight line' y=f(x)
            return B*(x-9.7) + A
        
        xi_sampled = data_temp['mass_BEAGLE_stellar']
        eta_sampled = data_temp['sfr_BEAGLE_instant']
        id_sampled = np.arange(len(xi_sampled))
            
        outliers = 1
        
        xi_outliers = []
        eta_outliers = []
        id_outliers = []
        
        count = 0
        print('')
        print('sigma clipping')
        while outliers > 0:
        
            popt, pcov = curve_fit(straight_line, xi_sampled, eta_sampled)
            beta = popt[1]
            alphaN = popt[0]
        
            # list of sfrs according to input values
            eta_from_relation = beta*(xi_sampled-9.7) + alphaN
            eta_residuals = eta_sampled - eta_from_relation
            eta_sigma = np.std(eta_residuals)
            eta_idx = (abs(eta_residuals)<3.0*eta_sigma)
        
        
            count+=1
            print(count, len(xi_outliers))
            plt.scatter(xi_sampled, eta_sampled)
            plt.scatter(xi_outliers, eta_outliers)
            plt.show() 
        
        
            outliers = sum(~eta_idx)
            
            for oo in range(outliers):
                xi_outliers.append(xi_sampled[~eta_idx][oo])
                eta_outliers.append(eta_sampled[~eta_idx][oo])
                id_outliers.append(id_sampled[~eta_idx][oo])
                
            xi_sampled = xi_sampled[eta_idx]
            eta_sampled = eta_sampled[eta_idx]
            id_sampled = id_sampled[eta_idx]
        
        id_outliers = np.array(id_outliers)

        print(id_outliers)

        idx_sigma = np.full(len(data_temp['mass_BEAGLE_stellar']), False)
        if len(id_outliers) > 0:
            idx_sigma[id_outliers] = True
        
        plt.scatter(data_temp['mass_BEAGLE_stellar'][~idx_sigma], data_temp['sfr_BEAGLE_instant'][~idx_sigma])
        plt.scatter(data_temp['mass_BEAGLE_stellar'][idx_sigma], data_temp['sfr_BEAGLE_instant'][idx_sigma])
        plt.title('scenario{} {} {}'.format(scenarioA, fields[s], z_bins[s]))
        plt.xlabel('stellar mass')
        plt.ylabel('instant sfr')
        plt.legend()
        plt.show()
        
        ### idx3 log error testing for 0 and -ve fluxes
#        t = (-2.5*np.log10(1.0*1e-6) + 8.90 < 27.5) # H band cut
#        print(t)# TRUE
#        t = (-2.5*np.log10(0.000001*1e-6) + 8.90 < 27.5) # H band cut
#        print(t)# FALSE
#        t = (-2.5*np.log10(0.0*1e-6) + 8.90 < 27.5) # H band cut
#        print(t)# FALSE
#        t = (-2.5*np.log10(0.000001*1e-6) + 8.90 < 27.5) # H band cut
#        print(t)# FALSE
        
        ### idx4 error testing - data['mass_BEAGLE_stellar'] contains some nans
        ### idx6 error testing - data['mass_BEAGLE_stellar'] contains some nans

        

        data_new = copy.deepcopy(data_temp)
        for key in data_new.keys():
            data_new[key] = data_new[key][~idx_sigma]


        # =============================================================================
        # BACK TO REMOVING OBJECTS and then PEAKS    
        # =============================================================================

        ### REMOVE OBJECTS WITH LARGE SFR STD DEV

#        z_med_hp = (z_lower[s]+z_upper[s])/2.0
#        z_med_hp_gap = (z_lower[s]+z_upper[s])/2.0 - z_lower[s]

        z_med_hp = np.empty(len(data_new['id_AD']))
        z_med_hp_gap = np.empty(len(data_new['id_AD']))

        for i in range(len(data_new['id_AD'])):
            if (data_new['redshift_BEAGLE'][i] > 1.25) & (data_new['redshift_BEAGLE'][i] < 2.0):
                z_med_hp[i] = (1.25 + 2.0)/2.0
                z_med_hp_gap[i] = (1.25 + 2.0)/2.0 - 1.25
            elif (data_new['redshift_BEAGLE'][i] > 2.0) & (data_new['redshift_BEAGLE'][i] < 3.0):
                z_med_hp[i] = (2.0 + 3.0)/2.0
                z_med_hp_gap[i] = (2.0 + 3.0)/2.0 - 2.0
            elif (data_new['redshift_BEAGLE'][i] > 3.0) & (data_new['redshift_BEAGLE'][i] < 4.0):
                z_med_hp[i] = (3.0 + 4.0)/2.0
                z_med_hp_gap[i] = (3.0 + 4.0)/2.0 - 3.0
            elif (data_new['redshift_BEAGLE'][i] > 4.0) & (data_new['redshift_BEAGLE'][i] < 5.0):
                z_med_hp[i] = (4.0 + 5.0)/2.0
                z_med_hp_gap[i] = (4.0 + 5.0)/2.0 - 4.0
            elif (data_new['redshift_BEAGLE'][i] > 5.0) & (data_new['redshift_BEAGLE'][i] < 6.0):
                z_med_hp[i] = (5.0 + 6.0)/2.0
                z_med_hp_gap[i] = (5.0 + 6.0)/2.0 - 5.0


        x_hp = np.array([])
        y_hp = np.array([])
        z_hp = np.array([])

        n_hp = 3000 # number of samples to take from GMM in total

        fig, ax = plt.subplots(1, 1, figsize=(15, 15))

        sfr_std = []

        count_blue = 1

        for i in range(len(data_new['id_AD'])):
            if i <100000:

                x_temp = np.array([])
                y_temp = np.array([])
                z_temp = np.array([])

                for G in range(3):

                    mean = np.array([data_new['x_GMM_3d'][i,G],data_new['y_GMM_3d'][i,G],data_new['z_GMM_3d'][i,G]])
                    cov = np.array([[np.power(data_new['xsig_GMM_3d'][i,G],2), data_new['xycov_GMM_3d'][i,G], data_new['xzcov_GMM_3d'][i,G]],[data_new['xycov_GMM_3d'][i,G], np.power(data_new['ysig_GMM_3d'][i,G],2), data_new['yzcov_GMM_3d'][i,G]],[data_new['xzcov_GMM_3d'][i,G], data_new['yzcov_GMM_3d'][i,G], np.power(data_new['zsig_GMM_3d'][i,G],2)]])

                    xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*data_new['amp_GMM_3d'][i,G]))


                    z_idx = (abs(xyz[:,2] - z_med_hp[i]) < z_med_hp_gap[i])

                    if True: # only keep samples within z bin
                        x_temp = np.concatenate((x_temp,xyz[:,0][z_idx]))
                        y_temp = np.concatenate((y_temp,xyz[:,1][z_idx]))
                        z_temp = np.concatenate((z_temp,xyz[:,2][z_idx]))

                        x_hp = np.concatenate((x_hp,xyz[:,0][z_idx]))
                        y_hp = np.concatenate((y_hp,xyz[:,1][z_idx]))
                        z_hp = np.concatenate((z_hp,xyz[:,2][z_idx]))

                    else:
                        x_hp = np.concatenate((x_hp,xyz[:,0]))
                        y_hp = np.concatenate((y_hp,xyz[:,1]))
                        z_hp = np.concatenate((z_hp,xyz[:,2]))

                        x_temp = np.concatenate((x_temp,xyz[:,0]))
                        y_temp = np.concatenate((y_temp,xyz[:,1]))
                        z_temp = np.concatenate((z_temp,xyz[:,2]))

                # ax.scatter(x_temp, y_temp, s=0.5)
                # confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='firebrick', linewidth=1)

                sfr_std.append(np.std(y_temp))

                sfr_limit = 2.5
                mass_limit = 9.5

                if np.std(y_temp) > 2:
                    confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='blue', linewidth=1)
                    print('MMMMMMM', count_blue)
                    count_blue += 1

                elif (data_new['x_GMM_3d'][i,0] > mass_limit and data_new['y_GMM_3d'][i,0] > sfr_limit) or (data_new['x_GMM_3d'][i,1] > mass_limit and data_new['y_GMM_3d'][i,1] > sfr_limit) or (data_new['x_GMM_3d'][i,2] > mass_limit and data_new['y_GMM_3d'][i,2] > sfr_limit):

                    confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='green', linewidth=1)
                    print(data_new['field_AD'][i], data_new['id_AD'][i], data_new['id_BEAGLE'][i])
                    print(data_new['x_GMM_3d'][i], data_new['y_GMM_3d'][i])
                else:
                    confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='firebrick', linewidth=1)




        ax.set_xlim(6, 11)
        # ax.set_ylim(-3, 4)
        ax.set_ylim(-20, 4)

        # ax.set_xlim(8, 11)
        # ax.set_ylim(0, 4)

        # ax.set_xlim(6.5, 10.5)
        # ax.set_ylim(-5.5, 2.5)

        plt.show()


        plt.hist(sfr_std, bins=30)
        plt.show()

        sfr_std = np.array(sfr_std)
        print(len(sfr_std))
        print(len(sfr_std[sfr_std>2]))


        plt.scatter(data_new['mass_BEAGLE_stellar'][sfr_std<2], data_new['sfr_BEAGLE_instant'][sfr_std<2], label=len(data_new['mass_BEAGLE_stellar'][sfr_std<2]))
        plt.scatter(data_temp['mass_BEAGLE_stellar'][idx_sigma], data_temp['sfr_BEAGLE_instant'][idx_sigma], label=len(data_temp['mass_BEAGLE_stellar'][idx_sigma]))
        plt.scatter(data_new['mass_BEAGLE_stellar'][sfr_std>2], data_new['sfr_BEAGLE_instant'][sfr_std>2], label=len(data_new['mass_BEAGLE_stellar'][sfr_std>2]))
        plt.title('scenario{} {} {}'.format(scenarioA, fields[s], z_bins[s]))
        plt.xlabel('stellar mass')
        plt.ylabel('instant sfr')
        plt.legend()
        plt.show()


        data_new_sfr = copy.deepcopy(data_new)
        for key in data_new_sfr.keys():
            data_new_sfr[key] = data_new_sfr[key][sfr_std<2]



        #%%
        # =============================================================================
        # remove 'bad' gmm peaks
        # =============================================================================

        plt.scatter(data_new_sfr['y_GMM_3d'][:,0], data_new_sfr['ysig_GMM_3d'][:,0])
        plt.scatter(data_new_sfr['y_GMM_3d'][:,1], data_new_sfr['ysig_GMM_3d'][:,1])
        plt.scatter(data_new_sfr['y_GMM_3d'][:,2], data_new_sfr['ysig_GMM_3d'][:,2])
        x_tmp = np.array([-40.0, 5.0])
        plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
        plt.show()


        idx_sort = np.argsort(data_new_sfr['amp_GMM_3d'].flatten())

        #MASS
        plt.title('mass, each GMM, sigma vs mean')
        plt.scatter(data_new_sfr['x_GMM_3d'].flatten()[idx_sort], data_new_sfr['xsig_GMM_3d'].flatten()[idx_sort], c=data_new_sfr['amp_GMM_3d'].flatten()[idx_sort])
        # x_tmp = np.array([-40.0, 5.0])
        # plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
        plt.colorbar()
        plt.show()

        cmap = mpl.cm.get_cmap("tab20", 20)
        #SFR
        plt.title('SFR, each GMM, sigma vs mean')
        plt.scatter(data_new_sfr['y_GMM_3d'].flatten()[idx_sort], data_new_sfr['ysig_GMM_3d'].flatten()[idx_sort], c=data_new_sfr['amp_GMM_3d'].flatten()[idx_sort], cmap=cmap)
        x_tmp = np.array([-40.0, 5.0])
        plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
        # plt.xlim(-5,5)
        # plt.ylim(0,6)
        plt.colorbar()
        plt.show()

        #Z
        plt.title('redshift, each GMM, sigma vs mean')
        plt.scatter(data_new_sfr['z_GMM_3d'].flatten()[idx_sort], data_new_sfr['zsig_GMM_3d'].flatten()[idx_sort], c=data_new_sfr['amp_GMM_3d'].flatten()[idx_sort])
        # x_tmp = np.array([-40.0, 5.0])
        # plt.plot(x_tmp, (-9.0/26.0)*x_tmp - (9.0/13.0))
        # plt.xlim(0.5,3)
        # plt.ylim(0,0.5)
        plt.colorbar()
        plt.show()



        #%% for objects where I remove a peak, show full GMM weighting + mean SFR + sig SFR

        print('bad gmm peaks')
        for i in range(len(data_new_sfr['amp_GMM_3d'])):
            if data_new_sfr['ysig_GMM_3d'][i,0] < (-9.0/26.0)*data_new_sfr['y_GMM_3d'][i,0] - (9.0/13.0) or data_new_sfr['ysig_GMM_3d'][i,1] < (-9.0/26.0)*data_new_sfr['y_GMM_3d'][i,1] - (9.0/13.0) or data_new_sfr['ysig_GMM_3d'][i,2] < (-9.0/26.0)*data_new_sfr['y_GMM_3d'][i,2] - (9.0/13.0):
                print(data_new_sfr['amp_GMM_3d'][i])
                print(data_new_sfr['y_GMM_3d'][i])
                print(data_new_sfr['ysig_GMM_3d'][i])
                print('')


        #%% actually remove bad peaks
        data_new_sfr['amp_GMM_3d'] = np.where(data_new_sfr['ysig_GMM_3d'] < (-9.0/26.0)*data_new_sfr['y_GMM_3d'] - (9.0/13.0), 0.0, data_new_sfr['amp_GMM_3d']) # set bad gmms to p==0
        data_new_sfr['amp_GMM_3d'] = data_new_sfr['amp_GMM_3d'] / np.array([np.sum(data_new_sfr['amp_GMM_3d'], axis=1)]*3).T # renormalise - GENIUS

        print('FINAL: {}'.format(len(data_new_sfr[key])))

#        pickle.dump(data_new_sfr, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, fields[s], z_bins[s]),'w'))
#        outputTable = Table(data_new_sfr)
#        outputTable.write('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.fits'.format(scenarioA, fields[s], z_bins[s]), overwrite=True)
#
#




