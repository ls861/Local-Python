#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 18:09:05 2020

@author: lester
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits

#mpl.rcdefaults()
#plt.style.use('lester_style')
mpl.rc('figure', figsize=(4,4))


param = 'DE'
revision = '108'
summary = '004'

# _a is astrodeep
# _b is beagle fitted

id_b = np.array([3, 5, 14, 17, 24, 25, 26, 27, 28, 30, 33, 35, 37, 38, 43, 44, 46, 52, 63, 64, 69, 74, 75, 77, 88, 100, 102, 104, 108, 109, 113, 117, 121, 142, 144]) - 1  # all remaining after run 1

id_b = np.array([3, 14, 17, 24, 26, 28, 30, 33, 35, 37, 38, 44, 46, 52, 69, 74, 75, 77, 100, 102, 104, 108, 117, 121, 142, 144]) - 1  # completed during run 2

id_b = np.array([5, 25, 27, 43, 63, 64, 88, 109, 113]) - 1 # all remaining after run 2


# =============================================================================
# OUTPUT - get fitted BEAGLE parameters 
# =============================================================================

fileName = '/Users/lester/Documents/PhD/param_100/fit_{}_{}/pyp-beagle/data_{}/BEAGLE_summary_catalogue.fits'.format(revision, param, summary)
data_fits = fits.open(fileName)

id_b = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1

#id_b = np.array([3, 5, 14, 17, 24, 25, 26, 27, 28, 30, 33, 35, 37, 38, 43, 44, 46, 52, 63, 64, 69, 74, 75, 77, 88, 100, 102, 104, 108, 109, 113, 117, 121, 142, 144]) - 1  # all remaining after run 1
#id_b = np.array([3, 5, 14, 17, 24, 25, 26, 27, 28, 30, 33, 35, 37, 38, 43, 44, 46, 52, 63, 69, 74, 75, 77, 88, 100, 102, 104, 108, 109, 113, 117, 121, 142, 144]) - 1  # all remaining after run 1 minus 64
#id_b = np.array([3, 14, 17, 24, 26, 28, 30, 33, 35, 37, 38, 44, 46, 52, 69, 74, 75, 77, 100, 102, 104, 108, 117, 121, 142, 144]) - 1  # completed during run 2
#id_b = np.array([5, 25, 27, 43, 63, 64, 88, 109, 113]) - 1 # all remaining after run 2

sfr_b = np.log10(data_fits['STAR FORMATION'].data['SFR_mean'])
sfr_68_b = np.log10(data_fits['STAR FORMATION'].data['SFR_68.00'])

ssfr_b = data_fits['STAR FORMATION'].data['sSFR_mean']
ssfr_68_b = data_fits['STAR FORMATION'].data['sSFR_68.00']

mass_b = data_fits['POSTERIOR PDF'].data['mass_mean']
msa_b = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_mean']
tauV_eff_b = data_fits['POSTERIOR PDF'].data['tauv_eff_mean']
metallicity_b = data_fits['POSTERIOR PDF'].data['metallicity_mean']
tau_b = 10**data_fits['POSTERIOR PDF'].data['tau_mean']

mass_68_b = data_fits['POSTERIOR PDF'].data['mass_68.00']
msa_68_b = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_68.00']
tauV_eff_68_b = data_fits['POSTERIOR PDF'].data['tauv_eff_68.00']
metallicity_68_b = data_fits['POSTERIOR PDF'].data['metallicity_68.00']
tau_68_b = 10**data_fits['POSTERIOR PDF'].data['tau_68.00']

if revision in ['999']:
    nebular_logU_b = np.full(len(id_b), -2.5)
    nebular_logU_68_b = np.full((len(id_b),2), -2.5)
else:
    nebular_logU_b = data_fits['POSTERIOR PDF'].data['nebular_logu_mean']        
    nebular_logU_68_b = data_fits['POSTERIOR PDF'].data['nebular_logu_68.00']      

if revision in ['108']:
    nebular_xi_b = np.full(len(id_b), 0.3)
    nebular_xi_68_b = np.full((len(id_b),2), 0.3)
else:
    nebular_xi_b = data_fits['POSTERIOR PDF'].data['nebular_xi_mean']      
    nebular_xi_68_b = data_fits['POSTERIOR PDF'].data['nebular_xi_68.00']
    
redshift_b = data_fits['POSTERIOR PDF'].data['redshift_mean']
redshift_68_b = data_fits['POSTERIOR PDF'].data['redshift_68.00']

H160_b = 10**((23.9 - data_fits['APPARENT MAGNITUDES'].data['HST_WFC3_IR_F160W_APP_mean']) / 2.5) #uJy
H160_68_b = 10**((23.9 - data_fits['APPARENT MAGNITUDES'].data['HST_WFC3_IR_F160W_APP_68.00']) / 2.5) #uJy


data_fits.close()

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/ASTRODEEP/{}_{}_chi2_{}.fits'.format(revision, param, summary)
data_fits = fits.open(fileName)
chi2 = data_fits[1].data['chi2']
data_fits.close()


# =============================================================================
# INPUT from ASTRODEEP subset, ordered by redshift
# =============================================================================

# dtype=[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]

D1_new = np.load('D1_new.npy')

redshift_a = D1_new['ZBEST'][id_b]
H160_a = D1_new['b_H160'][id_b]
H160_68_a = D1_new['b_errH160'][id_b]

#print(H160_a) # few small -ves, mainly between 1e-4 and 1e2, uJy
#print(H160_68_a) # between 1e-3 and 1e-1 with some -99s, uJy, absolute error
#print(H160_b)# uJy
#print(H160_68_b) # uJy, values

redshift_68_a = D1_new['ZBEST_SIQR'][id_b] 

# dtype=[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]

# TURNS OUT THIS WASN'T NEEDED BECAUSE ZBEST_SIQR IS IN D1_new anyway, amended above
#D = np.load('astrodeep_rawfile.npy')
#
#redshift_68_a = D[D['field']==1]
#redshift_68_a = redshift_68_a[np.isin(redshift_68_a['ID'], D1_new['ID'])]
#redshift_68_a = np.sort(redshift_68_a, order='ZBEST')['ZBEST_SIQR'][id_b]
#
#plt.scatter(redshift_68_a, test_redshift_68_a)
#plt.show()

# =============================================================================
# Timing Tests
# =============================================================================

runtime = np.empty(len(id_b))

for i, ID in enumerate(id_b + 1):

    data_fits = fits.open('/Users/lester/Documents/PhD/param_100/fit_{}_{}/{}_BEAGLE.fits.gz'.format(revision, param, ID))
    #print(data_fits.info())

    runtime[i] = data_fits[0].header['HIERARCH RUN TIME'] / 3600
    data_fits.close()


# =============================================================================
# ERRORS FOR PLOTS and renaming variables
# =============================================================================

redshift_68_a = redshift_68_a
redshift_68_b = [redshift_b - redshift_68_b[:, 0], redshift_68_b[:, 1] - redshift_b]
H160_68_a = H160_68_a + (0.02*abs(H160_a)) # assume abs is necessary here, due to -ve fluxes, also 8x H160_68_a are -99!!
H160_68_b = [H160_68_b[:, 0]-H160_b, H160_b-H160_68_b[:, 1]]


# =============================================================================
# PLOTS
# =============================================================================

plt.title('runtime')
plt.hist(runtime, bins=20)
plt.show()

plt.title('astrodeep Z vs runtime')
plt.scatter(runtime, redshift_a, marker='x', c=chi2, norm=colors.LogNorm())
plt.errorbar(runtime, redshift_a, yerr=redshift_68_a, linestyle="None", elinewidth=0.5, color='k', zorder=0)
plt.colorbar(label='chi2')
plt.show()

plt.title('beagle Z vs runtime')
plt.scatter(runtime, redshift_b, marker='x', c=chi2, norm=colors.LogNorm())
plt.errorbar(runtime, redshift_b, yerr=redshift_68_b, linestyle="None", elinewidth=0.5, color='k', zorder=0)
plt.colorbar(label='chi2')
plt.show()

plt.title('beagle Z vs astrodeep Z')
plt.scatter(redshift_a, redshift_b, marker='x', c=chi2, norm=colors.LogNorm(), zorder=2)
plt.errorbar(redshift_a, redshift_b, xerr=redshift_68_a, yerr=redshift_68_b, linestyle="None", elinewidth=0.5, color='k', zorder=1)
plt.plot((0, 15), (0, 15), color='r', zorder=0)
plt.xlim(0, 15)
plt.ylim(0, 15)
plt.colorbar(label='chi2')
plt.show()

plt.title('log astrodeep H160 vs runtime')
plt.scatter(runtime[H160_a>=0], np.log10(abs(H160_a[H160_a>=0])+1e-6), marker='x')
plt.scatter(runtime[H160_a<0], np.log10(abs(H160_a[H160_a<0])+1e-6), color='r', marker='x')
plt.errorbar(runtime, np.log10(abs(H160_a)+1e-6), yerr=np.log10(H160_68_a), linestyle="None", elinewidth=0.5, color='k', zorder=0)
plt.show()

# log 
plt.title('log beagle H160 vs runtime')
plt.scatter(runtime, np.log10(H160_b), marker='x', c=chi2, norm=colors.LogNorm())
plt.errorbar(runtime, np.log10(H160_b), yerr=np.log10(H160_68_b), linestyle="None", elinewidth=0.5, color='k', zorder=0)
plt.colorbar(label='chi2')
plt.show()

plt.title('log beagle H160 vs log astrodeep H160')
plt.scatter(np.log10(abs(H160_a)+1e-6), np.log10(H160_b), marker='x', c=chi2, norm=colors.LogNorm())
plt.errorbar(np.log10(abs(H160_a)+1e-6), np.log10(H160_b), xerr=np.log10(H160_68_a), yerr=np.log10(H160_68_b), linestyle="None", elinewidth=0.5, color='k', zorder=0)
plt.plot((-9, 3), (-9, 3))
plt.xlim(-9, 3)
plt.ylim(-9, 3)
plt.colorbar(label='chi2')
plt.show()

plt.title('log astrodeep H160 vs astrodeep Z')
plt.scatter(redshift_a, np.log10(abs(H160_a)+1e-6), marker='x', c=runtime)
plt.errorbar(redshift_a, np.log10(abs(H160_a)+1e-6), xerr=redshift_68_a, yerr=np.log10(H160_68_a), linestyle="None", elinewidth=0.5, color='k', zorder=0)
#plt.plot((-9, 3), (-9, 3))
plt.xlim(0, 11)
plt.ylim(-7, 3)
plt.colorbar(label='runtime (hours)')
plt.show()



print(str(len(id_b)) + ' / ' + str(len(D1_new)))

mpl.rcdefaults()

#np.save('DE_108_004_runtime', runtime)
#np.save('DE_108_004_log_H160_a', np.log10(abs(H160_a)+1e-6))


