
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import cosmolopy.distance as cd
import cosmolopy.constants as cc

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_004/astrodeep_002/mock_catalogue_004_002.fits'
fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_004/astrodeep_003/mock_catalogue_004_003.fits'


# =============================================================================
# get SFH parameters
# =============================================================================

data_fits = fits.open(fileName)

redshift = data_fits['STAR FORMATION'].data['redshift']
tau = data_fits['STAR FORMATION BINS'].data['bin_tau']
msa = data_fits['STAR FORMATION'].data['max_stellar_age']

lookback_age = data_fits['STAR FORMATION HISTORIES'].data['lookback_age']
SFR = data_fits['STAR FORMATION HISTORIES'].data['SFR']

data_fits.close()

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)
age_z    = cd.age(redshift, **cosmo)/cc.yr_s

beta = np.array([1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6., 11., 1., 6.])

alpha = np.array([1., 1., 1., 6., 6., 6., 11., 11., 11., 1., 1., 1., 6., 6., 6., 11., 11., 11., 1., 1., 1., 6., 6., 6., 11., 11., 11., 1., 1., 1., 6., 6., 6., 11., 11., 11., 1., 1., 1., 6., 6., 6., 11., 11., 11., 1., 1., 1., 6., 6.])

# =============================================================================
# get SED parameters
# =============================================================================

data_fits = fits.open(fileName)
wl = data_fits['FULL SED WL'].data['wl'][0]
flux = data_fits['FULL SED'].data
data_fits.close()

fig, ax = plt.subplots(figsize=(10, 5))
for i in range(len(flux)):
    ax.plot(wl, flux[i])

ax.set_xlim(0, 4E4)
ax.set_ylim(1E-27, 1E-15)
#ax.set_xscale('log')
ax.set_yscale('log')    
plt.show()


# =============================================================================
# SFHs when max_stellar_age has been set
# =============================================================================

#i=0
#
#for i in range(len(SFR)):
##for i in [16, 20, 25, 30, 40, 45]:
#    x_arr = np.linspace(0, 1.5E10, 1001)
#    y_arr = (  (((x_arr - age_z[i] + msa[i])/tau[i])**alpha[i]) + (((x_arr - age_z[i] + msa[i])/tau[i])**-beta[i])  )**-1
#    
#    norm = max(age_z[i] - lookback_age[i])
#    ratio = max(y_arr) * SFR[i][0] /  ( (  (((norm - age_z[i] + msa[i])/tau[i])**alpha[i]) + (((norm - age_z[i] + msa[i])/tau[i])**-beta[i])  )**-1 )
#    
#    plt.xlim(0, 1.4E10)
#    plt.ylim(bottom=0)
#    plt.plot(x_arr, y_arr / max(y_arr), color='#1f77b4')
##    plt.show()
#    
#    plt.xlim(0, 1.4E10)
#    plt.plot(age_z[i] - lookback_age[i], SFR[i] / ratio, color='#ff7f0e')
##    plt.xlim(0, 1.4E10)
##    plt.ylim(0, 1E5)
#    plt.show()


# =============================================================================
# SFHs when max_stellar_age has NOT been set
# =============================================================================

i=0

for i in range(len(SFR)):
#for i in [16, 20, 25, 30, 40, 45]:
    x_arr = np.linspace(0, 1.5E10, 1001)
    
#    s = cd.age(30, **cosmo)/cc.yr_s
    s = 0
    
    
    y_arr = (  (((x_arr - s)/tau[i])**alpha[i]) + (((x_arr - s)/tau[i])**-beta[i])  )**-1
    
    norm = max(age_z[i] - lookback_age[i])
    ratio = max(y_arr) * SFR[i][0] /  ( (  (((norm - s)/tau[i])**alpha[i]) + (((norm - s)/tau[i])**-beta[i])  )**-1 )
    
    plt.xlim(0, 1.4E10)
    plt.ylim(bottom=0)
    plt.plot(x_arr, y_arr / max(y_arr), color='#1f77b4')
#    plt.show()
    
    plt.xlim(0, 1.4E10)
    plt.plot(age_z[i] - lookback_age[i], SFR[i] / ratio, color='#ff7f0e')
#    plt.xlim(0, 1.4E10)
#    plt.ylim(0, 1E5)
    plt.show()
    
    print(alpha[i], beta[i], tau[i], msa[i], age_z[i])

