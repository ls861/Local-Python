
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_004/astrodeep_001/mock_catalogue_004_001.fits'

# =============================================================================
# get SFH parameters
# =============================================================================

data_fits = fits.open(fileName)
lookback_age = data_fits['STAR FORMATION HISTORIES'].data['lookback_age']
SFR = data_fits['STAR FORMATION HISTORIES'].data['SFR']
data_fits.close()

fig, ax = plt.subplots(figsize=(10, 5))
for i in range(len(SFR)):
    ax.plot(lookback_age[i], SFR[i]/max(SFR[i]))
ax.set_xlim(2E9, 0)
plt.show()

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

ax.set_ylim(1E-27, 1E-15)
ax.set_xscale('log')    
ax.set_yscale('log')    
plt.show()











