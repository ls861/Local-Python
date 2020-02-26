
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


# =============================================================================
# get astrodeep redshifts
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/ascii_to_fits/astrodeep_catalogue_001.fits'
data_fits = fits.open(fileName)
id_AD = data_fits[1].data['ID']
z_AD = data_fits[1].data['ZBEST']
H160 = data_fits[1].data['b_H160']
errH160 = data_fits[1].data['b_errH160']
data_fits.close()

# =============================================================================
# get beagle redshifts
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_001/pyp-beagle/data/BEAGLE_summary_catalogue.fits'
data_fits = fits.open(fileName)
id_B = np.asarray(data_fits[1].data['ID'], dtype=int)
z_B = data_fits[1].data['redshift_median']
z_Berr = data_fits[1].data['redshift_68.00']
#z_B = data_fits[1].data['redshift_mean']
data_fits.close()

z_Berr_m = z_B - z_Berr[:, 0]
z_Berr_p = z_Berr[:, 1] - z_B

# =============================================================================
# Resizing Arrays
# =============================================================================

z_AD = z_AD[id_B-1]
H160 = H160[id_B-1]
errH160 = errH160[id_B-1]
SN = abs(H160/errH160)

ind = (SN > 5) & (SN < 20)

z_B = z_B[ind]
z_Berr_m = z_Berr_m[ind]
z_Berr_p = z_Berr_p[ind]
z_AD = z_AD[ind]
H160 = H160[ind]
errH160 = errH160[ind]
SN = SN[ind]


# =============================================================================
# PLOT
# =============================================================================

pt = 20

#fig, ax = plt.subplots(figsize=np.array([6.4, 4.8]))
fig, ax = plt.subplots(figsize=np.array([17.6, 4.8]))
#fig.suptitle('BEAGLE vs Astrodeep redshifts', size=20, ha='center')
ax.set_title('BEAGLE vs Astrodeep redshifts', size=pt, ha='center')

h = ax.scatter(z_AD, z_B, c=SN, cmap='plasma', zorder=1)
ax.errorbar(z_AD, z_B, yerr=[z_Berr_m, z_Berr_p], linestyle="None", elinewidth=1, color='k', zorder=0)

x_lin = np.linspace(0, 20, 2)
ax.plot(x_lin, x_lin)

ax.set_xlim(0, 11)
ax.set_ylim(0, 13)

ax.set_xlabel(r'Astrodeep Redshift', size=pt)
ax.set_ylabel(r'BEAGLE Redshift', size=pt)

cbar = plt.colorbar(h)
cbar.set_label(label='SNR', size=pt)
cbar.ax.tick_params(labelsize=pt)
plt.show()

fig, ax = plt.subplots(figsize=(10, 5))
ax.hist(SN, bins=50)


