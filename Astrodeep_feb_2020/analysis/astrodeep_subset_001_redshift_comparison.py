
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# =============================================================================
# get beagle redshifts (len 222)
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_001/astrodeep_004/pyp-beagle/data/BEAGLE_summary_catalogue.fits'
data_fits = fits.open(fileName)
id_B = np.asarray(data_fits[1].data['ID'], dtype=int)
z_B = data_fits[1].data['redshift_median']
z_Berr = data_fits[1].data['redshift_68.00']
#z_B = data_fits[1].data['redshift_mean']
data_fits.close()

z_Berr_m = z_B - z_Berr[:, 0]
z_Berr_p = z_Berr[:, 1] - z_B

# =============================================================================
# get astrodeep redshifts (len 289 converted to len 222 to match BEAGLE fits)
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/ascii_to_fits/astrodeep_catalogue_002.fits'
data_fits = fits.open(fileName)
id_AD = data_fits[1].data['ID'][id_B-1]
field_AD = data_fits[1].data['field'][id_B-1]
ido_AD = data_fits[1].data['ID_original'][id_B-1]
z_AD = data_fits[1].data['ZBEST'][id_B-1]
H160_AD = data_fits[1].data['b_H160'][id_B-1]
H160_ADerr = data_fits[1].data['b_errH160'][id_B-1]
#print(data_fits[1].header)
data_fits.close()

# =============================================================================
# get astrodeep redshift errors from original file (len 29290, 222)
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/astrodeep_rawfile.npy'
cat = np.load(fileName)
#print(cat.dtype.names)

z_ADerr = np.empty(len(id_B))

for i in range(len(id_B)):
    z_ADerr[i] = cat['ZBEST_SIQR'][(cat['ID'] == ido_AD[i]) & (cat['field'] == field_AD[i])]

z_ADerr_m = z_AD - z_ADerr
z_ADerr_p = z_ADerr - z_AD


# =============================================================================
# PLOT
# =============================================================================

SN = abs(H160_AD/H160_ADerr)
ind = (SN > 5) & (SN < 50)

pt = 20

#fig, ax = plt.subplots(figsize=np.array([6.4, 4.8]))
fig, ax = plt.subplots(figsize=np.array([17, 5]))
#fig.suptitle('BEAGLE vs Astrodeep redshifts', size=20, ha='center')
ax.set_title('BEAGLE vs Astrodeep redshifts', size=pt, ha='center')

h = ax.scatter(z_AD[ind], z_B[ind], c=SN[ind], cmap='plasma', zorder=1)
ax.errorbar(z_AD[ind], z_B[ind], xerr=[z_ADerr[ind], z_ADerr[ind]], yerr=[z_Berr_m[ind], z_Berr_p[ind]], linestyle="None", elinewidth=1, color='k', zorder=0)

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
ax.hist(SN[ind], bins=50)
plt.show()


# =============================================================================
# investigating the galaxies which agree/disagree the most
# =============================================================================

distance = abs(z_AD-z_B)<0.02
print(len(distance))

'''
Need to plot input photometric points plus errors
plus fitted photometric points plus errors
plus best fit SED
'''

_arr = np.empty(len(distance))

for i, j in enumerate(id_B[distance]): # j is galaxy ID up to 289 etc (AD and B id)
    print(i, j)
    fileName = '/Users/lester/Documents/param_001/astrodeep_004/{}_BEAGLE.fits'.format(j)
    data_fits = fits.open(fileName)
    chi2 = data_fits['POSTERIOR PDF'].data['chi_square']    
    wl = data_fits['FULL SED WL'].data[0][0]
    flux = data_fits['FULL SED'].data[np.argmin(chi2)]
    data_fits.close()
    
    z = z_B[distance][i]
#    plt.plot(wl, flux)
    plt.plot(wl * (1+z), flux / (1+z))
    plt.xlim(5000, 15000)
#    plt.ylim(0, 1E-20)
#    plt.yscale('log')
    plt.show()

    



# NEED TO SOMEHOW GET THE ACTUAL REDSHIFTS I'M TRYING TO SORT OUT HERE.......

print(z_B[distance])
print(z_AD[distance])

















