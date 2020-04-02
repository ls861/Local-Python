
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

distance = abs(z_AD-z_B) < 0.02
print(len(distance))

filter_label = np.array(['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2'])

filters = ['HST_ACS_WFC_F435W_APP', 'HST_ACS_WFC_F606W_APP', 'HST_ACS_WFC_F814W_APP', 'HST_WFC3_IR_F105W_APP', 'HST_WFC3_IR_F125W_APP', 'HST_WFC3_IR_F140W_APP', 'HST_WFC3_IR_F160W_APP', 'Paranal_HAWKI_Ks_APP', 'Spitzer_IRAC_I1_APP', 'Spitzer_IRAC_I2_APP']

filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])

filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])


for i, j in enumerate(id_B[distance]): # j is galaxy ID up to 289 etc (AD and B id)
    
    z = z_B[distance][i]
    
    fileName = '/Users/lester/Documents/PhD/param_001/astrodeep_004/{}_BEAGLE.fits'.format(j)
    data_fits = fits.open(fileName)
    
    # get the best fit BEAGLE SED
    chi2 = data_fits['POSTERIOR PDF'].data['chi_square']    
    wl = data_fits['FULL SED WL'].data[0][0]
    flux = data_fits['FULL SED'].data[np.argmin(chi2)]
    
    # get the best fit BEAGLE fluxes
    m = []
    for k in range(len(filters)):
        m.append(np.array(data_fits['APPARENT MAGNITUDES'].data[filters[k]][np.argmin(chi2)]))
    m = np.array(m)
    
    f = []
    for k in range(len(filters)):
        f.append(10**( (23.9 - m[k]) / 2.5 ))
    f = np.array(f)
    
    data_fits.close()
    
    
    # get the input astrodeep fluxes (microJy)
    fileName_AD = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/ascii_to_fits/astrodeep_catalogue_002.fits'
    data_fits = fits.open(fileName_AD)

    f_AD = []
    for k in range(len(filter_label)):
        f_AD.append(data_fits[1].data[filter_label[k]][j])
    f_AD = np.array(f_AD)
    data_fits.close()

    

    plt.scatter(filter_fwhm_centre, f/1E19, marker='x')
    plt.scatter(filter_fwhm_centre, f_AD/1E20, marker='x')
    plt.plot(wl * (1+z), flux / (1+z))
    plt.xlim(0000, 55000)
#    plt.ylim(1E-22, 1E-20)
    plt.yscale('log')
    plt.show()
    print(i, j, z, int(1216*(1+z)), int(4000*(1+z)))
    













