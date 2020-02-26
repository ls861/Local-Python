
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sfr_calc import sfr_calc

size = 15
fsize = 7

# =============================================================================
# get "real" parameters
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/fit/mock_catalogue_005_010.fits'
data_fits = fits.open(fileName)

mtot_r = np.log10(data_fits['GALAXY PROPERTIES'].data['m_tot'])                    # mtot is the value we selected as prior
#r_mstar = data_fits['GALAXY PROPERTIES'].data['m_star']
sfr_r = data_fits['STAR FORMATION'].data['SFR']

data_fits.close()

plt.figure(figsize=(fsize, fsize))
plt.title('Plot showing SFR vs Mass', size=size)
plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.scatter(mtot_r, np.log10(sfr_r))
plt.show()

# =============================================================================
# get BEAGLE parameters
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_006/astrodeep_002/pyp-beagle/data/BEAGLE_summary_catalogue.fits'

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_007/astrodeep_001/pyp-beagle/data/BEAGLE_summary_catalogue.fits'

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_008/astrodeep_001/pyp-beagle/data/BEAGLE_summary_catalogue.fits'

data_fits = fits.open(fileName)

id_b = np.asarray(data_fits[1].data['ID'], dtype=int) - 1
z_b = data_fits[1].data['redshift_mean']
z_med_b = data_fits[1].data['redshift_median']
mtot_b = data_fits[1].data['mass_mean']
msa_b = data_fits[1].data['max_stellar_age_mean']
tau_b = data_fits[1].data['tau_mean']

data_fits.close()


# =============================================================================
# plot Redshifts
# =============================================================================

plt.figure(figsize=(fsize, fsize))
plt.title('DE - Histogram of BEAGLE fitted redshifts', size=size)
plt.xlabel(r'Redshift', size=size)
plt.ylabel(r'Count', size=size)
plt.hist(z_b, bins=50, histtype=u'step', label='mean')
plt.hist(z_med_b, bins=50, histtype=u'step', label='median')
plt.legend()
plt.show()

# =============================================================================
# plot main sequence
# =============================================================================

sfr_b, err = sfr_calc('DE', mtot_b, msa_b, tau_b, 0, 0, 0, 0)

plt.figure(figsize=(fsize, fsize))
plt.title('DE - Plot showing BEAGLE fitted SFR vs Mass', size=size)
plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.scatter(mtot_r[id_b], np.log10(sfr_r[id_b]), label='input')
plt.scatter(mtot_b, np.log10(sfr_b), label='BEAGLE output')
plt.scatter(np.delete(mtot_r,id_b), np.log10(np.delete(sfr_r,id_b)), label='not fitted')
plt.legend()
plt.show()

plt.figure(figsize=(fsize, fsize))
plt.title('DE - Plot showing BEAGLE fitted SFR vs Mass', size=size)
plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
for i in range(len(sfr_b)):
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    plt.plot( (mtot_r[id_b][i],mtot_b[i]), (np.log10(sfr_r[id_b][i]),np.log10(sfr_b[i])), '-', color=color)
    plt.plot( mtot_r[id_b][i], np.log10(sfr_r[id_b][i]), 'o', color=color, markersize=3)
plt.show()   
    

# =============================================================================
# plot main sequence as a 2d histogram
# =============================================================================

# input values
plt.figure(figsize=(1.2*fsize, fsize))
plt.hist2d(mtot_r, np.log10(sfr_r), range=[[7.5, 11], [-1, 3.5]], bins=100)
plt.colorbar()
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()

massh = np.empty(0)
sfrh = np.empty(0)

#for i in range(len(id_b)):
for i in range(10):

    beagleData = fits.open('/Users/lester/Documents/param_006/astrodeep_002/{}_BEAGLE.fits'.format(id_b[i]+1))
    
    #needs float64 to provide precision needed for the random.choice weights
    temp_probs = np.float64(beagleData['POSTERIOR PDF'].data['probability'])
    temp_probs = temp_probs/np.sum(temp_probs)

    #here's the key line - take weighted samples from the multinest output!
    idx = np.random.choice(len(temp_probs), size=10000, p=temp_probs)
    massh = np.append(massh, np.log10(beagleData['GALAXY PROPERTIES'].data['M_tot'][idx]))
    sfrh = np.append(sfrh, np.log10(beagleData['STAR FORMATION'].data['sfr'][idx]))

plt.figure(figsize=(fsize, fsize))
plt.scatter(massh, sfrh)
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()

plt.figure(figsize=(1.2*fsize, fsize))
plt.hist2d(massh, sfrh, range=[[7.5, 11], [-1, 3.5]], bins=100)
plt.colorbar()
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()

plt.figure(figsize=(1.2*fsize, fsize))
plt.hist2d(massh, sfrh, range=[[7.5, 11], [-1, 3.5]], bins=100, normed=True)
plt.colorbar()
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()

import matplotlib.colors as mcolors

plt.figure(figsize=(1.2*fsize, fsize))
plt.hist2d(massh, sfrh, range=[[7.5, 11], [-1, 3.5]], bins=100, norm=mcolors.LogNorm())
plt.colorbar()
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()



plt.figure(figsize=(1.2*fsize, fsize))
plt.hist2d(massh, sfrh, range=[[7.5, 11], [-1, 3.5]], bins=100, norm=mcolors.PowerNorm(gamma))
plt.colorbar()
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.show()






'''
# =============================================================================
# Resizing Arrays
# =============================================================================

z_AD = z_AD[id_B-1]
H160 = H160[id_B-1]
errH160 = errH160[id_B-1]
SN = abs(H160/errH160)

ind = (SN > 0) & (SN < 20)

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


'''