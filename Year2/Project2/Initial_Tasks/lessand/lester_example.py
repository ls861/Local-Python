import legac_utils as legac
import matplotlib.pyplot as plt

import os
import re
import warnings

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning

from astropy import constants, table, units, wcs

from time import perf_counter as clock

import ppxf
from ppxf import ppxf_util, miles_util

myspec = legac.spec1d.from_filename('legac_M32_v3.11_spec1d_54311.fits')
myspec.ppxf_fit(plot=True, clean=False) # clean=True is better but takes a long time.



wave_overlap=(3540., 7409.)*units.Angstrom
wave_overlap = (myspec.wave>wave_overlap[0]) & (myspec.wave<wave_overlap[1])



myspec.lspec  = myspec.spec[wave_overlap].value  # Dimensionless
myspec.lwave  = myspec.wave[wave_overlap].value  # Dimensionless
myspec.lnoise = myspec.noise[wave_overlap].value # Dimensionless
mask = myspec.mask[wave_overlap]
myspec.spec_norm = np.nanmedian(myspec.lspec[~mask])
myspec.lspec /= myspec.spec_norm
myspec.lnoise /= myspec.spec_norm
myspec.spec_norm *= myspec.spec.unit # Save units of spec_norm.

#%%
plt.scatter(myspec.wave[~wave_overlap], np.full(len(myspec.wave[~wave_overlap]),1), s=1, marker='.')
plt.scatter(myspec.lwave, myspec.lspec, s=1, marker='.')
plt.scatter(myspec.lwave[np.isnan(myspec.lspec)], np.full(len(myspec.lwave[np.isnan(myspec.lspec)]),0), s=1, marker='.')
# plt.scatter(myspec.lwave, myspec.lnoise, s=1, marker='.')
plt.show()



# plt.scatter(myspec.wave[~wave_overlap], np.full(len(myspec.wave[~wave_overlap]),1), s=1, marker='.')
plt.scatter(myspec.lwave, myspec.lspec, s=1, marker='.')
plt.scatter(myspec.lwave[np.isnan(myspec.lspec)], np.full(len(myspec.lwave[np.isnan(myspec.lspec)]),0), s=1, marker='.')
# plt.scatter(myspec.lwave, myspec.lnoise, s=1, marker='.')
plt.show()


# Log-rebin the wavelength, spectrum and noise vectors. pPXF works in
# log-wave space (i.e. velocity space).
c = constants.c.to('km/s').value

_, _, myspec.velscale = ppxf_util.log_rebin(
    myspec.lwave[[0, -1]], myspec.lspec) # Pixel scale in km/s

myspec.lspec, loglam, _ = ppxf_util.log_rebin(
    myspec.lwave[[0, -1]], myspec.lspec, velscale=myspec.velscale)

myspec.lnoise, _, _ = ppxf_util.log_rebin(
    myspec.lwave[[0, -1]], myspec.lnoise**2, velscale=myspec.velscale)

myspec.lnoise = np.sqrt(myspec.lnoise)


# log-rebinning of the mask is conservative: mask any pixel which had
# non-zero contribution from a masked pixel.
mask, _, _ = ppxf_util.log_rebin(
    myspec.lwave[[0, -1]], mask.astype(float), velscale=myspec.velscale)
mask = np.array([True if m>0 else False for m in mask])
mask = mask | (~np.isfinite(myspec.lspec*myspec.lnoise)) | (myspec.lnoise<=0.)

myspec.lnoise[mask] = 1.e10 # Large value. These pixels are not used (masked)
myspec.lspec[mask] = 0. # Again, these pixels are not fit.
goodpixels = np.arange(len(myspec.lspec), dtype=int)
goodpixels = goodpixels[~mask]
myspec.lwave = loglam


# plt.scatter(myspec.wave[~wave_overlap], np.full(len(myspec.wave[~wave_overlap]),1), s=1, marker='.')
plt.scatter(myspec.lwave, myspec.lspec, s=1, marker='.')
# plt.scatter(myspec.lwave[np.isnan(myspec.lspec)], np.full(len(myspec.lwave[np.isnan(myspec.lspec)]),0), s=1, marker='.')
# plt.scatter(myspec.lwave, myspec.lnoise, s=1, marker='.')
plt.show()







#%%

# LEGA-C assumes uniform spectral resolution equal to 0.9 Angstrom
# (sigma of a Gaussian line-spread function).
FWHM_gal = 0.9 * (2. * np.sqrt(2. * np.log(2.))) # AA.

#------------------- Setup templates -----------------------

# The templates are not normalized.
ppxf_dir = os.path.dirname(os.path.realpath(ppxf.__file__))
pathname = os.path.join(ppxf_dir, 'miles_models/Mun1.30*.fits')
templlib = miles_util.miles(pathname, myspec.velscale, FWHM_gal)

# The stellar templates are reshaped below into a 2-dim array with each
# spectrum as a column, however we save the original array dimensions,
# which are needed to specify the regularization dimensions
#
stars_templates = templlib.templates.reshape(templlib.templates.shape[0], -1)

# Construct a set of Gaussian emission line templates.
# The `emission_lines` function defines the most common lines, but additional
# lines can be included by editing the function in the file ppxf_ppxf_util.py.
myspec.gas_templates, gas_names, line_wave = ppxf_util.emission_lines(
    templlib.log_lam_temp, np.exp(myspec.lwave[[0, -1]]), FWHM_gal,
    tie_balmer=False, limit_doublets=False)

# Combines the stellar and gaseous templates into a single array.
# During the PPXF fit they will be assigned a different kinematic
# COMPONENT value
#
templates = np.column_stack([stars_templates, myspec.gas_templates])

#-----------------------------------------------------------

# The galaxy and the template spectra do not have the same starting wavelength.
# For this reason an extra velocity shift DV has to be applied to the template
# to fit the galaxy spectrum. We remove this artificial shift by using the
# keyword VSYST in the call to PPXF below, so that all velocities are
# measured with respect to DV. This assume the redshift is negligible.
# In the case of a high-redshift galaxy one should de-redshift its
# wavelength to the rest frame before using the line below as described
# in PPXF_EXAMPLE_KINEMATICS_SAURON and Sec.2.4 of Cappellari (2017)
#
dv = c*(templlib.log_lam_temp[0] - myspec.lwave[0])  # eq.(8) of Cappellari (2017)
start = [0., 100.]     # (km/s), starting guess for [V, sigma]

n_temps = stars_templates.shape[1]
n_forbidden = np.sum(["[" in a for a in gas_names])  # forbidden lines contain "[*]"
n_balmer = len(gas_names) - n_forbidden

# Assign component=0 to the stellar templates, component=1 to the Balmer
# gas emission lines templates and component=2 to the forbidden lines.
component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
gas_component = np.array(component) > 0  # gas_component=True for gas templates

# Fit (V, sig, h3, h4) moments=4 for the stars
# and (V, sig) moments=2 for the two gas kinematic components
moments = [2, 4, 4]

# Adopt the same starting value for the stars and the two gas components
start = [start, start, start]

# If the Balmer lines are tied one should allow for gas reddeining.
# The gas_reddening can be different from the stellar one, if both are fitted.
#gas_reddening = 0 if tie_balmer else None


print(len(myspec.lspec), len(myspec.lnoise), len(np.trim_zeros(myspec.lspec)), len(np.trim_zeros(myspec.lnoise)))



#%%

degree=-1
mdegree=10
clean=False



#%%
t = clock()
myspec.pp = ppxf.ppxf.ppxf(
    templates, myspec.lspec, myspec.lnoise, myspec.velscale, start,
    moments=moments, degree=degree, mdegree=mdegree,
    vsyst=dv, clean=clean,
    lam=np.exp(myspec.lwave), goodpixels=goodpixels,
    component=component, gas_component=gas_component,
    gas_names=gas_names, plot=True)#, gas_reddening=gas_reddening)



#%%

t1 = 0
t2 = 940



myspec.pp = ppxf.ppxf.ppxf(
    templates, myspec.lspec[t1:-t2], myspec.lnoise[t1:-t2], myspec.velscale, start,
    moments=moments, degree=degree, mdegree=mdegree,
    vsyst=dv, clean=clean,
    lam=np.exp(myspec.lwave[t1:-t2]), goodpixels=goodpixels[(goodpixels<len(myspec.lspec)-t2)],
    component=component, gas_component=gas_component,
    gas_names=gas_names, plot=True)#, gas_reddening=gas_reddening)


print(goodpixels)

print(len(goodpixels), len(goodpixels[goodpixels<len(myspec.lspec)-t2]))


#%%
# When the two Delta Chi^2 below are the same, the solution
# is the smoothest consistent with the observed spectrum.
#
print('Desired Delta Chi^2: %#.4g' % np.sqrt(2*myspec.lspec.size))
print('Current Delta Chi^2: %#.4g' % ((myspec.pp.chi2 - 1)*myspec.lspec.size))
print('Elapsed time in PPXF: %.2f s' % (clock() - t))

weights = myspec.pp.weights[~gas_component]  # Exclude weights of the gas templates

# pPXF gas fluxes are in normalised units * spectral pixel. We want
# Physical units.
print('Gas fluxes in relevant units are:\n')
print('Name flux [flux density * Angstrom]\n')
for name, flux, templ in zip(
    myspec.pp.gas_names, myspec.pp.gas_flux, myspec.gas_templates.T):

    # Get pixel size in Angstrom around the wavelength of the line.
    wavelength = np.argmax(templ)
    wavelength = np.exp(templlib.log_lam_temp[wavelength])
    wavelength = np.argmin(np.abs(myspec.pp.lam-wavelength))
    dwave = myspec.pp.lam[wavelength] - myspec.pp.lam[wavelength-1]
    dwave *= myspec.wave.unit

    print(name, flux*dwave*myspec.spec_norm)
    
















#%%
'''
myspec.continuum_subtraction()

fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True)

ax0.plot(myspec.wave, myspec.spec-myspec.gas_model, 'k-', alpha=0.5)
ax0.plot(myspec.wave, myspec.cont, 'r-', alpha=0.5)
ax1.plot(myspec.wave, myspec.contsub_spec, 'k-', alpha=0.5)
ax1.plot(myspec.wave, myspec.gas_model, 'r-', alpha=0.5)

ax1.scatter(6563,50)
ax1.scatter(4861,50)
ax1.scatter(4340,50)
ax1.scatter(4102,50)
ax1.scatter(3726,30)
ax1.scatter(3729,30)
ax1.scatter(5007,30)

ax1.scatter()



plt.show()


plt.plot(myspec.pp.mpoly)
plt.show()'''