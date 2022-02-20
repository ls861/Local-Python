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

__LEGAC_CATALOGUE__ = 'legac_team_DR3.fits'

def get_id_mask(filename):
    mask = re.findall(r'M[0-9]{1,3}', filename)[0]
    mask = int(mask[1:])
    legacid = re.findall(r'spec1d_[0-9]{4,6}', filename)[0]
    legacid = int(legacid[7:])

    return legacid, mask



def get_redshift(legacid, mask):

    catalogue = table.Table.read(__LEGAC_CATALOGUE__)

    match = np.where((catalogue['id']==legacid) & (catalogue['mask']==mask))[0]

    if len(match)==1:
        return catalogue[match]['z_spec'][0]
    else:
        raise ValueError(
            f'{legacid} M{mask} must match 1 file, but found {match}')

    

class spec1d():

    def __init__(self, wave, spec, noise, mask, redshift, legacid, obsmask):
        """
            wave : float 1-d array
                Rest-frame wavelength (air, Angstrom).
            spec : float 1-d array
                LEGA-C spectrum.
            noise : float 1-d array
                LEGA-C noise spectrum.
            mask : float 1-d array
                LEGA-C bad pixel map (True = bad pixel).
            redshift : float
                LEGA-C spectroscopic redshift.
            legacid : int
                LEGA-C ID (4-6 digits)
            obsmask: int
                LEGA-C observing mask (1-3 digits)
        """
        self.wave     = wave
        self.spec     = spec
        self.noise    = noise
        self.mask     = mask
        self.redshift = redshift
        self.legacid  = legacid
        self.obsmask  = obsmask
 


    def __repr__(self):
        return (f'<{self.__module__}.{self.__class__.__name__} instance '
            f'for LEGA-C galaxy {self.legacid} M{self.obsmask} at {hex(id(self))}>'
            )


    @classmethod
    def from_filename(cls, filename):

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', AstropyWarning)

            with fits.open(filename) as hdu:
                spec = hdu[0].data
                spec_wcs = wcs.WCS(hdu[0].header)
                wave = spec_wcs.all_pix2world(np.arange(
                    spec_wcs._naxis[0]), 0, 0)[0]
        
            wht_filename = filename.replace('spec1d', 'wht1d')
            with fits.open(wht_filename) as hdu:
                wht = hdu[0].data

        legacid, obsmask = get_id_mask(filename)
        z = get_redshift(legacid, obsmask)

        wave /= (1. + z)

        print(z)        

        mask = (wht==0.) | (~np.isfinite(wht))

        with np.errstate(divide='ignore'):
            noise = 1./np.sqrt(wht)

        wave  *= units.Angstrom
        spec  = units.Quantity(spec, '1e-18 erg/(s cm2 Angstrom)')
        noise = units.Quantity(noise, '1e-18 erg/(s cm2 Angstrom)')

        return cls(wave, spec, noise, mask, z, legacid, obsmask)


    def ppxf_fit(
        self, wave_overlap=(3540., 7409.)*units.Angstrom,
        degree=-1, mdegree=10, clean=True, **ppxf_kwargs
        ):
        """
        wave_overlap : 2-element array [Angstrom]
            wavelength limits of the template library (default is for MILES).
        degree : int, optional
            degree of additive Legendre polynomials for pPXF fit.
            -1=no polynomials.
        mdegree : int, optional
            degree of multiplicative Legendre polynomials for pPXF fit.
            0=no polynomials.
        clean : bool, optional
            if True, perform iterative sigma clipping in pPXF
        **ppxf_kwargs: optional
            any other valid pPXF keyword that is not included already in the
            pPXF call.
        """
        
        wave_overlap = (self.wave>wave_overlap[0]) & (self.wave<wave_overlap[1])
        self.lspec  = self.spec[wave_overlap].value  # Dimensionless
        self.lwave  = self.wave[wave_overlap].value  # Dimensionless
        
        print(self.wave, self.lwave)
        
        self.lnoise = self.noise[wave_overlap].value # Dimensionless
        mask = self.mask[wave_overlap]
        self.spec_norm = np.nanmedian(self.lspec[~mask])
        self.lspec /= self.spec_norm
        self.lnoise /= self.spec_norm
        self.spec_norm *= self.spec.unit # Save units of spec_norm.

        # Log-rebin the wavelength, spectrum and noise vectors. pPXF works in
        # log-wave space (i.e. velocity space).
        c = constants.c.to('km/s').value
        _, _, self.velscale = ppxf_util.log_rebin(
            self.lwave[[0, -1]], self.lspec) # Pixel scale in km/s
        self.lspec, loglam, _ = ppxf_util.log_rebin(
            self.lwave[[0, -1]], self.lspec, velscale=self.velscale)
        self.lnoise, _, _ = ppxf_util.log_rebin(
            self.lwave[[0, -1]], self.lnoise**2, velscale=self.velscale)
        self.lnoise = np.sqrt(self.lnoise)
        # log-rebinning of the mask is conservative: mask any pixel which had
        # non-zero contribution from a masked pixel.
        mask, _, _ = ppxf_util.log_rebin(
            self.lwave[[0, -1]], mask.astype(float), velscale=self.velscale)
        mask = np.array([True if m>0 else False for m in mask])
        mask = mask | (~np.isfinite(self.lspec*self.lnoise)) | (self.lnoise<=0.)

        self.lnoise[mask] = 1.e10 # Large value. These pixels are not used (masked)
        self.lspec[mask] = 0. # Again, these pixels are not fit.
        goodpixels = np.arange(len(self.lspec), dtype=int)
        goodpixels = goodpixels[~mask]
        self.lwave = loglam

        # LEGA-C assumes uniform spectral resolution equal to 0.9 Angstrom
        # (sigma of a Gaussian line-spread function).
        FWHM_gal = 0.9 * (2. * np.sqrt(2. * np.log(2.))) # AA.

        #------------------- Setup templates -----------------------
    
        # The templates are not normalized.
        ppxf_dir = os.path.dirname(os.path.realpath(ppxf.__file__))
        pathname = os.path.join(ppxf_dir, 'miles_models/Mun1.30*.fits')
        templlib = miles_util.miles(pathname, self.velscale, FWHM_gal)

        # The stellar templates are reshaped below into a 2-dim array with each
        # spectrum as a column, however we save the original array dimensions,
        # which are needed to specify the regularization dimensions
        #
        stars_templates = templlib.templates.reshape(templlib.templates.shape[0], -1)

        # Construct a set of Gaussian emission line templates.
        # The `emission_lines` function defines the most common lines, but additional
        # lines can be included by editing the function in the file ppxf_ppxf_util.py.
        self.gas_templates, gas_names, line_wave = ppxf_util.emission_lines(
            templlib.log_lam_temp, np.exp(self.lwave[[0, -1]]), FWHM_gal,
            tie_balmer=False, limit_doublets=False)

        # Combines the stellar and gaseous templates into a single array.
        # During the PPXF fit they will be assigned a different kinematic
        # COMPONENT value
        #
        templates = np.column_stack([stars_templates, self.gas_templates])

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
        dv = c*(templlib.log_lam_temp[0] - self.lwave[0])  # eq.(8) of Cappellari (2017)
        start = [0., 100.]     # (km/s), starting guess for [V, sigma]
    
        n_temps = stars_templates.shape[1]
        n_forbidden = np.sum(["[" in a for a in gas_names])  # forbidden lines contain "[*]"
        n_balmer = len(gas_names) - n_forbidden
    
        # Assign component=0 to the stellar templates, component=1 to the Balmer
        # gas emission lines templates and component=2 to the forbidden lines.
        component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
        print(component)
        gas_component = np.array(component) > 0  # gas_component=True for gas templates

        # Fit (V, sig, h3, h4) moments=4 for the stars
        # and (V, sig) moments=2 for the two gas kinematic components
        moments = [2, 4, 4]
    
        # Adopt the same starting value for the stars and the two gas components
        start = [start, start, start]
    
        # If the Balmer lines are tied one should allow for gas reddeining.
        # The gas_reddening can be different from the stellar one, if both are fitted.
        #gas_reddening = 0 if tie_balmer else None
        
        
        t = clock()
        self.pp = ppxf.ppxf.ppxf(
            templates, self.lspec, self.lnoise, self.velscale, start,
            moments=moments, degree=degree, mdegree=mdegree,
            vsyst=dv, clean=clean,
            lam=np.exp(self.lwave), goodpixels=goodpixels,
            component=component, gas_component=gas_component,
            gas_names=gas_names, **ppxf_kwargs)#, gas_reddening=gas_reddening)
        
        # When the two Delta Chi^2 below are the same, the solution
        # is the smoothest consistent with the observed spectrum.
        #
        print('Desired Delta Chi^2: %#.4g' % np.sqrt(2*self.lspec.size))
        print('Current Delta Chi^2: %#.4g' % ((self.pp.chi2 - 1)*self.lspec.size))
        print('Elapsed time in PPXF: %.2f s' % (clock() - t))
    
        weights = self.pp.weights[~gas_component]  # Exclude weights of the gas templates

        # pPXF gas fluxes are in normalised units * spectral pixel. We want
        # Physical units.
        print('Gas fluxes in relevant units are:\n')
        print('Name flux [flux density * Angstrom]\n')
        for name, flux, templ in zip(
            self.pp.gas_names, self.pp.gas_flux, self.gas_templates.T):

            # Get pixel size in Angstrom around the wavelength of the line.
            wavelength = np.argmax(templ)
            wavelength = np.exp(templlib.log_lam_temp[wavelength])
            wavelength = np.argmin(np.abs(self.pp.lam-wavelength))
            dwave = self.pp.lam[wavelength] - self.pp.lam[wavelength-1]
            dwave *= self.wave.unit

            print(name, flux*dwave*self.spec_norm)


    def continuum_subtraction(self):

        if not hasattr(self, 'pp'):
            raise ValueError('Run ppxf first')

        self.cont = np.interp(
            self.wave, np.exp(self.lwave)*units.Angstrom,
            self.pp.bestfit-self.pp.gas_bestfit) * self.spec_norm
        self.gas_model = np.interp(
            self.wave, np.exp(self.lwave)*units.Angstrom,
            self.pp.gas_bestfit) * self.spec_norm
        self.contsub_spec = self.spec - self.cont

    def display_galaxy(self,
        website_path=os.path.join('/export/data/fdeugenio/',
            'legac_team_share/legac_website/galaxyfiles')
        ):

        display_galaxy(
            self.legacid, obsmask=self.obsmask, website_path=website_path)



def display_galaxy(
    legacid, obsmask=None,
    website_path=os.path.join('/export/data/fdeugenio/',
        'legac_team_share/legac_website/galaxyfiles')
    ):

    obsmask = '*' if obsmask is None else obsmask

    filenames = os.path.join(website_path, f'LEGAC_M{obsmask}_{legacid}.html')
    filenames = glob.glob(filenames)

    if len(filenames)==1:
        filename = os.path.join('file://', filenames[0])
        print(f'Opening file {filename}...')
        subprocess.Popen(['firefox', '-new-window', filename])
    else:
        raise ValueError(
            f'ID {legacid}/Mask {obsmask} must match exactly one file.\n'
            f'(found {filenames}')
    return fspec
