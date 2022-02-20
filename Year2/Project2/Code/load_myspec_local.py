import legac_utils as legac
import numpy as np
import pickle
import matplotlib.pyplot as plt
import traceback
import datetime
from os import getcwd, listdir, path

# with open('.pickle', 'rb') as input_file:
#     x = pickle.load(input_file)

# with open('.pickle', 'wb') as output_file:
#     pickle.dump(x, output_file)

# from LEGAC UTILS
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

#%%

# =============================================================================
# LOAD PICKLE FILES
# =============================================================================

directory = './myspec/'

filenames_myspec = []
for filename in listdir(directory):
    if 'myspec' in filename and '.pickle' in filename:
        filenames_myspec.append(filename)

for filename in filenames_myspec:
    with open(directory+filename, 'rb') as input_file:
        myspec = pickle.load(input_file)

    # myspec.load(plot=True, clean=False) # clean=True is better but takes a long time.

    # LEGA-C assumes uniform spectral resolution equal to 0.9 Angstrom
    # (sigma of a Gaussian line-spread function).
    FWHM_gal = 0.9 * (2. * np.sqrt(2. * np.log(2.))) # AA.
    
    #------------------- Setup templates -----------------------
    
    # The templates are not normalized.
    ppxf_dir = os.path.dirname(os.path.realpath(ppxf.__file__))
    pathname = os.path.join(ppxf_dir, 'miles_models/Mun1.30*.fits')
    templlib = miles_util.miles(pathname, myspec.velscale, FWHM_gal)

    # pPXF gas fluxes are in normalised units * spectral pixel. We want
    # Physical units.

    print('\n')
    print(filename)
    print('Name flux [flux density * Angstrom]')
    
    for name, flux, templ in zip(
        myspec.pp.gas_names, myspec.pp.gas_flux, myspec.gas_templates.T):
    
        # Get pixel size in Angstrom around the wavelength of the line.
        wavelength = np.argmax(templ)
        wavelength = np.exp(templlib.log_lam_temp[wavelength])
        wavelength = np.argmin(np.abs(myspec.pp.lam-wavelength))
        dwave = myspec.pp.lam[wavelength] - myspec.pp.lam[wavelength-1]
        dwave *= myspec.wave.unit
    
        print(name, flux*dwave*myspec.spec_norm)
        

            
        
        
    
    





