#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 14:42:28 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Jan_20_1_004/mock_catalogue.fits'
data_fits = fits.open(fileName)

#info = data_fits.info()
#header = data_fits['GALAXY PROPERTIES'].header


wl_spec = data_fits[6].data[0][0]
redshift = data_fits[1].data

sl  = 3
z   = data_fits['GALAXY PROPERTIES'].data['redshift'][::sl]
wl  = data_fits['FULL SED WL'].data['wl'][0]
f   = data_fits['FULL SED'].data[::sl]    # erg s-1 cm-2 A-1
tau = np.arange(7.5, 10, 0.2)[::sl]

plt.figure(figsize=(10, 10))
plt.xlim(0, 20000)
#plt.ylim(bottom=0, top=5E-20)

for i in range(len(z)):
    plt.plot(wl, np.log10(f[i]), label=r'$\tau$ = {%.1f}' % (tau[i]) )


plt.legend()
plt.show()



data_fits.close()

print(len(f))