
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

import matplotlib
matplotlib.rcParams.update({'font.size': 18})


filters = ['HST_ACS_WFC_F435W', 'HST_ACS_WFC_F606W', 'HST_ACS_WFC_F814W', 'HST_WFC3_IR_F105W', 'HST_WFC3_IR_F125W', 'HST_WFC3_IR_F140W', 'HST_WFC3_IR_F160W', 'Paranal_HAWKI_Ks', 'Spitzer_IRAC_I1', 'Spitzer_IRAC_I2']

filter_label = ['F435W', 'F606W', 'F814W', 'F105W', 'F125W', 'F140W', 'F160W', 'Ks', 'IRAC 3.6 $\mu m$', 'IRAC 4.5 $\mu m$']

fileName = '/Users/lester/BEAGLE/Filter_Files/Ascii_Astrodeep/astrodeep_filters.fits'

wl = []
f = []

data_fits = fits.open(fileName)


for i in range(len(filters)):
    wl.append(data_fits['TRANSMISSION'].data[filters[i]][0][0])
    f.append(data_fits['TRANSMISSION'].data[filters[i]][0][1])

data_fits.close()


### PLOT ###

fig, ax = plt.subplots(figsize=(16, 6))
#fig.suptitle('Astrodeep Filters')

for i in range(len(filters)):
    ax.fill_between(wl[i], f[i], alpha=0.7, label=filter_label[i])

ax.set_xlim(0,55000)
ax.set_ylim(0,1)
ax.set_xlabel(r'Wavelength / $\AA$')
ax.set_ylabel('Transmission')
plt.legend()
plt.show()

### ###  ###

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

