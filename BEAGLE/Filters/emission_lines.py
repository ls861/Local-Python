#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 09:50:07 2021

@author: lester
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 10})

### FILTERS ###
# filter_label = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])
# filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])
# filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])

filter_label = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'H160', 'Ks', 'CH1', 'CH2'])
filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 15433.07, 21440.35, 35465.62, 45024.31])
filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 2874.18, 3249.92, 7431.71, 10096.82])

left_edges = filter_fwhm_centre - (filter_fwhm / 2.0)
right_edges = filter_fwhm_centre + (filter_fwhm / 2.0)

### EMISSION LINES ###
line_wl = np.array([912.0, 1025.73, 1215.67, 1548.19, 1908.73, 3646.0, 3729.23, 4101.74, 4340.47, 4861.33, 5006.84, 6548.05, 6562.82, 6583.45, 6730.82])
line_label = np.array(['LyBREAK', 'HLyB@1026', 'HLyA@1216','C4@1548','C3@1910','BaBREAK', 'O2@3729', 'HBaD@4102', 'HBaG@4340', 'HBaB@4861', 'O3@5007', 'N2@6548', 'HBaA@6563', 'N2@6584', 'S2@6731'])

z_lower = 0
z_upper = 6
z_tmp = np.linspace(z_lower, z_upper, 1000)

### ### ###
fig, ax = plt.subplots(figsize=(16,5))
plt.suptitle('Astrodeep Filters')
ax.set_xlim(0, 52000)
ax.set_ylim(z_lower, z_upper)
ax.set_xlabel(r'Observed Wavelength / $\AA$')
ax.set_ylabel(r'Redshift')

cmap = plt.get_cmap("tab10")
patterns = ['', '', '', '.', '', '', '', '', '']



for i in range(len(left_edges)):
    ax.axvspan(left_edges[i], right_edges[i], alpha=0.3, color=cmap(i), hatch=patterns[i])

for i in range(len(line_wl)):
    if line_label[i] in ['LyBREAK', 'BaBREAK']:
        ax.plot((1+z_tmp)*line_wl[i], z_tmp, color='r', linestyle='dashed')
        
    elif line_label[i] in ['HLyB@1026', 'HLyA@1216', 'HBaD@4102', 'HBaG@4340', 'HBaB@4861', 'HBaA@6563']:
        ax.plot((1+z_tmp)*line_wl[i], z_tmp, color='r')  
        
    elif line_label[i] in ['O2@3729', 'O3@5007']:
        ax.plot((1+z_tmp)*line_wl[i], z_tmp, color='blue')          

    elif line_label[i] in ['C4@1548', 'C3@1910']:
        ax.plot((1+z_tmp)*line_wl[i], z_tmp, color='green')   
        
    elif line_label[i] in ['N2@6548', 'N2@6584']:
        ax.plot((1+z_tmp)*line_wl[i], z_tmp, color='orange')   

    elif line_label[i] in ['S2@6731']:
        ax.plot((1+z_tmp)*line_wl[i], z_tmp, color='k')   
        
    else:
        ax.plot((1+z_tmp)*line_wl[i], z_tmp, color='gray')

for i in range(len(line_wl)):

    if line_label[i] == 'O2@3729':
        ax.text((1+z_upper)*line_wl[i]-1500, z_upper+0.2, 'O2@3729', color='blue')
    elif line_label[i] == 'O3@5007':
        ax.text((1+z_upper)*line_wl[i]-1500, z_upper+0.2, 'O3@5007', color='blue')
    elif line_label[i] == 'C4@1548':
        ax.text((1+z_upper)*line_wl[i]-1500, z_upper+0.2, 'C4@1548', color='green')
    elif line_label[i] == 'C3@1910':
        ax.text((1+z_upper)*line_wl[i]-1500, z_upper+0.5, 'C3@1910', color='green')

    elif line_label[i] == 'HLyB@1026':
        ax.text((1+z_upper)*line_wl[i]-1500, z_lower-1.1, 'HLyB@1026', color='red')
    elif line_label[i] == 'HLyA@1216':
        ax.text((1+z_upper)*line_wl[i]-1500, z_lower-1.4, 'HLyA@1216', color='red')
    elif line_label[i] == 'HBaD@4102':
        ax.text((1+z_upper)*line_wl[i]-1500, z_lower-1.1, 'HBaD@4102', color='red')
    elif line_label[i] == 'HBaG@4340':
        ax.text((1+z_upper)*line_wl[i]-1500, z_lower-1.4, 'HBaG@4340', color='red')
    elif line_label[i] == 'HBaB@4861':
        ax.text((1+z_upper)*line_wl[i]-1500, z_lower-1.7, 'HBaB@4861', color='red')
    elif line_label[i] == 'HBaA@6563':
        ax.text((1+z_upper)*line_wl[i]-1500, z_lower-1.1, 'HBaA@6563', color='red')
        
    elif line_label[i] == 'N2@6548':
        ax.text((1+z_upper)*line_wl[i]-1500, z_upper+0.2, 'N2@6548', color='orange')   
    elif line_label[i] == 'N2@6584':
        ax.text((1+z_upper)*line_wl[i]-1500, z_upper+0.5, 'N2@6584', color='orange')   
    elif line_label[i] == 'S2@6731':
        ax.text((1+z_upper)*line_wl[i]-1500, z_lower-1.4, 'S2@6731', color='k')


plt.grid(axis = 'y')
plt.show()


### ### ###
fig, ax = plt.subplots(figsize=(16,5))
plt.suptitle('Astrodeep Filters')
ax.set_xlim(0, 10000)
ax.set_ylim(z_lower, z_upper)
ax.set_xlabel(r'Restframe Wavelength / $\AA$')
ax.set_ylabel(r'Redshift')

for i in range(len(left_edges)):
    ax.fill_betweenx(z_tmp, left_edges[i]/(1+z_tmp), right_edges[i]/(1+z_tmp), alpha=0.3, color=cmap(i), hatch=patterns[i])
    
for i in range(len(line_wl)):
    if line_label[i] in ['LyBREAK', 'BaBREAK']:
        ax.plot((line_wl[i], line_wl[i]), (z_lower, z_upper), color='r', linestyle='dashed')
        
    elif line_label[i] in ['HLyB@1026', 'HLyA@1216', 'HBaD@4102', 'HBaG@4340', 'HBaB@4861', 'HBaA@6563']:
        ax.plot((line_wl[i], line_wl[i]), (z_lower, z_upper), color='r')  
    elif line_label[i] in ['O2@3729', 'O3@5007']:
        ax.plot((line_wl[i], line_wl[i]), (z_lower, z_upper), color='blue')          
    elif line_label[i] in ['C4@1548', 'C3@1910']:
        ax.plot((line_wl[i], line_wl[i]), (z_lower, z_upper), color='green')   
    elif line_label[i] in ['N2@6548', 'N2@6584']:
        ax.plot((line_wl[i], line_wl[i]), (z_lower, z_upper), color='orange')   
    elif line_label[i] in ['S2@6731']:
        ax.plot((line_wl[i], line_wl[i]), (z_lower, z_upper), color='k')   
    else:
        ax.plot((line_wl[i], line_wl[i]), (z_lower, z_upper), color='gray')

    if line_label[i] == 'O2@3729':
        ax.text(line_wl[i]-270, z_upper+0.2, 'O2@3729', color='blue')
    elif line_label[i] == 'O3@5007':
        ax.text(line_wl[i]-270, z_upper+0.2, 'O3@5007', color='blue')
    elif line_label[i] == 'C4@1548':
        ax.text(line_wl[i]-270, z_upper+0.2, 'C4@1548', color='green')
    elif line_label[i] == 'C3@1910':
        ax.text(line_wl[i]-270, z_upper+0.5, 'C3@1910', color='green')

    elif line_label[i] == 'HLyB@1026':
        ax.text(line_wl[i]-270, z_lower-1.1, 'HLyB@1026', color='red')
    elif line_label[i] == 'HLyA@1216':
        ax.text(line_wl[i]-270, z_lower-1.4, 'HLyA@1216', color='red')
    elif line_label[i] == 'HBaD@4102':
        ax.text(line_wl[i]-270, z_lower-1.1, 'HBaD@4102', color='red')
    elif line_label[i] == 'HBaG@4340':
        ax.text(line_wl[i]-270, z_lower-1.4, 'HBaG@4340', color='red')
    elif line_label[i] == 'HBaB@4861':
        ax.text(line_wl[i]-270, z_lower-1.7, 'HBaB@4861', color='red')
    elif line_label[i] == 'HBaA@6563':
        ax.text(line_wl[i]-270, z_lower-1.1, 'HBaA@6563', color='red')
        
    elif line_label[i] == 'N2@6548':
        ax.text(line_wl[i]-270, z_upper+0.2, 'N2@6548', color='orange')   
    elif line_label[i] == 'N2@6584':
        ax.text(line_wl[i]-270, z_upper+0.5, 'N2@6584', color='orange')   
    elif line_label[i] == 'S2@6731':
        ax.text(line_wl[i]-270, z_lower-1.4, 'S2@6731', color='k')


#plt.legend()
plt.grid(axis = 'y')
plt.show()










import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
    

    
#%%









































