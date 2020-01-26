#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 23:07:30 2019

@author: lester
"""


import numpy as np
import matplotlib.pyplot as plt


### FILTERS ###

filter_label = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])

filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])

filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])

height = np.array([1, 2, 3, 1, 2, 3, 4, 1, 2, 3])


### ### ###


plt.figure(figsize=(16,5)) 
plt.xlim(0, 60000)
plt.ylim(-0.5, 5)
ax = plt.axes()
ax.axes.get_yaxis().set_visible(False)


    
for i in range(len(height)):
    plt.errorbar(filter_fwhm_centre[i], height[i], xerr=filter_fwhm[i]/2, zorder = 0, elinewidth=10)
    plt.text(filter_fwhm_centre[i]-800, height[i]-0.3, filter_label[i])
    

plt.text(filter_fwhm_centre[0]-800, 0, 'HST ACS')
plt.text(filter_fwhm_centre[3]-800, 0, 'HST WFC3')
plt.text(filter_fwhm_centre[7]-800, 0, 'VLT')
plt.text(filter_fwhm_centre[8]+3000, 0, 'Spitzer')

plt.legend()
plt.show()






