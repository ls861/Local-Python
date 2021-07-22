# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 13:27:04 2021

@author: LSand
"""

import numpy as np
from astropy.io import fits
from os import listdir, path
import matplotlib.pyplot as plt


runtimes = []

folder = '/Users/LSand/Documents/GitHub/Local-Python/Year2/NSDC2/'
# folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/NSDC2/fit_R100/'

for filename in listdir(folder):
    if '.gz' in filename:
        filepath = folder+filename
        if path.getsize(filepath) > 100:
            data = fits.open(filepath)
            print(data[0].header['HIERARCH RUN TIME'])
            runtimes.append(data[0].header['HIERARCH RUN TIME'])
        
# np.save('runtimes.npy', runtimes)





#%%

runtimes = np.array(np.load(folder+'runtimes_cluster.npy')) / 3600

print(len(runtimes))
plt.hist(runtimes, bins=20)
plt.xlabel('runtime (hours)')
plt.show()



