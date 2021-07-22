#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 19:13:18 2020

@author: lester
"""




'''
THIS SCRIPT WAS JUST BITS AND BOBS FOR TESTING, NEVER ACTUALLY RAN

'''

import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import datetime
import time


folder = '/Users/lester/Documents/PhD/param_DPL_backup/fit_012_003/'
folder = '/Users/lester/BEAGLE/BEAGLE-general/results/astrodeep_002/'


for item in os.listdir(folder):
    if '.fits.gz' in item and os.stat(folder+item).st_size > 0:
        
        with fits.open(folder+item) as data:
#            print(data.info())
            print(item)
            dt = data['PRIMARY'].header['DATE']
            t = datetime.datetime(int(dt[:4]), int(dt[5:7]), int(dt[8:10]), int(dt[11:13]), int(dt[14:16]), int(dt[17:19]))
            t_unix = time.mktime(t.timetuple())
            
            
            print(os.stat(folder+item).st_mtime - t_unix)
            
            
            
            
            print()
            print(os.stat(folder+item).st_mtime)
            print(os.stat(folder+item).st_atime)
            print(os.stat(folder+item).st_ctime)
            print(t_unix)
            print(os.stat(folder+item).st_mtime - t_unix)




from astropy.io import fits
data = fits.open('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/danger_delayed_uniform_logtau/2/fit_001/1028_BEAGLE.fits.gz')

from astropy.io import fits
data = fits.open('/Users/lester/Documents/PhD/param_100/fit_108_DE/203_BEAGLE.fits.gz')
d = data['PRIMARY'].header['HIERARCH RUN TIME'] 
print(d)

