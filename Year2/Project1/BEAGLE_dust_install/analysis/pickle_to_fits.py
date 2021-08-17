#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 01:37:11 2021

@author: lester
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
import pickle

outputDict = {}

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_17t_data_z0p5.p'

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_18_subset_zgt5.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_19_subset_zgt3p5_lt5.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_20_data_z0p5.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_21_data_z0p5.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_21_no_chi2_data_z0p5.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_22_data_z0p5.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_22_no_chi2_data_z0p5.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_22_gt_chi2_data_z0p5.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/data.p'

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/astrodeep_pickle.p'
#
#data = pickle.load(open(fileName,'r'))
#
#outputTable = Table(data)
#outputTable.write(fileName.replace('.p','.fits'), overwrite=True)


#fileName2 = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/recalc_chi2/astrodeep_pickle.p'
#
#data = pickle.load(open(fileName2,'r'))
#
#outputTable = Table(data)
#outputTable.write(fileName.replace('.p','_original.fits'), overwrite=True)


#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_24_vis1.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_24_13.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_24_14.p'
fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/astrodeep_pickle.p'
fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/data.p'
#fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_24_vis2.p'
fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_24_vis3.p'
fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_25_vis3_check_clusters.p'
fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_25_vis3_check_parallels.p'



fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_26_clusters_z1p25-6p0.p'


fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_28_vis4.p'
fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_28_vis4_x35.p'

fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/vis5_low_slope_issue.p'

data = pickle.load(open(fileName,'r'))

outputTable = Table(data)
outputTable.write(fileName.replace('.p','.fits'), overwrite=True)




#%%

fields = ['clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels']
z_bins = ['z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0', 'z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0']

for s in range(len(fields)):
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_27_{}_{}.p'.format(fields[s], z_bins[s])
    
    data = pickle.load(open(fileName,'r'))
    
    outputTable = Table(data)
    outputTable.write(fileName.replace('.p','.fits'), overwrite=True)




#%%

fields = ['clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels']
z_bins = ['z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0', 'z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0']

fields = ['clusters']
z_bins = ['z3p0-6p0']

for s in range(len(fields)):
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_{}_{}.p'.format(fields[s], z_bins[s])
    
    data = pickle.load(open(fileName,'r'))
    
    outputTable = Table(data)
    outputTable.write(fileName.replace('.p','.fits'), overwrite=True)


#%%

fields = ['clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels', 'clusters+parallels']
z_bins = ['z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0', 'z1p25-6p0', 'z1p25-2p0', 'z2p0-3p0', 'z3p0-4p0', 'z4p0-5p0', 'z5p0-6p0']

#fields = ['clusters']
#z_bins = ['z3p0-6p0']

for s in range(len(fields)):
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_31_{}_{}.p'.format(fields[s], z_bins[s])
    
    data = pickle.load(open(fileName,'r'))
    
    outputTable = Table(data)
    outputTable.write(fileName.replace('.p','.fits'), overwrite=True)
    
    
    
    
    
    
    


