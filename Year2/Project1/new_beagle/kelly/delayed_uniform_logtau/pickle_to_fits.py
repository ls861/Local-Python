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
fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.p'

fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/data.p'

data = pickle.load(open(fileName,'r'))

outputTable = Table(data)
outputTable.write(fileName.replace('.p','.fits'), overwrite=True)


























