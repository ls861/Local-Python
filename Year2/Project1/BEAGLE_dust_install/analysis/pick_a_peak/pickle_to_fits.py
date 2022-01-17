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


fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/vis5_low_slope_issue.p'

data = pickle.load(open(fileName,'r'))

outputTable = Table(data)
outputTable.write(fileName.replace('.p','.fits'), overwrite=True)




    
    
    
    
    
    


