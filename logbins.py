#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 17:39:47 2019

@author: lester
"""

import numpy as np


logbins = np.geomspace(0.00000001, 25000, 100)
logbins = np.concatenate((-logbins[::-1], np.array([0]), logbins))

