#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:57:57 2019

@author: lester
"""

#import sys
import numpy as np
import matplotlib.pyplot as plt

#np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(threshold=25, edgeitems=10)

cat = np.load('cat.npy') # -ve z removed (z=-1), ID > 100000 removed (inc z=99)


###  ###  ###

title = ['Total', 'Clusters', 'Parallels']
tcp = [cat, cat[np.mod(cat['field'], 2) == 0], cat[np.mod(cat['field'], 2) == 1]]



### Finding filter which includes rest frame 1500 UV, yet excludes 1215.
    # will determine 'in' filter by fwhm opposed to min and max. Can recode if necessary.
    
### FILTERS ###

filter_label = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])

filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])

filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])

filter_hwhm = filter_fwhm / 2

leftedge = [ 3879.15,   4765. ,    7047.65 ,  9072.355, 10993.11,  12005.69,  13995.98, 19815.39,  31749.765 ,39975.9  ]

#leftedge   = [ 3879.15 4765.   7047.65  9072.355    10993.11  12005.69  13995.98 19815.39  31749.765 39975.9]

#rightedge  = [ 4818.15 7087.94 8903.65  11989.385   13998.31  15946.57  16870.16 23065.31  39181.475 50072.72 ]


### ### ###
    
z = np.array([5], dtype=float)


uv = (1+z) * 1500        # UV rest frame per galaxy
la = (1+z) * 1216        # lyman alpha rest frame per galaxy

filter_firstguess = np.zeros(len(uv), dtype=int)  # if UV in only 1 filter, this is that filter, if UV in > 1 filter, this is the nearest  to central filter, and if UV in 0 filters, this is next reddest filter. 

for i in range(len(uv)):
    
    in_filter = np.zeros(len(filter_fwhm_centre), dtype=float) 
    
    
    # determine which filters UV is IN
    for j in range(len(filter_fwhm_centre)):
        
        
        # UV in filter?
        if np.abs(uv[i] - filter_fwhm_centre[j]) < filter_hwhm[j]: 
            in_filter[j] = 1
        else:
            in_filter[j] = 0
            
    # UV in only 1 filter, select that filter 
    if sum(in_filter) == 1:
        
        filter_firstguess[i] = np.nonzero(in_filter)[0]
        
        
    # UV in > 1 filter, select nearest filter
    elif sum(in_filter) > 1:
        
        # select filter based on UV inside and nearest to centre of filter
        for j in range(len(filter_fwhm_centre)):
            
            in_filter[j] *= np.abs(uv[i] - filter_fwhm_centre[j])


        # this is finding the nearest central wavelength to UV[i]
        for j in np.nonzero(in_filter)[0]:
            
            if in_filter[j] <= in_filter[np.nonzero(in_filter)[0][0]]:
                
                filter_firstguess[i] = j
            

    # UV in NO filters, select nearest filter    
    elif sum(in_filter) == 0:

        # select filter nearest to centre of filter
        for j in range(len(filter_fwhm_centre)):
            
            # distance to each filter
            in_filter[j] = np.abs(uv[i] - filter_fwhm_centre[j])


        # this is finding the nearest central wavelength to UV[i]
        min_in_filter = in_filter[0]
        
        for j in range(len(in_filter)):
            
            if in_filter[j] <= min_in_filter:
                
                min_in_filter = in_filter[j]
                filter_firstguess[i] = j
        
    #### comparing with LA 1216, if LA in firstguess, then try redder filter and retest. 
    
    while filter_firstguess[i] <= len(filter_fwhm_centre):
        
        #print('while', i)
        #print(filter_firstguess[i])
        
        print('lester')
        print(filter_firstguess[i])
        print(filter_fwhm_centre[filter_firstguess[i]])
        print(filter_hwhm[filter_firstguess[i]])
        print(filter_fwhm_centre[filter_firstguess[i]] - filter_hwhm[filter_firstguess[i]])
        print(la[i])
        
        if filter_firstguess[i] == len(filter_fwhm_centre):
            print("there aren't enough red filters you muppet")
            break
           
            

        
        elif filter_fwhm_centre[filter_firstguess[i]] - filter_hwhm[filter_firstguess[i]]  > la[i]:
            # first guess is suitable
            print('lester')
            print(filter_firstguess[i])
            print(filter_fwhm_centre[filter_firstguess[i]])
            print(filter_hwhm[filter_firstguess[i]])
            print(filter_fwhm_centre[filter_firstguess[i]] - filter_hwhm[filter_firstguess[i]])
            print(la[i])
            break
        
        else:
            # must try redder filter
            filter_firstguess[i] += 1
            
            
            
            
            
print(filter_firstguess)
            
            
            
            
            
            
            
            
            
            
            
            
            
            