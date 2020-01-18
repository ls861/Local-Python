#!/usr/bin/env python
# -*- coding: utf-8 -*-

### EFFECTIVE WIDTH H ALPHA ###

import matplotlib.pyplot as plt
from astropy.io import fits




### Wavelength Values ### 

h_lower = [19659, 32765, 45871, 58977] #6553

h_centre = [19689, 32815, 45941, 59067] #6563

h_upper = [19719, 32865, 46011, 59157] #6573



### SPECTRA ###     #5/8 is 10^7, #6/7 is 10^8, #9 is 10^6
    
fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Lester_z_shift_grid_neb_9/my_first_test.fits'
data_fits = fits.open(fileName)

wl_spec = data_fits[6].data[0][0]
redshift = data_fits[1].data


for i in range(len(data_fits[7].data)):
    
    ### PLOT ###

    plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
    ax1=plt.gca()
    #ax1.set_yscale('log')
    ax1.set_ylabel('Flux $F_\lambda$ ($erg\:s^{-1}\:cm^{-2}\:\AA^{-1}$)', fontsize=14)
    
    z = redshift[i][0]
        
    data_spec = data_fits[7].data[i]        # [erg s^-1 cm^-2 A^-1]
    
    
    lambda_ = (1 + z) * wl_spec
    f_lambda = data_spec / (1+z)            # [erg s^-1 cm^-2 A^-1]


    lambda_temp = []
    f_lambda_temp = []
    
    for j in range(len(lambda_)):
        
        if lambda_[j] > h_lower[i] and lambda_[j] < h_upper[i]:
            
            lambda_temp.append(lambda_[j])
            f_lambda_temp.append(f_lambda[j])
            
    plt.plot(lambda_temp, f_lambda_temp)
    plt.show()
    print(lambda_temp) 
    print(f_lambda_temp)
    
    # area of line
    
    area = 0
    
    for k in range(len(lambda_temp)):
        
        area += (f_lambda_temp[k] - f_lambda_temp[0]) * (lambda_temp[1] - lambda_temp[0])

    print('Area under line: ', area)
    print('Equivalent Width: ', area / f_lambda_temp[0])
























