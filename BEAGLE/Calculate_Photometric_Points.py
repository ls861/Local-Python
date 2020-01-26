#!/usr/bin/env python
# -*- coding: utf-8 -*-

### Calculating Photometric Points compared with M_AB ###

import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy import interp
import numpy as np


c = 299792458 # m s^-2

### PLOT ###

plt.figure(figsize=(16,12)) 
plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)

ax1=plt.gca()
ax1.set_xlim(3000,55000)
ax1.set_ylim(10**-17, 10**-13)

ax1.set_yscale('log')
ax1.set_ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)


### SPECTRA ###     #8 is 10^7, #7 is 10^8, #9 is 10^6
    
fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Lester_z_shift_grid_neb_9/my_first_test.fits'
data_fits = fits.open(fileName)

wl_spec = data_fits[6].data[0][0]
redshift = data_fits[1].data



### PLOT SEDS ###
for i in range(len(redshift)):
        
    z = redshift[i][0]
    data_spec = data_fits[7].data[i]
    # [erg s^-1 cm^-2 A^-1]
    
    wl_spec_z = (1 + z) * wl_spec
    
    l = wl_spec_z
    lfl = wl_spec * data_spec
    
    #plt.plot(l, lfl, label = 'z = %i' % z, zorder = 2)

#plt.legend(loc = 'upper left')




### FILTERS ###
### JADES ### same as UVUDF plus Spitzer

fileName = '/Users/lester/BEAGLE/BEAGLE-general/filters/JADES_mock_filters_fixed.fits'
filter_fits = fits.open(fileName)

#print(catalogue.info())

filters = ['UVUDF-HST_ACS_WFC_F435W',
           'UVUDF-HST_ACS_WFC_F606W',
           'UVUDF-HST_ACS_WFC_F814W',
           'UVUDF-HST_WFC3_IR_F105W',
           'UVUDF-HST_WFC3_IR_F125W',
           'UVUDF-HST_WFC3_IR_F140W',
           'UVUDF-HST_WFC3_IR_F160W',
           'CANDELS_GS_filter12_IRACch1',
           'CANDELS_GS_filter12_IRACch2']

filter_label = ['F435W',
                'F606W',
                'F814W',
                'F105W',
                'F125W',
                'F140W',
                'F160W',
                'IRAC 3.6 $\mu m$',
                'IRAC 4.5 $\mu m$']

filter_centre = [4350, 6060, 8140, 10500, 12500, 14000, 16000, 36000, 45000]    #these are just related to the names of the filters, potentially the mean?

filter_fwhm_centre = [4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 35465.62, 45024.31]

filter_fwhm = [939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 7431.71, 10096.82]

filter_Ew = [823.16, 1771.39, 1886.72, 2371.86, 2674.41, 3569.84, 2750.2, 6836.16, 8649.93]



### M_AB Photometric Points ###

data_photo_AB = data_fits[8].data

for i in range(len(redshift)): # number of redshifts
    
    
    colour = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728']
    
    # lambda f_lambda

    data_photo = (c / np.array(filter_fwhm_centre)) * (10 ** (-(np.array(data_photo_AB[i])+23.6)/2.5))
    
    plt.scatter(filter_fwhm_centre, data_photo, marker = 'x', color = 'k', s = 100, zorder=3)
    plt.errorbar(filter_fwhm_centre, data_photo, xerr=np.array(filter_fwhm)/2, linestyle='None', color = colour[i], zorder = 1)




### Calculated Photometric Points using filter wavelengths ###

for j in range(len(redshift)): # number of redshifts
    
    z = redshift[j][0]
    data_spec = data_fits[7].data[j]
    # [erg s^-1 cm^-2 A^-1]
 
    # just easier for me to get my head around
    lambda_ = (1 + z) * wl_spec
    f_lambda = data_spec / (1+z)            # [erg s^-1 cm^-2 A^-1]
    
    

        
    # Interpolation

    f = interp1d(lambda_, f_lambda)

    
    photo_F_array = []
    
    for k in range(len(filters)): # number of filters
        
        wl_filter = filter_fits[1].data[filters[k]][0][0] # evenly spaced
        data_filter = filter_fits[1].data[filters[k]][0][1]
        
        f_lambda_interp = f(wl_filter)

        wl_spacing = wl_filter[1] - wl_filter[0]

        # photometric flux
        
        photo_F = sum( f_lambda_interp * data_filter * wl_spacing ) / sum( data_filter * wl_spacing ) 
        # JUST F_LAMBDA, need to multiply by lambda to plot on same as before.

        photo_F_array.append(photo_F * filter_fwhm_centre[k])  # this gives lambda f lambda
        
    plt.scatter(filter_fwhm_centre, photo_F_array, marker = 'o', color = 'y', s = 100, zorder=2)
    
    print(filter_fwhm_centre)
    print(photo_F_array)   



### Calculated Photometric Points using SED wavelengths ###

for j in range(len(redshift)): # number of redshifts
    
    z = redshift[j][0]
    data_spec = data_fits[7].data[j]
    # [erg s^-1 cm^-2 A^-1]
 
    # just easier for me to get my head around
    lambda_ = (1 + z) * wl_spec
    f_lambda = data_spec / (1+z)            # [erg s^-1 cm^-2 A^-1]
         
    
    
    # need to create an array of wavelength spacings to use  
    
    wl_spacing = np.empty(len(lambda_))
    
    for l in range(len(lambda_)):
        
        if l == 0:
            
            wl_spacing[l] = lambda_[1] - lambda_[0]
            
        elif l == len(lambda_) - 1:
            
            wl_spacing[l] = lambda_[len(lambda_)-1] - lambda_[len(lambda_)-2]
            
        else:
            
            wl_spacing[l] = 0.5 * (lambda_[l+1] - lambda_[l-1])


    photo_F_array = []
    
    for k in range(len(filters)): # number of filters
        
        
        wl_filter = filter_fits[1].data[filters[k]][0][0] # evenly spaced
        data_filter = filter_fits[1].data[filters[k]][0][1]      
        
        # Interpolation, NUMPY THIS TIME to INCLUDE EXTRAPOLATION
        
        
        data_filter_interp = np.interp(lambda_, wl_filter, data_filter, left=0, right=0) 
        

        # photometric flux
        
        photo_F = sum( f_lambda * data_filter_interp * wl_spacing ) / sum( data_filter_interp * wl_spacing ) 
        #photo_F = sum( f_lambda_interp * data_filter * wl_spacing ) / sum( data_filter * wl_spacing ) 
        # JUST F_LAMBDA, need to multiply by lambda to plot on same as before.





        photo_F_array.append(photo_F * filter_fwhm_centre[k])  # this gives lambda f lambda
        
    plt.scatter(filter_fwhm_centre, photo_F_array, marker = '+', color = 'b', s = 100, zorder=3)
    
    print(filter_fwhm_centre)
    print(photo_F_array)   


plt.show()



### just to get an idea of what wavelengths there are and the spread.
wl_filter = filter_fits[1].data[filters[8]][0][0] # evenly spaced

print(lambda_[(lambda_ >45000) & (lambda_ < 45100)])
print(wl_filter[(wl_filter >45000) & (wl_filter < 45100)])

wl_filter = filter_fits[1].data[filters[3]][0][0] # evenly spaced

print(lambda_[(lambda_ >10000) & (lambda_ < 10100)])
print(wl_filter[(wl_filter >10000) & (wl_filter < 10100)])

'''
ax1=plt.gca()
ax1.set_xlim(300,55000)
ax1.set_ylim(10**-23, 10**-16)
ax1.set_yscale('log')
plt.plot(lambda_, f_lambda)
plt.plot(wl_filter, f_lambda_interp)
plt.show()
'''     
        

### OLD METHOD ###
'''
photo = []

for i in range(len(filters)):
    wl_filter = filter_fits[1].data[filters[i]][0][0]
    data_filter = filter_fits[1].data[filters[i]][0][1]
    num = len(wl_filter)
    
    wl_spec = data_fits[4].data[0][0]
    data_spec = data_fits[5].data[0]

    total = 0
    
    for j in range(len(wl_filter)):
        
        if data_filter[j] > 0.2:
            
            ### the if statement was the only way I could get my head around the normalisation issues ###
            total += data_filter[j] * data_spec[index(wl_filter[j], wl_spec)]
    
    avg = total / (num * norm)
    photo.append(avg)
    
print(photo)
'''
##################


    





























