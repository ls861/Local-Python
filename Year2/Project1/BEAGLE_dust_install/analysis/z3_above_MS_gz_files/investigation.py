#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:42:49 2021

@author: lester
"""

#http://astrodeep.u-strasbg.fr/ff/?ffid=FF_M0717CL&id=1315&img=o_H160&cm=cubehelix
#http://astrodeep.u-strasbg.fr/ff/data/SED/MACS0717cl/SED_1315.png


#http://astrodeep.u-strasbg.fr/ff/?ffid=FF_A2744CL&id=15&img=o_H160&cm=cubehelix


import urllib
from IPython.display import display
# from wand.image import Image as WImage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from astropy.io import fits
from astropy.wcs import WCS
import webbrowser
import numpy as np
from astropy.nddata import Cutout2D

#filename = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_18_subset_zgt5.fits'
filename = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/z3_above_MS_gz_files/z3_above_MS_detailed.fits'


fields = ['0A2744C','1A2744P','2M0416C','3M0416P','4M0717C','5M0717P','6M1149C','7M1149P']
fields_FF = ['A2744CL', 'A2744PAR', 'M0416CL', 'M0416PAR', 'MACS0717cl', 'M0717PAR', 'MACS1149cl', 'M1149PAR']
fields_FF_URL = ['A2744CL', 'A2744PAR', 'M0416CL', 'M0416PAR', 'M0717CL', 'M0717PAR', 'M1149CL', 'M1149PAR']
fields_for_images = ['abell2744', 'abell2744', 'macs0416', 'macs0416', 'macs0717', 'macs0717', 'macs1149', 'macs1149']
filters = ['435', '606', '814', '105', '125', '140', '160']

d = fits.open(filename)
#print(d.info())
#print(d[1].columns)

field_AD = d[1].data['field_AD']

id_AD = d[1].data['id_AD']
z_AD = d[1].data['redshift_AD']
mass_AD = d[1].data['mass_AD_neb']

id_BEAGLE = d[1].data['id_BEAGLE']
z_BEAGLE = d[1].data['redshift_BEAGLE']
mass_BEAGLE = d[1].data['mass_BEAGLE_stellar']

RA = d[1].data['RA_AD']
DEC = d[1].data['DEC_AD']

chrome_path = 'open -a /Applications/Google\ Chrome.app %s'
#%%
# =============================================================================
# download the SEDs from FF
# =============================================================================

'''
for i in range(len(id_AD)):

    path = 'http://astrodeep.u-strasbg.fr/ff/data/SED/{}/SED_{}.png'.format(fields_FF[int(d[1].data['field_AD'][i])], int(id_AD[i]))
    output ='/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/z3_above_MS_gz_files/{}/AD_SED/SED_F-{}_AD-{}_B-{}.png'.format(int(d[1].data['field_AD'][i]), int(d[1].data['field_AD'][i]), int(id_AD[i]), int(id_BEAGLE[i]))

    print(path)
    urllib.urlretrieve(path, output)    

d.close()
'''

# =============================================================================
# basic plots to start categorising stuff
# =============================================================================
#%%

idx1 = z_BEAGLE > 3.*z_AD
idx2 = (z_BEAGLE > 13.*z_AD/6.)&(z_BEAGLE < 3.*z_AD)
idx3 = (z_BEAGLE > 12.*z_AD/8.)&(z_BEAGLE < 13.*z_AD/6.)
idx4 = (z_BEAGLE > 7.*z_AD/10.)&(z_BEAGLE < 12.*z_AD/8.)
idx5 = z_BEAGLE < 7.*z_AD/10.

'''
x = np.array((0, 15))
plt.scatter(z_AD, z_BEAGLE)
plt.scatter(z_AD[idx1], z_BEAGLE[idx1])
plt.scatter(z_AD[idx2], z_BEAGLE[idx2])
plt.scatter(z_AD[idx3], z_BEAGLE[idx3])
plt.scatter(z_AD[idx4], z_BEAGLE[idx4])
plt.scatter(z_AD[idx5], z_BEAGLE[idx5])

plt.scatter(z_AD[idx1][z_BEAGLE[idx1]>3.5], z_BEAGLE[idx1][z_BEAGLE[idx1]>3.5], marker='x', color='k')
plt.scatter(z_AD[idx2][z_BEAGLE[idx2]>3.5], z_BEAGLE[idx2][z_BEAGLE[idx2]>3.5], marker='x', color='k')
plt.scatter(z_AD[idx3][z_BEAGLE[idx3]>3.5], z_BEAGLE[idx3][z_BEAGLE[idx3]>3.5], marker='x', color='k')
plt.scatter(z_AD[idx4][z_BEAGLE[idx4]>3.5], z_BEAGLE[idx4][z_BEAGLE[idx4]>3.5], marker='x', color='k')
plt.scatter(z_AD[idx5][z_BEAGLE[idx5]>3.5], z_BEAGLE[idx5][z_BEAGLE[idx5]>3.5], marker='x', color='k')

plt.plot(x, 3*x, color='k')
plt.plot(x, (13.*x/6.), color='k')
plt.plot(x, (12.*x/8.), color='k')
plt.plot(x, (7.*x/10.), color='k')
plt.plot(x, x, color='red', linestyle='dashed')
plt.plot(x, (3.5, 3.5), color='red', linestyle='dashed')
plt.xlim(0, 5)
plt.ylim(0, 5)
plt.xlabel('z AD')
plt.ylabel('z BEAGLE')
plt.show()

print(sum(idx1), sum(idx2), sum(idx3), sum(idx4), sum(idx5))
print(len(z_AD))
'''
#%%
# =============================================================================
# show marginals and SEDs
# =============================================================================

idx = np.full(len(id_AD), True)
# idx = idx4
#idx = np.logical_and(idx, d[1].data['field_AD']==0)
print(len(id_AD[idx]))
# =============================================================================
# print numbers out for excel
# =============================================================================
'''
for i in range(len(id_AD[idx])):
    print(int(id_BEAGLE[idx][i]))

for i in range(len(id_AD[idx])):
    FIELD = int(d[1].data['field_AD'][idx][i])
    print(fields_FF[FIELD])

for i in range(len(id_AD[idx])):
    print(int(id_AD[idx][i]))

for i in range(len(id_AD[idx])):
    print(int(field_AD[idx][i]))


'''
#%%
'''
# =============================================================================
# the real stuff
# =============================================================================


for i in range(len(id_AD[idx])):
#for i in range(150):
#for i in range(150, len(id_AD[idx])):
#for i in range(30):
#for i in [10]:  
    
    FIELD = int(d[1].data['field_AD'][idx][i])
    SED_B = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/z3_above_MS_gz_files/{}/gz_files/pyp-beagle/plot/{}_BEAGLE_marginal_SED_phot.png'.format(FIELD, int(id_BEAGLE[idx][i]))
#    SED_AD = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_4/{}/AD_SED/SED_F-{}_AD-{}_B-{}.png'.format(FIELD, FIELD, int(id_AD[idx][i]), int(id_BEAGLE[idx][i]))
    img_B = mpimg.imread(SED_B)  
#    img_AD = mpimg.imread(SED_AD)
    print(' \n \n{}\nASTRODEEP: id{} z{:.2f} m{:.2f}\nBEAGLE: id{} z{:.2f} m{:.2f}'.format(fields_FF[FIELD], int(id_AD[idx][i]), z_AD[idx][i], mass_AD[idx][i], int(id_BEAGLE[idx][i]), z_BEAGLE[idx][i], mass_BEAGLE[idx][i]))
    
    
    
    # get image per filter
    if FIELD == 0 or FIELD == 2:
        epochs = ['','','','','','','']
    elif FIELD == 4:
        epochs = ['-epoch1','-epoch1','-epoch1','-epoch1','-epoch1','-epoch1','-epoch1']
    elif FIELD == 6:
        epochs = ['-epoch2','-epoch2','-epoch1','-epoch1','-epoch1','-epoch1','-epoch1']        
    fi = 6
    camera = 'wfc3'
    image_file = '/Users/lester/BEAGLE/FF_images/plotting_task/{}/hlsp_frontier_hst_{}-30mas_{}_f{}w_v1.0{}_drz.fits'.format(FIELD, camera, fields_for_images[FIELD], filters[fi], epochs[fi])
    image_data = fits.getdata(image_file)
    w = WCS(image_file)
    pix = w.wcs_world2pix(RA[idx][i], DEC[idx][i], 0)
    position = (pix[0], pix[1])
    s = 1500
    s_zoom = 200
    size = (0.8*s, s)     # pixels
    cutout = Cutout2D(image_data, position, size)



    # PLOT SEDs
    f, axarr = plt.subplots(1,3,figsize=(26,20))
#    f.suptitle("{}\n ASTRODEEP: id{} z{:.2f} m{:.2f}\n BEAGLE: id{} z{:.2f} m{:.2f}".format(fields_FF[FIELD], int(id_AD[idx][i]), z_AD[idx][i], mass_AD[idx][i], int(id_BEAGLE[idx][i]), z_BEAGLE[idx][i], mass_BEAGLE[idx][i]), fontsize=15)
    axarr[0].imshow(img_B)
    axarr[1].imshow(img_B)
    axarr[2].imshow(cutout.data, origin='lower', vmin=0, vmax=np.percentile(cutout.data.flatten(), 95))
    axarr[2].scatter(s/2., (0.8)*s/2., marker='x', color='r')
    axarr[2].plot(np.array((s/2. - s_zoom/2., s/2. - s_zoom/2., s/2. + s_zoom/2., s/2. + s_zoom/2., s/2. - s_zoom/2.)), (0.8)*np.array((s/2. - s_zoom/2., s/2. + s_zoom/2., s/2. + s_zoom/2., s/2. - s_zoom/2., s/2. - s_zoom/2.)), color='w')    
    f.tight_layout()
    plt.show()
    


    # PLOT individual filters
    f, axarr = plt.subplots(1,len(filters),figsize=(26,15))
    for fi in range(len(filters)):
        if fi <= 2:
            camera = 'acs'
        else:
            camera = 'wfc3'
        image_file = '/Users/lester/BEAGLE/FF_images/plotting_task/{}/hlsp_frontier_hst_{}-30mas_{}_f{}w_v1.0{}_drz.fits'.format(FIELD, camera, fields_for_images[FIELD], filters[fi], epochs[fi])
        image_data = fits.getdata(image_file)
#        print(image_file)
        w = WCS(image_file)
        pix = w.wcs_world2pix(RA[idx][i], DEC[idx][i], 0)
#        print(pix)
        position = (pix[0], pix[1])
        size = (s_zoom, s_zoom)     # pixels
        cutout = Cutout2D(image_data, position, size)
        axarr[fi].title.set_text('{}'.format(filters[fi]))
        axarr[fi].imshow(cutout.data, origin='lower', vmin=0, vmax=np.percentile(cutout.data.flatten(), 99.5))
#        plt.colorbar()
        axarr[fi].scatter(s_zoom/2., s_zoom/2., marker='o', color='r', s=500, facecolors='none')
        axarr[fi].imshow(cutout.data, origin='lower', vmin=0, vmax=np.percentile(cutout.data.flatten(), 95))
        axarr[fi].axis('off')
    f.tight_layout()
    plt.show()
    
#for i in range(len(id_AD[idx])):    
#    url = 'http://astrodeep.u-strasbg.fr/ff/?ffid=FF_{}&id={}&img=o_H160&cm=cubehelix'.format(fields_FF_URL[FIELD], int(id_AD[idx][i]))
##    webbrowser.get(chrome_path).open(url)

print(sum(idx))

'''



#%%

# =============================================================================
# for notion
# =============================================================================

# FILTERS

for i in range(len(id_AD[idx])):
    
    FIELD = int(d[1].data['field_AD'][idx][i])
    # SED_B = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/z3_above_MS_gz_files/{}/gz_files/pyp-beagle/plot/{}_BEAGLE_marginal_SED_phot.png'.format(FIELD, int(id_BEAGLE[idx][i]))
#    SED_AD = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_4/{}/AD_SED/SED_F-{}_AD-{}_B-{}.png'.format(FIELD, FIELD, int(id_AD[idx][i]), int(id_BEAGLE[idx][i]))
    # img_B = mpimg.imread(SED_B)  
#    img_AD = mpimg.imread(SED_AD)
    print(' \n \n{}\nASTRODEEP: id{} z{:.2f} m{:.2f}\nBEAGLE: id{} z{:.2f} m{:.2f}'.format(fields_FF[FIELD], int(id_AD[idx][i]), z_AD[idx][i], mass_AD[idx][i], int(id_BEAGLE[idx][i]), z_BEAGLE[idx][i], mass_BEAGLE[idx][i]))
    
    
    
    # get image per filter
    if FIELD == 0 or FIELD == 2:
        epochs = ['','','','','','','']
    elif FIELD == 4:
        epochs = ['-epoch1','-epoch1','-epoch1','-epoch1','-epoch1','-epoch1','-epoch1']
    elif FIELD == 6:
        epochs = ['-epoch2','-epoch2','-epoch1','-epoch1','-epoch1','-epoch1','-epoch1']        
    fi = 6
    camera = 'wfc3'
    image_file = '/Users/lester/BEAGLE/FF_images/plotting_task/{}/hlsp_frontier_hst_{}-30mas_{}_f{}w_v1.0{}_drz.fits'.format(FIELD, camera, fields_for_images[FIELD], filters[fi], epochs[fi])
    image_data = fits.getdata(image_file)
    w = WCS(image_file)
    pix = w.wcs_world2pix(RA[idx][i], DEC[idx][i], 0)
    position = (pix[0], pix[1])
    s = 1500
    s_zoom = 200
    size = (0.8*s, s)     # pixels
    cutout = Cutout2D(image_data, position, size)


    # PLOT individual filters
    f, axarr = plt.subplots(len(filters), 1,figsize=(26,15))
    for fi in range(len(filters)):
        if fi <= 2:
            camera = 'acs'
        else:
            camera = 'wfc3'
        image_file = '/Users/lester/BEAGLE/FF_images/plotting_task/{}/hlsp_frontier_hst_{}-30mas_{}_f{}w_v1.0{}_drz.fits'.format(FIELD, camera, fields_for_images[FIELD], filters[fi], epochs[fi])
        image_data = fits.getdata(image_file)
#        print(image_file)
        w = WCS(image_file)
        pix = w.wcs_world2pix(RA[idx][i], DEC[idx][i], 0)
#        print(pix)
        position = (pix[0], pix[1])
        size = (s_zoom, s_zoom)     # pixels
        cutout = Cutout2D(image_data, position, size)
        if fi == 0:
            axarr[fi].title.set_text('{} {}'.format(int(field_AD[idx][i]), int(id_AD[idx][i])))
        else:
            axarr[fi].title.set_text('{}'.format(filters[fi]))
        
        axarr[fi].imshow(cutout.data, origin='lower', vmin=0, vmax=np.percentile(cutout.data.flatten(), 99.5))
#        plt.colorbar()
        axarr[fi].scatter(s_zoom/2., s_zoom/2., marker='o', color='r', s=500, facecolors='none')
        axarr[fi].imshow(cutout.data, origin='lower', vmin=0, vmax=np.percentile(cutout.data.flatten(), 95))
        axarr[fi].axis('off')
    f.tight_layout()

    plt.show()





#for i in range(len(id_AD[idx])):    
#    url = 'http://astrodeep.u-strasbg.fr/ff/?ffid=FF_{}&id={}&img=o_H160&cm=cubehelix'.format(fields_FF_URL[FIELD], int(id_AD[idx][i]))
##    webbrowser.get(chrome_path).open(url)

























#%%
## https://learn.astropy.org/FITS-images.html
#from astropy.io import fits
#image_file = '/Users/lester/BEAGLE/FF_images/0A2744C/hlsp_frontier_hst_acs-30mas_abell2744_f435w_v1.0_drz.fits'
#hdu_list = fits.open(image_file)
#hdu_list.info()
#image_data = hdu_list[0].data
#print(type(image_data))
#print(image_data.shape)
#hdu_list.close()
#
#image_data = fits.getdata(image_file)
#print(type(image_data))
#print(image_data.shape)
#
#plt.imshow(image_data, cmap='gray')
#plt.colorbar()
#
#print('Min:', np.min(image_data))
#print('Max:', np.max(image_data))
#print('Mean:', np.mean(image_data))
#print('Stdev:', np.std(image_data))
#
#print(type(image_data.flatten()))
#
#NBINS = 1000
#histogram = plt.hist(image_data.flatten(), NBINS)

#from matplotlib.colors import LogNorm
#plt.imshow(image_data, cmap='gray', norm=LogNorm())
#
## I chose the tick marks based on the histogram above
#cbar = plt.colorbar(ticks=[5.e3,1.e4,2.e4])
#cbar.ax.set_yticklabels(['5,000','10,000','20,000'])



## https://docs.astropy.org/en/stable/wcs/
#from astropy.io import fits
#from astropy.wcs import WCS
#from astropy.utils.data import get_pkg_data_filename
#w = WCS(image_file)
#pix = w.wcs_world2pix(3.600050, -30.389715, 0)
#print(pix[0], pix[1])
#
## need +- 400

## https://docs.astropy.org/en/stable/nddata/utils.html
#from astropy.nddata import Cutout2D
#from astropy import units as u
#position = (pix[0], pix[1])
#size = (800, 800)     # pixels
#cutout = Cutout2D(image_data, position, size)
#plt.imshow(cutout.data, origin='lower', vmin=0, vmax=np.percentile(cutout.data.flatten(), 99.5))
#plt.colorbar()
#plt.show()
#
#plt.hist(cutout.data.flatten(), 100)
#print(np.percentile(cutout.data.flatten(), 10))
#































