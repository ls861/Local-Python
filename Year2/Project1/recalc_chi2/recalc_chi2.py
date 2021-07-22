#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 22:39:19 2021

@author: lester
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
import os


MRE = [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.05, 0.10, 0.10]
#MRE = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.10, 0.10]

fields = ['A2744_c', 'A2744_p', 'M0416_c', 'M0416_p', 'M0717_c', 'M0717_p', 'M1149_c', 'M1149_p']
fields2 = ['0A2744C', '1A2744P', '2M0416C', '3M0416P', '4M0717C', '5M0717P', '6M1149C', '7M1149P']
input_filters =  ['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2']
input_err_filters =  ['b_errB435', 'b_errV606', 'b_errI814', 'b_errY105', 'b_errJ125', 'b_errJH140', 'b_errH160', 'b_errKs', 'b_errCH1', 'b_errCH2']

output_filters = ['HST_ACS_WFC_F435W_APP', 'HST_ACS_WFC_F606W_APP', 'HST_ACS_WFC_F814W_APP', 'HST_WFC3_IR_F105W_APP', 'HST_WFC3_IR_F125W_APP', 'HST_WFC3_IR_F140W_APP', 'HST_WFC3_IR_F160W_APP', 'Paranal_HAWKI_Ks_APP', 'Spitzer_IRAC_I1_APP', 'Spitzer_IRAC_I2_APP']

for field_id, field in enumerate(fields):

#    inputFolder = '/Users/lester/BEAGLE/BEAGLE-general/data/astrodeep_{}_subset_RF1_001.fits'.format(field)
    inputFolder = '/home/ls861/BEAGLE-general/data/astrodeep_{}_subset_RF1_001.fits'.format(field)
    input_fits = fits.open(inputFolder)
    #print(input_fits[1].header)
#    outputFolder = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_3/test/'
    outputFolder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/fit_001/'.format(fields2[field_id])

    id_arr = []
    old_chi2_arr = []
    new_chi2_arr = []

    # go through .gz files
    fileList = os.listdir(outputFolder)
    for file in fileList:
        if '.fits.gz' in file:
            #print(file)
            new_chi2 = 0.0
            output_fits = fits.open(outputFolder+file)
            id_BEAGLE = int(file.replace('_BEAGLE.fits.gz',''))

            # get gz row with min chi2
            chi2_arr_temp = output_fits['POSTERIOR PDF'].data['chi_square']
            chi2_idx = np.isin(chi2_arr_temp, np.amin(chi2_arr_temp))
#            chi2_idx = np.where(chi2_arr_temp==np.amin(chi2_arr_temp))[0][0]
            old_chi2 = float(chi2_arr_temp[chi2_idx])

            # go to input file and get fluxes for bands with +ve error
            for f, input_filter in enumerate(input_filters):
                id_idx = np.isin(input_fits[1].data['ID'], id_BEAGLE)
                input_err_flux_temp = float(input_fits[1].data[input_err_filters[f]][id_idx])

                if input_err_flux_temp > 0:
                    input_flux_temp = float(input_fits[1].data[input_filter][id_idx])
                    MRE_err_flux_temp = ((input_err_flux_temp**2)+((MRE[f]*input_flux_temp)**2))**0.5

                    # find corresponding apparent magnitudes (for bands which had +ve error) and convert to fluxes
                    #print(output_fits['APPARENT MAGNITUDES'].header)
                    output_mag_temp = float(output_fits['APPARENT MAGNITUDES'].data[output_filters[f]][chi2_idx])
                    output_flux_temp = (10.0 ** ( (output_mag_temp - 8.90) / (-2.5) )) * 1e6 # uJy
                    #print(float(input_flux_temp), float(output_flux_temp))

                    # calculate new chi2
                    new_chi2 += ((input_flux_temp - output_flux_temp)**2) / (MRE_err_flux_temp**2)

            output_fits.close()

            id_arr.append(id_BEAGLE)
            old_chi2_arr.append(old_chi2)
            new_chi2_arr.append(new_chi2)

    input_fits.close()

    sort_idx = np.argsort(id_arr)
    id_arr = np.array(id_arr)[sort_idx]
    old_chi2_arr = np.array(old_chi2_arr)[sort_idx]
    new_chi2_arr = np.array(new_chi2_arr)[sort_idx]

#    print(id_arr, old_chi2_arr, new_chi2_arr)
    outputDict = {}
    outputDict['id']                = id_arr
    outputDict['old_chi2']          = old_chi2_arr
    outputDict['new_chi2']          = new_chi2_arr

    outputTable = Table(outputDict)
    #os.system("mkdir ./chi2") # should already exist for the original chi2
    outputTable.write("./chi2/{}_new_chi2.fits".format(field_id), overwrite=True)
