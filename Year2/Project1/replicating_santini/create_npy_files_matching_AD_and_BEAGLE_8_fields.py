#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""
# =============================================================================
# NOTES
# =============================================================================

'''
ALL THE FOLLOWING NPY FILES SHOULD BE EQUAL IN LENGTH:
    
field_AD
id_AD
mag_AD
redshift_AD
mass_AD         MAGNIFICATION ADJUSTED
mass_AD_neb     MAGNIFICATION ADJUSTED
sfr_AD          MAGNIFICATION ADJUSTED
sfr_AD_neb      MAGNIFICATION ADJUSTED
sfr_SAN         MAGNIFICATION ADJUSTED

id_BEAGLE
mass_BEAGLE_tot         NOT LOG
mass_BEAGLE_stellar     NOT LOG
sfr_BEAGLE_instant
redshift_BEAGLE

NOTE sfr_SAN has -99.0 when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)

AND BEAGLE has -101.0 when the object was not a BEAGLE input and -102.0 if the object was a BEAGLE input, but was NOT fitted

AND AD masses (normal + neb) set to 0 are now -40 in log space
AND AD sfr and redshift set to 0 are now -40 in log space
'''




import numpy as np
from astropy.io import fits
import pickle

AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'
sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy'

AD = np.load(AD_location)
sfr_SAN = np.load(sfr_SAN_location)
#print(AD.dtype.names)
#print(len(AD))

# =============================================================================
# get some key columns
# =============================================================================
field_AD = AD['field']
id_AD = AD['ID']
mag_AD = AD['MAGNIF']
redshift_AD = AD['ZBEST']

with np.errstate(divide='ignore'):
    mass_AD = np.log10(AD['MSTAR']*1e9) - np.log10(mag_AD)
    mass_AD_neb = np.log10(AD['MASTAR_NEB']*1e9) - np.log10(mag_AD)
    sfr_AD = np.log10(AD['SFR']) - np.log10(mag_AD)
    sfr_AD_neb = np.log10(AD['SFR_NEB']) - np.log10(mag_AD)    
    
mass_AD[np.isneginf(mass_AD)]=-40
mass_AD_neb[np.isneginf(mass_AD_neb)]=-40
sfr_AD[np.isneginf(sfr_AD)]=-40
sfr_AD_neb[np.isneginf(sfr_AD_neb)]=-40

sfr_SAN = sfr_SAN - np.log10(mag_AD)

# NO BEAGLE OUTPUTS USED YET
field_AD = np.array(field_AD)
id_AD = np.array(id_AD)
mag_AD = np.array(mag_AD)
redshift_AD = np.array(redshift_AD)
mass_AD = np.array(mass_AD)
mass_AD_neb = np.array(mass_AD_neb)
sfr_AD = np.array(sfr_AD)
sfr_AD_neb = np.array(sfr_AD_neb)
sfr_SAN = np.array(sfr_SAN)

# =============================================================================
# BEAGLE OUTPUTS
# =============================================================================
id_BEAGLE = []
mass_BEAGLE_tot = []
mass_BEAGLE_stellar = []
sfr_BEAGLE_instant = []
redshift_BEAGLE = []

for i, field in enumerate(field_AD):

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/report_BEAGLE_inputs/BEAGLE_input_{}.fits'.format(int(field))
    data_fits = fits.open(fileName)
#    print(data_fits.info())
#    print(data_fits[1].header)
    id_input = np.asarray(data_fits[1].data['ID'], dtype=int)
    field_original = np.asarray(data_fits[1].data['field'], dtype=int)
    id_original = np.asarray(data_fits[1].data['ID_original'], dtype=int)
    data_fits.close()

#    GET BEAGLE ID FROM ORIGINAL

    if len(np.where(id_original==id_AD[i])[0]) == 0:
        print('{} {} AD OBJECT WAS NOT A BEAGLE INPUT'.format(int(field), id_AD[i]))
        id_BEAGLE.append(-101)
        mass_BEAGLE_tot.append(-101.0)
        mass_BEAGLE_stellar.append(-101.0)
        redshift_BEAGLE.append(-101.0)     
        sfr_BEAGLE_instant.append(-101.0)
    else:
        id_BEAGLE.append(id_input[np.where(id_original==id_AD[i])[0][0]])
        
        fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/report_BEAGLE_summary_catalogues/BEAGLE_summary_catalogue_{}.fits'.format(int(field))
        data_fits = fits.open(fileName)
    #    print(data_fits.info())
    #    print(data_fits[2].header)
        id_input = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
        
    #    GET BEAGLE MASS, SFR and Z
        
        if len(np.where(id_input==id_BEAGLE[i])[0]) == 0:
            print('{} {} AD OBJECT WAS NOT FITTED BY BEAGLE'.format(int(field), id_AD[i]))
            mass_BEAGLE_tot.append(-102.0)
            mass_BEAGLE_stellar.append(-102.0)
            redshift_BEAGLE.append(-102.0)     
            sfr_BEAGLE_instant.append(-102.0)
            data_fits.close()
        else:  
            idx = np.where(id_input==id_BEAGLE[i])[0][0]
            mass_BEAGLE_tot.append(data_fits['GALAXY PROPERTIES'].data['M_tot_median'][idx])
            mass_BEAGLE_stellar.append(data_fits['GALAXY PROPERTIES'].data['M_star_median'][idx])
            redshift_BEAGLE.append(data_fits['GALAXY PROPERTIES'].data['redshift_median'][idx])
            data_fits.close()
        
            data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/report_BEAGLE_instant_SFRs/{}_instant_sfr_medians.p'.format(int(field)),'r'))
        #    print(data.keys())
            idx_sort = np.argsort(np.asarray(data['ID'], float))
            idx = np.where(np.asarray(data['ID'], float)[idx_sort]==id_BEAGLE[i])[0][0]
            sfr_BEAGLE_instant.append(data['log_instant_sfr_median.npy'][idx_sort][idx])
        
#id_BEAGLE = np.array(id_BEAGLE)
#mass_BEAGLE_tot = np.array(mass_BEAGLE_tot)
#mass_BEAGLE_stellar = np.array(mass_BEAGLE_stellar)
#sfr_BEAGLE_instant = np.array(sfr_BEAGLE_instant)
#redshift_BEAGLE = np.array(redshift_BEAGLE)
#
#np.save('id_BEAGLE',id_BEAGLE)
#np.save('mass_BEAGLE_tot',mass_BEAGLE_tot)
#np.save('mass_BEAGLE_stellar',mass_BEAGLE_stellar)
#np.save('sfr_BEAGLE_instant',sfr_BEAGLE_instant)
#np.save('redshift_BEAGLE',redshift_BEAGLE)
#
#np.save('field_AD',field_AD)
#np.save('id_AD',id_AD)
#np.save('mag_AD',mag_AD)
#np.save('redshift_AD',redshift_AD)
#np.save('mass_AD',mass_AD)
#np.save('mass_AD_neb',mass_AD_neb)
#np.save('sfr_AD',sfr_AD)
#np.save('sfr_AD_neb',sfr_AD_neb)
#np.save('sfr_SAN',sfr_SAN)
            





