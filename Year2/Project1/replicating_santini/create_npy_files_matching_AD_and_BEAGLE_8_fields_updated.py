#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""

import numpy as np
from astropy.io import fits
import pickle

AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
#print(AD.dtype.names)
#print(len(AD))

field_AD = np.array(AD['field'])
id_AD = np.array(AD['ID'])

# =============================================================================
# BEAGLE OUTPUTS
# =============================================================================
id_BEAGLE = []
mass_BEAGLE_tot = []
mass_BEAGLE_stellar = []
sfr_BEAGLE_instant = []
redshift_BEAGLE = []

tau_BEAGLE = []
tauv_BEAGLE = []
msa_BEAGLE = []
metallicity_BEAGLE = []

min_chi2_BEAGLE = []

# =============================================================================
# GMM outputs
# =============================================================================

id_GMM = []
x_GMM = []
y_GMM = []
xsig_GMM = []
ysig_GMM = []
xycov_GMM = []
amp_GMM = []


#for i, field in enumerate(field_AD[:10]):
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
        
        tau_BEAGLE.append(-101.0)
        tauv_BEAGLE.append(-101.0)
        msa_BEAGLE.append(-101.0)
        metallicity_BEAGLE.append(-101.0)
        
        min_chi2_BEAGLE.append(-101.0)
        
        # GMM
        id_GMM.append(-101.0)
        x_GMM.append(np.array([-101.0,-101.0,-101.0]))
        y_GMM.append(np.array([-101.0,-101.0,-101.0]))
        xsig_GMM.append(np.array([-101.0,-101.0,-101.0]))
        ysig_GMM.append(np.array([-101.0,-101.0,-101.0]))
        xycov_GMM.append(np.array([-101.0,-101.0,-101.0]))
        amp_GMM.append(np.array([-101.0,-101.0,-101.0]))
     
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
            data_fits.close()
            
            mass_BEAGLE_tot.append(-102.0)
            mass_BEAGLE_stellar.append(-102.0)
            redshift_BEAGLE.append(-102.0)     
            sfr_BEAGLE_instant.append(-102.0)

            tau_BEAGLE.append(-102.0)
            tauv_BEAGLE.append(-102.0)
            msa_BEAGLE.append(-102.0)
            metallicity_BEAGLE.append(-102.0)
            
            min_chi2_BEAGLE.append(-102.0)
            
            # GMM
            id_GMM.append(-102.0)
            x_GMM.append(np.array([-102.0,-102.0,-102.0]))
            y_GMM.append(np.array([-102.0,-102.0,-102.0]))
            xsig_GMM.append(np.array([-102.0,-102.0,-102.0]))
            ysig_GMM.append(np.array([-102.0,-102.0,-102.0]))
            xycov_GMM.append(np.array([-102.0,-102.0,-102.0]))
            amp_GMM.append(np.array([-102.0,-102.0,-102.0]))
        
        else:  
            idx = np.where(id_input==id_BEAGLE[i])[0][0]
            mass_BEAGLE_tot.append(data_fits['GALAXY PROPERTIES'].data['M_tot_median'][idx])
            mass_BEAGLE_stellar.append(data_fits['GALAXY PROPERTIES'].data['M_star_median'][idx])
            redshift_BEAGLE.append(data_fits['GALAXY PROPERTIES'].data['redshift_median'][idx])
            
            tau_BEAGLE.append(data_fits['POSTERIOR PDF'].data['tau_median'][idx])
            tauv_BEAGLE.append(data_fits['POSTERIOR PDF'].data['tauv_eff_median'][idx])
            msa_BEAGLE.append(data_fits['POSTERIOR PDF'].data['max_stellar_age_median'][idx])
            metallicity_BEAGLE.append(data_fits['POSTERIOR PDF'].data['metallicity_median'][idx])
        
            data_fits.close()
            
            # CHI2
            fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/from_cluster_chi2/{}_chi2.fits'.format(int(field))
            data_fits = fits.open(fileName)
            #print(data_fits.info())
            #print(data_fits[1].header)
            #id_chi2 = np.asarray(data_fits[1].data['id'], dtype=int)
            chi2 = data_fits[1].data['chi2']
            data_fits.close()
            
            min_chi2_BEAGLE.append(chi2[idx])
            
            # INSTANT SFR
            data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/report_BEAGLE_instant_SFRs/{}_instant_sfr_medians.p'.format(int(field)),'r'))
        #    print(data.keys())
            idx_sort = np.argsort(np.asarray(data['ID'], float))
            idx = np.where(np.asarray(data['ID'], float)[idx_sort]==id_BEAGLE[i])[0][0]
            sfr_BEAGLE_instant.append(data['log_instant_sfr_median.npy'][idx_sort][idx])
            
            # GMM
            data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/GMM/kelly_2d_GMM_inputs/{}_mStar_delayed_GMM_2d_fits.p'.format(int(field)),'r'))
        #    print(data.keys())
            idx_sort = np.argsort(np.asarray(data['id'], float))
            idx = np.where(np.asarray(data['id'], float)[idx_sort]==id_BEAGLE[i])[0][0]

            id_GMM.append(float(data['id'][idx_sort][idx]))
            x_GMM.append(data['x'][idx_sort][idx])
            y_GMM.append(data['y'][idx_sort][idx])
            xsig_GMM.append(data['xsig'][idx_sort][idx])
            ysig_GMM.append(data['ysig'][idx_sort][idx])
            xycov_GMM.append(data['xycov'][idx_sort][idx])
            amp_GMM.append(data['amp'][idx_sort][idx])
        
id_BEAGLE = np.array(id_BEAGLE)
mass_BEAGLE_tot = np.array(mass_BEAGLE_tot)
mass_BEAGLE_stellar = np.array(mass_BEAGLE_stellar)
sfr_BEAGLE_instant = np.array(sfr_BEAGLE_instant)
redshift_BEAGLE = np.array(redshift_BEAGLE)

tau_BEAGLE = np.array(tau_BEAGLE)
tauv_BEAGLE = np.array(tauv_BEAGLE)
msa_BEAGLE = np.array(msa_BEAGLE)
metallicity_BEAGLE = np.array(metallicity_BEAGLE)

min_chi2_BEAGLE = np.array(min_chi2_BEAGLE)

# GMM
id_GMM = np.array(id_GMM)
x_GMM = np.array(x_GMM)
y_GMM = np.array(y_GMM)
xsig_GMM = np.array(xsig_GMM)
ysig_GMM = np.array(ysig_GMM)
xycov_GMM = np.array(xycov_GMM)
amp_GMM = np.array(amp_GMM)

# =============================================================================
# SAVE
# =============================================================================
#np.save('id_BEAGLE',id_BEAGLE)
#np.save('mass_BEAGLE_tot',mass_BEAGLE_tot)
#np.save('mass_BEAGLE_stellar',mass_BEAGLE_stellar)
#np.save('sfr_BEAGLE_instant',sfr_BEAGLE_instant)
#np.save('redshift_BEAGLE',redshift_BEAGLE)
#
#np.save('tau_BEAGLE',tau_BEAGLE)
#np.save('tauv_BEAGLE',tauv_BEAGLE)
#np.save('msa_BEAGLE',msa_BEAGLE)
#np.save('metallicity_BEAGLE',metallicity_BEAGLE)
#
#np.save('min_chi2_BEAGLE',min_chi2_BEAGLE)
#
###GMM
#np.save('id_GMM',id_GMM)
#np.save('x_GMM',x_GMM)
#np.save('y_GMM',y_GMM)
#np.save('xsig_GMM',xsig_GMM)
#np.save('ysig_GMM',ysig_GMM)
#np.save('xycov_GMM',xycov_GMM)
#np.save('amp_GMM',amp_GMM)



#print(id_GMM)
#print(x_GMM)
#print(y_GMM)
#print(xsig_GMM)
#print(ysig_GMM)
#print(xycov_GMM)
#print(amp_GMM)
#
#
#
#print(id_BEAGLE)
#print(mass_BEAGLE_stellar)


