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
redshift_BEAGLE_mean = []

tau_BEAGLE = []
tauv_BEAGLE = []
msa_BEAGLE = []
metallicity_BEAGLE = []

min_chi2_BEAGLE = []
new_min_chi2_BEAGLE = []

Ks = []
CH1 = []
CH2 = []

ch1_beagle_mag_median = []
ch1_beagle_mag_lower = []
ch1_beagle_mag_upper = []

ch2_beagle_mag_median = []
ch2_beagle_mag_lower = []
ch2_beagle_mag_upper = []


# =============================================================================
# GMM outputs
# =============================================================================

id_GMM_2d = []
x_GMM_2d = []
y_GMM_2d = []
xsig_GMM_2d = []
ysig_GMM_2d = []
xycov_GMM_2d = []
amp_GMM_2d = []

id_GMM_3d = []
x_GMM_3d = []
y_GMM_3d = []
z_GMM_3d = []
xsig_GMM_3d = []
ysig_GMM_3d = []
zsig_GMM_3d = []
xycov_GMM_3d = []
xzcov_GMM_3d = []
yzcov_GMM_3d = []
amp_GMM_3d = []

cut = 0
cut_i = 0
#for i, field in enumerate(field_AD[:10]):
for i, field in enumerate(field_AD):

    if field > cut:
        cut_i = i

        id_BEAGLE = np.array(id_BEAGLE)
        mass_BEAGLE_tot = np.array(mass_BEAGLE_tot)
        mass_BEAGLE_stellar = np.array(mass_BEAGLE_stellar)
        sfr_BEAGLE_instant = np.array(sfr_BEAGLE_instant)
        redshift_BEAGLE = np.array(redshift_BEAGLE)
        redshift_BEAGLE_mean = np.array(redshift_BEAGLE_mean)
        
        tau_BEAGLE = np.array(tau_BEAGLE)
        tauv_BEAGLE = np.array(tauv_BEAGLE)
        msa_BEAGLE = np.array(msa_BEAGLE)
        metallicity_BEAGLE = np.array(metallicity_BEAGLE)

        min_chi2_BEAGLE = np.array(min_chi2_BEAGLE)
        new_min_chi2_BEAGLE = np.array(new_min_chi2_BEAGLE)

        Ks = np.array(Ks)
        CH1 = np.array(CH1)
        CH2 = np.array(CH2)

        ch1_beagle_mag_median = np.array(ch1_beagle_mag_median)
        ch1_beagle_mag_lower = np.array(ch1_beagle_mag_lower)
        ch1_beagle_mag_upper = np.array(ch1_beagle_mag_upper)
        
        ch2_beagle_mag_median = np.array(ch2_beagle_mag_median)
        ch2_beagle_mag_lower = np.array(ch2_beagle_mag_lower)
        ch2_beagle_mag_upper = np.array(ch2_beagle_mag_upper)

        # GMM
        id_GMM_2d = np.array(id_GMM_2d)
        x_GMM_2d = np.array(x_GMM_2d)
        y_GMM_2d = np.array(y_GMM_2d)
        xsig_GMM_2d = np.array(xsig_GMM_2d)
        ysig_GMM_2d = np.array(ysig_GMM_2d)
        xycov_GMM_2d = np.array(xycov_GMM_2d)
        amp_GMM_2d = np.array(amp_GMM_2d)

        id_GMM_3d = np.array(id_GMM_3d)
        x_GMM_3d = np.array(x_GMM_3d)
        y_GMM_3d = np.array(y_GMM_3d)
        z_GMM_3d = np.array(z_GMM_3d)
        xsig_GMM_3d = np.array(xsig_GMM_3d)
        ysig_GMM_3d = np.array(ysig_GMM_3d)
        zsig_GMM_3d = np.array(zsig_GMM_3d)
        xycov_GMM_3d = np.array(xycov_GMM_3d)
        xzcov_GMM_3d = np.array(xzcov_GMM_3d)
        yzcov_GMM_3d = np.array(yzcov_GMM_3d)
        amp_GMM_3d = np.array(amp_GMM_3d)

        # =============================================================================
        # SAVE
        # =============================================================================
        np.save('id_BEAGLE_{}'.format(cut),id_BEAGLE)
        np.save('mass_BEAGLE_tot_{}'.format(cut),mass_BEAGLE_tot)
        np.save('mass_BEAGLE_stellar_{}'.format(cut),mass_BEAGLE_stellar)
        np.save('sfr_BEAGLE_instant_{}'.format(cut),sfr_BEAGLE_instant)
        np.save('redshift_BEAGLE_{}'.format(cut),redshift_BEAGLE)
        np.save('redshift_BEAGLE_mean_{}'.format(cut),redshift_BEAGLE_mean)
        
        np.save('tau_BEAGLE_{}'.format(cut),tau_BEAGLE)
        np.save('tauv_BEAGLE_{}'.format(cut),tauv_BEAGLE)
        np.save('msa_BEAGLE_{}'.format(cut),msa_BEAGLE)
        np.save('metallicity_BEAGLE_{}'.format(cut),metallicity_BEAGLE)

        np.save('min_chi2_BEAGLE_{}'.format(cut),min_chi2_BEAGLE)
        np.save('new_min_chi2_BEAGLE_{}'.format(cut),new_min_chi2_BEAGLE)

        np.save('Ks_{}'.format(cut),Ks)
        np.save('CH1_{}'.format(cut),CH1)
        np.save('CH2_{}'.format(cut),CH2)

        np.save('ch1_beagle_mag_median_{}'.format(cut),ch1_beagle_mag_median)
        np.save('ch1_beagle_mag_lower_{}'.format(cut),ch1_beagle_mag_lower)
        np.save('ch1_beagle_mag_upper_{}'.format(cut),ch1_beagle_mag_upper)

        np.save('ch2_beagle_mag_median_{}'.format(cut),ch2_beagle_mag_median)
        np.save('ch2_beagle_mag_lower_{}'.format(cut),ch2_beagle_mag_lower)
        np.save('ch2_beagle_mag_upper_{}'.format(cut),ch2_beagle_mag_upper)

        ##GMM
        np.save('id_GMM_2d_{}'.format(cut),id_GMM_2d)
        np.save('x_GMM_2d_{}'.format(cut),x_GMM_2d)
        np.save('y_GMM_2d_{}'.format(cut),y_GMM_2d)
        np.save('xsig_GMM_2d_{}'.format(cut),xsig_GMM_2d)
        np.save('ysig_GMM_2d_{}'.format(cut),ysig_GMM_2d)
        np.save('xycov_GMM_2d_{}'.format(cut),xycov_GMM_2d)
        np.save('amp_GMM_2d_{}'.format(cut),amp_GMM_2d)

        np.save('id_GMM_3d_{}'.format(cut),id_GMM_3d)
        np.save('x_GMM_3d_{}'.format(cut),x_GMM_3d)
        np.save('y_GMM_3d_{}'.format(cut),y_GMM_3d)
        np.save('z_GMM_3d_{}'.format(cut),z_GMM_3d)
        np.save('xsig_GMM_3d_{}'.format(cut),xsig_GMM_3d)
        np.save('ysig_GMM_3d_{}'.format(cut),ysig_GMM_3d)
        np.save('zsig_GMM_3d_{}'.format(cut),zsig_GMM_3d)
        np.save('xycov_GMM_3d_{}'.format(cut),xycov_GMM_3d)
        np.save('xzcov_GMM_3d_{}'.format(cut),xzcov_GMM_3d)
        np.save('yzcov_GMM_3d_{}'.format(cut),yzcov_GMM_3d)
        np.save('amp_GMM_3d_{}'.format(cut),amp_GMM_3d)

        # =============================================================================
        # BEAGLE OUTPUTS
        # =============================================================================

        id_BEAGLE = []
        mass_BEAGLE_tot = []
        mass_BEAGLE_stellar = []
        sfr_BEAGLE_instant = []
        redshift_BEAGLE = []
        redshift_BEAGLE_mean = []
        
        tau_BEAGLE = []
        tauv_BEAGLE = []
        msa_BEAGLE = []
        metallicity_BEAGLE = []

        min_chi2_BEAGLE = []
        new_min_chi2_BEAGLE = []
        
        Ks = []
        CH1 = []
        CH2 = []
        
        ch1_beagle_mag_median = []
        ch1_beagle_mag_lower = []
        ch1_beagle_mag_upper = []
        
        ch2_beagle_mag_median = []
        ch2_beagle_mag_lower = []
        ch2_beagle_mag_upper = []

        # =============================================================================
        # GMM outputs
        # =============================================================================

        id_GMM_2d = []
        x_GMM_2d = []
        y_GMM_2d = []
        xsig_GMM_2d = []
        ysig_GMM_2d = []
        xycov_GMM_2d = []
        amp_GMM_2d = []

        id_GMM_3d = []
        x_GMM_3d = []
        y_GMM_3d = []
        z_GMM_3d = []
        xsig_GMM_3d = []
        ysig_GMM_3d = []
        zsig_GMM_3d = []
        xycov_GMM_3d = []
        xzcov_GMM_3d = []
        yzcov_GMM_3d = []
        amp_GMM_3d = []

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/report_BEAGLE_inputs/astrodeep_field_{}_subset_RF1_001.fits'.format(int(field))
    data_fits = fits.open(fileName)
#    print(data_fits.info())
#    print(data_fits[1].header)
    id_input = np.asarray(data_fits[1].data['ID'], dtype=int)
    field_original = np.asarray(data_fits[1].data['field'], dtype=int)
    id_original = np.asarray(data_fits[1].data['ID_original'], dtype=int)
    b_Ks = np.asarray(data_fits[1].data['b_Ks'], dtype=float)
    b_CH1 = np.asarray(data_fits[1].data['b_CH1'], dtype=float)
    b_CH2 = np.asarray(data_fits[1].data['b_CH2'], dtype=float)
    data_fits.close()

#    GET BEAGLE ID FROM ORIGINAL

    if len(np.where(id_original==id_AD[i])[0]) == 0:
        print('{} {} AD OBJECT WAS NOT A BEAGLE INPUT'.format(int(field), id_AD[i]))
        id_BEAGLE.append(-101)
        mass_BEAGLE_tot.append(-101.0)
        mass_BEAGLE_stellar.append(-101.0)
        redshift_BEAGLE.append(-101.0)
        redshift_BEAGLE_mean.append(-101.0)
        sfr_BEAGLE_instant.append(-101.0)

        tau_BEAGLE.append(-101.0)
        tauv_BEAGLE.append(-101.0)
        msa_BEAGLE.append(-101.0)
        metallicity_BEAGLE.append(-101.0)

        min_chi2_BEAGLE.append(-101.0)
        new_min_chi2_BEAGLE.append(-101.0)
        
        Ks.append(-101.0)
        CH1.append(-101.0)
        CH2.append(-101.0)
        
        ch1_beagle_mag_median.append(-101.0)
        ch1_beagle_mag_lower.append(-101.0)
        ch1_beagle_mag_upper.append(-101.0)
        
        ch2_beagle_mag_median.append(-101.0)
        ch2_beagle_mag_lower.append(-101.0)
        ch2_beagle_mag_upper.append(-101.0)

        # GMM
        id_GMM_2d.append(-101.0)
        x_GMM_2d.append(np.array([-101.0,-101.0,-101.0]))
        y_GMM_2d.append(np.array([-101.0,-101.0,-101.0]))
        xsig_GMM_2d.append(np.array([-101.0,-101.0,-101.0]))
        ysig_GMM_2d.append(np.array([-101.0,-101.0,-101.0]))
        xycov_GMM_2d.append(np.array([-101.0,-101.0,-101.0]))
        amp_GMM_2d.append(np.array([-101.0,-101.0,-101.0]))

        id_GMM_3d.append(-101.0)
        x_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))
        y_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))
        z_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))
        xsig_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))
        ysig_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))
        zsig_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))
        xycov_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))
        xzcov_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))
        yzcov_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))
        amp_GMM_3d.append(np.array([-101.0,-101.0,-101.0]))

    else:
        id_BEAGLE.append(id_input[np.where(id_original==id_AD[i])[0][0]])
        Ks.append(b_Ks[np.where(id_original==id_AD[i])[0][0]])
        CH1.append(b_CH1[np.where(id_original==id_AD[i])[0][0]])
        CH2.append(b_CH2[np.where(id_original==id_AD[i])[0][0]])

        fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/from_cluster/summary_catalogues/{}_BEAGLE_summary_catalogue.fits'.format(int(field))

        data_fits = fits.open(fileName)
    #    print(data_fits.info())
    #    print(data_fits[2].header)
        id_input = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)

    #    GET BEAGLE MASS, SFR and Z

        if len(np.where(id_input==id_BEAGLE[i-cut_i])[0]) == 0:
            print('{} {} AD OBJECT WAS NOT FITTED BY BEAGLE'.format(int(field), id_AD[i]))
            data_fits.close()

            mass_BEAGLE_tot.append(-102.0)
            mass_BEAGLE_stellar.append(-102.0)
            redshift_BEAGLE.append(-102.0)
            redshift_BEAGLE_mean.append(-102.0)
            sfr_BEAGLE_instant.append(-102.0)

            tau_BEAGLE.append(-102.0)
            tauv_BEAGLE.append(-102.0)
            msa_BEAGLE.append(-102.0)
            metallicity_BEAGLE.append(-102.0)

            min_chi2_BEAGLE.append(-102.0)
            new_min_chi2_BEAGLE.append(-102.0)
            
            ch1_beagle_mag_median.append(-102.0)
            ch1_beagle_mag_lower.append(-102.0)
            ch1_beagle_mag_upper.append(-102.0)
            
            ch2_beagle_mag_median.append(-102.0)
            ch2_beagle_mag_lower.append(-102.0)
            ch2_beagle_mag_upper.append(-102.0)

            # GMM
            id_GMM_2d.append(-102.0)
            x_GMM_2d.append(np.array([-102.0,-102.0,-102.0]))
            y_GMM_2d.append(np.array([-102.0,-102.0,-102.0]))
            xsig_GMM_2d.append(np.array([-102.0,-102.0,-102.0]))
            ysig_GMM_2d.append(np.array([-102.0,-102.0,-102.0]))
            xycov_GMM_2d.append(np.array([-102.0,-102.0,-102.0]))
            amp_GMM_2d.append(np.array([-102.0,-102.0,-102.0]))

            id_GMM_3d.append(-102.0)
            x_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))
            y_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))
            z_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))
            xsig_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))
            ysig_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))
            zsig_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))
            xycov_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))
            xzcov_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))
            yzcov_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))
            amp_GMM_3d.append(np.array([-102.0,-102.0,-102.0]))

        else:
            idx = np.where(id_input==id_BEAGLE[i-cut_i])[0][0]
            mass_BEAGLE_tot.append(data_fits['GALAXY PROPERTIES'].data['M_tot_median'][idx])
            mass_BEAGLE_stellar.append(data_fits['GALAXY PROPERTIES'].data['M_star_median'][idx])
            redshift_BEAGLE.append(data_fits['GALAXY PROPERTIES'].data['redshift_median'][idx])
            redshift_BEAGLE_mean.append(data_fits['GALAXY PROPERTIES'].data['redshift_mean'][idx])
            
            tau_BEAGLE.append(data_fits['POSTERIOR PDF'].data['tau_median'][idx])
            tauv_BEAGLE.append(data_fits['POSTERIOR PDF'].data['tauv_eff_median'][idx])
            msa_BEAGLE.append(data_fits['POSTERIOR PDF'].data['max_stellar_age_median'][idx])
            metallicity_BEAGLE.append(data_fits['POSTERIOR PDF'].data['metallicity_median'][idx])

            ch1_beagle_mag_median.append(data_fits['APPARENT MAGNITUDES'].data['Spitzer_IRAC_I1_APP_median'][idx])
            ch1_beagle_mag_lower.append(data_fits['APPARENT MAGNITUDES'].data['Spitzer_IRAC_I1_APP_68.00'][:,1][idx])
            ch1_beagle_mag_upper.append(data_fits['APPARENT MAGNITUDES'].data['Spitzer_IRAC_I1_APP_68.00'][:,0][idx])
            
            ch2_beagle_mag_median.append(data_fits['APPARENT MAGNITUDES'].data['Spitzer_IRAC_I2_APP_median'][idx])
            ch2_beagle_mag_lower.append(data_fits['APPARENT MAGNITUDES'].data['Spitzer_IRAC_I2_APP_68.00'][:,1][idx])
            ch2_beagle_mag_upper.append(data_fits['APPARENT MAGNITUDES'].data['Spitzer_IRAC_I2_APP_68.00'][:,0][idx])
        
            data_fits.close()

            # CHI2
            fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/for_cluster/slurm/BEAGLE_dust_install/chi2/{}_chi2.fits'.format(int(field))

            data_fits = fits.open(fileName)
            #print(data_fits.info())
            #print(data_fits[1].header)
            #id_chi2 = np.asarray(data_fits[1].data['id'], dtype=int)
            chi2 = data_fits[1].data['chi2']
            data_fits.close()

            min_chi2_BEAGLE.append(chi2[idx])

            # NEW CHI2
            fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/for_cluster/slurm/BEAGLE_dust_install/chi2/{}_chi2.fits'.format(int(field))
            data_fits = fits.open(fileName)
            #print(data_fits.info())
            #print(data_fits[1].header)
            #id_chi2 = np.asarray(data_fits[1].data['id'], dtype=int)
            new_chi2 = data_fits[1].data['chi2']
            data_fits.close()

            new_min_chi2_BEAGLE.append(new_chi2[idx])

            # INSTANT SFR
            data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/for_cluster/slurm/BEAGLE_dust_install/instant_sfr_medians/{}_instant_sfr_medians.p'.format(int(field)),'r'))
        #    print(data.keys())
            idx_sort = np.argsort(np.asarray(data['ID'], float))
            idx = np.where(np.asarray(data['ID'], float)[idx_sort]==id_BEAGLE[i-cut_i])[0][0]
            sfr_BEAGLE_instant.append(data['log_instant_sfr_median.npy'][idx_sort][idx])

            # GMM 2d
            data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/for_cluster/slurm/BEAGLE_dust_install/kelly_2d_GMM_inputs/{}_mStar_delayed_GMM_2d_fits.p'.format(int(field)),'r'))
        #    print(data.keys())
            idx_sort = np.argsort(np.asarray(data['id'], float))

            if len(np.where(np.asarray(data['id'], float)[idx_sort]==id_BEAGLE[i-cut_i])[0]) == 0:
                print('{} {} AD OBJECT WAS NOT FITTED BY GMM_2d'.format(int(field), id_AD[i]))
                data_fits.close()
                id_GMM_2d.append(-103.0)
                x_GMM_2d.append(np.array([-103.0,-103.0,-103.0]))
                y_GMM_2d.append(np.array([-103.0,-103.0,-103.0]))
                xsig_GMM_2d.append(np.array([-103.0,-103.0,-103.0]))
                ysig_GMM_2d.append(np.array([-103.0,-103.0,-103.0]))
                xycov_GMM_2d.append(np.array([-103.0,-103.0,-103.0]))
                amp_GMM_2d.append(np.array([-103.0,-103.0,-103.0]))

            else:
                idx = np.where(np.asarray(data['id'], float)[idx_sort]==id_BEAGLE[i-cut_i])[0][0]
                id_GMM_2d.append(float(data['id'][idx_sort][idx]))
                x_GMM_2d.append(data['x'][idx_sort][idx])
                y_GMM_2d.append(data['y'][idx_sort][idx])
                xsig_GMM_2d.append(data['xsig'][idx_sort][idx])
                ysig_GMM_2d.append(data['ysig'][idx_sort][idx])
                xycov_GMM_2d.append(data['xycov'][idx_sort][idx])
                amp_GMM_2d.append(data['amp'][idx_sort][idx])
                data_fits.close()

            # GMM 3d
            data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/for_cluster/slurm/BEAGLE_dust_install/kelly_3d_GMM_inputs/{}_mStar_delayed_GMM_3d_fits.p'.format(int(field)),'r'))
        #    print(data.keys())
            idx_sort = np.argsort(np.asarray(data['id'], float))

            if len(np.where(np.asarray(data['id'], float)[idx_sort]==id_BEAGLE[i-cut_i])[0]) == 0:
                print('{} {} AD OBJECT WAS NOT FITTED BY GMM_3d'.format(int(field), id_AD[i]))
                data_fits.close()
                id_GMM_3d.append(-103.0)
                x_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))
                y_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))
                z_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))
                xsig_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))
                ysig_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))
                zsig_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))
                xycov_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))
                xzcov_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))
                yzcov_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))
                amp_GMM_3d.append(np.array([-103.0,-103.0,-103.0]))

            else:
                idx = np.where(np.asarray(data['id'], float)[idx_sort]==id_BEAGLE[i-cut_i])[0][0]
                id_GMM_3d.append(float(data['id'][idx_sort][idx]))
                x_GMM_3d.append(data['x'][idx_sort][idx])
                y_GMM_3d.append(data['y'][idx_sort][idx])
                z_GMM_3d.append(data['z'][idx_sort][idx])
                xsig_GMM_3d.append(data['xsig'][idx_sort][idx])
                ysig_GMM_3d.append(data['ysig'][idx_sort][idx])
                zsig_GMM_3d.append(data['zsig'][idx_sort][idx])
                xycov_GMM_3d.append(data['xycov'][idx_sort][idx])
                xzcov_GMM_3d.append(data['xzcov'][idx_sort][idx])
                yzcov_GMM_3d.append(data['yzcov'][idx_sort][idx])
                amp_GMM_3d.append(data['amp'][idx_sort][idx])
                data_fits.close()

    cut = int(field)

# =============================================================================
# this is just for final field
# =============================================================================

id_BEAGLE = np.array(id_BEAGLE)
mass_BEAGLE_tot = np.array(mass_BEAGLE_tot)
mass_BEAGLE_stellar = np.array(mass_BEAGLE_stellar)
sfr_BEAGLE_instant = np.array(sfr_BEAGLE_instant)
redshift_BEAGLE = np.array(redshift_BEAGLE)
redshift_BEAGLE_mean = np.array(redshift_BEAGLE_mean)

tau_BEAGLE = np.array(tau_BEAGLE)
tauv_BEAGLE = np.array(tauv_BEAGLE)
msa_BEAGLE = np.array(msa_BEAGLE)
metallicity_BEAGLE = np.array(metallicity_BEAGLE)

min_chi2_BEAGLE = np.array(min_chi2_BEAGLE)
new_min_chi2_BEAGLE = np.array(new_min_chi2_BEAGLE)

Ks = np.array(Ks)
CH1 = np.array(CH1)
CH2 = np.array(CH2)

ch1_beagle_mag_median = np.array(ch1_beagle_mag_median)
ch1_beagle_mag_lower = np.array(ch1_beagle_mag_lower)
ch1_beagle_mag_upper = np.array(ch1_beagle_mag_upper)

ch2_beagle_mag_median = np.array(ch2_beagle_mag_median)
ch2_beagle_mag_lower = np.array(ch2_beagle_mag_lower)
ch2_beagle_mag_upper = np.array(ch2_beagle_mag_upper)
        
# GMM
id_GMM_2d = np.array(id_GMM_2d)
x_GMM_2d = np.array(x_GMM_2d)
y_GMM_2d = np.array(y_GMM_2d)
xsig_GMM_2d = np.array(xsig_GMM_2d)
ysig_GMM_2d = np.array(ysig_GMM_2d)
xycov_GMM_2d = np.array(xycov_GMM_2d)
amp_GMM_2d = np.array(amp_GMM_2d)

id_GMM_3d = np.array(id_GMM_3d)
x_GMM_3d = np.array(x_GMM_3d)
y_GMM_3d = np.array(y_GMM_3d)
z_GMM_3d = np.array(z_GMM_3d)
xsig_GMM_3d = np.array(xsig_GMM_3d)
ysig_GMM_3d = np.array(ysig_GMM_3d)
zsig_GMM_3d = np.array(zsig_GMM_3d)
xycov_GMM_3d = np.array(xycov_GMM_3d)
xzcov_GMM_3d = np.array(xzcov_GMM_3d)
yzcov_GMM_3d = np.array(yzcov_GMM_3d)
amp_GMM_3d = np.array(amp_GMM_3d)

# =============================================================================
# SAVE
# =============================================================================
np.save('id_BEAGLE_{}'.format(cut),id_BEAGLE)
np.save('mass_BEAGLE_tot_{}'.format(cut),mass_BEAGLE_tot)
np.save('mass_BEAGLE_stellar_{}'.format(cut),mass_BEAGLE_stellar)
np.save('sfr_BEAGLE_instant_{}'.format(cut),sfr_BEAGLE_instant)
np.save('redshift_BEAGLE_{}'.format(cut),redshift_BEAGLE)
np.save('redshift_BEAGLE_mean_{}'.format(cut),redshift_BEAGLE_mean)

np.save('tau_BEAGLE_{}'.format(cut),tau_BEAGLE)
np.save('tauv_BEAGLE_{}'.format(cut),tauv_BEAGLE)
np.save('msa_BEAGLE_{}'.format(cut),msa_BEAGLE)
np.save('metallicity_BEAGLE_{}'.format(cut),metallicity_BEAGLE)

np.save('min_chi2_BEAGLE_{}'.format(cut),min_chi2_BEAGLE)
np.save('new_min_chi2_BEAGLE_{}'.format(cut),new_min_chi2_BEAGLE)

np.save('Ks_{}'.format(cut),Ks)
np.save('CH1_{}'.format(cut),CH1)
np.save('CH2_{}'.format(cut),CH2)

np.save('ch1_beagle_mag_median_{}'.format(cut),ch1_beagle_mag_median)
np.save('ch1_beagle_mag_lower_{}'.format(cut),ch1_beagle_mag_lower)
np.save('ch1_beagle_mag_upper_{}'.format(cut),ch1_beagle_mag_upper)

np.save('ch2_beagle_mag_median_{}'.format(cut),ch2_beagle_mag_median)
np.save('ch2_beagle_mag_lower_{}'.format(cut),ch2_beagle_mag_lower)
np.save('ch2_beagle_mag_upper_{}'.format(cut),ch2_beagle_mag_upper)
        
##GMM
np.save('id_GMM_2d_{}'.format(cut),id_GMM_2d)
np.save('x_GMM_2d_{}'.format(cut),x_GMM_2d)
np.save('y_GMM_2d_{}'.format(cut),y_GMM_2d)
np.save('xsig_GMM_2d_{}'.format(cut),xsig_GMM_2d)
np.save('ysig_GMM_2d_{}'.format(cut),ysig_GMM_2d)
np.save('xycov_GMM_2d_{}'.format(cut),xycov_GMM_2d)
np.save('amp_GMM_2d_{}'.format(cut),amp_GMM_2d)

np.save('id_GMM_3d_{}'.format(cut),id_GMM_3d)
np.save('x_GMM_3d_{}'.format(cut),x_GMM_3d)
np.save('y_GMM_3d_{}'.format(cut),y_GMM_3d)
np.save('z_GMM_3d_{}'.format(cut),z_GMM_3d)
np.save('xsig_GMM_3d_{}'.format(cut),xsig_GMM_3d)
np.save('ysig_GMM_3d_{}'.format(cut),ysig_GMM_3d)
np.save('zsig_GMM_3d_{}'.format(cut),zsig_GMM_3d)
np.save('xycov_GMM_3d_{}'.format(cut),xycov_GMM_3d)
np.save('xzcov_GMM_3d_{}'.format(cut),xzcov_GMM_3d)
np.save('yzcov_GMM_3d_{}'.format(cut),yzcov_GMM_3d)
np.save('amp_GMM_3d_{}'.format(cut),amp_GMM_3d)


#print(id_GMM_2d)
#print(x_GMM_2d)
#print(y_GMM_2d)
#print(xsig_GMM_2d)
#print(ysig_GMM_2d)
#print(xycov_GMM_2d)
#print(amp_GMM_2d)

#print(id_GMM_3d)
#print(x_GMM_3d)
#print(y_GMM_3d)
#print(z_GMM_3d)
#print(xsig_GMM_3d)
#print(ysig_GMM_3d)
#print(zsig_GMM_3d)
#print(xycov_GMM_3d)
#print(xzcov_GMM_3d)
#print(yzcov_GMM_3d)
#print(amp_GMM_3d)

#print(id_BEAGLE)
#print(mass_BEAGLE_stellar)
