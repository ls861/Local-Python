
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 13:00:38 2020

@author: lester
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import cosmolopy.distance as cd
import cosmolopy.constants as cc

fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P']
runs = ['001']

fsize = 5
size = 8

for field in fields:    
    
    for run in runs:
        
        # =============================================================================
        # OUTPUT - get BEAGLE parameters
        # =============================================================================
        
        fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jun_2020/from_cluster/{}/data/BEAGLE_summary_catalogue.fits'.format(field[0])
        
        
        data_fits = fits.open(fileName)
        
        id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    
        sfr_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_median'])
        sfr_68_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_68.00'])
        
        ssfr_b1 = data_fits['STAR FORMATION'].data['sSFR_median']
        ssfr_68_b1 = data_fits['STAR FORMATION'].data['sSFR_68.00']
        
        redshift_b1 = data_fits['POSTERIOR PDF'].data['redshift_median']
        mass_b1 = data_fits['POSTERIOR PDF'].data['mass_median']
        msa_b1 = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_median']
        tauV_eff_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_median']
        metallicity_b1 = data_fits['POSTERIOR PDF'].data['metallicity_median']
        tau_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_median']
 
        redshift_68_b1 = data_fits['POSTERIOR PDF'].data['redshift_68.00']       
        mass_68_b1 = data_fits['POSTERIOR PDF'].data['mass_68.00']
        msa_68_b1 = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_68.00']
        tauV_eff_68_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_68.00']
        metallicity_68_b1 = data_fits['POSTERIOR PDF'].data['metallicity_68.00']
        tau_68_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_68.00']

        nebular_xi_b1 = np.full(len(id_b1), 0.3)
        nebular_xi_68_b1 = np.full((len(id_b1),2), 0.3)

        data_fits.close()
        
        # =============================================================================
        # calculate fitted output gradient for rising or falling
        # =============================================================================

        cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
        cosmo = cd.set_omega_k_0(cosmo)
        ageUniv2 = cd.age(2.0, **cosmo)/cc.yr_s # not actually used
        ageUniv999 = cd.age(999.0, **cosmo)/cc.yr_s # not actually used

        grad_out = np.empty(len(mass_b1))
        
        for i in range(len(mass_b1)):
    
            sfr_at_msa = msa_b1[i]*np.exp(-msa_b1[i]/tau_b1[i])
            sfr_at_msa_plus1 = (1.01*msa_b1[i])*np.exp(-(1.01*msa_b1[i])/tau_b1[i])
            grad_out[i] = sfr_at_msa_plus1 - sfr_at_msa
        
        idx_r1_out = grad_out >= 0
        idx_f1_out = grad_out < 0
        
        # =============================================================================
        # PLOT MAIN SEQUENCE
        # =============================================================================
            
        plt.scatter(mass_b1[idx_r1_out], sfr_b1[idx_r1_out])
        plt.scatter(mass_b1[idx_f1_out], sfr_b1[idx_f1_out])       
        plt.show()

        # =============================================================================
        # PLOT MAIN SEQUENCE for z=2 to 3
        # =============================================================================
            
        idx = abs(redshift_b1 - 2.5) < 0.5 
        plt.scatter(mass_b1[idx], sfr_b1[idx])      
        plt.show()
        
        
        

        # =============================================================================
        # plot Z fitted with BEAGLE vs ZBEST from Astrodeep
        # =============================================================================

        fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jun_2020/from_cluster/{}/astrodeep_{}_{}_subset_RF1_001.fits'.format(field[0], field[1:-1], field[-1].lower())   
        
        data_fits = fits.open(fileName)
        # print(data_fits.info())
        # print(data_fits[1].header)

        id_input = data_fits[1].data['ID'][id_b1-1]
        zbest = data_fits[1].data['ZBEST'][id_b1-1]

        data_fits.close()



        # Calculate the point density
        from scipy.stats import gaussian_kde
        x = zbest
        y = redshift_b1
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        
        plt.title('fitted redshift vs zbest')
        plt.scatter(x, y, c=z)
        plt.show()





























        # =============================================================================
        # create fits file of average chi2 - TAKES AGES TO RUN COMMENT OUT FOR NOW
        # =============================================================================
        
   
        ''' UNCOMMENT WHEN NEEDED
        
        samples = 10
        IDs = id_b1
        chi2_fit_arr_allID = [] 
        for ID in IDs:
            
#            title = '{}-{} ID{} {} samples'.format(param, revision.replace('_', '-'), str(ID), str(samples))
            
            # =============================================================================
            # Get OUTPUT SEDs
            # =============================================================================
        
            data_fits = fits.open('/Users/lester/Documents/PhD/{}/fit_{}/{}_BEAGLE.fits.gz'.format(field, run, ID))
            
            # print(data_fits.info())
            # print(data_fits['POSTERIOR PDF'].header)
            
            chi2_fit_total = data_fits['POSTERIOR PDF'].data['chi_square']
            
            #needs float64 to provide precision needed for the random.choice weights
            temp_probs = np.float64(data_fits['POSTERIOR PDF'].data['probability'])
            temp_probs = temp_probs/np.sum(temp_probs)
              
            chi2_fit_arr = []
            
            for i in range(samples):
                
                #here's the key line - take weighted samples from the multinest output!
                idx = np.random.choice(len(temp_probs), size=1, p=temp_probs)
                
                # CHI SQUARED
                chi2_fit_arr.append(data_fits['POSTERIOR PDF'].data['chi_square'][idx][0])
        
            chi2_fit_arr_allID.append(np.average(chi2_fit_arr))

        
        print(len(id_b1.astype(str)))
        print(len(chi2_fit_arr_allID))
        
        
        outputDict = {}
        outputDict['id']                = id_b1.astype(str)
        outputDict['chi2']              = chi2_fit_arr_allID
        outputTable = Table(outputDict)
        
#        outputTable.write('{}_{}_chi2.fits'.format(field, run), overwrite=True)
        print('LESTER')
        plt.hist(outputTable['chi2'], bins=100)
        plt.show()
                
        plt.hist(outputTable['chi2'], bins=100)
        plt.ylim(0, 10)
        plt.show()
        
        plt.hist(outputTable['chi2'], bins=100, range=(0, 200))
        plt.show()
        print('LESTER')
        print(len(outputTable['chi2']))
        print(len(outputTable['chi2'][outputTable['chi2']<=25]))
        print(len(outputTable['chi2'][outputTable['chi2']>25]))
        
        
        UNCOMMENT WHEN NEEDED '''
              
        # =============================================================================
        # PLOT SPECIFIC SFHs GIVEN ID (eg for diagnostics given poor chi2)        
        # =============================================================================
                
#        fsize=2        
#        ageUniv = cd.age(redshift_b1, **cosmo)/cc.yr_s
#        
#        xlin = np.linspace(1, 1.1e10, 100000)
#        
#        IDs = ([20])
#        IDs = id_b1[-3:]
#        IDs=[21, 115, 3030]
#        
#        for ID in IDs:
#            
#            i = (np.abs(id_b1 - ID)).argmin()
#            
#            plt.figure(figsize=(4*fsize, fsize))
#            plt.xlim(0, 1.1e10)
#            plt.ylim(0.0, 1.1)
#            
#            sfr_out = 1 * (xlin-(ageUniv[i]-msa_b1[i]))*np.exp(-(xlin-(ageUniv[i]-msa_b1[i]))/tau_b1[i])
#            plt.plot(xlin, sfr_out/max(sfr_out), label='OUTPUT SFH')
#            
#            plt.plot((ageUniv[i], ageUniv[i]), (0, 1))
#            plt.legend()
#            plt.show()
#            
#            print(10**mass_b1[i], msa_b1[i], tau_b1[i])
        

            
        # =============================================================================
        # ASTRODEEP filter info
        # =============================================================================
        
#        # column name from apparent mag table in BEAGLE mock (same as left column in config + _APP)
#        filters = ['HST_ACS_WFC_F435W_APP', 'HST_ACS_WFC_F606W_APP', 'HST_ACS_WFC_F814W_APP', 'HST_WFC3_IR_F105W_APP', 'HST_WFC3_IR_F125W_APP', 'HST_WFC3_IR_F140W_APP', 'HST_WFC3_IR_F160W_APP', 'Paranal_HAWKI_Ks_APP', 'Spitzer_IRAC_I1_APP', 'Spitzer_IRAC_I2_APP']
#        
#        # column name from ASTRODEEP config file
#        filter_label = ['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2'] 
#        
#        # for taking errors of the perturbed input fluxes
#        filter_err = ['b_errB435', 'b_errV606', 'b_errI814', 'b_errY105', 'b_errJ125', 'b_errJH140', 'b_errH160', 'b_errKs', 'b_errCH1', 'b_errCH2'] 
#        
#        filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])
#        filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])
        
        # =============================================================================
        # get INPUT perturbed fluxes and errors
        # =============================================================================
        
#        fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jun_2020/from_cluster/{}/astrodeep_{}_{}_subset_RF1_001.fits'.format(field[0], field[1:-1], field[-1].lower())   
#        
#        data_fits = fits.open(fileName)
#        # print(data_fits.info())
#        # print(data_fits[1].header)
#        
#        # PHOTOMETRY 
#        ptbf_phot_mock = np.zeros(len(filter_label))
#        for i in range(len(filters)):
#            ptbf_phot_mock[i] = data_fits[1].data[filter_label[i]][ID-1] # uJy
#
#        ptbferr_phot_mock = np.zeros(len(filter_label))
#        for i in range(len(filters)):
#            ptbferr_phot_mock[i] = data_fits[1].data[filter_err[i]][ID-1] # uJy    
#            
#        # adding min rel error to errors: 
#        min_rel_error = np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.05, 0.1, 0.1])
#        ptbferr_phot_mock = ( (ptbferr_phot_mock**2) + ((min_rel_error*ptbf_phot_mock)**2) ) ** 0.5
#        
#        ptblfl_phot_mock = (ptbf_phot_mock/filter_fwhm_centre)*(1e-10 / 3.34) # erg cm-2 s-1
#        ptblflerr_phot_mock = (ptbferr_phot_mock/filter_fwhm_centre)*(1e-10 / 3.34) # erg cm-2 s-1
#        
#        data_fits.close()
        
















