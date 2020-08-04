
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
from scipy.stats import norm

fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P']
#fields = ['0A2744C']
runs = ['001']

fsize = 5
size = 8

for field in fields:    
    
    for run in runs:
        
        # =============================================================================
        # GET DATA
        # =============================================================================
        
        directory = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jun_2020/from_cluster/'
        
        # BEAGLE OUTPUT SUMMARY
        fileName = directory+'{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(field[0])
        data_fits = fits.open(fileName)
#        print(data_fits.info())
#        print(data_fits[1].header)
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

        # BEAGLE OUTPUT CHI SQUARED
        fileName = directory+'lester_run/{}_chi2.fits'.format(field)
        data_fits = fits.open(fileName)
#        print(data_fits.info())
#        print(data_fits[1].header)
        id_chi2 = np.asarray(data_fits[1].data['id'], dtype=int)
        chi2 = data_fits[1].data['chi2']
        data_fits.close()

        # BEAGLE INPUT FLUXES
        fileName = directory+'{}/astrodeep_{}_{}_subset_RF1_001.fits'.format(field[0], field[1:-1], field[-1].lower()) 
        data_fits = fits.open(fileName)
#        print(data_fits.info())
#        print(data_fits[1].header)
        id_input = np.asarray(data_fits[1].data['ID'][id_b1-1], dtype=int)
        field_original = np.asarray(data_fits[1].data['field'][id_b1-1], dtype=int)
        id_original = np.asarray(data_fits[1].data['ID_original'][id_b1-1], dtype=int)
        zbest = data_fits[1].data['ZBEST'][id_b1-1]
        data_fits.close()

        # ASTRODEEP CATALOG
        catalog = np.load(directory[:-13]+'astrodeep_rawfile_ABCZ.npy')
#        print(catalog.dtype.names)
        
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
        # PLOT MAIN SEQUENCE per redshift bin
        # =============================================================================
         
        
        z_med = np.linspace(1.5, 9.5, 9)
        z_gap = (z_med[1] - z_med[0]) / 2

        for z in z_med:
          idx = abs(redshift_b1 - z) < z_gap
          plt.scatter(mass_b1[idx], sfr_b1[idx], label='{}  - {}'.format(z-z_gap, z+z_gap))  
        plt.xlim(6, 11)
        plt.ylim(-5, 5)
        plt.legend(title='redshift')
        plt.show()
        
        
        # =============================================================================
        # PLOT adding in mass dependence somehow   
        # =============================================================================
        
        m_med = np.linspace(5.5, 11.5, 7)
        m_gap = (m_med[1] - m_med[0]) / 2

        for m in m_med:
          idx = abs(mass_b1 - m) < m_gap
          plt.scatter(mass_b1[idx], sfr_b1[idx], label='{}  - {}'.format(m-m_gap, m+m_gap))  
        plt.xlim(6, 11)
        plt.ylim(-5, 5)
        plt.legend(title='mass')
        plt.show()
        
        # =============================================================================
        # PLOT combining redshift and mass     
        # =============================================================================
        
#        for z in z_med:
#          for m in m_med:
#            z_idx = abs(redshift_b1 - z) < z_gap
#            m_idx = abs(mass_b1 - m) < m_gap
#            idx = z_idx & m_idx
#            
#            plt.scatter(mass_b1[idx], sfr_b1[idx], label='{}  - {}'.format(m-m_gap, m+m_gap))  
#          plt.title('{} redshift {}  - {}'.format(field, z-z_gap, z+z_gap))
#          plt.xlim(6, 11)
#          plt.ylim(-5, 5)
#          plt.legend(title='mass')
#          plt.show()          
            
        # =============================================================================
        # calculating alpha and beta for given redshifts (and sigma)
        # =============================================================================
           #%%   
           
#           
#        alpha = []
#        beta = []
#        sigma = []
#        
#        for z in z_med:
#          idx = abs(redshift_b1 - z) < z_gap
#          
#          # np.where is to sort the -inf SFR values at low redshift...
#          lin = np.polyfit(mass_b1[idx], np.where(sfr_b1[idx]<-50, -50, sfr_b1[idx]), 1)
#          
#          alpha.append(lin[1])
#          beta.append(lin[0])
#  
#
#          plt.xlim(6, 11)
#          plt.ylim(-5, 5)
#          plt.plot((6, 11), ((6*lin[0] + lin[1]),(11*lin[0] + lin[1])))
#          plt.scatter(mass_b1[idx], sfr_b1[idx], label='{}  - {}'.format(z-z_gap, z+z_gap))  
#          plt.legend(title='redshift')
#          plt.show()          
#          
#          # np.where is to sort the -inf SFR values at low redshift...
#          data = np.where(sfr_b1[idx]<-50, -50, sfr_b1[idx]) - (mass_b1[idx]*lin[0] + lin[1])
#          mean,std=norm.fit(data)
#          sigma.append(std)
#          
#          
#          testx = np.linspace(min(data), max(data), 1000)
#          plt.hist(data, density=True)
#          plt.plot(testx, norm.pdf(testx, mean, std))
#          plt.show()
#          
#      
#        plt.title('{} alpha vs redshift'.format(field))
#        plt.scatter(z_med, alpha)
#        plt.show()
#        
#        plt.title('{} beta vs redshift'.format(field))
#        plt.scatter(z_med, beta)
#        plt.show()
#        
#        plt.title('{} sigma vs redshift'.format(field))
#        plt.scatter(z_med, sigma)
#        plt.show()        
#        
          
          
     #%%     
          
#        # sigma vs redshift for given mass bin
#          
#        for m in m_med:     
#          sigma = []
#          for z in z_med:
#
#            z_idx = abs(redshift_b1 - z) < z_gap
#            m_idx = abs(mass_b1 - m) < m_gap
#            idx = z_idx & m_idx
#            
#
#            if len(mass_b1[idx]) == 0:
#              sigma.append(-1)
#            else:
#              # np.where is to sort the -inf SFR values at low redshift...
#              lin = np.polyfit(mass_b1[idx], np.where(sfr_b1[idx]<-50, -50, sfr_b1[idx]), 1)
#              data = np.where(sfr_b1[idx]<-50, -50, sfr_b1[idx]) - (mass_b1[idx]*lin[0] + lin[1])
#              mean,std=norm.fit(data)
#              sigma.append(std)          
#          
#          plt.title('{} sigma vs redshift: mass {}  - {}'.format(field, m-m_gap, m+m_gap))
#          plt.scatter(z_med, sigma)
#          plt.show()            
#          
#        # sigma vs mass for given redshift bin
#          
#        for z in z_med:     
#          sigma = []
#          for m in m_med:
#
#            z_idx = abs(redshift_b1 - z) < z_gap
#            m_idx = abs(mass_b1 - m) < m_gap
#            idx = z_idx & m_idx
#            
#
#            if len(mass_b1[idx]) == 0:
#              sigma.append(-1)
#            else:
#              # np.where is to sort the -inf SFR values at low redshift...
#              lin = np.polyfit(mass_b1[idx], np.where(sfr_b1[idx]<-50, -50, sfr_b1[idx]), 1)
#              data = np.where(sfr_b1[idx]<-50, -50, sfr_b1[idx]) - (mass_b1[idx]*lin[0] + lin[1])
#              mean,std=norm.fit(data)
#              sigma.append(std)          
#          
#          plt.title('{} sigma vs mass: redshift {}  - {}'.format(field, z-z_gap, z+z_gap))
#          plt.scatter(m_med, sigma)
#          plt.show()             
#          
          
          #%%
          

        # =============================================================================
        # plot Z fitted with BEAGLE vs ZBEST from Astrodeep
        # =============================================================================

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
        # MAIN SEQUENCE inc chi2 and magnification
        # =============================================================================
        
        # note id_b1 == id_chi2 == id_input, for image 0, length 3033, max 3043
        
        #id_b1
        #sfr_b1
        #mass_b1
        
        #id_chi2
        #chi2
        
        #field_original
        #id_original
        #id_input
        
        #catalog['MAGNIF']
        
        xsize = 6
        ysize = 6
        
        idx_field = (catalog['field']==float(field[0]))
        catalog_subset = catalog[idx_field] # take field subset of total catalog
        
        idx_id = np.isin(catalog_subset['ID'], id_original)
        catalog_subset = catalog_subset[idx_id] # take ID subset of total catalog
        
        idx_chi2 = (chi2 < 9.5)
        catalog_subset = catalog_subset[idx_chi2] # take chi2 subset of total catalog
        
        mag = np.log10(catalog_subset['MAGNIF'])
        
        
        # plot main sequence
        plt.figure(figsize=(xsize, ysize))
        plt.title('{} Main Sequence'.format(field))
        x = mass_b1
        y = sfr_b1
        xerr = mass_68_b1
        yerr = sfr_68_b1
        plt.xlim(6, 12)
        plt.ylim(-5, 5)
        plt.scatter(x, y, marker='x', zorder=2)
        plt.errorbar(x, y, xerr=[x-xerr[:,0],xerr[:,1]-x], yerr=[y-yerr[:,0],yerr[:,1]-y], linestyle="None", elinewidth=0.5, color='k', zorder=1)
        plt.show()
        
        # plot main sequence filtered by chi2
        plt.figure(figsize=(xsize, ysize))
        plt.title('{} Main Sequence, chi2 less than 9.5'.format(field))
        x = mass_b1[idx_chi2]
        y = sfr_b1[idx_chi2]
        xerr = mass_68_b1[idx_chi2]
        yerr = sfr_68_b1[idx_chi2]
        plt.xlim(6, 12)
        plt.ylim(-5, 5)
        plt.scatter(x, y, marker='x', zorder=2)
        plt.errorbar(x, y, xerr=[x-xerr[:,0],xerr[:,1]-x], yerr=[y-yerr[:,0],yerr[:,1]-y], linestyle="None", elinewidth=0.5, color='k', zorder=1)
        plt.show()
        
        # plot main sequence filtered by chi2 including magnification
        plt.figure(figsize=(xsize, ysize))
        plt.title('{} Main Sequence, chi2 less than 9.5, adjusted for magnification'.format(field))
        x = mass_b1[idx_chi2]-mag
        y = sfr_b1[idx_chi2]-mag
        xerr = mass_68_b1[idx_chi2]-np.array([mag, mag]).transpose()
        yerr = sfr_68_b1[idx_chi2]-np.array([mag, mag]).transpose()
        plt.xlim(6, 12)
        plt.ylim(-5, 5)
        plt.scatter(x, y, marker='x', zorder=2)
        plt.errorbar(x, y, xerr=[x-xerr[:,0],xerr[:,1]-x], yerr=[y-yerr[:,0],yerr[:,1]-y], linestyle="None", elinewidth=0.5, color='k', zorder=1)
        plt.show()
        
        
















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


        
















