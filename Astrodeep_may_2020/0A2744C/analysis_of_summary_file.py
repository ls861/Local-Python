
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 13:00:38 2020

@author: lester
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import cosmolopy.distance as cd
import cosmolopy.constants as cc
#from scipy.integrate import quad

#params = ['DE']
#revisions = ['100', '101', '102', '103', '104']
#revisions = ['106']

fields = ['0A2744C']
runs = ['001']



fsize = 5
size = 8

for field in fields:    
    
    for run in runs:
        
#        title1 = param1 + ' ' + revision1.replace('_', ' ')
        
#        # =============================================================================
#        # get INPUT params (len 100)
#        # =============================================================================
#        
#        fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/{}/mock_MS_parameters_100_{}.fits'.format(param1, param1)
#        data_fits = fits.open(fileName)
#    #    print(data_fits[1].header)
#    
#        id_i = data_fits[1].data['id']
#        
#        mass = data_fits[1].data['mass']
#        msa = 10**data_fits[1].data['max_stellar_age']
#        tauV_eff = data_fits[1].data['tauV_eff']
#        metallicity = data_fits[1].data['metallicity']
#        nebular_logU = data_fits[1].data['nebular_logU']
#        tau = 10**(data_fits[1].data['tau'])
#        nebular_xi = data_fits[1].data['nebular_xi']
#        redshift = np.full(len(id_i), 2.)
#        
#        A = data_fits[1].data['A']
#        sfr = np.log10(data_fits[1].data['sfr'])
#        
#        data_fits.close()
        
        # =============================================================================
        # calculate input gradient for rising or falling (len 100)
        # =============================================================================
        

#        grad_in = np.empty(len(A))
#        
#        # nice trick to find index in xlin which has value closest to ageUniv2
#        idx = (np.abs(xlin - ageUniv2)).argmin()
#        
#        for i in range(len(A)):
#            sfr_calc = A[i] * (xlin-(ageUniv2-msa[i]))*np.exp(-(xlin-(ageUniv2-msa[i]))/tau[i])
#            grad_in[i] = np.gradient(sfr_calc, xlin)[idx]
#        
#        print(msa[i])
#        idx_rising_in = grad_in >= 0
#        idx_falling_in = grad_in < 0
        
        # =============================================================================
        # OUTPUT - get BEAGLE parameters (<100)
        # =============================================================================
        
#        fileName = '/Users/lester/Documents/PhD/param_100/fit_{}_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(revision1, param1)
        fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_may_2020/{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(field)
        data_fits = fits.open(fileName)
        
        id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    
        sfr_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_mean'])
        sfr_68_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_68.00'])
        
        ssfr_b1 = data_fits['STAR FORMATION'].data['sSFR_mean']
        ssfr_68_b1 = data_fits['STAR FORMATION'].data['sSFR_68.00']
        
        redshift_b1 = data_fits['POSTERIOR PDF'].data['redshift_mean']
        mass_b1 = data_fits['POSTERIOR PDF'].data['mass_mean']
        msa_b1 = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_mean']
        tauV_eff_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_mean']
        metallicity_b1 = data_fits['POSTERIOR PDF'].data['metallicity_mean']
        tau_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_mean']
 
        redshift_68_b1 = data_fits['POSTERIOR PDF'].data['redshift_68.00']       
        mass_68_b1 = data_fits['POSTERIOR PDF'].data['mass_68.00']
        msa_68_b1 = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_68.00']
        tauV_eff_68_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_68.00']
        metallicity_68_b1 = data_fits['POSTERIOR PDF'].data['metallicity_68.00']
        tau_68_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_68.00']

        nebular_xi_b1 = np.full(len(id_b1), 0.3)
        nebular_xi_68_b1 = np.full((len(id_b1),2), 0.3)
        
        print(msa_b1)
        
#        if revision1 in ['101']:
#            nebular_logU_b1 = np.full(len(id_b1), -2.5)
#            nebular_logU_68_b1 = np.full((len(id_b1),2), -2.5)
#        else:
#            nebular_logU_b1 = data_fits['POSTERIOR PDF'].data['nebular_logu_mean']        
#            nebular_logU_68_b1 = data_fits['POSTERIOR PDF'].data['nebular_logu_68.00']      


        data_fits.close()
        
        # =============================================================================
        # calculate output gradient for rising or falling
        # =============================================================================


        cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
        cosmo = cd.set_omega_k_0(cosmo)
        ageUniv2 = cd.age(2.0, **cosmo)/cc.yr_s
        ageUniv999 = cd.age(999.0, **cosmo)/cc.yr_s
#
#    

        
        grad_out = np.empty(len(mass_b1))
#        A_b1 = np.empty(len(mass_b1))
#        integral = np.empty(len(mass_b1))
        

        
        
        for i in range(len(mass_b1)):
    
            sfr_at_msa = msa_b1[i]*np.exp(-msa_b1[i]/tau_b1[i])
            sfr_at_msa_plus1 = (1.01*msa_b1[i])*np.exp(-(1.01*msa_b1[i])/tau_b1[i])
            grad_out[i] = sfr_at_msa_plus1 - sfr_at_msa
        

            
        idx_r1_out = grad_out >= 0
        idx_f1_out = grad_out < 0
        
        
        # =============================================================================
        # PLOTS      
        # =============================================================================
            
        plt.scatter(mass_b1[idx_r1_out], sfr_b1[idx_r1_out])
        plt.scatter(mass_b1[idx_f1_out], sfr_b1[idx_f1_out])       
        
        
        
        
        #        idx_r1_out = idx_rising_out
        #        idx_f1_out = idx_falling_out
                
        #        idx_r1_in = idx_rising_in[id_b1]
        #        idx_f1_in = idx_falling_in[id_b1]
                
        #        idx_rr = np.logical_and(idx_r1_in, idx_r1_out)
        #        idx_ff = np.logical_and(idx_f1_in, idx_f1_out)
        #        idx_rf = np.logical_and(idx_r1_in, idx_f1_out)
        #        idx_fr = np.logical_and(idx_f1_in, idx_r1_out)
        #        
        #        sum_rr = np.sum(np.logical_and(idx_r1_in, idx_r1_out))
        #        sum_ff = np.sum(np.logical_and(idx_f1_in, idx_f1_out))
        #        sum_rf = np.sum(np.logical_and(idx_r1_in, idx_f1_out))
        #        sum_fr = np.sum(np.logical_and(idx_f1_in, idx_r1_out))
 

        '''       
        # =======================================================================f======
        # PLOT comparing individual parameters
        # =============================================================================
        
        params_names = ['mass', 'msa', 'tauVeff', 'metallicity', 'nebularlogU', 'tau', 'nebularxi', 'redshift']
        params = [mass, np.log10(msa), tauV_eff, metallicity, nebular_logU, np.log10(tau), nebular_xi, redshift]
        params_b1 = [mass_b1, np.log10(msa_b1), tauV_eff_b1, metallicity_b1, nebular_logU_b1, np.log10(tau_b1), nebular_xi_b1, redshift_b1]
        params_68_b1 = [mass_68_b1, np.log10(msa_68_b1), tauV_eff_68_b1, metallicity_68_b1, nebular_logU_68_b1, np.log10(tau_68_b1), nebular_xi_68_b1, redshift_68_b1]
        
        fig, axs = plt.subplots(2, 4, figsize=(3.2*fsize, 1.6*fsize))
        fig.suptitle(title1)
        for j in [0, 1]:
            for i in range(4):
                if j==1 and i==3 and revision1 not in ['105', '106']:
                    break
                
                axs[j,i].set_title(params_names[i+4*j])
                
                axs[j,i].scatter(params[i+4*j][id_b1][idx_rr], params_b1[i+4*j][idx_rr], s=10, zorder=2, color='r', label='RR {}'.format(sum_rr))
                axs[j,i].scatter(params[i+4*j][id_b1][idx_ff], params_b1[i+4*j][idx_ff], s=10, zorder=2, color='b', label='FF {}'.format(sum_ff))      
                
                axs[j,i].scatter(params[i+4*j][id_b1][idx_rf], params_b1[i+4*j][idx_rf], s=10, zorder=2, color='c', label='RF {}'.format(sum_rf), marker='o')
                axs[j,i].scatter(params[i+4*j][id_b1][idx_fr], params_b1[i+4*j][idx_fr], s=10, zorder=2, color='m', label='FR {}'.format(sum_fr), marker='o')   
                
                axs[j,i].errorbar(params[i+4*j][id_b1], params_b1[i+4*j], yerr=[params_b1[i+4*j] - params_68_b1[i+4*j][:, 0], params_68_b1[i+4*j][:, 1] - params_b1[i+4*j]], linestyle="None", elinewidth=0.5, color='k', zorder=1)
                
                min_ax = min(min(params[i+4*j][id_b1]), min(params_b1[i+4*j]))
                #            min_ax = 9.3
                max_ax = max(max(params[i+4*j][id_b1]), max(params_b1[i+4*j]))
                
                axs[j,i].plot((min_ax, max_ax), (min_ax, max_ax), color='k', zorder=0)
                
                axs[j,i].set_xlim(min_ax, max_ax)
                axs[j,i].set_ylim(min_ax, max_ax)
                axs[j,i].legend()
        plt.show()


  
        # =============================================================================
        # PLOT - input mass vs output mass 1
        # =============================================================================
        
        plt.figure(figsize=(fsize, fsize))
        plt.title('Input Mass (DE) vs Output Mass ({})'.format(title1), size=size)
        plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
        plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)
    
        plt.scatter(mass[id_b1][idx_rr], mass_b1[idx_rr], s=10, zorder=2, color='r', label='RR {}'.format(sum_rr))
        plt.scatter(mass[id_b1][idx_ff], mass_b1[idx_ff], s=10, zorder=2, color='b', label='FF {}'.format(sum_ff))      
        
        plt.scatter(mass[id_b1][idx_rf], mass_b1[idx_rf], s=10, zorder=2, color='c', label='RF {}'.format(sum_rf), marker='o')
        plt.scatter(mass[id_b1][idx_fr], mass_b1[idx_fr], s=10, zorder=2, color='m', label='FR {}'.format(sum_fr), marker='o')   
        
        plt.errorbar(mass[id_b1], mass_b1, yerr=[mass_b1 - mass_68_b1[:, 0], mass_68_b1[:, 1] - mass_b1], linestyle="None", elinewidth=0.5, color='k', zorder=1)
        
        min_ax = 7.5
        max_ax = 11.0
        
        plt.plot((min_ax, max_ax), (min_ax, max_ax), color='k', zorder=0)
        plt.xlim(min_ax, max_ax)
        plt.ylim(min_ax, max_ax)
        plt.legend()
        plt.show()
    
        # =============================================================================
        # PLOT - input sfr vs output sfr 1
        # =============================================================================
    
        plt.figure(figsize=(fsize, fsize))
        plt.title('Input SFR (DE) vs Output SFR ({})'.format(title1), size=size)
        plt.xlabel(r'$\text{Input - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
        plt.ylabel(r'$\text{Output - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
    
        plt.scatter(sfr[id_b1][idx_rr], sfr_b1[idx_rr], s=10, zorder=2, color='r', label='RR {}'.format(sum_rr))
        plt.scatter(sfr[id_b1][idx_ff], sfr_b1[idx_ff], s=10, zorder=2, color='b', label='FF {}'.format(sum_ff))      
        
        plt.scatter(sfr[id_b1][idx_rf], sfr_b1[idx_rf], s=10, zorder=2, color='c', label='RF {}'.format(sum_rf), marker='o')
        plt.scatter(sfr[id_b1][idx_fr], sfr_b1[idx_fr], s=10, zorder=2, color='m', label='FR {}'.format(sum_fr), marker='o')   
        
        plt.errorbar(sfr[id_b1], sfr_b1, yerr=[sfr_b1 - sfr_68_b1[:, 0], sfr_68_b1[:, 1] - sfr_b1], linestyle="None", elinewidth=0.5, color='k', zorder=1)
        
        min_ax = -1.0
        max_ax = 3.5
        
        plt.plot((min_ax, max_ax), (min_ax, max_ax), color='k', zorder=0)
        plt.xlim(min_ax, max_ax)
        plt.ylim(min_ax, max_ax)
        plt.legend()
        plt.show()
    
    
        # =============================================================================
        # PLOT - input sfr vs output ssfr 1
        # =============================================================================
    
        ssfr = sfr - mass # log space    
       
        plt.figure(figsize=(fsize, fsize))
        plt.title('Input SSFR (DE) vs Output SSFR ({})'.format(title1), size=size)
        plt.xlabel(r'$\text{Input - log}(\text{SSFR} / yr^{-1})$', size=size)
        plt.ylabel(r'$\text{Output - log}(SSFR / yr^{-1})$', size=size)
        
        plt.scatter(ssfr[id_b1][idx_rr], ssfr_b1[idx_rr], s=10, zorder=2, color='r', label='RR {}'.format(sum_rr))
        plt.scatter(ssfr[id_b1][idx_ff], ssfr_b1[idx_ff], s=10, zorder=2, color='b', label='FF {}'.format(sum_ff))      
        
        plt.scatter(ssfr[id_b1][idx_rf], ssfr_b1[idx_rf], s=10, zorder=2, color='c', label='RF {}'.format(sum_rf), marker='o')
        plt.scatter(ssfr[id_b1][idx_fr], ssfr_b1[idx_fr], s=10, zorder=2, color='m', label='FR {}'.format(sum_fr), marker='o')   
        
        plt.errorbar(ssfr[id_b1], ssfr_b1, yerr=[ssfr_b1 - ssfr_68_b1[:, 0], ssfr_68_b1[:, 1] - ssfr_b1], linestyle="None", elinewidth=0.5, color='k', zorder=1)
    
        min_ax = -10.0
        max_ax = -7.0
        
        plt.plot((min_ax, max_ax), (min_ax, max_ax), color='k', zorder=0)
        plt.xlim(min_ax, max_ax)
        plt.ylim(min_ax, max_ax)
        plt.legend()
        plt.show()    
    
    
        # =============================================================================
        # calculate distance between points to find ID of "good or "bad" ones
        # =============================================================================
        
        plt.title('Distance between input and output SSFR')
        plt.hist(abs(ssfr[id_b1]-ssfr_b1))
        plt.show()
        
        good_fit_idx = id_b1[abs(ssfr[id_b1]-ssfr_b1) < 0.35]
        print(good_fit_idx+1)
        
        bad_fit_idx = id_b1[abs(ssfr[id_b1]-ssfr_b1) >= 0.35]
        print(bad_fit_idx+1)
    
        # =============================================================================
        # just messing around
        # =============================================================================
    
        print(id_b1[idx_rr]+1)
        print(id_b1[idx_ff]+1)
        print(id_b1[idx_rf]+1)
        print(id_b1[idx_fr]+1)
    
    
    
    
        if revision1 in ['105', '106']:
            plt.hist(redshift_b1)
    
    


        '''



        # =============================================================================
        # create fits file of average chi2
        # =============================================================================
        
        
        #!/usr/bin/env python2
        # -*- coding: utf-8 -*-
        """
        Created on Sun Apr  5 16:22:23 2020
        
        @author: lester
        """
        
        #import numpy as np
        from astropy.table import Table
        #from astropy.io import fits
        #import matplotlib.pyplot as plt
        

        
        samples = 10
        
        # =============================================================================
        # OUTPUT - get BEAGLE parameters (<100)
        # =============================================================================
        
#        fileName = '/Users/lester/Documents/PhD/param_100/fit_{}_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(revision, param)
#        data_fits = fits.open(fileName)
        
#        id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
        
#        data_fits.close()
        
        chi2_fit_arr_allID = [] 
        
        IDs = id_b1
        print(IDs)
        
        ''' UNCOMMENT WHEN NEEDED
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
                
                
        #!/usr/bin/env python2
        # -*- coding: utf-8 -*-
        """
        Created on Fri Apr 10 09:09:20 2020
        
        @author: lester
        """
        
        #import numpy as np
        #import matplotlib.pyplot as plt
        #from astropy.io import fits
        #import cosmolopy.distance as cd
        #import cosmolopy.constants as cc
        
        #revision = 'DE'
        fsize=2
            
        #fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/{}_combined.fits'.format(revision)
        #data_fits = fits.open(fileName)
        ##    print(data_fits[1].header)
        #
        #ID = np.array(data_fits[1].data['id_1']).astype(float)
        #
        #mass_in = 10**(data_fits[1].data['mass'])
        #msa_in = 10**data_fits[1].data['max_stellar_age']
        #tau_in = 10**(data_fits[1].data['tau'])
        
        mass_out = 10**mass_b1
        msa_out = msa_b1
        tau_out = tau_b1
        
        
        #data_fits.close()
        
        ## calculated using python file 8
        #ID_rr_DE = np.array([ 4,  5,  8, 10, 11, 12, 13, 14, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 34, 35, 36, 37, 39, 41, 42, 43, 44, 45, 47, 49, 50, 51, 53, 55, 56, 58, 61, 62, 63, 64, 65, 66, 68, 70, 71, 73, 75, 78, 79, 81, 83, 84, 88, 89, 90, 92, 93, 94, 95, 96, 98, 99])
        #ID_ff_DE = np.array([ 3, 59])
        #ID_rf_DE = np.array([])
        #ID_fr_DE = np.array([2,   7,   9,  15,  18, 33,  38,  40,  52,  69,  72,  74,  76,  77,  80,  82,  85,  86,87,  91,  97, 100])
        
        
        
        
        #cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
        #cosmo = cd.set_omega_k_0(cosmo)
        #ageUniv2 = cd.age(2.0, **cosmo)/cc.yr_s
        #ageUniv999 = cd.age(999.0, **cosmo)/cc.yr_s
        
        ageUniv = cd.age(redshift_b1, **cosmo)/cc.yr_s
        
        
        xlin = np.linspace(1, 1.1e10, 100000)
        
        IDs = ([20])
        IDs = id_b1[-3:]
        
        for j in IDs:
            
            i = (np.abs(id_b1 - j)).argmin()
            
            plt.figure(figsize=(4*fsize, fsize))
        #    plt.title('{} ID {}'.format(revision.replace('_', ' '), j))
            plt.xlim(0, 1.1e10)
        #    plt.xlim(0.9*msa, 1.1*msa)
            plt.ylim(0.0, 1.1)
        
        #    sfr_in = 1 * (xlin-(ageUniv2-msa_in[i]))*np.exp(-(xlin-(ageUniv2-msa_in[i]))/tau_in[i])
        #    plt.plot(xlin, sfr_in/max(sfr_in), label='INPUT SFH')
            
            sfr_out = 1 * (xlin-(ageUniv[i]-msa_out[i]))*np.exp(-(xlin-(ageUniv[i]-msa_out[i]))/tau_out[i])
            plt.plot(xlin, sfr_out/max(sfr_out), label='OUTPUT SFH')
            
            plt.plot((ageUniv[i], ageUniv[i]), (0, 1))
            plt.legend()
            plt.show()
            
        #    print(mass_in[i], msa_in[i], tau_in[i])
            print(mass_out[i], msa_out[i], tau_out[i])
        
            
        
        
        
                
# =============================================================================
#         SEDs
# =============================================================================
        
        
        #!/usr/bin/env python2
        # -*- coding: utf-8 -*-
        """
        Created on Sun Apr  5 16:22:23 2020
        
        @author: lester
        """
        
        #import numpy as np
        #import matplotlib.pyplot as plt
        #from astropy.io import fits
        
        #param = 'DE'
        #revision = '100'
        
        samples = 10
        
        IDs=[22, 61, 63]
        IDs = id_b1[:10]
        
        fsize=10
        c = 299792458 # m s^-2
        
        
        for ID in IDs:
            
            title = '{}-{} ID{} {} samples'.format(field, run, str(ID), str(samples))
            
            # =============================================================================
            # ASTRODEEP filter info
            # =============================================================================
            
            # column name from apparent mag table in BEAGLE mock (same as left column in config + _APP)
            filters = ['HST_ACS_WFC_F435W_APP', 'HST_ACS_WFC_F606W_APP', 'HST_ACS_WFC_F814W_APP', 'HST_WFC3_IR_F105W_APP', 'HST_WFC3_IR_F125W_APP', 'HST_WFC3_IR_F140W_APP', 'HST_WFC3_IR_F160W_APP', 'Paranal_HAWKI_Ks_APP', 'Spitzer_IRAC_I1_APP', 'Spitzer_IRAC_I2_APP']
            
            # column name from ASTRODEEP config file
            filter_label = ['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2'] 
            
            # for taking errors of the perturbed input fluxes
            filter_err = ['b_errB435', 'b_errV606', 'b_errI814', 'b_errY105', 'b_errJ125', 'b_errJH140', 'b_errH160', 'b_errKs', 'b_errCH1', 'b_errCH2'] 
            
            filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])
            filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])
            
            # =============================================================================
            # from INPUT fits file get SEDs and apparent ABMAG
            # =============================================================================
        
            
#            fileName = '/Users/lester/Documents/PhD/param_100/mock_100_{}/mock_catalogue_100_{}.fits'.format(param, param)
#            data_fits = fits.open(fileName)
            #print(data_fits.info())
            
#            z_mock = data_fits['GALAXY PROPERTIES'].data['redshift'][ID-1]
#            wl_spec_mock = data_fits['FULL SED WL'].data['wl'][0]*(1+z_mock)
#            f_spec_mock = data_fits['FULL SED'].data[ID-1]/(1+z_mock)
#            lfl_spec_mock = f_spec_mock*wl_spec_mock
            
#            # PHOTOMETRY 
#            fileName = '/Users/lester/Documents/PhD/param_100/mock_100_{}/mock_catalogue_100_{}.fits'.format(param, param)
#            data_fits = fits.open(fileName)
#            print(data_fits.info())
#
#
#
#            appmag_phot_mock = np.zeros(len(filters))
#            for i in range(len(filters)):
#                appmag_phot_mock[i] = data_fits['APPARENT MAGNITUDES'].data[filters[i]][ID-1]
#                
#            data_fits.close()
#            
#            f_phot_mock_mJy = (10**( (23.9 - appmag_phot_mock) / 2.5 )) # not used [mJy]
#            lfl_test = (3631 * 10**( -appmag_phot_mock / 2.5 ) )   /    ( 3.34e4 * filter_fwhm_centre) # not used [erg cm-2 s-1]
#            lfl_phot_mock = (c / filter_fwhm_centre) * (10 ** (-(appmag_phot_mock+23.6)/2.5)) # [erg cm-2 s-1]
#            
            
            # =============================================================================
            # get INPUT perturbed fluxes and errors
            # =============================================================================
            
            fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_may_2020/0A2744C/input/astrodeep_A2744_c_subset_RF1_001.fits'
            
            data_fits = fits.open(fileName)
            # print(data_fits.info())
            # print(data_fits[1].header)
            
            # PHOTOMETRY 
            ptbf_phot_mock = np.zeros(len(filter_label))
            for i in range(len(filters)):
                ptbf_phot_mock[i] = data_fits[1].data[filter_label[i]][ID-1] # uJy
                
            ptbferr_phot_mock = np.zeros(len(filter_label))
            for i in range(len(filters)):
                ptbferr_phot_mock[i] = data_fits[1].data[filter_err[i]][ID-1] # uJy    
                
            # adding min rel error to errors: 
            min_rel_error = np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.05, 0.1, 0.1])
            ptbferr_phot_mock = ( (ptbferr_phot_mock**2) + ((min_rel_error*ptbf_phot_mock)**2) ) ** 0.5
            
            ptblfl_phot_mock = (ptbf_phot_mock/filter_fwhm_centre)*(1e-10 / 3.34) # erg cm-2 s-1
            ptblflerr_phot_mock = (ptbferr_phot_mock/filter_fwhm_centre)*(1e-10 / 3.34) # erg cm-2 s-1
            
            data_fits.close()
            
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
              
            z_fit_arr = []
            wl_spec_fit_arr = []
            lfl_spec_fit_arr = []
            
            lfl_phot_fit_arr = []
            
            chi2_fit_arr = []
            
            for i in range(samples):
                
                #here's the key line - take weighted samples from the multinest output!
                idx = np.random.choice(len(temp_probs), size=1, p=temp_probs)
            
                z_fit_arr.append(data_fits['GALAXY PROPERTIES'].data['redshift'][idx])
                wl_spec_fit_arr.append(data_fits['FULL SED WL'].data['wl'][0]*(1+z_fit_arr[i]))
                lfl_spec_fit_arr.append(data_fits['FULL SED'].data[idx][0]/(1+z_fit_arr[i]) * wl_spec_fit_arr[i])
                
            
                # PHOTOMETRY 
                appmag_phot_fit = np.zeros(len(filters))
                for i in range(len(filters)):
                    appmag_phot_fit[i] = data_fits['APPARENT MAGNITUDES'].data[filters[i]][idx]
                    
                    
                data_fits.close()
                
                lfl_phot_fit = (c / filter_fwhm_centre) * (10 ** (-(appmag_phot_fit+23.6)/2.5)) # [erg cm-2 s-1]
                lfl_phot_fit_arr.append(lfl_phot_fit)
                
                # CHI SQUARED
                chi2_fit_arr.append(data_fits['POSTERIOR PDF'].data['chi_square'][idx][0])
            
            # =============================================================================
            # averaging the samples to assist in plotting
            # =============================================================================
            
            # spec
            lfl_spec_fit_min = lfl_spec_fit_arr[0]
            lfl_spec_fit_max = lfl_spec_fit_arr[0]
            
            # phot
            lfl_phot_fit_min = lfl_phot_fit_arr[0]
            lfl_phot_fit_max = lfl_phot_fit_arr[0]
            
            
            for i in range(samples-1):
                lfl_spec_fit_min = np.minimum(lfl_spec_fit_min, lfl_spec_fit_arr[i+1])
                lfl_spec_fit_max = np.maximum(lfl_spec_fit_max, lfl_spec_fit_arr[i+1])
            
                lfl_phot_fit_min = np.minimum(lfl_phot_fit_min, lfl_phot_fit_arr[i+1])
                lfl_phot_fit_max = np.maximum(lfl_phot_fit_max, lfl_phot_fit_arr[i+1])
            
            
            # =============================================================================
            # PLOT with fitted SEDs shaded
            # =============================================================================
            
            plt.figure(figsize=(fsize, 0.5*fsize))
            plt.title(title + ', chi2 = {0:.3g}'.format(np.average(chi2_fit_arr)), fontsize=14)
            plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
            plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)
        
#            print(wl_spec_fit_arr)
#            print(wl_spec_mock)
#        
            wl_spec_mock = wl_spec_fit_arr[0]
        
        
        
            
#            plt.plot(wl_spec_mock, lfl_spec_mock, linewidth=0.5, zorder=2, label='Real SED')
            
            plt.fill_between(wl_spec_mock, lfl_spec_fit_min, lfl_spec_fit_max, alpha=0.4, color='k', linewidth=0, zorder=1, label='Output SED samples') 
            
#            plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='c', zorder=3, label='Real Phot')
            
            plt.errorbar(filter_fwhm_centre, lfl_phot_fit_min, yerr=[np.zeros(len(lfl_phot_fit_min)), lfl_phot_fit_max - lfl_phot_fit_min], linestyle="None", linewidth=10, color='r', zorder=4, label='Output Phot samples')
            
            plt.scatter(filter_fwhm_centre, ptblfl_phot_mock, color='k', marker='x', zorder=5, label='Perturbed Input Fluxes')
            plt.errorbar(filter_fwhm_centre, ptblfl_phot_mock, yerr= ptblflerr_phot_mock, linestyle="None", color='k', zorder=5)
        
#            idx_lim1 = (np.abs(wl_spec_mock - 45000)).argmin()
#            y_lim1 = lfl_spec_mock[idx_lim1]
            
#            idx_lim2 = (np.abs(wl_spec_mock - 12500)).argmin()
#            y_lim2 = lfl_spec_mock[idx_lim2]
            
            plt.xlim(0, 50000)
        #    plt.xlim(10000, 15000)
        #    plt.ylim(0.9*y_lim1, 1.5*y_lim2)
            plt.ylim(0.8*min(ptblfl_phot_mock), 1.5*max(ptblfl_phot_mock))
            
            plt.yscale('log')
            plt.legend()
            plt.show()
            
        
            # =============================================================================
            # chi_squared random comparisons
            # =============================================================================
            
        #    chi2_mock = np.sum(((ptblfl_phot_mock - lfl_phot_mock)**2) / (ptblflerr_phot_mock**2) )
        #    
        #    #print(chi2_fit_total)
        #    plt.figure(figsize=(0.2*fsize, 0.2*fsize))
        #    plt.hist(chi2_fit_total, bins=20)
        #    plt.show()
        #    
        #    #print(chi2_fit_arr)
        #    plt.figure(figsize=(0.2*fsize, 0.2*fsize))
        #    plt.hist(chi2_fit_arr, bins=20)
        #    plt.show()
        #    
        #    # just choose a number from the samples
        #    i=0
        #    
        #    test = np.sum(((lfl_phot_fit_arr[i] - ptblfl_phot_mock)**2) / (ptblflerr_phot_mock**2) )
        #    print('lester')
        #    print(test)
        #    print(chi2_fit_arr[i])
            
        
        
        
        
        
        
        
        # =============================================================================
        # PLOT mock
        # =============================================================================
        
        #plt.figure(figsize=(2*fsize, fsize))
        #plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
        #plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)
        #plt.xlim(0, 50000)
        #plt.ylim(1e-15, 1e-14)
        #plt.yscale('log')
        #plt.plot(wl_spec_mock, f_spec_mock*wl_spec_mock, linewidth=0.5, zorder=0)
        ##plt.scatter(filter_fwhm_centre, lfl_phot_mock, marker='x', color='r', zorder=2)
        #plt.scatter(filter_fwhm_centre, lfl_phot_fit, marker='x', color='k', zorder=2)
        #plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='r', zorder=1)
        ##plt.legend()
        #plt.show()
        
        
        
        # =============================================================================
        # PLOT a few SEDs
        # =============================================================================
        
        #plt.figure(figsize=(6*fsize, 3*fsize))
        #plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
        #plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)
        #plt.xlim(0, 50000)
        #plt.ylim(1e-15, 1e-14)
        #plt.yscale('log')
        #plt.plot(wl_spec_mock, f_spec_mock*wl_spec_mock, linewidth=0.5, zorder=0)
        #for i in range(samples):
        #    plt.plot(wl_spec_fit_arr[i], f_spec_fit_arr[i]*wl_spec_fit_arr[i], linewidth=0.5, zorder=1)
        #    
        ##plt.scatter(filter_fwhm_centre, lfl_phot_mock, marker='x', color='r', zorder=3)
        #plt.scatter(filter_fwhm_centre, lfl_phot_fit, marker='x', color='k', zorder=3)
        #plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='r', zorder=2)
        ##plt.legend()
        #plt.show()
        
        # =============================================================================
        # PLOT with fitted SEDs shaded
        # =============================================================================
        
        #plt.figure(figsize=(6*fsize, 2*fsize))
        #plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
        #plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)
        #plt.xlim(0, 50000)
        #plt.ylim(1e-15, 1e-14)
        #
        #plt.xlim(15000, 20000)
        #plt.ylim(2.5e-15, 3.5e-15)
        #
        #plt.plot(wl_spec_mock, f_spec_mock*wl_spec_mock, linewidth=0.5, zorder=1)
        #
        #
        ##plt.fill_between(wl_spec_mock, f_spec_fit_min*wl_spec_mock, f_spec_fit_max*wl_spec_mock) 
        #plt.plot(wl_spec_mock, f_spec_fit_min*wl_spec_mock, alpha=0.3, color='k', zorder=0) 
        #plt.plot(wl_spec_mock, f_spec_fit_max*wl_spec_mock, alpha=0.3, color='k', zorder=0) 
        #
        ##plt.scatter(filter_fwhm_centre, lfl_phot_mock, marker='x', color='r', zorder=3)
        ##plt.scatter(filter_fwhm_centre, lfl_phot_fit, marker='x', color='k', zorder=3)
        ##plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='r', zorder=2)
        ##plt.legend()
        #
        #for i in range(samples):
        #    plt.plot(wl_spec_fit_arr[i], f_spec_fit_arr[i]*wl_spec_fit_arr[i], linewidth=0.5, zorder=1)
        #    
        #plt.yscale('log')
        #plt.show()
        
        
        # =============================================================================
        # PLOT with fitted SEDs shaded
        # =============================================================================
        
        #plt.figure(figsize=(6*fsize, 2*fsize))
        #plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)
        #plt.ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)
        #plt.xlim(0, 50000)
        #plt.ylim(1e-15, 1e-14)
        #
        ##plt.xlim(10000, 20000)
        ##plt.ylim(2.5e-15, 3.5e-15)
        #
        #plt.plot(wl_spec_mock, f_spec_mock*wl_spec_mock, linewidth=0.5, zorder=1)
        #
        #
        ##plt.fill_between(wl_spec_mock, f_spec_fit_min*wl_spec_mock, f_spec_fit_max*wl_spec_mock) 
        #plt.plot(wl_spec_mock, f_spec_fit_min*wl_spec_mock, alpha=0.3, color='k', zorder=0) 
        #plt.plot(wl_spec_mock, f_spec_fit_max*wl_spec_mock, alpha=0.3, color='k', zorder=0) 
        #
        ##plt.scatter(filter_fwhm_centre, lfl_phot_mock, marker='x', color='r', zorder=3)
        ##plt.scatter(filter_fwhm_centre, lfl_phot_fit, marker='x', color='k', zorder=3)
        #plt.errorbar(filter_fwhm_centre, lfl_phot_mock, xerr=filter_fwhm/2, linestyle="None", color='c', zorder=2)
        #plt.errorbar(filter_fwhm_centre, lfl_phot_fit_min, yerr=[np.zeros(len(lfl_phot_fit_min)), lfl_phot_fit_max - lfl_phot_fit_min], linestyle="None", linewidth=10, color='r', zorder=5)
        #
        #    
        #
        ##plt.legend()
        #
        #plt.yscale('log')
        #plt.show()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        


























