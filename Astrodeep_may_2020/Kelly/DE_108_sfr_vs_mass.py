#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:26:47 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from astropy.io import fits

title = 'DE'
param = '006'
revision = '108'

size = 15
fsize = 10

# =============================================================================
# INPUT - get "real" parameters
# =============================================================================

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/param_{}/astrodeep_{}/mock_catalogue_005_010.fits'.format(param, revision)
#data_fits = fits.open(fileName)
#
#mtot_r = np.log10(data_fits['GALAXY PROPERTIES'].data['m_tot'])
#sfr_r = data_fits['STAR FORMATION'].data['SFR']
#msa_r = np.log10(data_fits['STAR FORMATION'].data['max_stellar_age'])
#tau_r = np.log10(data_fits['STAR FORMATION BINS'].data['bin_tau'])
#
#data_fits.close()

# =============================================================================
# OUTPUT - get BEAGLE parameters
# =============================================================================

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/param_{}/astrodeep_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(param, revision)

fileName = "/Users/lester/Documents/PhD/param_100/fit_{}_DE/pyp-beagle/data/BEAGLE_summary_catalogue.fits".format(revision)
data_fits = fits.open(fileName)

id_b = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1
mtot_b = data_fits['POSTERIOR PDF'].data['mass_mean']
mtot_68_b = data_fits['POSTERIOR PDF'].data['mass_68.00']
sfr_b = data_fits['STAR FORMATION'].data['SFR_mean']
sfr_68_b = data_fits['STAR FORMATION'].data['SFR_68.00']

data_fits.close()

# =============================================================================
# PLOT - input mass vs output mass
# =============================================================================

#plt.figure(figsize=(fsize, fsize))
#plt.title('Input Mass (DE) vs Output Mass ({})'.format(title), size=size)
#plt.xlabel(r'$\text{Input - log}(m_{tot}/M_{\odot})$', size=size)
#plt.ylabel(r'$\text{Output - log}(m_{tot}/M_{\odot})$', size=size)
#plt.plot((7.5, 11), (7.5, 11))
#plt.scatter(mtot_r[id_b], mtot_b, s=10, c=msa_r[id_b]-tau_r[id_b], zorder=10)
#plt.errorbar(mtot_r[id_b], mtot_b, yerr=[mtot_b - mtot_68_b[:, 0], mtot_68_b[:, 1] - mtot_b], linestyle="None", elinewidth=0.5, color='k')
#
#plt.colorbar()
#
#plt.xlim(7.5, 11)
#plt.ylim(7.5, 11)
#plt.show()

# =============================================================================
# PLOT - input sfr vs output sfr
# =============================================================================

#plt.figure(figsize=(fsize, fsize))
#plt.title('Input SFR (DE) vs Output SFR ({})'.format(title), size=size)
#plt.xlabel(r'$\text{Input - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
#plt.ylabel(r'$\text{Output - log}(\Psi / M_{\odot} yr^{-1})$', size=size)
#plt.plot((-1, 3.5), (-1, 3.5))
#plt.scatter(np.log10(sfr_r[id_b]), np.log10(sfr_b), s=10)
#plt.errorbar(np.log10(sfr_r[id_b]), np.log10(sfr_b), yerr=[np.log10(sfr_b / sfr_68_b[:, 0]), np.log10(sfr_68_b[:, 1] / sfr_b)], linestyle="None", elinewidth=0.5, color='k')
#
#plt.xlim(-1, 3.5)
#plt.ylim(-1, 3.5)
#plt.show()

# =============================================================================
# PLOT - main sequence as a 2d histogram
# =============================================================================

# input values
#plt.figure(figsize=(1.2*fsize, fsize))
#plt.title('INPUT (DE) SFR vs Mass', size=size)
#plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
#plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
#plt.hist2d(mtot_r, np.log10(sfr_r), range=[[7.5, 11], [-1, 3.5]], bins=100)
#plt.colorbar()
#plt.xlim(7.5, 11)
#plt.ylim(-1, 3.5)
#plt.show()

massh = np.empty(0)
sfrh = np.empty(0)

#for i in range(len(id_b)):
#for i in range(1):
for i in [119]:
#for i in range(len(id_b))[120:121]:

#    beagleData = fits.open('/Users/lester/Documents/PhD/param_{}/astrodeep_{}/{}_BEAGLE.fits'.format(param, revision, id_b[i]+1))
    
    beagleData = fits.open('/Users/lester/Documents/PhD/param_100/fit_{}_DE/{}_BEAGLE.fits.gz'.format(revision, id_b[i]+1))
        

    print(id_b[i]+1)
    
    #needs float64 to provide precision needed for the random.choice weights
    temp_probs = np.float64(beagleData['POSTERIOR PDF'].data['probability'])
    temp_probs = temp_probs/np.sum(temp_probs)

    #here's the key line - take weighted samples from the multinest output!
    idx = np.random.choice(len(temp_probs), size=1000, p=temp_probs)
    massh = np.append(massh, np.log10(beagleData['GALAXY PROPERTIES'].data['M_tot'][idx]))
    sfrh = np.append(sfrh, np.log10(beagleData['STAR FORMATION'].data['sfr'][idx]))


plt.figure(figsize=(1.2*fsize, fsize))
plt.title('FITTED ({}) SFR vs Mass'.format(title), size=size)
plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.hist2d(massh, sfrh, range=[[6.5, 12], [-2, 4.5]], bins=100, norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
plt.colorbar()

cmap = cm.get_cmap('viridis')
rgba = cmap(0)
ax = plt.axes()
#ax.set_facecolor(rgba)

#plt.xlim(7.5, 11)
#plt.ylim(-1, 3.5)
plt.show()
















# =============================================================================
# SEDs
# =============================================================================


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

fields = ['0A2744C']
runs = ['001']

fsize = 5
size = 8

for field in fields:    
    
    for run in runs:
        
        # =============================================================================
        # OUTPUT - get BEAGLE parameters
        # =============================================================================
        
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

        data_fits.close()
        
        # =============================================================================
        # calculate output gradient for rising or falling
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
                
        fsize=2        
        ageUniv = cd.age(redshift_b1, **cosmo)/cc.yr_s
        
        xlin = np.linspace(1, 1.1e10, 100000)
        
        IDs = ([20])
        IDs = id_b1[-3:]
        IDs=[21, 115, 3030]
        
        for ID in IDs:
            
            i = (np.abs(id_b1 - ID)).argmin()
            
            plt.figure(figsize=(4*fsize, fsize))
            plt.xlim(0, 1.1e10)
            plt.ylim(0.0, 1.1)
            
            sfr_out = 1 * (xlin-(ageUniv[i]-msa_b1[i]))*np.exp(-(xlin-(ageUniv[i]-msa_b1[i]))/tau_b1[i])
            plt.plot(xlin, sfr_out/max(sfr_out), label='OUTPUT SFH')
            
            plt.plot((ageUniv[i], ageUniv[i]), (0, 1))
            plt.legend()
            plt.show()
            
            print(10**mass_b1[i], msa_b1[i], tau_b1[i])
        

        # =============================================================================
        # SEDs for given ID
        # =============================================================================
        
        samples = 10
        
        IDs=[22, 61, 63]
        IDs = id_b1[:10]
        IDs=[120]
        
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
            # get INPUT perturbed fluxes and errors
            # =============================================================================
            
            fileName = "/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/ASTRODEEP/astrodeep_A2744_p_subset_001.fits"
            
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
            
            data_fits = fits.open('/Users/lester/Documents/PhD/param_100/fit_108_DE/{}_BEAGLE.fits.gz'.format(ID))
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

            wl_spec_mock = wl_spec_fit_arr[0]
            
            plt.fill_between(wl_spec_mock, lfl_spec_fit_min, lfl_spec_fit_max, alpha=0.4, color='k', linewidth=0, zorder=1, label='Output SED samples')
            
            plt.errorbar(filter_fwhm_centre, lfl_phot_fit_min, yerr=[np.zeros(len(lfl_phot_fit_min)), lfl_phot_fit_max - lfl_phot_fit_min], linestyle="None", linewidth=10, color='r', zorder=4, label='Output Phot samples')
            
            plt.scatter(filter_fwhm_centre[ptblfl_phot_mock>0], ptblfl_phot_mock[ptblfl_phot_mock>0], color='k', marker='x', zorder=5, label='Input Fluxes')
            plt.scatter(filter_fwhm_centre[(ptblfl_phot_mock<0)&(ptbf_phot_mock>-65)], abs(ptblfl_phot_mock[(ptblfl_phot_mock<0)&(ptbf_phot_mock>-65)]), color='lime', marker='x', zorder=5, label='Negative Input Fluxes')
            
            
            plt.errorbar(filter_fwhm_centre[ptbf_phot_mock>-65], abs(ptblfl_phot_mock)[ptbf_phot_mock>-65], yerr= ptblflerr_phot_mock[ptbf_phot_mock>-65], linestyle="None", color='k', zorder=5)
            
            plt.xlim(0, 50000)
            plt.ylim(0.8*min(abs(ptblfl_phot_mock)[ptbf_phot_mock>-65]), 1.5*max(abs(ptblfl_phot_mock)[ptbf_phot_mock>-65]))
            
            plt.yscale('log')
            plt.legend()
            plt.show()
            
            print(ptblfl_phot_mock)
            print(ptblflerr_phot_mock)


















































