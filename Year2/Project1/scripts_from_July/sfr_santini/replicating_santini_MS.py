#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""

# Santini approx 1711 sources 1.3 < z < 6.0

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'

AD = np.load(AD_location)
print(AD.dtype.names)
print(len(AD))

# =============================================================================
# H band mag cut of < 27.5
# =============================================================================

AD = AD[AD['H160']<27.5]
print(len(AD))

# =============================================================================
# Photo-z from median of 6 methods
# =============================================================================

zLow = 1.3
zHigh = 2.0

wLim = (zHigh - zLow) / 2.0
zLim = zLow + wLim

AD = AD[abs(AD['ZBEST']-zLim) < wLim]
print(len(AD))


# =============================================================================
# Need sample to be complete above given mass (magnification NOT already included)
# =============================================================================

#AD = AD[np.log10(AD['MSTAR']*1e9) > 8.3]
#print(len(AD))

#AD = AD[(np.log10(AD['MASTAR_NEB']*1e9)+np.log10(AD['MAGNIF'])) > 8.3]
AD = AD[(np.log10(AD['MASTAR_NEB']*1e9)) > 8.4]
print(len(AD))

# =============================================================================
# Adding RELFLAG just to be sure...
# =============================================================================

AD = AD[AD['RELFLAG']==1.0]
print('RELFLAG==1', len(AD))

# =============================================================================
# Use high mass cutoff according to Tomczak (between 10.2 up to 10.8 increasing with z)
# =============================================================================

#AD = AD[np.log10(AD['MSTAR']*1e9) < 10.2]
#print(len(AD))

AD = AD[np.log10(AD['MASTAR_NEB']*1e9) - np.log10(AD['MAGNIF']) < 10.2]
print(len(AD))

# =============================================================================
# get some key columns
# =============================================================================

field_AD = AD['field']
id_AD = AD['ID']
mag_AD = AD['MAGNIF']

mass_AD = np.log10(AD['MSTAR']*1e9) - np.log10(mag_AD)
mass_AD_neb = np.log10(AD['MASTAR_NEB']*1e9) - np.log10(mag_AD)
sfr_AD = np.log10(AD['SFR']) - np.log10(mag_AD)
sfr_AD_neb = np.log10(AD['SFR_NEB']) - np.log10(mag_AD)


print(field_AD)
print(id_AD)



# =============================================================================
# compare original and BEAGLE IDs for next section
# =============================================================================

fields = ['A2744_c', 'A2744_p', 'M0416_c', 'M0416_p', 'M0717_c', 'M0717_p', 'M1149_c', 'M1149_p']
#print([fields[1]])
#fields = ['A2744_c', 'M0416_c', 'M0717_c', 'M1149_c']
#fields = ['A2744_p', 'M0416_p', 'M0717_p', 'M1149_p']
#fields = fields[0]
#fields = ['A2744_c']
#fields = ['A2744_p']
#fields = ['M0416_c']
#fields = ['M0416_p']
#fields = ['M0717_c']
#fields = ['M0717_p']
#fields = ['M1149_c']
#fields = ['M1149_p']


for field in fields:
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/data/astrodeep_{}_subset_RF1_001.fits'.format(field)
    data_fits = fits.open(fileName)
    id_input = np.asarray(data_fits[1].data['ID'], dtype=int)
    field_original = np.asarray(data_fits[1].data['field'], dtype=int)
    id_original = np.asarray(data_fits[1].data['ID_original'], dtype=int)
    data_fits.close()
    
    if field == fields[0]:
        print('yes')
        id_input_arr = id_input
        field_original_arr = field_original
        id_original_arr = id_original
        
    else:
        id_input_arr = np.concatenate((id_input_arr, id_input))
        field_original_arr = np.concatenate((field_original_arr, field_original))
        id_original_arr = np.concatenate((id_original_arr, id_original))

    
# =============================================================================
# Santini SFRs (and BEAGLE masses)
# =============================================================================
        
folder_santini = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/scripts/sfr_santini/results/'

mass_BEAGLE_arr = []
sfr_santini_arr = []

for j in range(len(id_AD)):
    
    
    # id_AD and field_AD are AD catalog
    # BEAGLE IDs are needed from id_input_arr, which is same length as id_original and field_original
    # ie match id_AD and field_AD, to id and field ORIGINAL, to then get BEAGLE input id...
    # to do this I need array the size of ORIGINAL / id_input
    

    
    # find BEAGLE input ID, note id_AD and field_AD is new subset
    # field_original_arr only includes selected fields in variable "fields"
    if field_AD[j] in field_original_arr:
        
        arg_field = np.isin(field_original_arr, field_AD[j])
        arg_ID = np.isin(id_original_arr[arg_field], id_AD[j])      

        if sum(arg_ID) != 1:
            print('arg_ID ERROR', sum(arg_ID), field_AD[j], id_AD[j])
            
#        print(field_original_arr[arg_field][arg_ID], id_original_arr[arg_field][arg_ID], id_input_arr[arg_field][arg_ID])
        
        
        
        '''
    #    arg_field = np.where(field_AD[j]==field_original_arr)[0]
    #    arg_field = np.isin(field_AD[j], field_original_arr)
        arg_field = np.isin(field_original_arr, field_AD[1100])
        print(arg_field)
        print(len(arg_field))  
        print(len(arg_field[arg_field==False]))  
            
        plt.hist(field_original_arr)
        
        
    #    arg_ID = np.where(id_AD[j]==id_original_arr[arg_field])[0]
    #    arg_ID = np.isin(id_original_arr[arg_field], id_AD[j]) 
        
        
        print(len(id_original_arr), len(arg_field))
        print(len(arg_field), sum(arg_field))    
        print(id_AD[j])
        print(id_original_arr[arg_field])    
        print(arg_ID)    
       ''' 
       
   
   
   
        ID_BEAGLE_input = id_input_arr[arg_field][arg_ID]
        
        ID_santini = np.load(folder_santini+'{}_sfr_ID_santini.npy'.format(str(int(field_AD[j]))))
        mass_BEAGLE = np.load(folder_santini+'{}_mass_santini.npy'.format(str(int(field_AD[j]))))
        sfr_santini = np.load(folder_santini+'{}_sfr_santini.npy'.format(str(int(field_AD[j]))))
        sfr_err_santini = np.load(folder_santini+'{}_sfr_err_santini.npy'.format(str(int(field_AD[j]))))
    

        idx = np.isin(ID_santini, ID_BEAGLE_input)
        
        if sum(idx) == 0:
            mass_BEAGLE_arr.append(-11.0)
            sfr_santini_arr.append(-11.0)
            print('1 ERROR')
        elif sum(idx) == 1:
            mass_BEAGLE_arr.append(mass_BEAGLE[idx][0] - np.log10(mag_AD[j]))
            sfr_santini_arr.append(sfr_santini[idx][0] - np.log10(mag_AD[j]))
#            sfr_santini_arr.append(sfr_santini[idx][0] - 0)
        else:
            print('2 ERROR')
            pass
            
    else:
        mass_BEAGLE_arr.append(-12.0)
        sfr_santini_arr.append(-12.0)
#        print('0 ERROR')       

        
##    idx = np.where(ID_BEAGLE_input==ID_santini)[0]
#    idx = np.isin(ID_santini, ID_BEAGLE_input)
#    if sum(idx) == 0:
#        mass_BEAGLE_arr.append(-10.0)
#        sfr_santini_arr.append(-10.0)
#        print('0 ERROR')
#    elif sum(idx) == 1:
#        mass_BEAGLE_arr.append(mass_BEAGLE[idx][0] - np.log10(mag_AD[j]))
#        sfr_santini_arr.append(sfr_santini[idx][0] - np.log10(mag_AD[j]))
#    else:
#        print('ERROR')




print(len(id_AD), len(mass_BEAGLE_arr), len(ID_santini))


'''
Summary so far...
AD: mass_AD, mass_AD_neb, sfr_AD, sfr_AD_neb
BEAGLE: mass_BEAGLE_arr
Santini: sfr_santini_arr
'''


field_AD = np.array(field_AD)
id_AD = np.array(id_AD)
mag_AD = np.array(mag_AD)

mass_AD = np.array(mass_AD)
mass_AD_neb = np.array(mass_AD_neb)
mass_BEAGLE_arr = np.array(mass_BEAGLE_arr)
sfr_AD = np.array(sfr_AD)
sfr_AD_neb = np.array(sfr_AD_neb)
#sfr_santini_arr = np.array(sfr_santini_arr)
sfr_santini_arr = np.array(sfr_santini_arr) - np.log10(mag_AD)





# if using santini, creates santini subset
idx = (mass_BEAGLE_arr>0.0) # these are the ones which santini ID didn't match AD ID, because of...


# if not using santini, keeps full set, but need correct fields
#idx = (field_AD%2 == 0)
#idx = (field_AD%2 == 1)


print(len(mass_AD), len(mass_AD_neb), len(mass_BEAGLE_arr), len(sfr_AD), len(sfr_santini_arr), sum(idx))

#plt.hist(sfr_santini_arr)
#plt.show()
#plt.hist(mass_AD_neb)
#plt.show()


print('SANTINI CUT', len(field_AD), len(field_AD[idx]))


field_AD = field_AD[idx]
id_AD = id_AD[idx]
mag_AD = mag_AD[idx]

mass_AD = mass_AD[idx]
mass_AD_neb = mass_AD_neb[idx]
mass_BEAGLE_arr = mass_BEAGLE_arr[idx]
sfr_AD = sfr_AD[idx]
sfr_AD_neb = sfr_AD_neb[idx]
sfr_santini_arr = sfr_santini_arr[idx]












#plt.hist(sfr_santini_arr)
#plt.show()
#plt.hist(mass_AD_neb)
#plt.show()

print(min(sfr_santini_arr))

# =============================================================================
# some plots
# =============================================================================

#plt.scatter(mass_AD, sfr_AD)
#plt.show()
#
#plt.scatter(mass_BEAGLE_arr, sfr_AD)
#plt.show()
#
#plt.scatter(mass_AD, sfr_santini_arr)
#plt.show()
#
#plt.scatter(mass_BEAGLE_arr, sfr_santini_arr)
#plt.show()
#
#plt.scatter(mass_AD, mass_BEAGLE_arr)
#plt.show()
#
#plt.scatter(sfr_AD, sfr_santini_arr)
#plt.show()

# =============================================================================
# 2 sigma clipping to remove outliers
# =============================================================================

#mass = mass_BEAGLE_arr
#sfr = sfr_santini_arr

#mass = mass_AD
#sfr = sfr_AD

mass = mass_AD_neb
sfr = sfr_santini_arr


#mass = mass_AD_neb
#sfr = sfr_AD_neb



outliers = 1

mass1 = mass
sfr1 = sfr
fit1 = np.polyfit(mass1, sfr1, 1)
sfr_residuals1= sfr - (fit1[0]*mass + fit1[1])
sigma1 = np.std(sfr_residuals1)

while outliers > 0:
    
    fit = np.polyfit(mass, sfr, 1)


    sfr_residuals= sfr - (fit[0]*mass + fit[1])
  
    sigma = np.std(sfr_residuals)
    
    
    print(fit, sigma)   
    
    idx = (abs(sfr_residuals)<2.0*sigma)
    
    outliers = len(mass) - sum(idx)

    mass = mass[idx]
    sfr = sfr[idx]
    
    plt.scatter(mass, sfr)
    plt.show()
    print('LEN MASS', len(mass))


# =============================================================================
# Santini
# =============================================================================

# logSFR = alpha log(M / M_9.7) + beta

z_med_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))

alpha_san = B_san - 9.7*A_san
beta_san = A_san

alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
beta_err_san = A_err_san

x = np.linspace(7, 10.5)

plt.figure(figsize=(10, 10))
plt.scatter(mass1, sfr1)
plt.scatter(mass, sfr)
plt.plot(x, beta_san[0]*x + alpha_san[0], color='g', label='santini', linewidth=4)
plt.plot(x, fit[0]*x + fit[1], color='k', label='clipping', linewidth=4)
plt.plot(x, fit1[0]*x + fit1[1], color='r', label='no clipping', linewidth=4)
#plt.title(field.replace('_',''))
plt.xlim(7, 11)
plt.ylim(-2, 3)
plt.legend()
plt.show()

print(fit[0], fit[1])


print(beta_san[0], alpha_san[0])


#print(id_AD[66])
#print(mass_AD_neb[66])
#print(mass_AD_neb_orig[66])
#
#print(mass_AD[66])
#print((10**mass_AD[66])/1e9)
#print(field_AD[66])
#
#print(np.log10(4.075*1e9))
#
#
#print(mag_AD[66])


xlow=6.8
xhigh=11.8
ylow=-2.5
yhigh=3.6

plt.figure(figsize=(10, 10))
plt.hist2d(mass1,sfr1, bins=((xhigh-xlow)*3, (yhigh-ylow)*8))
plt.plot(x, beta_san[0]*x + alpha_san[0], color='w', label='santini', linewidth=4)
plt.plot(x, fit[0]*x + fit[1], color='k', label='clipping', linewidth=4)
plt.plot(x, fit1[0]*x + fit1[1], color='r', label='no clipping', linewidth=4)
plt.xlim(xlow, xhigh)
plt.ylim(ylow, yhigh)
plt.colorbar()
plt.legend()
plt.show()



# =============================================================================
# FIRST YEAR REPORT
# =============================================================================

#%%

plt.figure(figsize=(10, 10))
plt.xlabel(r'$M_\mathrm{tot}$')
plt.ylabel(r'$\Psi$')
plt.scatter(mass1, sfr1, marker='x')
plt.scatter(mass, sfr, marker='x')
plt.plot(x, fit1[0]*x + fit1[1], color='#1f77b4', label='Without clipping', linewidth=4)
plt.plot(x, fit[0]*x + fit[1], color='r', label='With clipping', linewidth=4)
plt.plot(x, beta_san[0]*x + alpha_san[0], color='#2ca02c', label='Santini+17', linewidth=4)
#plt.title(field.replace('_',''))
plt.xlim(7, 10.5)
plt.ylim(-2, 3)
plt.legend()
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_santini.png')
plt.show()


# =============================================================================
# bootstrapping
# =============================================================================
#%%

mass = mass_AD_neb
sfr = sfr_santini_arr

print(len(mass_AD_neb))

iterations = 10000 # 10000 for report
samples = 1000      # len(mass)== 1310

alpha = []
beta = []

for i in range(iterations):
    idx_random = np.random.choice(range(len(mass)), samples, replace=True)
#    print(np.sort(idx_random))
    mass_bs = mass[idx_random]
    sfr_bs = sfr[idx_random]
    
    outliers = 1
    
    mass1 = mass_bs
    sfr1 = sfr_bs
    fit1 = np.polyfit(mass1, sfr1, 1)
    sfr_residuals1= sfr_bs - (fit1[0]*mass_bs + fit1[1])
    sigma1 = np.std(sfr_residuals1)
    
    while outliers > 0:
        
        fit = np.polyfit(mass_bs, sfr_bs, 1)
        sfr_residuals= sfr_bs - (fit[0]*mass_bs + fit[1])  
        sigma = np.std(sfr_residuals)
#        print(fit, sigma)       
        idx = (abs(sfr_residuals)<2.0*sigma)    
        outliers = len(mass_bs) - sum(idx)
        mass_bs = mass_bs[idx]
        sfr_bs = sfr_bs[idx]    
#        plt.scatter(mass, sfr)
#        plt.show()
    # logSFR = alpha log(M / M_9.7) + beta
    
    z_med_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
    A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
    A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
    B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
    B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))
    
    alpha_san = B_san - 9.7*A_san
    beta_san = A_san
    
    alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
    beta_err_san = A_err_san
    
    x = np.linspace(min(mass_bs), max(mass_bs))
    
#    plt.figure(figsize=(10, 10))
#    plt.scatter(mass1, sfr1)
#    plt.scatter(mass_bs, sfr_bs)
#    plt.plot(x, beta_san[0]*x + alpha_san[0], color='g', label='santini', linewidth=4)
#    plt.plot(x, fit[0]*x + fit[1], color='k', label='clipping', linewidth=4)
#    plt.plot(x, fit1[0]*x + fit1[1], color='r', label='no clipping', linewidth=4)
#    #plt.title(field.replace('_',''))
#    plt.xlim(7, 11)
#    plt.ylim(-2, 3)
#    plt.legend()
#    plt.show()
    
#    print(fit[0], fit[1])
    
    
    alpha.append(fit[1])
    beta.append(fit[0])
    
#plt.title('alpha - intercept')
plt.xlabel(r'$\alpha$')
plt.ylabel('Count')
y, x, _ = plt.hist(alpha, bins=20)
plt.plot((alpha_san[0], alpha_san[0]), (0, y.max()), label='Santini+17', linewidth=4)
plt.plot((np.median(alpha), np.median(alpha)), (0, y.max()), label='Median', linewidth=4)
plt.legend()
plt.show()



#plt.title('beta - slope')
plt.xlabel(r'$\beta$')
plt.ylabel('Count')
y, x, _ = plt.hist(beta, bins=20)
plt.plot((beta_san[0], beta_san[0]), (0, y.max()), label='Santini+17', linewidth=4)
plt.plot((np.median(beta), np.median(beta)), (0, y.max()), label='Median', linewidth=4)
plt.legend()
plt.show()


# =============================================================================
# FIRST YEAR REPORT
# =============================================================================

#%%



#plt.title('alpha - intercept')
plt.xlabel(r'$\alpha$')
plt.ylabel('Count')
y, x, _ = plt.hist(alpha, bins=20)

plt.plot((np.median(alpha), np.median(alpha)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
plt.plot((alpha_san[0], alpha_san[0]), (0, y.max()), label='Santini+17', color='#2ca02c', linewidth=4)
plt.legend()
plt.tight_layout()
plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_alpha_bootstrap.png')
plt.show()

#plt.title('beta - slope')
plt.xlabel(r'$\beta$')
plt.ylabel('Count')
y, x, _ = plt.hist(beta, bins=20)

plt.plot((np.median(beta), np.median(beta)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
plt.plot((beta_san[0], beta_san[0]), (0, y.max()), label='Santini+17', color='#2ca02c', linewidth=4)
plt.legend()
plt.tight_layout()
plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_beta_bootstrap.png')
plt.show()



