#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 23:08:27 2021

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner
from astropy.io import fits
import matplotlib as mpl
from pylab import cm
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit

import scipy.stats as stats
import math

def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

filenames = ['scenario_27_clusters_z1p25-2p0_4x5000', # usechains  0 and 1, not 2 and 3
             'scenario_27_clusters_z2p0-3p0_4x5000',        
             'scenario_27_clusters_z3p0-4p0_4x5000',         
             'scenario_27_clusters_z4p0-5p0_4x5000',     
         
             'scenario_27_clusters_z1p25-6p0_4x5000',
             'scenario_27_clusters_z1p25-6p0_4x5000_linear_alphaN',
         
             'scenario_27_clusters_z1p25-2p0_4x5000_pbad0',
             'scenario_27_clusters_z2p0-3p0_4x5000_pbad0',        
             'scenario_27_clusters_z3p0-4p0_4x5000_pbad0',         
             'scenario_27_clusters_z4p0-5p0_4x5000_pbad0',     

             'scenario_27_clusters_z1p25-6p0_4x5000_linear_alphaN_pbad0',
    
             'scenario_27_clusters_z1p25-6p0_4x5000_test4',
         
             'scenario_27_clusters_z1p25-2p0_4x5000_fixed_mass_distribution',
             'scenario_27_clusters_z2p0-3p0_4x5000_fixed_mass_distribution',         
             'scenario_27_clusters_z3p0-4p0_4x5000_fixed_mass_distribution',         
             'scenario_27_clusters_z4p0-5p0_4x5000_fixed_mass_distribution',         
             'scenario_27_clusters_z1p25-6p0_4x5000_linear_alphaN_fixed_mass_distribution' 
         
         ]   

filenames = ['scenario_27_clusters_z1p25-6p0_4x5000_test4_z1']

filenames = ['scenario_27_clusters_z1p25-6p0_4x5000_test4']

filenames = ['scenario_27_clusters_z4p0-5p0_4x5000']


#filenames = ['scenario_27_clusters_z1p25-2p0_4x5000_fixed_mass_distribution',
#             'scenario_27_clusters_z2p0-3p0_4x5000_fixed_mass_distribution',         
#             'scenario_27_clusters_z3p0-4p0_4x5000_fixed_mass_distribution',         
#             'scenario_27_clusters_z4p0-5p0_4x5000_fixed_mass_distribution',         
#             'scenario_27_clusters_z1p25-6p0_4x5000_linear_alphaN_fixed_mass_distribution'    ]

import numpy as np
print(np.log10(0))



filenames = ['scenario_27_clusters_z1p25-6p0_4x2000_test4_z4_run2']
filenames = ['scenario_27_clusters_z1p25-6p0_4x2000_test4_z4_run4']


filenames = ['scenario_27_clusters_z1p25-2p0_1x16x8x8_k_fitted']
filenames = ['scenario_27_clusters_z1p25-2p0_4x5000x500x50_k_fitted']



filenames = ['scenario_29_clusters_z1p25-2p0_4x5000', 
             'scenario_29_clusters_z2p0-3p0_4x5000',        
             'scenario_29_clusters_z3p0-4p0_4x5000',         
             'scenario_29_clusters_z4p0-5p0_4x5000',     
             'scenario_29_clusters_z5p0-6p0_4x5000', 

             'scenario_29_clusters+parallels_z1p25-2p0_4x5000', 
             'scenario_29_clusters+parallels_z2p0-3p0_4x5000',        
             'scenario_29_clusters+parallels_z3p0-4p0_4x5000',         
             'scenario_29_clusters+parallels_z4p0-5p0_4x5000',     
             'scenario_29_clusters+parallels_z5p0-6p0_4x5000',
             ]  

filenames = ['scenario_29_clusters_z1p25-6p0_4x5000_ssfr_alpha', 
             'scenario_29_clusters+parallels_z1p25-6p0_4x5000_ssfr_alpha',
             ]  


filenames = ['scenario_29_clusters_z1p25-6p0_4x5000_z1-2', 
             'scenario_29_clusters_z1p25-6p0_4x5000_z2-3',        
             'scenario_29_clusters_z1p25-6p0_4x5000_z3-4',         
             'scenario_29_clusters_z1p25-6p0_4x5000_z4-5',     
             'scenario_29_clusters_z1p25-2p0_4x5000_k_fitted'
             ]

filenames = ['scenario_29_clusters_z1p25-6p0_4x5000_z3-4']


filenames = ['scenario_29_clusters_z1p25-2p0_4x5000', 
             'scenario_29_clusters_z2p0-3p0_4x5000',        
             'scenario_29_clusters_z3p0-4p0_4x5000',         
             'scenario_29_clusters_z4p0-5p0_4x5000',     
             'scenario_29_clusters_z5p0-6p0_4x5000', 
             ]  

#filenames = ['scenario_29_clusters_z1p25-6p0_1x160_ssfr_alpha_fixed_pbad']
#filenames = ['scenario_29_clusters+parallels_z1p25-6p0_1x16_ssfr_alpha_fixed_mass_distribution']

filenames = ['scenario_29_clusters_z1p25-6p0_4x5000_ssfr_alpha_fixed_pbad']

filenames = ['scenario_29_clusters_z1p25-6p0_4x5000_ssfr_alpha']

filenames = ['scenario_29_clusters_z3p0-6p0_4x5000_z4-5']
filenames = ['scenario_29_clusters_z1p25-6p0_4x5000_z4-5']

filenames = ['scenario_mock_002_4x80_z4-5_001']
filenames = ['scenario_mock_003_4x80_z4-5_002']
filenames = ['scenario_mock_003_4x200_z4-5_002']


filenames = ['scenario_29_clusters_z1p25-2p0_4x10000']
filenames = ['scenario_29_clusters_z1p25-2p0_4x10000_k_fitted']

filenames = ['scenario_29_clusters_z1p25-2p0_4x50000']


filenames = ['scenario_29_clusters_z1p25-2p0_4x32_true_linear_alpha']


# checking 1-2

#filenames = ['scenario_29_clusters_z1p25-2p0_4x5000']
#filenames = ['scenario_29_clusters_z1p25-2p0_4x10000']
#filenames = ['scenario_29_clusters_z1p25-2p0_4x50000']
filenames = ['scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha']

#filenames = ['scenario_29_clusters_z1p25-2p0_4x5000_k_fitted']
#filenames = ['scenario_29_clusters_z1p25-2p0_4x10000_k_fitted']

filenames = ['scenario_29_clusters_z1p25-2p0_4x32_true_constant_alpha_update_1d_norm']

filenames = ['scenario_29_clusters_z1p25-2p0_4x50000_true_constant_alpha_update_1d_norm']


filenames = ['scenario_29_clusters_z1p25-2p0_4x20000_true_constant_alpha_mass_sfr_z_correction']

filenames = ['scenario_29_clusters_z1p25-2p0_4x2000_true_linear_alpha_mass_sfr_z_correction_alpha_beta_covprop']

filenames = ['scenario_29_clusters_z1p25-2p0_4x2000_true_linear_alpha_mass_sfr_z_correction_alpha_beta_covprop_propscale']

filenames = ['scenario_29_clusters_z1p25-2p0_4x200_true_linear_alpha_mass_sfr_z_correction_alpha_beta_covprop_propscale']

filenames = ['scenario_29_clusters_z1p25-2p0_4x20_true_linear_alpha_mass_sfr_z_correction_alpha_beta_covprop_propscale']


filenames = ['scenario_29_clusters_z4p0-5p0_4x2000_true_linear_alpha_mass_sfr_z_correction_alpha_beta_covprop_propscale']


filenames = ['scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_mass_sfr_z_correction_alpha_beta_covprop_propscale']

filenames = ['scenario_29_clusters_z4p0-5p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_gmm']



filenames = ['scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all', 
             'scenario_29_clusters_z2p0-3p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all', 
             'scenario_29_clusters_z3p0-4p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all',
             'scenario_29_clusters_z4p0-5p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all',
             'scenario_29_clusters_z5p0-6p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all']


filenames = ['scenario_29_clusters_z1p25-2p0_4x5000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_mock']


filenames = ['scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_lm8p5']

filenames = ['scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_lm9p0']

#%%
for f, filename in enumerate(filenames):

    chain_original = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_{}.p'.format(filename),'r'))
    names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']
#    names = ['pi', 'mu', 'tausqr']
    names = chain_original.dtype.names
    nChains=4
    minIter=len(chain_original)/nChains
    
    burn=int(0.5*minIter)
#    burn=0
    cap=minIter
    #cap=1900
    chain_arr = []
    for i in range(nChains):
        print(i)
        start = minIter*i+burn
        finish = (minIter*(i+1))-(minIter-cap)


        if filename == 'scenario_27_clusters_z1p25-2p0_4x5000':
            if i <= 1:
                chain_arr.append(chain_original[start:finish])      

        elif filename == 'scenario_27_clusters_z1p25-6p0_4x5000_test4_z1':
            if i != 1:
                chain_arr.append(chain_original[start:finish])    

        elif filename == 'scenario_27_clusters_z1p25-2p0_4x5000x500x50_k_fitted':
            if i != 3:
                chain_arr.append(chain_original[start:finish])   

        elif filename == 'scenario_29_clusters_z1p25-6p0_4x5000_z3-4':
            if i == 0 or i == 3:
                chain_arr.append(chain_original[start:finish])

        elif filename == 'scenario_29_clusters_z1p25-2p0_4x10000_k_fitted':
            if i == 1 or i == 2:
                chain_arr.append(chain_original[start:finish])
                
        else:
            if i < 20:
                chain_arr.append(chain_original[start:finish])                
            
            
    chain = np.concatenate(chain_arr)

    for name in names:            
        plt.figure(figsize=(16, 4))
        plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
        plt.plot(chain[name])
        plt.show()
    
    print(' ')
    dic = {}    
    for name in names:
        dic[name+'_16'] = []
        dic[name] = []
        dic[name+'_84'] = []
    
    for name in names:            
        dic[name+'_16'].append(float('{:.3f}'.format(np.percentile(chain[name], 16))))
        dic[name].append(float('{:.3f}'.format(np.median(chain[name]))))
        dic[name+'_84'].append(float('{:.3f}'.format(np.percentile(chain[name], 84))))                
    
    for name in names: 
        print(dic[name]) 


    pickle.dump(chain, open('/Users/lester/Documents/linmix_files/lm_chain_{}.p'.format(filename),'w'))


plt.plot(chain['zeta'][:,0])
#print(chain['zeta'][:,0])

#%%

# =============================================================================
# propscales
# =============================================================================

for f, filename in enumerate(filenames):

    chain_original = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_{}_{}.p'.format(filename, 'propscales'),'r'))
#    names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']
#    names = ['pi', 'mu', 'tausqr']
    names = chain_original.dtype.names
    print(names)

    for name in names:            
        plt.figure(figsize=(16, 4))
        plt.title('{} {}'.format(filename, name).replace('_',' '))
        plt.plot(chain_original[name])
        plt.show()



#%%
        
print(chain['pi'][0, :, 0])
print(len(chain['pi'][0, :, 0]))
plt.plot(chain['pi'][:, 3, 0])
plt.show()


for i in range(len(chain['pi'][0, :, 0])):
    pass

for i in range(100):
    plt.plot(chain['pi'][:, i, 0], color='blue')
    plt.plot(chain['pi'][:, i, 1], color='orange')
    plt.plot(chain['pi'][:, i, 2], color='green')
plt.show()

#%%
# =============================================================================
# SLOPE
# =============================================================================

plots = np.array([(np.linspace(1.25,2.0,100), 0.847, 0.0),
                  (np.linspace(2.0,3.0,100), 0.904, 0.0),
                  (np.linspace(3.0,4.0,100), 0.794, 0.0),
                  (np.linspace(4.0,5.0,100), 0.911, 0.0),
                  (np.linspace(1.25,6.0,100), 0.782, 0.0),
                  (np.linspace(1.25,6.0,100), 0.762, 0.0),
                  
                  (np.linspace(1.25,2.0,100), 0.961, 0.0),
                  (np.linspace(2.0,3.0,100), 0.862, 0.0),
                  (np.linspace(3.0,4.0,100), 0.729, 0.0),
                  (np.linspace(4.0,5.0,100), 0.953, 0.0),          
                  (np.linspace(1.25,6.0,100), 0.835, 0.0),
                  
                  (np.linspace(4.0,5.0,100), 0.606, 0.0),
                  
                  (np.linspace(1.25,2.0,100), 0.826, 0.0),
                  (np.linspace(2.0,3.0,100), 0.871, 0.0),
                  (np.linspace(3.0,4.0,100), 0.682, 0.0),
                  (np.linspace(4.0,5.0,100), 0.828, 0.0),          
                  (np.linspace(1.25,6.0,100), 0.725, 0.0), 
                  
                  (np.linspace(1.25,2.0,100), 0.835, 0.0),                
                  
                  (np.linspace(1.25,2.0,100), 0.808, 0.0),
                  (np.linspace(2.0,3.0,100), 0.922, 0.0),
                  (np.linspace(3.0,4.0,100), 0.843, 0.0),
                  (np.linspace(4.0,5.0,100), 0.84, 0.0),              
                  (np.linspace(5.0,6.0,100), 1.289, 0.0)                     
                  
])

for p in range(len(plots)):
    
    if p <= 4:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='k')
        
    elif p == 5:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='b')
        
    elif p <= 9:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='r')
    
    elif p == 10:    
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='r')

    elif p == 11 or p == 17:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='g')
        
    elif p <= 16:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='k', linestyle='dashed')
               
    elif p <= 21:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='orange')       
        
        

      
plt.show()
        
          
        
        

    
# =============================================================================
# NORMALISATION 
# =============================================================================

plots = np.array([(np.linspace(1.25,2.0,100), 1.805, 0.0),
                  (np.linspace(2.0,3.0,100), 2.228, 0.0),
                  (np.linspace(3.0,4.0,100), 2.62, 0.0),
                  (np.linspace(4.0,5.0,100), 9.344, 0.0),
                  (np.linspace(1.25,6.0,100), 0.138, 2.194),
                  (np.linspace(1.25,6.0,100), 0.273, 0.286),
                  
                  (np.linspace(1.25,2.0,100), 3.477, 0.0),
                  (np.linspace(2.0,3.0,100), 2.926, 0.0),
                  (np.linspace(3.0,4.0,100), 4.964, 0.0),
                  (np.linspace(4.0,5.0,100), 10.108, 0.0),              
                  (np.linspace(1.25,6.0,100), 0.514, 0.298),
                  
                  (np.linspace(4.0,5.0,100), 9.21, 0.0),

                  (np.linspace(1.25,2.0,100), 1.772, 0.0),
                  (np.linspace(2.0,3.0,100), 2.053, 0.0),
                  (np.linspace(3.0,4.0,100), 2.353, 0.0),
                  (np.linspace(4.0,5.0,100), 9.71, 0.0),              
                  (np.linspace(1.25,6.0,100), 0.248, 0.282),

                  (np.linspace(1.25,2.0,100), 1.699, 0.0),
                  
                  (np.linspace(1.25,2.0,100), 0.869, 0.0),
                  (np.linspace(2.0,3.0,100), 1.048, 0.0),
                  (np.linspace(3.0,4.0,100), 1.151, 0.0),
                  (np.linspace(4.0,5.0,100), 1.606, 0.0),              
                  (np.linspace(5.0,6.0,100), 2.42, 0.0)                  

])

for p in range(len(plots)):
    
    if p <= 4:
        plt.plot(plots[p][0], np.log10(plots[p][1] * ((1+plots[p][0])**plots[p][2])) + 0.7, color='k')
        
    elif p == 5:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='b')

    elif p <= 9:
        plt.plot(plots[p][0], np.log10(plots[p][1] * ((1+plots[p][0])**plots[p][2])) + 0.7, color='r')
        
    elif p == 10:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='r')
        
    elif p == 11 or p == 17:
        plt.plot(plots[p][0], np.log10(plots[p][1] * ((1+plots[p][0])**plots[p][2])) + 0.7, color='g')

    elif p <= 15:
        plt.plot(plots[p][0], np.log10(plots[p][1] * ((1+plots[p][0])**plots[p][2])) + 0.7, color='k', linestyle='dashed')

    elif p == 16:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='k', linestyle='dashed')
        
        
    elif p <= 21:
        plt.plot(plots[p][0], plots[p][1] + (plots[p][0]*plots[p][2]), color='orange')       
        
        
plt.show()
        
        


#%%
# =============================================================================
# fixing pbad from redshift bins
# =============================================================================

z_means = np.array((1.25+2.0, 2.0+3.0, 3.0+4.0, 4.0+5.0, 5.0+6.0))/2.0
pbads = np.array((0.189, 0.237, 0.269, 0.022, 0.109))

def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

popt, pcov = curve_fit(f, z_means[:3], pbads[:3]) # your data x, y to fit

popt[0], popt[1]


plt.scatter(z_means, pbads)
plt.plot((1.25,4.0),(1.25*popt[0]+popt[1],4.0*popt[0]+popt[1]), color='k')
plt.plot((4.0,6.0),(0.0,0.0), color='k')
plt.show

print('slope:{} intercept:{}'.format(popt[0], popt[1]))
        
# slope:0.0424142011813 intercept:0.123863905324




z = 1.25 + 4.75*np.random.rand(1000)
pbad = np.where(z>4, 0.0, z)
pbad = np.where(z<=4, 0.0424142011813*z + 0.123863905324, pbad)

plt.hist(z)
plt.show()

plt.hist(pbad, bins=20)
plt.show()


#%%
# =============================================================================
# fixing mass distributions from redshift bins
# =============================================================================

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z5p0-6p0_4x5000.p','r'))

names = ['pi', 'mu', 'tausqr']
#names = ['pi']


for name in names:            
    plt.figure(figsize=(16, 4))
    plt.plot(np.sort(chain[name]))
#    plt.ylim(0, 1)
    plt.show()
    

dic = {}    
for name in names:
    print(np.median(np.sort(chain[name]), axis=0))

#%%

'''
z1.25 to z2.0
[0.13548469 0.33323309 0.52685976]
[8.13009999 8.53596604 8.78357908]
[0.03308468 0.17993438 0.40268562]

z2.0 to z3.0
[0.0993328  0.32233578 0.55496863]
[8.3466861  8.56228457 8.92190759]
[0.09806075 0.18239447 0.37332221]

z3.0 to z4.0
[0.0910511  0.28266587 0.60341715]
[8.62193205 8.91168967 9.28067668]
[0.09491375 0.18239957 0.41371418]

z4.0 to z5.0
[0.08804127 0.27108384 0.61970209]
[8.72618459 8.94678562 9.22033123]
[0.08542093 0.18509107 0.49231697]

z5.0 to z6.0
[0.08279497 0.26805673 0.626205  ]
[8.86066061 8.93750707 9.01095146]
[0.00545111 0.01416271 0.04848477]
'''



import scipy.stats as stats
import math

pi = np.array(((0.13548469, 0.33323309, 0.52685976), 
              (0.0993328,  0.32233578, 0.55496863), 
              (0.0910511,  0.28266587, 0.60341715), 
              (0.08804127, 0.27108384, 0.61970209), 
              (0.08279497, 0.26805673, 0.626205)))

mu = np.array(((8.13009999, 8.53596604, 8.78357908), 
              (8.3466861,  8.56228457, 8.92190759), 
              (8.62193205, 8.91168967, 9.28067668), 
              (8.72618459, 8.94678562, 9.22033123), 
              (8.86066061, 8.93750707, 9.01095146)))

tausqr = np.array(((0.03308468, 0.17993438, 0.40268562), 
              (0.09806075, 0.18239447, 0.37332221), 
              (0.09491375, 0.18239957, 0.41371418), 
              (0.08542093, 0.18509107, 0.49231697), 
              (0.00545111, 0.01416271, 0.04848477)))


for z in range(5):
    for i in range(len(pi[z])):
        x = np.linspace(5, 12, 100)
        plt.plot(x, pi[z][i]*stats.norm.pdf(x, mu[z][i], math.sqrt(tausqr[z][i])))
    
    plt.plot(x, pi[z][0]*stats.norm.pdf(x, mu[z][0], math.sqrt(tausqr[z][0])) + pi[z][1]*stats.norm.pdf(x, mu[z][1], math.sqrt(tausqr[z][1])) + pi[z][2]*stats.norm.pdf(x, mu[z][2], math.sqrt(tausqr[z][2])))
    plt.show()

plt.figure(figsize=(10,10))
for z in range(5):
    x = np.linspace(5, 12, 100)
    plt.plot(x, pi[z][0]*stats.norm.pdf(x, mu[z][0], math.sqrt(tausqr[z][0])) + pi[z][1]*stats.norm.pdf(x, mu[z][1], math.sqrt(tausqr[z][1])) + pi[z][2]*stats.norm.pdf(x, mu[z][2], math.sqrt(tausqr[z][2])), label='z{}-{}'.format(z+1, z+2))
plt.legend()
plt.show()





pi = np.array(((0.13548469, 0.33323309, 0.52685976), 
              (0.0993328,  0.32233578, 0.55496863), 
              (0.0910511,  0.28266587, 0.60341715), 
              (0.08804127, 0.27108384, 0.61970209), 
              (0.08279497, 0.26805673, 0.626205)))

mu = np.array(((8.13009999, 8.53596604, 8.78357908), 
              (8.3466861,  8.56228457, 8.92190759), 
              (8.62193205, 8.91168967, 9.28067668), 
              (8.72618459, 8.94678562, 9.22033123), 
              (8.86066061, 8.93750707, 9.01095146)))

tausqr = np.array(((0.03308468, 0.17993438, 0.40268562), 
              (0.09806075, 0.18239447, 0.37332221), 
              (0.09491375, 0.18239957, 0.41371418), 
              (0.08542093, 0.18509107, 0.49231697), 
              (0.00545111, 0.01416271, 0.04848477)))



z = 1.25 + 4.75*np.random.rand(10)

pi = np.empty((len(z),3))
mu = np.empty((len(z),3))
tausqr = np.empty((len(z),3))


pi[:,0] = np.where(z<=2, 0.13548469, pi[:,0])
pi[:,0] = np.where((z>2)&(z<=3), 0.0993328, pi[:,0])
pi[:,0] = np.where(z>3, 0.0910511, pi[:,0])

pi[:,1] = np.where(z<=2, 0.33323309, pi[:,1])
pi[:,1] = np.where((z>2)&(z<=3), 0.32233578, pi[:,1])
pi[:,1] = np.where(z>3, 0.28266587, pi[:,1])

pi[:,2] = np.where(z<=2, 0.52685976, pi[:,2])
pi[:,2] = np.where((z>2)&(z<=3), 0.55496863, pi[:,2])
pi[:,2] = np.where(z>3, 0.60341715, pi[:,2])

pi = pi/np.sum(pi, axis=1)[:, np.newaxis]


mu[:,0] = np.where(z<=2, 8.13009999, mu[:,0])
mu[:,0] = np.where((z>2)&(z<=3), 8.3466861, mu[:,0])
mu[:,0] = np.where(z>3, 8.62193205, mu[:,0])

mu[:,1] = np.where(z<=2, 8.53596604, mu[:,1])
mu[:,1] = np.where((z>2)&(z<=3), 8.56228457, mu[:,1])
mu[:,1] = np.where(z>3, 8.91168967, mu[:,1])

mu[:,2] = np.where(z<=2, 8.78357908, mu[:,2])
mu[:,2] = np.where((z>2)&(z<=3), 8.92190759, mu[:,2])
mu[:,2] = np.where(z>3, 9.28067668, mu[:,2])

mu = mu/np.sum(mu, axis=1)[:, np.newaxis]


tausqr[:,0] = np.where(z<=2, 0.03308468, tausqr[:,0])
tausqr[:,0] = np.where((z>2)&(z<=3), 0.09806075, tausqr[:,0])
tausqr[:,0] = np.where(z>3, 0.09491375, tausqr[:,0])

tausqr[:,1] = np.where(z<=2, 0.17993438, tausqr[:,1])
tausqr[:,1] = np.where((z>2)&(z<=3), 0.18239447, tausqr[:,1])
tausqr[:,1] = np.where(z>3, 0.18239957, tausqr[:,1])

tausqr[:,2] = np.where(z<=2, 0.40268562, tausqr[:,2])
tausqr[:,2] = np.where((z>2)&(z<=3), 0.37332221, tausqr[:,2])
tausqr[:,2] = np.where(z>3, 0.41371418, tausqr[:,2])

tausqr = tausqr/np.sum(tausqr, axis=1)[:, np.newaxis]


print(pi)
print(mu)
print(tausqr)

#%%
## =============================================================================
## Corner Plots
## =============================================================================

        

filenames = ['scenario_27_clusters_z1p25-2p0_4x5000x500x50_k_fitted']
filenames = ['scenario_29_clusters_z1p25-2p0_4x50000']
filenames = ['scenario_29_clusters_z1p25-2p0_4x50000_true_constant_alpha']
filenames = ['scenario_29_clusters_z1p25-2p0_4x20000_true_constant_alpha_mass_sfr_z_correction']
filenames = ['scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_mass_sfr_z_correction_alpha_beta_covprop_propscale']

for f, filename in enumerate(filenames):
    chain_corner = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_{}.p'.format(filename),'r'))

    if filename in ['scenario_27_clusters_z1p25-2p0_4x5000', # usechains  0 and 1, not 2 and 3
                    'scenario_27_clusters_z2p0-3p0_4x5000',        
                    'scenario_27_clusters_z3p0-4p0_4x5000',         
                    'scenario_27_clusters_z4p0-5p0_4x5000', 
                    'scenario_27_clusters_z1p25-2p0_4x5000_fixed_mass_distribution',
                    'scenario_27_clusters_z2p0-3p0_4x5000_fixed_mass_distribution',         
                    'scenario_27_clusters_z3p0-4p0_4x5000_fixed_mass_distribution',         
                    'scenario_27_clusters_z4p0-5p0_4x5000_fixed_mass_distribution',
                    'scenario_27_clusters_z1p25-6p0_4x5000_test4',
                    'scenario_29_clusters_z1p25-2p0_4x50000',
                    'scenario_29_clusters_z1p25-2p0_4x50000_true_constant_alpha',
                    'scenario_29_clusters_z1p25-2p0_4x20000_true_linear_alpha_mass_sfr_z_correction',
                    'scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_mass_sfr_z_correction_alpha_beta_covprop_propscale']:
        
        names_plot = ['ssfr a', 'beta a', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
        data_corner = np.array([chain_corner['alphaN_a'][:], chain_corner['beta_a'][:], chain_corner['sig0'][:], chain_corner['pbad'][:], chain_corner['outlier_mean'][:], chain_corner['outlier_sigma'][:]]).T
    
    elif filename in ['scenario_27_clusters_z1p25-2p0_4x5000_pbad0',
                    'scenario_27_clusters_z2p0-3p0_4x5000_pbad0',        
                    'scenario_27_clusters_z3p0-4p0_4x5000_pbad0',         
                    'scenario_27_clusters_z4p0-5p0_4x5000_pbad0']:
        
        names_plot = ['ssfr a', 'beta a', 'sig0', 'outlier mean', 'outlier sigma']
        data_corner = np.array([chain_corner['alphaN_a'], chain_corner['beta_a'], chain_corner['sig0'], chain_corner['outlier_mean'], chain_corner['outlier_sigma']]).T
    
    elif filename in ['scenario_27_clusters_z1p25-6p0_4x5000',
                    'scenario_27_clusters_z1p25-6p0_4x5000_linear_alphaN',
                    'scenario_27_clusters_z1p25-6p0_4x5000_linear_alphaN_fixed_mass_distribution']:
        
        names_plot = ['ssfr a', 'ssfr b', 'beta a', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
        data_corner = np.array([chain_corner['alphaN_a'], chain_corner['alphaN_b'], chain_corner['beta_a'], chain_corner['sig0'], chain_corner['pbad'], chain_corner['outlier_mean'], chain_corner['outlier_sigma']]).T
    
    elif filename in ['scenario_27_clusters_z1p25-6p0_4x5000_linear_alphaN_pbad0']:
        
        names_plot = ['ssfr a', 'ssfr b', 'beta a', 'sig0', 'outlier mean', 'outlier sigma']
        data_corner = np.array([chain_corner['alphaN_a'], chain_corner['alphaN_b'], chain_corner['beta_a'], chain_corner['sig0'], chain_corner['outlier_mean'], chain_corner['outlier_sigma']]).T
        
    elif filename in ['scenario_27_clusters_z1p25-2p0_4x5000x500x50_k_fitted']:
        
        names_plot = ['ssfr a', 'beta a', 'sig0', 'k', 'outlier mean', 'outlier sigma']
        data_corner = np.array([chain_corner['alphaN_a'], chain_corner['beta_a'], chain_corner['sig0'], chain_corner['k'], chain_corner['outlier_mean'], chain_corner['outlier_sigma']]).T
        
        
    
    figure = corner.corner(data_corner, labels=names_plot,
                           quantiles=[0.16, 0.5, 0.84],
                           show_titles=True, title_kwargs={"fontsize": 20})






#%%
# =============================================================================
# heatplots
# =============================================================================
        
if False:

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_27_clusters_z1p25-2p0.fits'
    s27_clusters_z1p25_2p0 = fits.open(fileName)[1].data
            
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_27_clusters_z2p0-3p0.fits'
    s27_clusters_z2p0_3p0 = fits.open(fileName)[1].data

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_27_clusters_z3p0-4p0.fits'
    s27_clusters_z3p0_4p0 = fits.open(fileName)[1].data

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_27_clusters_z4p0-5p0.fits'
    s27_clusters_z4p0_5p0 = fits.open(fileName)[1].data

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_27_clusters_z1p25-6p0.fits'
    s27_clusters_z1p25_6p0 = fits.open(fileName)[1].data            



#%%
# =============================================================================
# s26 loop for parallels and clusters
# =============================================================================


# Edit the font, font size, and axes width
mpl.rcParams['font.family'] = 'Arial'
mpl.rc('image', cmap='jet')
cmap = mpl.cm.get_cmap('jet')
plt.rcParams['font.size'] = 24
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['text.usetex'] = True
fontsize_legend = 12
fontsize_axes = 20
figuresize = 7

s=s27_clusters_z3p0_4p0
#s=s27_clusters_z1p25_6p0
z_med_hp_low = 3.0
z_med_hp_high = 4.0

z_med_hp = (z_med_hp_low+z_med_hp_high)/2.0
z_med_hp_gap = (z_med_hp_low+z_med_hp_high)/2.0 - z_med_hp_low

idx_rdm = np.arange(len(s))[:] 
print(len(idx_rdm))

x_hp = np.array([])
y_hp = np.array([])
z_hp = np.array([])

n_hp = 300 # number of samples to take from GMM in total

for i in idx_rdm:

    x_temp = np.array([])
    y_temp = np.array([])
    z_temp = np.array([])

    for G in range(3):
        
        mean = np.array([s['x_GMM_3d'][i,G],s['y_GMM_3d'][i,G],s['z_GMM_3d'][i,G]])
        cov = np.array([[np.power(s['xsig_GMM_3d'][i,G],2), s['xycov_GMM_3d'][i,G], s['xzcov_GMM_3d'][i,G]],[s['xycov_GMM_3d'][i,G], np.power(s['ysig_GMM_3d'][i,G],2), s['yzcov_GMM_3d'][i,G]],[s['xzcov_GMM_3d'][i,G], s['yzcov_GMM_3d'][i,G], np.power(s['zsig_GMM_3d'][i,G],2)]])

        xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s['amp_GMM_3d'][i,G]))

        x_hp = np.concatenate((x_hp,xyz[:,0]))
        y_hp = np.concatenate((y_hp,xyz[:,1]))
        z_hp = np.concatenate((z_hp,xyz[:,2]))

        x_temp = np.concatenate((x_temp,xyz[:,0]))
        y_temp = np.concatenate((y_temp,xyz[:,1]))
        z_temp = np.concatenate((z_temp,xyz[:,2]))
        
    plt.scatter(x_temp,y_temp,c=z_temp)
    plt.title(i)
    plt.colorbar()
    plt.show()
        
# only keep GMM samples within the redshift bin
x_hp = x_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]
y_hp = y_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]
z_hp = z_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]

y_temp = y_hp 

# only keep GMM samples with sfr between -3 and 3
#x_hp = x_hp[abs(y_temp) < 3.0]
#y_hp = y_hp[abs(y_temp) < 3.0]
#z_hp = z_hp[abs(y_temp) < 3.0]


fig = plt.figure(figsize=(0.7*figuresize, 0.7*figuresize))

ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]
#    ax1.set_title('{} - {}'.format(z_med_hp_low, z_med_hp_high))

xlow = 5.5
xhigh = 11.5
ylow = -3.5
yhigh = 3.5

ximin = 8.5
ximax = 10.0

z_s = np.linspace(z_med_hp_low,z_med_hp_high,1000)

h = ax1.hist2d(x_hp, y_hp, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)

ax1.plot((9.7,9.7), (ylow, yhigh), color='w')


# fit straight line 
popt, pcov = curve_fit(f, x_hp, y_hp) # your data x, y to fit
ax1.plot((xlow,xhigh), (popt[1] + popt[0]*xlow, popt[1] + popt[0]*xhigh), color='w')
ax1.title.set_text('slope:{:.2f} intercept:{:.2f}\n{}-{}'.format(popt[0], popt[1]+9.7*popt[0], z_med_hp_low, z_med_hp_high))

#8.5 and 10
ax1.set_xlim(xlow, xhigh)
ax1.set_xlabel(r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(ylow, yhigh)
ax1.set_ylabel(r'$\log(\psi \, / \, \mathrm{M_{\odot} \, yr^{-1}})$', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(h[3], cax=cbaxes)
cb.set_ticks([3, 5, 10, 20])
cb.set_ticklabels([3, 5, 10, 20])
#cbaxes.set_ylabel(r'TEST', rotation=270, labelpad=30)
#cbaxes.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=1))
cbaxes.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
cbaxes.yaxis.set_tick_params(which='minor', size=0)
cbaxes.yaxis.set_tick_params(labelsize=fontsize_axes)

cmap = cm.get_cmap('viridis')
rgba = cmap(0)
ax1.set_facecolor(rgba)




#%%
        
plt.hist(x_hp)
plt.show()
plt.hist(y_hp[y_hp<-3.0])
plt.show()
plt.hist(y_hp[abs(y_hp)<3.0])
plt.show()
plt.hist(y_hp[y_hp>3.0])
plt.show()





# =============================================================================
# investigation which objects scatter into z4-5 from full
# =============================================================================

scenarioA = '27'
field = 'clusters'
z_bin = 'z1p25-6p0'
data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, field, z_bin),'r'))
#%%
for f, filename in enumerate(filenames):

    chain_original = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_{}.p'.format(filename),'r'))
#    names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']
#    names = ['pi', 'mu', 'tausqr']
    names = ['zeta']
#    names = chain_original.dtype.names
    nChains=4
    minIter=len(chain_original)/nChains
    
    burn=int(0.5*minIter)
#        burn=0
    cap=minIter
    #cap=1900
    chain_arr = []
    for i in range(nChains):
        print(i)
        start = minIter*i+burn
        finish = (minIter*(i+1))-(minIter-cap)

        if filename == 'scenario_27_clusters_z1p25-2p0_4x5000':
            if i <= 1:
                chain_arr.append(chain_original[start:finish])      

        elif filename == 'scenario_27_clusters_z1p25-6p0_4x5000_test4_z1':
            if i != 1:
                chain_arr.append(chain_original[start:finish])    
                
        else:
            if i < 20:
                chain_arr.append(chain_original[start:finish])                
            
            
    chain = np.concatenate(chain_arr)



#%%

    chain_start = 2400
    chain_end = 520
    object_start = 0
    object_end = 979 #979 max

    for name in ['zeta']:            
        plt.figure(figsize=(12, 12))
        plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
        for i in range(len(data['redshift_BEAGLE'][object_start:object_end])):
            if data['redshift_BEAGLE'][i] > 4.0 and data['redshift_BEAGLE'][i] < 5.0:
                plt.plot(chain[name][chain_start:chain_end,i], color='green')
            elif data['redshift_BEAGLE'][i] > 3.5 and data['redshift_BEAGLE'][i] < 5.5:
                plt.plot(chain[name][chain_start:chain_end,i], color='orange')
            else:
                plt.plot(chain[name][chain_start:chain_end,i], color='red')
        plt.ylim(1.5, 7.5)
        plt.show()
   
    for name in ['xi']:            
        plt.figure(figsize=(12, 12))
        plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
        for i in range(len(data['redshift_BEAGLE'][object_start:object_end])):
            if data['redshift_BEAGLE'][i] > 4.0 and data['redshift_BEAGLE'][i] < 5.0:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='green')
            elif data['redshift_BEAGLE'][i] > 3.5 and data['redshift_BEAGLE'][i] < 5.5:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='orange')
            elif data['redshift_BEAGLE'][i] > 5.5:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='black')
            else:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='red')
        plt.ylim(7, 11)
        plt.show()
    
    for name in ['eta']:            
        plt.figure(figsize=(12, 12))
        plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
        for i in range(len(data['redshift_BEAGLE'][object_start:object_end])):
            if data['redshift_BEAGLE'][i] > 4.0 and data['redshift_BEAGLE'][i] < 5.0:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='green')
            elif data['redshift_BEAGLE'][i] > 3.5 and data['redshift_BEAGLE'][i] < 5.5:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='orange')
            elif data['redshift_BEAGLE'][i] > 5.5:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='black')
            else:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='red')
        plt.ylim(-1, 3)
        plt.show()    
    
    for name in ['zeta']:            
        plt.figure(figsize=(12, 12))
        plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
        for i in range(len(data['redshift_BEAGLE'][object_start:object_end])):
            if data['redshift_BEAGLE'][i] > 4.0 and data['redshift_BEAGLE'][i] < 5.0:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='green')
            elif data['redshift_BEAGLE'][i] > 3.5 and data['redshift_BEAGLE'][i] < 5.5:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='orange')
            elif data['redshift_BEAGLE'][i] > 5.5:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='black')
            else:
                plt.plot(np.where((chain['zeta'][chain_start:chain_end,i]<=4)|(chain['zeta'][chain_start:chain_end,i]>=5), float('nan'), chain[name][chain_start:chain_end,i]), color='red')
        plt.ylim(3.8, 5.2)
        plt.show()    

#%%
    # HISTOGRAM
    for name in ['zeta']:
        chain_start = 7200          
        plt.figure(figsize=(12, 12))
        plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
        plt.hist(data['redshift_BEAGLE'][object_start:object_end][(chain['zeta'][chain_start,object_start:object_end]>=4)&(chain['zeta'][chain_start,object_start:object_end]<=5)], bins=20)
        plt.show()    



















#%%
# =============================================================================
# new investigation for objects entering z4-5 bin when sampling from z1.25-6.0
# =============================================================================

data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z1p25-6p0.p','r'))
chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_z4-5.p','r'))
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_ssfr_alpha.p','r'))




#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_z3-4.p','r'))

#data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z3p0-6p0.p','r'))
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z3p0-6p0_4x5000_z4-5.p','r'))

#%%

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors


xi = chain['xi'].flatten()
eta = chain['eta'].flatten()
zeta = chain['zeta'].flatten()

idx = np.random.choice(len(xi), 3000000, replace=False)

xi = xi[idx]
eta = eta[idx]
zeta = zeta[idx]

z = np.full(np.shape(chain['zeta']), data['redshift_BEAGLE']).flatten()[idx]


field_AD = np.full(np.shape(chain['zeta']), data['field_AD']).flatten()[idx]
id_AD = np.full(np.shape(chain['zeta']), data['id_AD']).flatten()[idx]
id_BEAGLE = np.full(np.shape(chain['zeta']), data['id_BEAGLE']).flatten()[idx]

idx = np.array((zeta>4.0) & (zeta<5.0))
#idx = np.array((zeta>3.0) & (zeta<4.0))

field_AD = field_AD[idx]
id_AD = id_AD[idx]
id_BEAGLE = id_BEAGLE[idx]

plt.figure(figsize=(12, 12))
cmap = cm.get_cmap('tab10', 6)    # 11 discrete colors
norm = matplotlib.colors.BoundaryNorm([0,1.25,2,3,4,5,6], 6)
c = np.linspace(0,6)
plt.scatter(xi[idx], eta[idx], c=z[idx], cmap=cmap, norm=norm)
plt.colorbar(spacing="uniform")
#plt.colorbar(spacing="proportional")

x = np.array((6,12))
y = np.array((-4,4))

plt.xlim(x)
plt.ylim(y)

#z3
#plt.plot(x,0.7224322491319064*(x-9.7) + 1.098423703873932, color='k')
#plt.plot(x,0.5833882583655428*(x-9.7) + 1.0383493410505018, color='k', linestyle='dashed')

#z4
plt.plot(x,0.9252070312634476*(x-9.7) + 1.6797732943151367, color='r', label='sampling z4-5')
plt.plot(x,0.6180518580111605*(x-9.7) + 1.7351855940565848, color='r', linestyle='dashed', label='sampling z1-6')
plt.plot(x,0.7831723848638836*(x-9.7) + 1.7292793235350037, color='r', linestyle='dotted', label='sampling z3-6')

plt.title('z4-5')
plt.legend()

specific_idx = (data['field_AD']==4)&(data['id_AD']==1284)
idx_chain_z1to2 = (chain['zeta'][:,specific_idx][::15]>1)&(chain['zeta'][:,specific_idx][::15]<2)
idx_chain_z4to5 = (chain['zeta'][:,specific_idx][::15]>4)&(chain['zeta'][:,specific_idx][::15]<5)


plt.scatter(chain['xi'][:,specific_idx][::15][idx_chain_z4to5], chain['eta'][:,specific_idx][::15][idx_chain_z4to5], marker='x', color='k')
plt.scatter(chain['xi'][:,specific_idx][::15][idx_chain_z1to2], chain['eta'][:,specific_idx][::15][idx_chain_z1to2], marker='.', color='k')

plt.show()


#print(id_AD)

'''
                'id_GMM_3d':            astrodeep_pickle['id_GMM_3d'],
                'x_GMM_3d':             astrodeep_pickle['x_GMM_3d'] - mag_GMM,
                'y_GMM_3d':             astrodeep_pickle['y_GMM_3d'] - mag_GMM,
                'z_GMM_3d':             astrodeep_pickle['z_GMM_3d'],
                'xsig_GMM_3d':          astrodeep_pickle['xsig_GMM_3d'],
                'ysig_GMM_3d':          astrodeep_pickle['ysig_GMM_3d'],
                'zsig_GMM_3d':          astrodeep_pickle['zsig_GMM_3d'],
                'xycov_GMM_3d':         astrodeep_pickle['xycov_GMM_3d'],
                'xzcov_GMM_3d':         astrodeep_pickle['xzcov_GMM_3d'],
                'yzcov_GMM_3d':         astrodeep_pickle['yzcov_GMM_3d'],
                'amp_GMM_3d':           astrodeep_pickle['amp_GMM_3d'],
'''


print(len(id_AD))
field_AD = field_AD[z[idx]<2.0]
id_AD = id_AD[z[idx]<2.0]
id_BEAGLE = id_BEAGLE[z[idx]<2.0]
print(len(id_AD))

#print(field_AD)
#print(id_AD)
#%%

plt.scatter(xi[idx][z[idx]<2.0], eta[idx][z[idx]<2.0], c=z[idx][z[idx]<2.0], cmap=cmap, norm=norm)
plt.show()

#%%

from collections import Counter
c = []
for i in range(len(id_AD)):
    c.append(str(int(field_AD[i]))+'_'+str(int(id_AD[i]))+'_'+str(int(id_BEAGLE[i])))
    
count = Counter(c)
print(count.keys())
print(count.values())

print(np.array(count.keys())[np.argsort(count.values())])
print(np.array(count.values())[np.argsort(count.values())])

c_k = np.array(count.keys())[np.argsort(count.values())]
c_v = np.array(count.values())[np.argsort(count.values())]



plt.hist(field_AD)
plt.show()


split = np.char.split(c_k, sep='_')


#%%
# =============================================================================
# create fits file 
# =============================================================================
from astropy.table import Table
from astropy.io import fits
outputDict = {}

print(len(split))
print(c_v)


outputDict['field_AD'] = []
outputDict['id_AD'] = []
outputDict['id_BEAGLE'] = []
outputDict['count'] = []

for i in range(len(split)):
    outputDict['field_AD'].append(split[i][0])
    outputDict['id_AD'].append(split[i][1])
    outputDict['id_BEAGLE'].append(split[i][2])
    outputDict['count'].append(c_v[i])

outputTable = Table(outputDict)
outputTable.write("low_slope_issue.fits", overwrite=True)

print(outputDict['id_BEAGLE'])



#%%
# =============================================================================
# redshifts
# =============================================================================
for g in range(len(split)):
    
    print(split[g][0], split[g][1])
    
    gidx = (data['field_AD']==float(split[g][0]))&(data['id_AD']==float(split[g][1]))
    
#    print(data['field_AD'][gidx])
#    print(data['id_AD'][gidx])
#    print(data['redshift_AD'][gidx])
#    print(data['redshift_BEAGLE'][gidx])
#    print(data['z_GMM_3d'][gidx][0])
#    print(data['zsig_GMM_3d'][gidx][0])
#    print(data['amp_GMM_3d'][gidx][0])
    
    x = np.linspace(0, 6, 1000)

    for i in range(3):
        mu = data['z_GMM_3d'][gidx][0][i]
        sigma = data['zsig_GMM_3d'][gidx][0][i]
        plt.plot(x, data['amp_GMM_3d'][gidx][0][i] * stats.norm.pdf(x, mu, sigma))
    

    plt.plot(x,     data['amp_GMM_3d'][gidx][0][0] * stats.norm.pdf(x, data['z_GMM_3d'][gidx][0][0], data['zsig_GMM_3d'][gidx][0][0]) +\
                    data['amp_GMM_3d'][gidx][0][1] * stats.norm.pdf(x, data['z_GMM_3d'][gidx][0][1], data['zsig_GMM_3d'][gidx][0][1]) +\
                    data['amp_GMM_3d'][gidx][0][2] * stats.norm.pdf(x, data['z_GMM_3d'][gidx][0][2], data['zsig_GMM_3d'][gidx][0][2]))
    
    plt.title(c_k[g].replace('_',' ') + ', count:' + str(c_v[g]))
    plt.show()   


#%%
# =============================================================================
# masses
# =============================================================================
for g in range(len(split)):
    
    print(split[g][0], split[g][1])
    
    gidx = (data['field_AD']==float(split[g][0]))&(data['id_AD']==float(split[g][1]))
    
#    print(data['field_AD'][gidx])
#    print(data['id_AD'][gidx])
#    print(data['redshift_AD'][gidx])
#    print(data['redshift_BEAGLE'][gidx])
#    print(data['z_GMM_3d'][gidx][0])
#    print(data['zsig_GMM_3d'][gidx][0])
#    print(data['amp_GMM_3d'][gidx][0])
    
    x = np.linspace(5, 13, 1000)

    for i in range(3):
        mu = data['x_GMM_3d'][gidx][0][i]
        sigma = data['xsig_GMM_3d'][gidx][0][i]
        plt.plot(x, data['amp_GMM_3d'][gidx][0][i] * stats.norm.pdf(x, mu, sigma))
    

    plt.plot(x,     data['amp_GMM_3d'][gidx][0][0] * stats.norm.pdf(x, data['x_GMM_3d'][gidx][0][0], data['xsig_GMM_3d'][gidx][0][0]) +\
                    data['amp_GMM_3d'][gidx][0][1] * stats.norm.pdf(x, data['x_GMM_3d'][gidx][0][1], data['xsig_GMM_3d'][gidx][0][1]) +\
                    data['amp_GMM_3d'][gidx][0][2] * stats.norm.pdf(x, data['x_GMM_3d'][gidx][0][2], data['xsig_GMM_3d'][gidx][0][2]))
    
    plt.title(c_k[g].replace('_',' ') + ', count:' + str(c_v[g]))
    plt.show()   

#%%
# =============================================================================
# sfrs
# =============================================================================
for g in range(len(split)):
    
    print(split[g][0], split[g][1])
    
    gidx = (data['field_AD']==float(split[g][0]))&(data['id_AD']==float(split[g][1]))
    
#    print(data['field_AD'][gidx])
#    print(data['id_AD'][gidx])
#    print(data['redshift_AD'][gidx])
#    print(data['redshift_BEAGLE'][gidx])
#    print(data['z_GMM_3d'][gidx][0])
#    print(data['zsig_GMM_3d'][gidx][0])
#    print(data['amp_GMM_3d'][gidx][0])
    
    x = np.linspace(-4, 4, 1000)

    for i in range(3):
        mu = data['y_GMM_3d'][gidx][0][i]
        sigma = data['ysig_GMM_3d'][gidx][0][i]
        plt.plot(x, data['amp_GMM_3d'][gidx][0][i] * stats.norm.pdf(x, mu, sigma))
    

    plt.plot(x,     data['amp_GMM_3d'][gidx][0][0] * stats.norm.pdf(x, data['y_GMM_3d'][gidx][0][0], data['ysig_GMM_3d'][gidx][0][0]) +\
                    data['amp_GMM_3d'][gidx][0][1] * stats.norm.pdf(x, data['y_GMM_3d'][gidx][0][1], data['ysig_GMM_3d'][gidx][0][1]) +\
                    data['amp_GMM_3d'][gidx][0][2] * stats.norm.pdf(x, data['y_GMM_3d'][gidx][0][2], data['ysig_GMM_3d'][gidx][0][2]))
    
    plt.title(c_k[g].replace('_',' ') + ', count:' + str(c_v[g]))
    plt.show()   


#%%
    

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_z4-5_no_burn.p','r'))


#%%
print(len(chain['zeta'][:,0]))


print(chain['zeta'][:,0])

specific_idx = (data['field_AD']==4)&(data['id_AD']==2831)


plt.plot(chain['xi'][:,specific_idx])
plt.show()

plt.plot(chain['eta'][:,specific_idx])
plt.show()

plt.plot(chain['zeta'][:,specific_idx])
plt.show()


print(len(specific_idx))
print(len(chain['pi']))
#%%
'''
g = 0

#for g in range(20):
for g in range(len(id_AD)):
    gidx = (data['field_AD']==field_AD[g])&(data['id_AD']==id_AD[g])
    
#    print(data['field_AD'][gidx])
    print(data['id_AD'][gidx])
#    print(data['redshift_AD'][gidx])
#    print(data['redshift_BEAGLE'][gidx])
#    print(data['z_GMM_3d'][gidx][0])
#    print(data['zsig_GMM_3d'][gidx][0])
#    print(data['amp_GMM_3d'][gidx][0])
    
    x = np.linspace(0, 6, 1000)

#    for i in range(3):
#        mu = data['z_GMM_3d'][gidx][0][i]
#        sigma = data['zsig_GMM_3d'][gidx][0][i]
#        plt.plot(x, data['amp_GMM_3d'][gidx][0][i] * stats.norm.pdf(x, mu, sigma))
#    plt.show()
    
    
        
#    plt.plot(x,     data['amp_GMM_3d'][gidx][0][0] * stats.norm.pdf(x, data['z_GMM_3d'][gidx][0][0], data['zsig_GMM_3d'][gidx][0][0]) +\
#                    data['amp_GMM_3d'][gidx][0][1] * stats.norm.pdf(x, data['z_GMM_3d'][gidx][0][1], data['zsig_GMM_3d'][gidx][0][1]) +\
#                    data['amp_GMM_3d'][gidx][0][2] * stats.norm.pdf(x, data['z_GMM_3d'][gidx][0][2], data['zsig_GMM_3d'][gidx][0][2]))
#    plt.show()
    '''
#%%
# =============================================================================
# taking sample from GMM and plotting on my redshift bin MS fits
# =============================================================================


#print(np.random.choice(3, 100))

filename = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/vis5_low_slope_issue.fits'

d = fits.open(filename)[1].data
#print(d.info())
#print(d[1].columns)


#print(d['id_AD'])
#


x_plot = np.array((6, 12))
y_plot = np.array((-4, 4))

#16 is 4 1284
#25 is 4 2831

for i in range(len(d['field_AD'])):
#for i in [25]:
    print(d['field_AD'][i], d['id_AD'][i])
    for j in range(1000):
        Gs = np.random.choice(3, p=d['amp_GMM_3d'][i]) 
        mean = [d['x_GMM_3d'][i][Gs], d['y_GMM_3d'][i][Gs], d['z_GMM_3d'][i][Gs]]
        cov = [[(d['xsig_GMM_3d'][i][Gs])**2,d['xycov_GMM_3d'][i][Gs],d['xzcov_GMM_3d'][i][Gs]],\
               [d['xycov_GMM_3d'][i][Gs],(d['ysig_GMM_3d'][i][Gs])**2,d['yzcov_GMM_3d'][i][Gs]],\
               [d['xzcov_GMM_3d'][i][Gs],d['yzcov_GMM_3d'][i][Gs],(d['zsig_GMM_3d'][i][Gs])**2]]
        
        if j==0:
            xi = np.random.multivariate_normal(mean, cov)[0]
            eta = np.random.multivariate_normal(mean, cov)[1]
            zeta = np.random.multivariate_normal(mean, cov)[2]
            
        else:
            xi = np.append(xi, np.random.multivariate_normal(mean, cov)[0])
            eta = np.append(eta, np.random.multivariate_normal(mean, cov)[1])
            zeta = np.append(zeta, np.random.multivariate_normal(mean, cov)[2])
    
    plt.figure(figsize=(10,10))
    plt.title('{} {}'.format(int(d['field_AD'][i]), int(d['id_AD'][i])))
    plt.xlim(x_plot)
    plt.ylim(y_plot)

    plt.scatter(xi, eta, c=zeta, s=5)

    alpha_a = [1.885, 2.068, 2.503, 9.545, 93.583, 0.13]
    alpha_b = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    beta_a = [0.86, 0.886, 0.722	, 0.925, 1.538, 0.768]
    beta_b = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    redshift = [1.625, 2.5, 3.5, 4.5, 5.5, 3.125]

#    alpha_a = [1.885, 2.068, 2.503, 9.545, 93.583, 0.13]
#    alpha_b = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#    beta_a = [0.86, 0.886, 0.722	, 0.925, 1.538, 0.768]
#    beta_b = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#    redshift = [1.625, 2.5, 3.5, 4.5, 5.5, 3.125]
    
#    for k in range(len(alpha_a)):
    for k in range(5):
        ssfr = np.log10(alpha_a[k]*(1.0+redshift[k])**alpha_b[k]) - 9.0
        alpha = ssfr + 9.7
        beta = beta_a[k] + redshift[k]*beta_b[k]
        plt.plot(x_plot, beta*(x_plot-9.7) + alpha, label=str(k+1))
    
    plt.legend()
    plt.colorbar()
    plt.show()
    
#    print(d['amp_GMM_3d'][i])
    print(d['mag_AD'][i])

#    # =============================================================================
#    # create fits file 
#    # =============================================================================
#    outputDict = {}
#    outputDict['xi']                = xi
#    outputDict['eta']               = eta
#    outputDict['zeta']               = zeta
#    outputDict['zeta']               = zeta
#    outputTable = Table(outputDict)
#    outputTable.write('{}_{}_GMM_samples.fits'.format(int(d['field_AD'][i]), int(d['id_AD'][i])), overwrite=True)








