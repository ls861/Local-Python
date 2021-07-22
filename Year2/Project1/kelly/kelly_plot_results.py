#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 13:11:43 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle


#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_new_001.p','r'))



#options = 'z1p3_test2' # 263 items 

options = 'z1p3'
#options = 'z2p0'
#options = 'z3p0'
#options = 'z4p0'
#options = 'z5p0'
'''

#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_new_{}.p'.format(options),'r'))
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_5_7_data_{}_k0.p'.format(options),'r')) # 4x1000, burn at 100
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_5_7_data_{}_k.p'.format(options),'r')) # 4x1000, burn at 100

#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_5_7_data_{}_k0_002.p'.format(options),'r')) # 4x2000, burn at 200
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_5_7_data_{}_k_002.p'.format(options),'r')) # 4x2000, burn at 200

#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k0_002.p'.format(options),'r')) # 4x2000, burn at 200
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k_002.p'.format(options),'r')) # 4x2000, burn at 200

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_z5p0_k0_convergence_test_001.p'.format(options),'r')) # 4x12000, burn at 1000

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_z1p3_k0_convergence_test_002.p'.format(options),'r')) # 4x1000, burn at 100

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_z1p3_k0_convergence_test_003.p'.format(options),'r')) # 4x1000, burn at 100, extended to 1600

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_z1p3_k0_convergence_test_004.p'.format(options),'r')) # 4x1000, burn at 100

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_z1p3_k_convergence_test_005.p'.format(options),'r')) # 4x80, burn at 8

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_z5p0_k_convergence_test_006.p'.format(options),'r')) # 4x80, burn at 8

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_z1p3_k_convergence_test_007.p'.format(options),'r')) # 4x2000, burn at 100

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_redshift_test_00q.p'.format(options),'r')) # 4 chains

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_redshift_test_001.p'.format(options),'r')) # 4 chains
#Extended Iteration:  10000
#Rhat values:
#alpha: 1.1882903236145372 1.1974383385165992 1.0946523512435733
#beta: 1.0462844808637575 1.214802992921653 1.2321855079852468
#sig_0: 1.0212075623045882
#k: 1.025943666588585
#ximean: 1.0007778381202068
#xivar: 1.000005744144727

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_redshift_test_002.p'.format(options),'r')) # 4 chains
#Extended Iteration:  10000
#Rhat values:
#alpha: 1.0172913707649116 1.0358361266104281 1.0687698349350978
#beta: 1.0169513253994102 1.0677362864924183 1.0689584599941822
#sig_0: 1.0098163528851192
#k: 1.0127237168177754
#ximean: 1.0002563656116887
#xivar: 1.000037014889825



chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_z5p0_k_convergence_test_008.p'.format(options),'r'))

chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_redshift_test_009.p'.format(options),'r'))


chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_redshift_test_010.p'.format(options),'r'))





chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_redshift_test_011.p'.format(options),'r'))

'''

#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_mock_004_fit_001.p'.format(options),'r')) # 4 chains
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_mock_004_fit_002.p'.format(options),'r')) # 4 chains
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_mock_005_fit_001.p'.format(options),'r')) # 4 chains
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_mock_005_fit_004.p'.format(options),'r')) # 4 chains
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_mock_006_fit_016.p'.format(options),'r')) # 4 chains


#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_redshift_test_031.p'.format(options),'r'))
#chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_bias_test_103.p'.format(options),'r'))





# hogg 3d onwards
chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_mock_hogg_redshift_109.p','r')) # 4 chains


print(type(chain))
#print(chain['alpha'])
print(len(chain))
print(chain.dtype.names)

nChains=4
minIter=len(chain)/nChains
burn=minIter/2
#burn=0



# =============================================================================
# combine chains and remove burn in
# =============================================================================

chain_arr = []
for i in range(nChains):
    start = minIter*i+burn
    finish = minIter*(i+1)
    chain_arr.append(chain[start:finish])
chain = np.concatenate(chain_arr)

# =============================================================================
# plot chains 
# =============================================================================

# ('alpha', 'alpha_a', 'alpha_b', 'alpha_c', 'beta', 'beta_a', 'beta_b', 'beta_c', 'sig0', 'k', 'xi_min', 'xi_max', 'xi', 'eta', 'zeta', 'pi', 'mu', 'tausqr', 'mu0', 'usqr', 'wsqr', 'ximean', 'xisig')

#names = ['alpha_a', 'alpha_b', 'beta_a', 'sig0', 'k', 'pi', 'mu', 'tausqr', 'mu0', 'usqr', 'wsqr', 'ximean', 'xisig']

names = ['alpha_a', 'alpha_b', 'beta_a', 'sig0', 'k', 'pi', 'mu', 'tausqr', 'mu0', 'usqr', 'wsqr', 'ximean', 'xisig', 'pbad', 'outlier_mean', 'outlier_sigma']

#for i in range(4):
#for i in range(len(chain.dtype.names)):
#    plt.title(chain.dtype.names[i].replace('_',''))
#    plt.plot(chain[chain.dtype.names[i]])
#    plt.show()

for name in names: 
    plt.title(name.replace('_',''))
    plt.plot(chain[name])
    plt.show()


# =============================================================================
# single object chains
# =============================================================================
'''
plt.title('xi[0]')
for i in range(1):
    plt.plot(chain['xi'][:,i])
plt.show()
print(len(chain['xi'][0]))

plt.title('eta[0]')
for i in range(1):
    plt.plot(chain['eta'][:,i])
plt.show()
print(len(chain['eta'][0]))

plt.title('zeta[0]')
for i in range(1):
    plt.plot(chain['zeta'][:,i])
plt.show()
print(len(chain['zeta'][0]))
'''
# =============================================================================
# histograms
# =============================================================================
    
#for i in range(8):
#    plt.title(chain.dtype.names[i].replace('_',''))
#    plt.hist(chain[chain.dtype.names[i]])
#    plt.show()
#    print('{:.3f} {:.3f} {:.3f}'.format(np.percentile(chain[chain.dtype.names[i]], 16), np.median(chain[chain.dtype.names[i]]), np.percentile(chain[chain.dtype.names[i]], 84)))

names = ['alpha_a', 'alpha_b', 'beta_a', 'sig0', 'k']
names = ['alpha_a', 'alpha_b', 'beta_a', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']

for name in names: 
    plt.title(name.replace('_',''))
    plt.hist(chain[name])
    plt.show()
    print('{:.3f} {:.3f} {:.3f}'.format(np.percentile(chain[name], 16), np.median(chain[name]), np.percentile(chain[name], 84)))


#plt.scatter(chain['pi'][:,0], chain['mu'][:,0], alpha=0.1, s=3)
#plt.scatter(chain['pi'][:,1], chain['mu'][:,1], alpha=0.1, s=3)
#plt.scatter(chain['pi'][:,2], chain['mu'][:,2], alpha=0.1, s=3)
#plt.show()
#
#
#print(chain['pi'])




# =============================================================================
# bias test histograms
# =============================================================================

tests = ['101', '102', '103', '104', '105', '106', '107', '108', '109', '110']
names = ['alpha_a', 'alpha_b', 'beta_a', 'sig0', 'k']
trues = [-9.0, 0.5, 1.0, 0.3, 1.0]

# HOGG 3d
tests = ['100', '101', '102', '103', '104', '105', '106', '107', '108', '109']
tests = ['100', '102', '103', '104', '105', '107']
names = ['alpha_a', 'alpha_b', 'beta_a', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']
trues = [-9.0, 0.5, 1.0, 0.3, 1.0, 0.3, 0.0, 2.0]

for j, name in enumerate(names): 
    plt.title(name.replace('_',''))
    
    for test in tests: 
#        chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_bias_test_{}.p'.format(test),'r'))
        chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_mock_hogg_redshift_{}.p'.format(test),'r'))
        nChains=4
        minIter=len(chain)/nChains
        burn=minIter/2
        #burn=0
        chain_arr = []
        for i in range(nChains):
            start = minIter*i+burn
            finish = minIter*(i+1)
            chain_arr.append(chain[start:finish])
        chain = np.concatenate(chain_arr)
        plt.hist(chain[name], histtype=u'step')
    plt.show()
    
    count= 1
    plt.plot((trues[j], trues[j]),(0,11),linestyle=':',color='g')
    for test in tests:
#        chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_bias_test_{}.p'.format(test),'r'))
        chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_mock_hogg_redshift_{}.p'.format(test),'r'))
        nChains=4
        minIter=len(chain)/nChains
        burn=minIter/2
        #burn=0
        chain_arr = []
        for i in range(nChains):
            start = minIter*i+burn
            finish = minIter*(i+1)
            chain_arr.append(chain[start:finish])
        chain = np.concatenate(chain_arr)
        
        plt.plot((np.percentile(chain[name], 16), np.percentile(chain[name], 84)), (count, count), color='k')
        plt.scatter(np.median(chain[name]), count, color='r', marker='x')
        count += 1
    plt.yticks([])        
    plt.show()    
    


    
#    print('{:.3f} {:.3f} {:.3f}'.format(np.percentile(chain[name], 16), np.median(chain[name]), np.percentile(chain[name], 84)))


























