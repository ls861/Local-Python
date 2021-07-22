import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from scipy import integrate

local = True 

if local == True:
    fields = ['1A2744P']
else:
    fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P', '4M0717C', '5M0717P', '6M1149C', '7M1149P']


for field in fields:
    
    if local == True:
        folder = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/fields/1/'
    else:
        folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/fit_001/'.format(field)


    fileList = os.listdir(folder)
    for f, file in enumerate(fileList):
        if '.fits.gz' in file:
            #check if the results for this individual object already stored
            print(field, file, f, len(fileList))
            data = fits.open(folder+file)

#            mTot = np.float64(data['GALAXY PROPERTIES'].data['M_tot'])
#            mStar = np.float64(data['GALAXY PROPERTIES'].data['M_star'])
            
            tau = np.power(10,np.float64(data['POSTERIOR PDF'].data['tau']))
            age = np.power(10,np.float64(data['POSTERIOR PDF'].data['max_stellar_age']))
            norm_denom = (-tau*np.exp(-age/tau)*(tau+age))  +np.power(tau,2)
 
#            tau = np.power(10,data['POSTERIOR PDF'].data['tau'])
#            age = np.power(10,data['POSTERIOR PDF'].data['max_stellar_age'])
#            norm_denom3 = (-tau*np.exp(-age/tau)*(tau+age))  +np.power(tau,2)
            
#            norm_denom1 = []
#            for i in range(len(tau)):
#                sfr_function = lambda t: t*np.exp(-t/tau[i])
#
#                integral = integrate.quad(sfr_function, 0, age[i])
#                norm_denom1.append(integral[0])
#
#                if integral[1]/integral[0] > 1e-5:
#                    print('integration error', i, file, integral[0], integral[1]/integral[0])
#            norm_denom1 = np.array(norm_denom1)

            log10_age = np.log10(age)
            log10_mass = np.float64(data['POSTERIOR PDF'].data['mass']) # mTot always needed for calculation
            log10_norm = log10_mass - np.log10(norm_denom)

            exp_term = -age/tau
            exp_idx = (exp_term <= -100.0)
            exp_term[exp_idx] = -100.0

            temp_sfr = log10_norm + log10_age + np.log10(np.exp(exp_term))
            log_instant_sfr = np.where(temp_sfr < -30.0, -30.0, temp_sfr)

#            for i in range(len(log10_norm)):
#                if exp_term[i] == -100.0:
#                    print(log10_norm[i], log10_age[i], np.log10(np.exp(exp_term))[i], temp_sfr[i] )

            ID = file.replace('_BEAGLE.fits.gz','')
            
            np.save(folder+'{}_log_instant_sfr.npy'.format(ID), log_instant_sfr)

#            print(min(log_instant_sfr), max(log_instant_sfr))
#            plt.hist(log_instant_sfr)
#            plt.show()


