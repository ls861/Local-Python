import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from scipy import integrate

subfolder = 'danger_delayed_uniform_logtau'
params = ['0', '1', '2', '3', '4', '5', '6', '7']
revision = '002'

for param in params:
    folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/{}/fit_{}/'.format(subfolder, param, revision)

    fileList = os.listdir(folder)
    for f, file in enumerate(fileList):
        if '.fits.gz' in file:
            print(param, file, f, len(fileList))
            data = fits.open(folder+file)

            tau = np.power(10,np.float64(data['POSTERIOR PDF'].data['tau']))
            age = np.power(10,np.float64(data['POSTERIOR PDF'].data['max_stellar_age']))
            norm_denom = (-tau*np.exp(-age/tau)*(tau+age))  +np.power(tau,2)

            log10_age = np.log10(age)
            log10_mass = np.float64(data['POSTERIOR PDF'].data['mass']) # mTot always needed for calculation
            log10_norm = log10_mass - np.log10(norm_denom)

            exp_term = -age/tau
            exp_idx = (exp_term <= -100.0)
            exp_term[exp_idx] = -100.0

            temp_sfr = log10_norm + log10_age + np.log10(np.exp(exp_term))
            log_instant_sfr = np.where(temp_sfr < -30.0, -30.0, temp_sfr)

            ID = file.replace('_BEAGLE.fits.gz','')

            np.save(folder+'{}_log_instant_sfr.npy'.format(ID), log_instant_sfr)

#            print(min(log_instant_sfr), max(log_instant_sfr))
#            plt.hist(log_instant_sfr)
#            plt.show()
