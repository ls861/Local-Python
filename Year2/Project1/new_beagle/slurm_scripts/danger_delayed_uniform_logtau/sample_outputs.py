import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import pickle

subfolder = 'danger_delayed_uniform_logtau'
params = ['0', '1', '2', '3', '4', '5', '6', '7']
revision = '002'

os.system("mkdir ./sample_outputs")

# mass_sfr_options = ['_mTot', '_mTot_delayed', '_mStar', '_mStar_delayed']
mass_sfr_option = '_mStar_delayed'
samples = 100

# =============================================================================
# BEAGLE samples for heatplots
# =============================================================================
'''
need BEAGLE output per object - mTot, mStar, sfr, instant_sfr
'''
for param in params:

    ID_arr = []
    log_mStar_arr = []
#    log_sfr_arr = []
    log_sfr_instant_arr = []

    folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/{}/fit_{}/'.format(subfolder, param, revision)

    fileList = os.listdir(folder)
    for file in fileList:
        if '_BEAGLE.fits.gz' in file:

            ID = file.replace('_BEAGLE.fits.gz','')
            data = fits.open(folder+ID+'_BEAGLE.fits.gz')

            log_mStar_gz = np.log10(np.array(data['GALAXY PROPERTIES'].data['M_star'], np.float64))
#            temp_sfr_gz = np.array(data['STAR FORMATION'].data['SFR'], np.float64)
#            log_sfr_gz = np.log10(np.where(temp_sfr_gz==0, 1e-30, temp_sfr_gz))
            log_instant_sfr_gz = np.load(folder+ID+'_log_instant_sfr.npy')

            probs_prop = np.array(data['POSTERIOR PDF'].data['probability'], np.float64)
            probs_prop = probs_prop/probs_prop.sum().astype(np.float64)

            data.close()

            idx = np.random.choice(len(probs_prop), size=samples, p=probs_prop)

            log_mStar = log_mStar_gz[idx]
#            log_sfr = log_sfr_gz[idx]
            log_sfr_instant = log_instant_sfr_gz[idx]

            ID_arr.append(ID)
            log_mStar_arr.append(log_mStar)
#            log_sfr_arr.append(log_sfr)
            log_sfr_instant_arr.append(log_sfr_instant)

    idx_sort = np.argsort(np.asarray(np.array(ID_arr), float))

    ID_arr = np.array(ID_arr)[idx_sort]
    log_mStar_arr = np.array(log_mStar_arr)[idx_sort]
#    log_sfr_arr = np.array(log_sfr_arr)[idx_sort]
    log_sfr_instant_arr = np.array(log_sfr_instant_arr)[idx_sort]

    pickle.dump({'ID':ID_arr, 'log_mStar_arr':log_mStar_arr, 'log_sfr_instant_arr':log_sfr_instant_arr}, open('./sample_outputs/{}{}_BEAGLE_samples.p'.format(param, mass_sfr_option),'w'))

# =============================================================================
# GMM samples for heatplots
# =============================================================================
'''
need GMM fit per object - mTot, mStar, sfr, instant_sfr
'''
for param in params:

    ID_arr = []
    hp_x_arr = []
    hp_y_arr = []

    GMM_2d = pickle.load(open('./kelly_2d_GMM_inputs/{}{}_GMM_2d_fits.p'.format(param, mass_sfr_option),'r'))

    IDs = GMM_2d['id']
    GMMx = GMM_2d['x']
    GMMy = GMM_2d['y']
    GMMxsig = GMM_2d['xsig']
    GMMysig = GMM_2d['ysig']
    GMMxycov = GMM_2d['xycov']
    pi_err = GMM_2d['amp']

    for obj, ID in enumerate(IDs):

        hp_x = []
        hp_y = []

        draws = np.random.choice([0, 1, 2], samples, p=pi_err[obj]/sum(pi_err[obj]))
        for draw in draws:

            mean = np.array([GMMx[obj,draw],GMMy[obj,draw]])
            cov = np.array([[np.square(GMMxsig[obj,draw]), GMMxycov[obj,draw]],[GMMxycov[obj,draw], np.square(GMMysig[obj,draw])]])

            hpx, hpy = np.random.multivariate_normal(mean, cov)
            hp_x.append(hpx)
            hp_y.append(hpy)

        ID_arr.append(ID)
        hp_x_arr.append(hp_x)
        hp_y_arr.append(hp_y)

    idx_sort = np.argsort(np.asarray(np.array(ID_arr), float))

    ID_arr = np.array(ID_arr)[idx_sort]
    hp_x_arr = np.array(hp_x_arr)[idx_sort]
    hp_y_arr = np.array(hp_y_arr)[idx_sort]

    pickle.dump({'ID':ID_arr, 'hp_x_arr':hp_x_arr, 'hp_y_arr':hp_y_arr}, open('./sample_outputs/{}{}_GMM_2d_samples.p'.format(param, mass_sfr_option),'w'))

for param in params:

    ID_arr = []
    hp_x_arr = []
    hp_y_arr = []
    hp_z_arr = []

    GMM_3d = pickle.load(open('./kelly_3d_GMM_inputs/{}{}_GMM_3d_fits.p'.format(param, mass_sfr_option),'r'))

    IDs = GMM_3d['id']
    GMMx = GMM_3d['x']
    GMMy = GMM_3d['y']
    GMMz = GMM_3d['z']
    GMMxsig = GMM_3d['xsig']
    GMMysig = GMM_3d['ysig']
    GMMzsig = GMM_3d['zsig']
    GMMxycov = GMM_3d['xycov']
    GMMxzcov = GMM_3d['xzcov']
    GMMyzcov = GMM_3d['yzcov']
    pi_err = GMM_3d['amp']

    for obj, ID in enumerate(IDs):

        hp_x = []
        hp_y = []
        hp_z = []
        draws = np.random.choice([0, 1, 2], samples, p=pi_err[obj]/sum(pi_err[obj]))
        for draw in draws:

            mean = np.array([GMMx[obj,draw],GMMy[obj,draw],GMMz[obj,draw]])
            cov = np.array([[np.square(GMMxsig[obj,draw]), GMMxycov[obj,draw], GMMxzcov[obj,draw]],[GMMxycov[obj,draw], np.square(GMMysig[obj,draw]), GMMyzcov[obj,draw]],[GMMxzcov[obj,draw], GMMyzcov[obj,draw], np.square(GMMzsig[obj,draw])]])

            hpx, hpy, hpz = np.random.multivariate_normal(mean, cov)
            hp_x.append(hpx)
            hp_y.append(hpy)
            hp_z.append(hpz)

        ID_arr.append(ID)
        hp_x_arr.append(hp_x)
        hp_y_arr.append(hp_y)
        hp_z_arr.append(hp_z)

    idx_sort = np.argsort(np.asarray(np.array(ID_arr), float))

    ID_arr = np.array(ID_arr)[idx_sort]
    hp_x_arr = np.array(hp_x_arr)[idx_sort]
    hp_y_arr = np.array(hp_y_arr)[idx_sort]
    hp_z_arr = np.array(hp_z_arr)[idx_sort]

    pickle.dump({'ID':ID_arr, 'hp_x_arr':hp_x_arr, 'hp_y_arr':hp_y_arr, 'hp_z_arr':hp_z_arr}, open('./sample_outputs/{}{}_GMM_3d_samples.p'.format(param, mass_sfr_option),'w'))
