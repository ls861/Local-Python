import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import pickle


'''

pickle.dump({'ID':np.array(ID_arr), 'log_instant_sfr_mean':np.array(mean_arr), 'log_instant_sfr_median.npy':np.array(median_arr), 'log_instant_sfr_median_68_00.npy':np.array(interval_arr)}, open(output+'{}_instant_sfr_medians.p'.format(field[0]),'w'))

pickle.dump({'mTot':log_mTot, 'mStar':log_mStar, 'sfr':log_sfr, 'sfr_instant':log_sfr_instant}, open(output+'{}_{}_BEAGLE_samples.p'.format(field[0], ID),'w'))

pickle.dump({'mass':hp_x, 'sfr':hp_y}, open(output+'{}_{}_GMM_samples{}.p'.format(field[0], ID, mass_sfr_option),'w'))

'''


def get1DInterval(param_values, probability, levels):

    """
    Compute several quantities from a 1D probability density function

    Parameters
    ----------
    param_values : numpy array
        Contains the values of the parameter.

    probability : numpy array
        Contains the probability associated with each value of the
        parameter, hence must have same dimension as 'param_values'.

    levels : numpy array or list containing float
        Contains the (percentage) levels used to compute the credible
        regions, e.g. levels=[68.,95.] will compute 68 % and 95 % (central)
        credible regions

    Returns
    -------
    mean : float
        Mean of the parameter, computed as
        sum(probability*param_values)/sum(probability)
    median : float
        Median of the parameter, computed from the cumulative integral of
        the PDF
    interval : list of float
        2-dimensional list containing the lower and upper value of the
        parameter corresponding to the different `levels`

    """

    sort_ = np.argsort(param_values)

    # ******************************************************************
    # Here you must simply use `cumsum`, and not `cumtrapz` as in
    # beagle_utils.prepare_violin_plot, since the output of MultiNest are a set
    # of weights (which sum up to 1) associated to each set of parameters (the
    # `p_j` of equation 9 of Feroz+2009), and not a probability density (as the
    # MultiNest README would suggest).
    # ******************************************************************
    cumul_pdf = np.cumsum(probability[sort_])
    cumul_pdf /= cumul_pdf[len(cumul_pdf)-1]

    # Get the interpolant of the cumulative probability
    f_interp = interp1d(cumul_pdf, param_values[sort_])

    # You shoud integrate rather than summing here
    mean = np.sum(probability * param_values) / np.sum(probability)

    median = f_interp(0.5)

    interval = list()
    for lev in levels:

        low, high = f_interp([0.5*(1.-lev/100.), 1.-0.5*(1.-lev/100.)])
        interval.append([low,high])

    return mean, median, interval


local = True

if local == True:
    fields = ['1A2744P']
else:
    fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P', '4M0717C', '5M0717P', '6M1149C', '7M1149P']

mass_sfr_options = ['_mTot', '_mTot_delayed', '_mStar', '_mStar_delayed']
samples = 100

if local == True:
    output = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/scripts/sfr_instant_and_m_star/instant_sfr_medians_and_sample_outputs/'
else:
    output = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/kelly/instant_sfr_medians_and_sample_outputs/'

# =============================================================================
# median instant sfr
# =============================================================================

for field in fields:

    if local == True:
        folder = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/fields/1/'
    else:
        folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/fit_001/'.format(field)

    fileList = os.listdir(folder)

    ID_arr = []
    mean_arr = []
    median_arr = []
    interval_arr = []
    
    for file in fileList:
        if '_log_instant_sfr.npy' in file:
        
            ID = file.replace('_log_instant_sfr.npy','')
            log_instant_sfr = np.load(folder+file)

            data = fits.open(folder+ID+'_BEAGLE.fits.gz')
            probs_prop = np.array(data['POSTERIOR PDF'].data['probability'], np.float64)
            probs_prop = probs_prop/probs_prop.sum().astype(np.float64)
            data.close()

            mean, median, interval = get1DInterval(log_instant_sfr, probs_prop, [68.00])

            ID_arr.append(ID)
            mean_arr.append(mean)
            median_arr.append(median)
            interval_arr.append(interval)
            
    ID_arr = np.array(ID_arr, int)  
    idx_sort = np.argsort(ID_arr)
    ID_arr = ID_arr[idx_sort]
    
    mean_arr = np.array(mean_arr)[idx_sort]
    median_arr = np.array(median_arr)[idx_sort]
    interval_arr = np.array(interval_arr)[idx_sort]

    pickle.dump({'ID':ID_arr, 'log_instant_sfr_mean':mean_arr, 'log_instant_sfr_median':median_arr, 'log_instant_sfr_median_68_00':interval_arr}, open(output+'{}_instant_sfr_medians.p'.format(field[0]),'w'))



# =============================================================================
# BEAGLE samples for heatplots
# =============================================================================

'''
need BEAGLE output per object - mTot, mStar, sfr, instant_sfr
'''

for field in fields:

    if local == True:
        folder = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/fields/1/'
    else:
        folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/fit_001/'.format(field)

    fileList = os.listdir(folder)
    for file in fileList:
        if '_log_instant_sfr.npy' in file:

            ID = file.replace('_log_instant_sfr.npy','')
            data = fits.open(folder+ID+'_BEAGLE.fits.gz')

            log_mTot_arr = np.log10(np.array(data['GALAXY PROPERTIES'].data['M_tot'], np.float64))
            log_mStar_arr = np.log10(np.array(data['GALAXY PROPERTIES'].data['M_star'], np.float64))
            temp_sfr_arr = np.array(data['STAR FORMATION'].data['SFR'], np.float64)
            log_sfr_arr = np.log10(np.where(temp_sfr_arr==0, 1e-30, temp_sfr_arr))
            log_instant_sfr_arr = np.load(folder+file)

            probs_prop = np.array(data['POSTERIOR PDF'].data['probability'], np.float64)
            probs_prop = probs_prop/probs_prop.sum().astype(np.float64)

            data.close()

            idx = np.random.choice(len(probs_prop), size=samples, p=probs_prop)

            log_mTot = log_mTot_arr[idx]
            log_mStar = log_mStar_arr[idx]
            log_sfr = log_sfr_arr[idx]
            log_sfr_instant = log_instant_sfr_arr[idx]

            pickle.dump({'mTot':log_mTot, 'mStar':log_mStar, 'sfr':log_sfr, 'sfr_instant':log_sfr_instant}, open(output+'{}_{}_BEAGLE_samples.p'.format(field[0], ID),'w'))

            '''
            print(ID)
            plt.hist2d(log_mTot, log_sfr, bins=50)
            plt.show()
            plt.hist2d(log_mStar, log_sfr_instant, bins=50)
            plt.show()
            '''

# =============================================================================
# GMM samples for heatplots
# =============================================================================

'''
need GMM fit per object - mTot, mStar, sfr, instant_sfr
'''

for field in fields:

    if local == True:
        folder = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/fields/1/'
    else:
        folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/fit_001/'.format(field)

    for mass_sfr_option in mass_sfr_options:

        if local == True:
            GMM_folder = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/linmix_inputs/linmix_inputs_GMM_2d_{}{}/1p0_'.format(field[0], mass_sfr_option)
        else:
            GMM_folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/kelly/linmix_inputs_GMM_2d_{}{}/1p0_'.format(field[0], mass_sfr_option)

        IDs = np.load(GMM_folder+'id.npy')
        GMMx = np.load(GMM_folder+'GMMx.npy')
        GMMy = np.load(GMM_folder+'GMMy.npy')
        GMMxsig = np.load(GMM_folder+'GMMxsig.npy')
        GMMysig = np.load(GMM_folder+'GMMysig.npy')
        GMMxycov = np.load(GMM_folder+'GMMxycov.npy')
        pi_err = np.load(GMM_folder+'pi_err.npy')

        for obj, ID in enumerate(IDs):

            hp_x = []
            hp_y = []

            draws = np.random.choice([0, 1, 2], samples, p=pi_err[obj]/sum(pi_err[obj]))
            for draw in draws:
                mean = (GMMx[obj,draw], GMMy[obj,draw])
                cov = ((GMMxsig[obj,draw],GMMxycov[obj,draw]),(GMMxycov[obj,draw],GMMysig[obj,draw]))
                hpx, hpy = np.random.multivariate_normal(mean, cov)
                hp_x.append(hpx)
                hp_y.append(hpy)

            pickle.dump({'mass':hp_x, 'sfr':hp_y}, open(output+'{}_{}_GMM_samples{}.p'.format(field[0], ID, mass_sfr_option),'w'))


#data = pickle.load(open(args.folder+"sklearn_GMM_"+fileStr+"/"+file,"r"))
