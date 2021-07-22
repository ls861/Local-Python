import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import pickle

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

    interval = list()
    for lev in levels:
        lower_lev = 0.5*(1.-lev/100.)
        upper_lev = 1.-0.5*(1.-lev/100.)
        #check if either limit is lower or higher than current samples
        if lower_lev < cumul_pdf[0]:
            interval.append([-99,-99])
            median = -99
        elif upper_lev > cumul_pdf[-1]:
            interval.append([-99,-99])
            median = -99
        else:
            low, high = f_interp([lower_lev, upper_lev])
            interval.append([low,high])
            median = f_interp(0.5)
    return mean, median, interval


subfolder = 'danger_delayed_uniform_logtau'
params = ['0', '1', '2', '3', '4', '5', '6', '7']
revision = '002'

os.system("mkdir ./instant_sfr_medians")

# =============================================================================
# median instant sfr
# =============================================================================

for param in params:

    folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/{}/fit_{}/'.format(subfolder, param, revision)

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

    idx_sort = np.argsort(np.asarray(np.array(ID_arr), float))

    ID_arr = np.array(ID_arr)[idx_sort]
    mean_arr = np.array(mean_arr)[idx_sort]
    median_arr = np.array(median_arr)[idx_sort]
    interval_arr = np.array(interval_arr)[idx_sort]

    pickle.dump({'ID':ID_arr, 'log_instant_sfr_mean':mean_arr, 'log_instant_sfr_median.npy':median_arr, 'log_instant_sfr_median_68_00.npy':interval_arr}, open('./instant_sfr_medians/{}_instant_sfr_medians.p'.format(param),'w'))
