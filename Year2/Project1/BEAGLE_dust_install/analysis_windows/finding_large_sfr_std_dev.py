# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 15:11:26 2021

@author: LSand
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

# =============================================================================
# EDIT THESE AS REQUIRED
# =============================================================================
scenarioA = '29'
fields = ['clusters']
z_bins = ['z1p25-2p0']
z_lower = 1.25
z_upper = 2.0

for field in fields:
    for z_bin in z_bins:

        with open('/Users/LSand/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis_windows/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, field, z_bin), 'rb') as f:
            data = pickle.load(f, encoding='latin1')

        s = data        
        z_med_hp = (z_lower+z_upper)/2.0
        z_med_hp_gap = (z_lower+z_upper)/2.0 - z_lower
        
        n_hp = 3000 # number of samples to take from GMM in total
        
        fig, ax = plt.subplots(1, 1, figsize=(15, 15))

        sfr_std = []
        
        print(len(s['id_AD']))
        
        for i in range(len(s['id_AD'])):
            
            x_temp = np.array([])
            y_temp = np.array([])
            z_temp = np.array([])
        
            for G in range(3):
                
                mean = np.array([s['x_GMM_3d'][i,G],s['y_GMM_3d'][i,G],s['z_GMM_3d'][i,G]])
                cov = np.array([[np.power(s['xsig_GMM_3d'][i,G],2), s['xycov_GMM_3d'][i,G], s['xzcov_GMM_3d'][i,G]],[s['xycov_GMM_3d'][i,G], np.power(s['ysig_GMM_3d'][i,G],2), s['yzcov_GMM_3d'][i,G]],[s['xzcov_GMM_3d'][i,G], s['yzcov_GMM_3d'][i,G], np.power(s['zsig_GMM_3d'][i,G],2)]])
        
                xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s['amp_GMM_3d'][i,G]))

                z_idx = (abs(xyz[:,2] - z_med_hp) < z_med_hp_gap)
        
                # x_temp = np.concatenate((x_temp,xyz[:,0][z_idx]))
                # y_temp = np.concatenate((y_temp,xyz[:,1][z_idx]))
                # z_temp = np.concatenate((z_temp,xyz[:,2][z_idx]))

                x_temp = np.concatenate((x_temp,xyz[:,0]))
                y_temp = np.concatenate((y_temp,xyz[:,1]))
                z_temp = np.concatenate((z_temp,xyz[:,2]))
                
            sfr_std.append(np.std(y_temp))
            
            # if i == 0:
            #     plt.hist(z_temp, bins=30)
            #     plt.show()
            
            if i == 0:
                if np.std(y_temp) > 2:
                    confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='blue', linewidth=1)
                    ax.scatter(x_temp, y_temp, s=0.5)
    
                else:
                    confidence_ellipse(x_temp, y_temp, ax, n_std=1, label=r'$1\sigma$', edgecolor='firebrick', linewidth=1)
                    ax.scatter(x_temp, y_temp, s=0.5)
                    pass
    
        ax.set_xlim(6, 11)
        ax.set_ylim(-3, 4)
        # ax.set_ylim(-20, 4)        
        plt.show()


        sfr_std = np.array(sfr_std)
        field_AD = s['field_AD'][sfr_std<2]
        id_AD = s['id_AD'][sfr_std<2]
        id_BEAGLE = s['id_BEAGLE'][sfr_std<2]
          

        print(len(field_AD))
