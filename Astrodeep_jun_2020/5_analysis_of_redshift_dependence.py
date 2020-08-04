#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 13:00:38 2020

@author: lester
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde

fsize = 5
size = 8
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# removed mass > 10.5, SFR < -2 and SFR > 2.1, first 4 fields 
cut_mass_high = 10.5
cut_sfr_low = -2.0
cut_sfr_high = 2.1

# options:
# 7p0 2p0 10p0 9p5 001
# 8p5 2p0 10p0 9p5 001
# above with mass from 7p0 to 9p5

massLim = '8p0'
zLim = '2p0'
wLim = '10p0' # either side
chi2Lim = '9p5'
comment = '001'

# =============================================================================
# SCRIPT OPTIONS
# =============================================================================

option_plot_MS                  = 0
option_plot_MS_per_z            = 0
option_plot_mass_dependence     = 0
option_z_and_mass               = 0

option_alpha_beta               = 0 # SCATTER - shows random cuts on MS
option_sigma                    = 0 # histograms plus fit
option_sig0_and_k               = 0 # MS plot of fits 
option_sig0_and_k_scatter       = 0 # SCATTER + 2D histogram

option_plot_MS_fits             = 0
option_santini                  = 0
option_santini_ff               = 1
option_speagle                  = 0
option_speagle_paper            = 0

option_alpha_beta_sigma_z       = 0
option_alpha_beta_sigma_t       = 0


# plotting
MS_mass_min                     = 6.0
MS_mass_max                     = 11.0
MS_sfr_min                      = -3.0
MS_sfr_max                      = 4.0

# redshift bins
#z_med = np.linspace(1.5, 9.5, 9)
z_med = np.linspace(1.5, 5.5, 5) # TO MAKE MORE SIMILAR TO SANTINI REDSHIFTS
z_gap = (z_med[1] - z_med[0]) / 2

#sig0 and k
xi_max = 10.0
xi_min = 8.5

# =============================================================================
# GET DATA
# =============================================================================

sbf = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jun_2020/linmix_npy_files_z{}_w{}_chi{}_{}/{}_'.format(zLim, wLim, chi2Lim, comment, massLim)

pi_err      = np.load(sbf+'pi_err.npy')         # 3x probability of each posterior gaussian
GMMx        = np.load(sbf+'GMMx.npy')           # 3x posterior means per mass
GMMy        = np.load(sbf+'GMMy.npy')           # 3x posterior means per sfr
GMMxsig     = np.load(sbf+'GMMxsig.npy')        # 3x posterior sigmas per mass
GMMysig     = np.load(sbf+'GMMysig.npy')        # 3x posterior sigmas per sfr
GMMxycov    = np.load(sbf+'GMMxycov.npy')       # 3x posterior covar per mass-sfr pair

nK          = np.load(sbf+'nK.npy')             # 3 #gaussians modelling xi
nGauss      = np.load(sbf+'nGauss.npy')         # 3 #gaussians modelling BEAGLE posterior

#    nChains     = np.load(sbf+'nChains.npy')        # 2
#    minIter     = np.load(sbf+'minIter.npy')        # 3000
#    maxIter     = np.load(sbf+'maxIter.npy')        # 3000

GMMz        = np.load(sbf+'GMMz.npy') # posterior median
GMMchi2     = np.load(sbf+'GMMchi2.npy') # calculated elsewhere by me
GMMmass     = np.load(sbf+'GMMmass.npy') # posterior median
GMMsfr      = np.load(sbf+'GMMsfr.npy') # posterior median
        
# =============================================================================
# PLOT MAIN SEQUENCE
# =============================================================================

if option_plot_MS == 1:
    plt.scatter(GMMmass, GMMsfr)       
    plt.show()

# =============================================================================
# PLOT MAIN SEQUENCE per redshift bin
# =============================================================================

if option_plot_MS_per_z == 1:
    print(z_med, z_gap)
    for z in z_med:
        idx = (abs(GMMz - z) < z_gap)
        plt.scatter(GMMmass[idx], GMMsfr[idx], label='{}  - {}'.format(z-z_gap, z+z_gap))  
    #plt.xlim(6, 11)
    #plt.ylim(-5, 5)
    plt.legend(title='redshift')
    plt.show()

#%%
# =============================================================================
# PLOT adding in mass dependence somehow   
# =============================================================================

m_med = np.linspace(5.75, 10.75, 11)
m_gap = (m_med[1] - m_med[0]) / 2

if option_plot_mass_dependence == 1:
    print(m_med, m_gap)
    for m in m_med:
        idx = abs(GMMmass - m) < m_gap
        plt.scatter(GMMmass[idx], GMMsfr[idx], label='{}  - {}'.format(m-m_gap, m+m_gap))  
    #plt.xlim(6, 11)
    #plt.ylim(-5, 5)
    plt.legend(title='mass')
    plt.show()

#%%
# =============================================================================
# PLOT combining redshift and mass     
# =============================================================================

if option_z_and_mass == 1:
    for z in z_med:
        for m in m_med:
            z_idx = abs(GMMz - z) < z_gap
            m_idx = abs(GMMmass - m) < m_gap
            idx = z_idx & m_idx
        
            plt.scatter(GMMmass[idx], GMMsfr[idx], label='{}  - {}'.format(m-m_gap, m+m_gap))  
            plt.title('redshift {}  - {}'.format(z-z_gap, z+z_gap))
    #        plt.xlim(6, 11)
    #        plt.ylim(-5, 5)
        plt.legend(title='mass')
        plt.show()     

#%%
# =============================================================================
# calculating alpha and beta for given redshifts (and sigma)
# =============================================================================
           
alpha = []
beta = []
sigma = []

for z in z_med:
    idx = abs(GMMz - z) < z_gap
    
    # np.where is to sort the -inf SFR values at low redshift... (and replace -inf with -50)
    lin = np.polyfit(GMMmass[idx], np.where(GMMsfr[idx]<-50, -50, GMMsfr[idx]), 1)
        
    alpha.append(lin[1])
    beta.append(lin[0])

    if option_alpha_beta == 1:
        x = np.array((MS_mass_min, MS_mass_max))
        y = np.array((MS_sfr_min, MS_sfr_max))
        plt.xlim(x)
        plt.ylim(y)
        plt.plot(x, x*lin[0] + lin[1])
       
        # Calculate the point density
        v = np.vstack([GMMmass[idx],GMMsfr[idx]])
        density = gaussian_kde(v)(v)
        
        # Sort the points by density, so that the densest points are plotted last
        c_idx = density.argsort()
        c_x, c_y, c_z = GMMmass[idx][c_idx], GMMsfr[idx][c_idx], density[c_idx]

        plt.scatter(c_x, c_y, c=c_z, s=10, label='{}  - {}'.format(z-z_gap, z+z_gap))  
        
        # cuts        
        plt.plot(x,(cut_sfr_low, cut_sfr_low),color='r')
        plt.plot(x,(cut_sfr_high, cut_sfr_high),color='r')
        plt.plot((cut_mass_high, cut_mass_high),y,color='r')
        
        plt.legend(title='redshift')
        plt.title(massLim)
        plt.show()
    
    # this subtracts relation from each point, then fits normal to residual SFRs to calculate average sigma
    # np.where is to sort the -inf SFR values at low redshift... (and replace -inf with -50)
    data = np.where(GMMsfr[idx]<-50, -50, GMMsfr[idx]) - (GMMmass[idx]*lin[0] + lin[1])
    mean,std=norm.fit(data)
    sigma.append(std)
    
    if option_sigma == 1:
    # plotting residuals (SFR - relation), and the fitted normal
        testx = np.linspace(min(data), max(data), 1000)
        plt.hist(data, density=True)
        plt.plot(testx, norm.pdf(testx, mean, std))
        plt.show()

#%%
# =============================================================================
# plotting linear fits all in one go
# =============================================================================

if option_plot_MS_fits == 1:
    for i in range(len(alpha)):
        x = np.array((MS_mass_min, MS_mass_max))
        y = np.array((MS_sfr_min, MS_sfr_max))
        plt.xlim(x)
        plt.ylim(y)
        plt.plot(x, x*beta[i] + alpha[i], label='z={:.1f}'.format(z_med[i]))
    plt.legend()
    plt.show()

#%%
# =============================================================================
# sig0 and k
# =============================================================================

sig0=[]
k = []

for i, z in enumerate(z_med):
    idx = abs(GMMz - z) < z_gap
    
    # np.where is to sort the -inf SFR values at low redshift... (and replace -inf with -50)
    lin = np.polyfit(GMMmass[idx], np.where(GMMsfr[idx]<-50, -50, GMMsfr[idx]), 1)
      
    # this subtracts relation from each point, then fits normal to residual SFRs to calculate average sigma
    # np.where is to sort the -inf SFR values at low redshift... (and replace -inf with -50)
    data = np.where(GMMsfr[idx]<-50, -50, GMMsfr[idx]) - (GMMmass[idx]*lin[0] + lin[1])    
    lin2 = np.polyfit(GMMmass[idx], abs(data), 1)
    
    sig0.append(lin2[0]*xi_max + lin2[1])
    k.append((lin2[0]*xi_min + lin2[1]) / sig0[i])
    
    x_lin2 = np.array((MS_mass_min, MS_mass_max))
    y_lin2 = lin2[0]*x_lin2 + lin2[1]

    if option_sig0_and_k == 1:
        plt.plot(x_lin2, y_lin2, label=z)

if option_sig0_and_k == 1:
    plt.title('fitted scatter on residuals per redshift')
    plt.xlim(MS_mass_min, MS_mass_max)
    plt.legend()
    plt.show()

#%%
# =============================================================================
# scatter plots for sig0 and k 
# =============================================================================

for i, z in enumerate(z_med):
    idx = abs(GMMz - z) < z_gap
    
    # np.where is to sort the -inf SFR values at low redshift... (and replace -inf with -50)
    lin = np.polyfit(GMMmass[idx], np.where(GMMsfr[idx]<-50, -50, GMMsfr[idx]), 1)
      
    # this subtracts relation from each point, then fits normal to residual SFRs to calculate average sigma
    # np.where is to sort the -inf SFR values at low redshift... (and replace -inf with -50)
    data = np.where(GMMsfr[idx]<-50, -50, GMMsfr[idx]) - (GMMmass[idx]*lin[0] + lin[1])    
    lin2 = np.polyfit(GMMmass[idx], abs(data), 1)
    
    x_lin2 = np.array((MS_mass_min, MS_mass_max))
    y_lin2 = lin2[0]*x_lin2 + lin2[1]
    
    if option_sig0_and_k_scatter == 1:
        
        # Calculate the point density
        # v = np.vstack([GMMmass[idx],abs(data)])
        # density = gaussian_kde(v)(v)
        
        plt.scatter(GMMmass[idx], abs(data), marker='x')
        plt.plot(x_lin2, y_lin2, label=z)
        plt.title('fitted scatter on residuals per redshift')
        plt.xlim(MS_mass_min, MS_mass_max)
        plt.ylim(bottom=0)
        plt.legend()
        plt.show()
        
        plt.title(massLim)
        plt.hist2d(GMMmass[idx], abs(data), bins=20)
        plt.show()

#%%
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

#for i in range(len(A_san)):
#    plt.plot(x, A_san[i]*(x-9.7) + B_san[i], linestyle='-', label='San17 z={}'.format(z_med_san[i]))
#plt.xlim(6, 11)
#plt.ylim(-3, 4)
#plt.legend()
#plt.show()

if option_santini == 1:
    # SAME AS ABOVE, JUST CHECKING I REARRANGED EQN OK
    x = np.array((MS_mass_min, MS_mass_max))
    y = np.array((MS_sfr_min, MS_sfr_max))
    for i in range(len(alpha_san)):
        plt.plot(x, x*beta_san[i] + alpha_san[i], linestyle='-', label='San17 z={}'.format(z_med_san[i]))
    plt.xlim(x)
    plt.ylim(y)
    plt.legend()
    plt.show()
    
# =============================================================================
# functional form fitting to Santini
# =============================================================================

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
    
if option_santini_ff == 1:

    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
    
#    xdata = np.linspace(0, 4, 50)
#    y = func(xdata, 2.5, 1.3, 0.5)
#    np.random.seed(1729)
#    y_noise = 0.2 * np.random.normal(size=xdata.size)
#    ydata = y + y_noise
#    plt.plot(xdata, ydata, 'b-', label='data')
#    
#    popt, pcov = curve_fit(func, xdata, ydata)
#    
#    plt.plot(xdata, func(xdata, *popt), 'r-',
#             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
#        
#    popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [3., 1., 0.5]))
#    
#    plt.plot(xdata, func(xdata, *popt), 'g--',
#             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
#        
#    plt.xlabel('x')
#    plt.ylabel('y')
#    plt.legend()
#    plt.show()
    
    def func2(x, a, b, c):
        return a + b*x + c*x**2
    
    def func3(x, a, b, c):
        return a + b*(1/x) + c*(1/x)**2

    def func4(x, a, b, c, d):
        return a + b*(1/(x-d)) + c*(1/(x-d))**2

    def func5(x, a, b, c, d):
        return a + b*x + c*x**2 + d*x**3   
    
    def func6(x, a, b, c, d, e):
        return a + b*x + c*x**2 + d*x**3 + e*x**4
    
    def func7(x, a, b, c):
        return a + b*(1/x) + c*(1/x)**2
    
    
    
    f = func3
 
    xdata = z_med_san
    ydata = alpha_san
    xdata_arr = np.linspace(min(xdata), max(xdata), 100)   
    plt.scatter(xdata, ydata, label='data')
    
    popt, pcov = curve_fit(f, xdata, ydata)
    plt.plot(xdata_arr, f(xdata_arr, *popt), 'r-')
         
    print(popt)
    plt.title('Santini alpha')
    plt.xlabel('redshift')
    plt.ylabel('alpha')
    plt.legend()
    plt.show()

    xdata = z_med
    ydata = alpha
    xdata_arr = np.linspace(min(xdata), max(xdata), 100)   
    plt.scatter(xdata, ydata, label='data')
    
    popt, pcov = curve_fit(f, xdata, ydata)
    plt.plot(xdata_arr, f(xdata_arr, *popt), 'r-')
          
    plt.title('Fitted alpha')
    plt.xlabel('redshift')
    plt.ylabel('alpha')
    plt.legend()
    plt.show()
             
    xdata = z_med_san
    ydata = beta_san
    xdata_arr = np.linspace(min(xdata), max(xdata), 100)   
    plt.scatter(xdata, ydata, label='data')
    
    popt, pcov = curve_fit(f, xdata, ydata)
    plt.plot(xdata_arr, f(xdata_arr, *popt), 'r-')
          
    print(popt)
    plt.title('Santini beta')
    plt.xlabel('redshift')
    plt.ylabel('beta')
    plt.legend()
    plt.show()
    


#%%
# =============================================================================
# Speagle
# =============================================================================

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z_med_spe = np.linspace(min(z_med), max(z_med), 100)

t_med_spe = cosmo.age(z_med_spe).value
u = cosmo.age(z_med_spe).unit

alpha_spe = 0.11*t_med_spe - 6.51
beta_spe = 0.84 - 0.026*t_med_spe

alpha_err_spe = ((0.03*t_med_spe)**2 + 0.24**2) ** 0.5
beta_err_spe = ((0.003*t_med_spe)**2 + 0.02**2) ** 0.5

if option_speagle == 1:

    z_med_spe_subset = z_med_san
    t_med_spe_subset = cosmo.age(z_med_spe_subset).value
    alpha_spe_subset = 0.11*t_med_spe_subset - 6.51
    beta_spe_subset = 0.84 - 0.026*t_med_spe_subset

    x = np.array((MS_mass_min, MS_mass_max))
    y = np.array((MS_sfr_min, MS_sfr_max))
    
    for i in range(len(alpha_spe_subset)):
        plt.plot(x, x*beta_spe_subset[i] + alpha_spe_subset[i], linestyle='-', label='Spe14 z={}'.format(z_med_spe_subset[i]))
        
    plt.xlim(x)
    plt.ylim(y)
    plt.legend()
    plt.show()

if option_speagle_paper == 1:
    # for plotting
    z_med_spe_subset = np.array((0, 0.25, 0.5, 1, 2, 4))
    t_med_spe_subset = cosmo.age(z_med_spe_subset).value
    alpha_spe_subset = 0.11*t_med_spe_subset - 6.51
    beta_spe_subset = 0.84 - 0.026*t_med_spe_subset
    
    for i in range(len(alpha_spe_subset)):
        plt.plot(x, x*beta_spe_subset[i] + alpha_spe_subset[i], linestyle='-', label='Spe14 z={}'.format(z_med_spe_subset[i]))
    plt.xlim(9.7, 11)
    plt.ylim(-0.5, 3)
    plt.legend()
    plt.show()

#%%
# =============================================================================
# PLOT alpha, beta, sigma vs redshift 
# =============================================================================

if option_alpha_beta_sigma_z == 1:
    
    plt.title('{} - alpha vs redshift'.format(massLim))
    plt.scatter(z_med, alpha, marker='x', label='This Work')
    plt.plot(z_med_san, alpha_san, color=colors[1], label='Santini')
    plt.plot(z_med_san, alpha_san+alpha_err_san, linestyle=':', color=colors[1])
    plt.plot(z_med_san, alpha_san-alpha_err_san, linestyle=':', color=colors[1])
    plt.plot(z_med_spe, alpha_spe, color=colors[2], label='Speagle')
    plt.plot(z_med_spe, alpha_spe+alpha_err_spe, linestyle=':', color=colors[2])
    plt.plot(z_med_spe, alpha_spe-alpha_err_spe, linestyle=':', color=colors[2])
    plt.legend()
    plt.show()
    
    plt.title('{} - beta vs redshift'.format(massLim))
    plt.scatter(z_med, beta, marker='x', label='This Work')
    plt.plot(z_med_san, beta_san, color=colors[1], label='Santini')
    plt.plot(z_med_san, beta_san+beta_err_san, linestyle=':', color=colors[1])
    plt.plot(z_med_san, beta_san-beta_err_san, linestyle=':', color=colors[1])
    plt.plot(z_med_spe, beta_spe, color=colors[2], label='Speagle')
    plt.plot(z_med_spe, beta_spe+beta_err_spe, linestyle=':', color=colors[2])
    plt.plot(z_med_spe, beta_spe-beta_err_spe, linestyle=':', color=colors[2])
    plt.legend()
    plt.show()
    
    plt.title('{} - sigma vs redshift'.format(massLim))
    plt.scatter(z_med, sigma, marker='x', label='This Work')
    plt.legend()
    plt.show()   
    
    plt.title('{} - sig0 vs redshift'.format(massLim))
    plt.scatter(z_med, sig0, marker='x', label='This Work')
    plt.legend()
    plt.show()     
    
    plt.title('{} - k vs redshift'.format(massLim))
    plt.scatter(z_med, k, marker='x', label='This Work')
    plt.legend()
    plt.show()          

#%%
# =============================================================================
# PLOT alpha, beta, sigma vs age of universe
# =============================================================================

if option_alpha_beta_sigma_t  == 1:
    
    t_med = cosmo.age(z_med).value
    t_med_san = cosmo.age(z_med_san).value
    
    plt.title('{} - alpha vs age of universe'.format(massLim))
    plt.scatter(t_med, alpha, marker='x', label='This Work')
    plt.plot(t_med_san, alpha_san, color=colors[1], label='Santini')
    plt.plot(t_med_san, alpha_san+alpha_err_san, linestyle=':', color=colors[1])
    plt.plot(t_med_san, alpha_san-alpha_err_san, linestyle=':', color=colors[1])
    plt.plot(t_med_spe, alpha_spe, color=colors[2], label='Speagle')
    plt.plot(t_med_spe, alpha_spe+alpha_err_spe, linestyle=':', color=colors[2])
    plt.plot(t_med_spe, alpha_spe-alpha_err_spe, linestyle=':', color=colors[2])
    plt.legend()
    plt.show()
    
    plt.title('{} - beta vs age of universe'.format(massLim))
    plt.scatter(t_med, beta, marker='x', label='This Work')
    plt.plot(t_med_san, beta_san, color=colors[1], label='Santini')
    plt.plot(t_med_san, beta_san+beta_err_san, linestyle=':', color=colors[1])
    plt.plot(t_med_san, beta_san-beta_err_san, linestyle=':', color=colors[1])
    plt.plot(t_med_spe, beta_spe, color=colors[2], label='Speagle')
    plt.plot(t_med_spe, beta_spe+beta_err_spe, linestyle=':', color=colors[2])
    plt.plot(t_med_spe, beta_spe-beta_err_spe, linestyle=':', color=colors[2])
    plt.legend()
    plt.show()
    
    plt.title('{} - sigma vs age of universe'.format(massLim))
    plt.scatter(t_med, sigma, marker='x', label='This Work')
    plt.legend()
    plt.show()        
    
    plt.title('{} - sig0 vs age of universe'.format(massLim))
    plt.scatter(t_med, sig0, marker='x', label='This Work')
    plt.legend()
    plt.show()     
    
    plt.title('{} - k vs age of universe'.format(massLim))
    #plt.ylim(0, 5)
    plt.scatter(t_med, k, marker='x', label='This Work')
    plt.legend()
    plt.show()     
    


#%%








