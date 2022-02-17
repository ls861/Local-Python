#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 21:35:05 2022

@author: lester
"""


import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# IMF corrections
# =============================================================================

# https://ned.ipac.caltech.edu/level5/March14/Madau/Madau3.html

s2c_M = np.log10(0.61) # -0.21 x 
k2c_M = np.log10(0.61/0.66) # -0.034 x 
 
s2c_SFR = np.log10(0.63) # -0.20 y 
k2c_SFR = np.log10(0.63/0.67) # -0.027 y 

# =============================================================================
# Reddy+12
# =============================================================================
# z_redd = np.array([2.3006, 3.0027]) # original data thief 
ssfr_redd = np.array([0.3663, 0.3472]) - 9.0 # original data thief from Santini paper
# ssfr_p_redd = np.array([0.7475, 0.7236]) - ssfr_redd - 9.0 # original data thief 
# ssfr_m_redd = ssfr_redd + 9.0 - np.array([-5.3603e-3, -0.0244]) # original data thief 

# from fig18 caption
z_redd = np.array([2.3, 3.0])
ssfr_p_redd = np.array([0.37, 0.37])
ssfr_m_redd = np.array([0.37, 0.37])
ssfr_redd += s2c_SFR - s2c_M

# =============================================================================
# de Barros+14
# =============================================================================
# z_barr = np.array([3.2294, 3.7121, 4.8237, 5.9062]) # original data thief 
z_barr = np.array([3.32, 3.88, 4.94, 6.00]) # from table 1
ssfr_barr = np.array([0.8046, 0.8046, 1.2382, 1.2335]) - 9.0 # original data thief 
ssfr_p_barr = np.array([1.3002, 1.3145, 1.3096, 1.2906]) - ssfr_barr - 9.0 # original data thief 
ssfr_m_barr = ssfr_barr + 9.0 - np.array([-0.0578, -0.0482, 0.0899, 0.0995]) # original data thief 
ssfr_barr += s2c_SFR - s2c_M

# =============================================================================
# Gonzalez+14
# =============================================================================
z_gonz = np.array([3.8072, 5.0212, 5.9135]) # original data thief 
ssfr_gonz = np.array([0.5426, 0.5331, 0.6808]) - 9.0 # original data thief 
ssfr_p_gonz = np.array([0.8285, 0.8237, 0.9762]) - ssfr_gonz - 9.0 # original data thief 
ssfr_m_gonz = ssfr_gonz + 9.0 - np.array([0.2376, 0.2328, 0.3806]) # original data thief 
ssfr_gonz += s2c_SFR - s2c_M

# =============================================================================
# Smit+14 FROM SMIT figure 6
# =============================================================================
z_smit = 6.8
ssfr_smit = 1.73
ssfr_p_smit = 2.02-ssfr_smit
ssfr_m_smit = ssfr_smit-1.06
ssfr_smit -= 9.0 
ssfr_smit += s2c_SFR - s2c_M

# =============================================================================
# Speagle+14 - errors by error propagation
# =============================================================================

# log SFR(M∗, t) = (0.84 ± 0.02 − 0.026 ± 0.003 × t ) logM∗−(6.51 ± 0.24 − 0.11 ± 0.03 × t ), where t is the age of the universe in Gyr.

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z_spea = np.linspace(0.5, 6.5, 100)
t_spea = cosmo.age(z_spea).value # Gyr

alpha_spea = -(6.51 - (0.11*t_spea) )
alpha_err_spea = ((0.24)**2+(0.03*t_spea)**2)**0.5

beta_spea = 0.84 - (0.026*t_spea)
beta_err_spea = ((0.02)**2+(0.003*t_spea)**2)**0.5

alphaN_spea = alpha_spea + ((9.7-k2c_M)*beta_spea) + k2c_SFR # santini normalised
alphaN_err_spea = ((alpha_err_spea)**2+(beta_err_spea*9.7)**2)**0.5

ssfr_spea = alphaN_spea - 9.7
ssfr_err_spea = alphaN_err_spea

# sigma_t is "true" scatter. σt (Δt ) = (0.20 ± 0.02) + (0.007 ± 0.006) Δt
# Δt = 0 for us, therefore σt = (0.20 ± 0.02)
# refer to eqn 27: σt ∼ 0.2dex

sig_spea = np.full(100, 0.20)
sig_err_spea = np.full(100, 0.02)

# =============================================================================
# Steinhardt+14
# =============================================================================
#log SFR(M yr−1) = α × (logM∗/M − 10.5) + β,
z_stei = np.array(((4.0 + 4.8)/2, (4.8 + 6.0)/2))

beta_stei = np.array((0.78, 0.78)) # α
beta_err_stei = np.array((0.02, 0.02)) # α

# alpha_stei = np.array((1.976, 2.110)) - (10.5*beta_stei) # error on 1.976 and 2.110 is 0.005 and 0.003
# alpha_err_stei = np.array((0.005, 0.003))
# alphaN_stei = alpha_stei + (9.7*beta_stei) # santini normalised
# alphaN_err_stei = ((alpha_err_stei)**2+(beta_err_stei*9.7)**2)**0.5

alphaN_stei = -0.8*beta_stei + np.array((1.976, 2.110)) # santini normalised
alphaN_err_stei = ((np.array((0.005, 0.003)))**2+(0.8*beta_err_stei)**2)**0.5

ssfr_stei = alphaN_stei - 9.7
ssfr_err_stei = alphaN_err_stei

# from a single sentence
sig_stei = np.array((0.24, 0.24))

# =============================================================================
# Whitaker+14 - CHABRIER IMF
# =============================================================================

# there is a full parameterisation, but I am using low mass regime from TABLE 3
# note the erratum does NOT account for this, so I assume it is fine
# mass < 10.2 relates to basically all of my points

# logSFR = a [logM - 10.2] + b = a*logM - 10.2a + b

z_whit = np.array([(0.5 + 1.0)/2, (1.0 + 1.5)/2, (1.5 + 2.0)/2, (2.0 + 2.5)/2])

a_low = np.array([0.94,0.99,1.04,0.91])
a_low_err = np.array([0.03,0.04,0.05,0.06])

b = np.array([1.11,1.31,1.49,1.62])
b_err = np.array([0.03,0.02,0.02,0.02])

# alpha_whit = b - 10.2*a_low
# alpha_err_whit = ((b_err)**2+(10.2*a_low_err)**2)**0.5

beta_whit = a_low
beta_err_whit = a_low_err

# alphaN_whit = alpha_whit + (9.7*beta_whit) # santini normalised
# alphaN_err_whit = ((alpha_err_whit)**2+(beta_err_whit*9.7)**2)**0.5

alphaN_whit = -0.5*a_low + b # santini normalised
alphaN_err_whit = ((b_err)**2+(0.5*a_low_err)**2)**0.5

ssfr_whit = alphaN_whit - 9.7
ssfr_err_whit = alphaN_err_whit

# This method allows us to reach far below the limits placed by the individual MIPS 24μm image depths, with the trade off that we are not able to measure the intrinsic scatter in the average log Ψ– log M relation in this study.

# =============================================================================
# Furlong+15 EAGLE SIM - not obvious, maybe just for SSFR?
# =============================================================================

# DATA FROM DATA THIEF, MIDDLE PANEL, figure 7

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z_furl_arr = np.linspace(0.0, 6.0, 100000)
t_furl_arr = cosmo.lookback_time(z_furl_arr).value # Gyr

t_furl = np.array((0.0788, 1.1823, 2.798, 4.0985, 5.3202, 5.9506, 6.6995, 7.803, 8.867, 10.4828, 10.9163, 12.0197, 12.4138, 12.6897)) # lookback time in Gyr
z_furl = np.zeros(len(t_furl))

for i in range(len(t_furl)):
    idx = np.abs(t_furl[i] - t_furl_arr).argmin()
    z_furl[i] = z_furl_arr[idx]
    # print(t_furl[i], t_furl_arr[idx])

ssfr_furl = np.array((-1.2778, -1.1993, -1.0725, -0.9398, -0.8191, -0.7528, -0.6623, -0.4935, -0.391, -0.1197, -0.0835, 0.1755, 0.3923, 0.597)) - 9.0

# =============================================================================
# Lee+15 TABLE 1 - CHABRIER IMF - ONLY PLOTTED SSFR WITHOUT ERROR 
# =============================================================================

# log SFR = s0 - log [1+ (10^M / 10^M0)**-g ]

# using median redshift, not centre of range

z_lee = np.array([0.36,0.55,0.70,0.85,0.99,1.19])
s0_lee = np.array([0.80,0.99,1.23,1.35,1.53,1.72])
s0_err_lee = np.array([0.019,0.015,0.016,0.014,0.017,0.024])
M0_lee = np.array([10.03,9.82,9.93,9.96,10.10,10.31])
M0_err_lee = np.array([0.042,0.031,0.031,0.025,0.029,0.043])
g_lee = np.array([0.92,1.13,1.11,1.28,1.26,1.07])
g_err_lee = np.array([0.017,0.033,0.025,0.034,0.032,0.028])

# low mass approx
beta_lee = g_lee
beta_err_lee = g_err_lee

alphaN_lee = s0_lee - np.log10(1+(((10**9.7)/(10**M0_lee))**-g_lee)) # errors seem complicated


alphaN_err_lee = np.zeros(len(alphaN_lee))

for i in range(len(s0_lee)):
    s0_lee_d = np.random.normal(s0_lee[i], s0_err_lee[i], int(1e5))
    M0_lee_d = np.random.normal(M0_lee[i], M0_err_lee[i], int(1e5))
    g_lee_d = np.random.normal(g_lee[i], g_err_lee[i], int(1e5))

    alphaN_err_lee[i] = np.std(s0_lee_d - np.log10(1+(((10**9.7)/(10**M0_lee_d))**-g_lee_d)))

ssfr_lee = alphaN_lee - 9.7
ssfr_err_lee = alphaN_err_lee
# print(ssfr_err_lee)

'''
# error on 10**M0_lee

a = 10.0
b = 1.0
A = M0_lee
sA = M0_err_lee
f = a**(b*A)
sf = abs(f) * abs(b*np.log(a)*sA)

# error on (10**9.7)/(10**M0_lee)

A = (10**9.7)
sA = 0.0
B = f
sB = sf
f = A/B
sf = abs(f) * (((sA/A)**2 + (sB/B)**2)**0.5)

# error on (((10**9.7)/(10**M0_lee))**-g_lee)

A = f
sA = sf
B = -g_lee
sB = g_err_lee
f = A**B
sf = abs(f) * (((B*sA/A)**2 + (sB*np.log(A))**2)**0.5)

# error on np.log10(1+(((10**9.7)/(10**M0_lee))**-g_lee))

a = 1.0
b = 1.0
A = 1.0 + f
sA = sf
f = a*np.log10(b*A)
sf = abs((a*sA) / (A*np.log(10)))

# error on s0_lee - np.log10(1+(((10**9.7)/(10**M0_lee))**-g_lee))

A = s0_lee
sA = s0_err_lee
B = f
sB = sf
f = A-B
sf = (sA**2 + sB**2)**0.5


alphaN_err_lee = sf

ssfr_lee = alphaN_lee - 9.7
ssfr_err_lee = alphaN_err_lee
print(ssfr_err_lee)
'''
# The standard deviation of the SFR in each stellar mass bin remains mostly constant at all masses and at all redshifts, with σ ∼ 0.36 dex in all bins. 

# sig_lee = np.full(len(z_lee), 0.36) # only plotting ssfr

# =============================================================================
# Salmon+15 - table 4 - SALPETER IMF, but explain how to convert to Chabrier (-0.25dex to mass and sfr)
# =============================================================================
#log(SFR/M yr−1) = a log(M/Msun) + b
z_salm = np.array((4.0, 5.0, 6.0))
'''
alpha_salm = np.array((-5.7, -4.4, -3.9))
alpha_err_salm = np.array((2.1, 2.6, 1.6))

beta_salm = np.array((0.7, 0.59, 0.54))
beta_err_salm = np.array((0.21, 0.26, 0.16))

alphaN_salm = alpha_salm + (9.7*beta_salm) # santini normalised
alphaN_err_salm = ((alpha_err_salm)**2+(beta_err_salm*9.7)**2)**0.5

ssfr_salm = alphaN_salm - 9.7
ssfr_err_salm = alphaN_err_salm
'''
# we use the median absolute deviation (MAD) to compute the equivalent standard deviation, σMAD, as the measure of scatter in given quantities (Beers et al. 1990), including the quoted scatter for redshift, stellar mass, and SFR. The σMAD is an analog for the 68% confidence, σ, if the error distribution were Gaussian and is therefore less sensitive to outliers 

# values also available for slope fixed at 1 in paper, but not used here
# errors not supplied
sig_salm = np.array([0.35, 0.41, 0.21])


### SLOPE a == 1
alpha_salm1 = np.array((-8.64, -8.49, -8.45))
alpha_err_salm1 = np.array((0.11, 0.14, 0.06))

beta_salm1 = np.array((1.0, 1.0, 1.0))

alphaN_salm1 = alpha_salm1 + ((9.7-s2c_M)*beta_salm1) + s2c_SFR # santini normalised
alphaN_err_salm1 = alpha_err_salm1

ssfr_salm1 = alphaN_salm1 - 9.7
ssfr_err_salm1 = alphaN_err_salm1

# =============================================================================
# Schreiber+15 - ignoring high mass - SALPETER IMF
# =============================================================================
#r ≡ log10(1 + z) and m ≡ log10(M∗/10**9 M):
#log10(SFRMS[M/yr]) = m − m0 + a0r − a1 max(0,m − m1 − a2 r)**2,
#with m0 = 0.5 ± 0.07, a0 = 1.5 ± 0.15, a1 = 0.3 ± 0.08, m1 = 0.36 ± 0.3 and a2 = 2.5 ± 0.6.

# errors available in paper
z_schr = np.linspace(0.5, 4.0, 100)
r_schr = np.log10(1+z_schr)

m0_schr = 0.5
a0_schr = 1.5

m0_err_schr = 0.07
a0_err_schr = 0.15

# m - m1 - a2r is usually < 0, except high mass, low redshift, IGNORED FOR NOW
# at least for m=9.7, for 0.5<z<4, this final term can be ignored
# a1 = 0.3
# m1 = 0.36
# a2 = 2.5
# print(np.log10((10**9.9)/(10**9.0)) - m1 - (a2*r_schr))
# print(z_schr)

# plt.plot(z_schr, 9.0+m1+(a2*r_schr))
# plt.xlabel('z')
# plt.ylabel('mass at which slope =/= 1')
# plt.show()

# alpha_schr = - (9.0 + m0_schr - (a0_schr*r_schr)) # at M=0
# alphaN_schr = alpha_schr + (9.7*beta_schr) # santini normalised
# alpha_err_schr = ((0.07)**2+(0.15**r_schr)**2)**0.5
# alphaN_err_schr = alpha_err_schr


beta_schr = np.linspace(1.0, 1.0, 100)
alphaN_schr = (9.7-s2c_M-9.0) - m0_schr + (a0_schr*r_schr) + s2c_SFR
ssfr_schr = alphaN_schr - 9.7

alphaN_err_schr = ((m0_err_schr)**2+(a0_err_schr*r_schr)**2)**0.5
ssfr_err_schr = alphaN_err_schr

# fig 12 mass 10.2, solid diamonds BY EYE

z_sig_schr = np.array([0.5, 1.0, 1.5])
sig_schr = np.array([0.28, 0.27, 0.24])

# =============================================================================
# Kurczynski+16 SALPETER IMF - z0.5 - 3 - TABLE 1, 
# =============================================================================
# also mentions purpose of rescaling which I haven't included here
# log SFR = a ´ log M* + b + N (0, sIS).
z_kurc = np.array(((0.5 + 1.0)/2, (1.0 + 1.5)/2, (1.5 + 2.0)/2, (2.0 + 2.5)/2, (2.5 + 3.0)/2))

alpha_kurc = np.array((-8.394, -7.474, -7.484, -7.513, -7.729))
alpha_err_kurc = np.array((0.011, 0.010, 0.011, 0.018, 0.015))
beta_kurc = np.array((0.919, 0.825, 0.867, 0.849, 0.899))
beta_err_kurc = np.array((0.017, 0.012, 0.013, 0.021, 0.017))
alphaN_kurc = alpha_kurc + ((9.7-s2c_M)*beta_kurc) + s2c_SFR # santini normalised
alphaN_err_kurc = ((alpha_err_kurc)**2+(beta_err_kurc*9.7)**2)**0.5

ssfr_kurc = alphaN_kurc - 9.7
ssfr_err_kurc = alphaN_err_kurc

# total scatter is available in table 2 (instead of intrinsic)
# scatter per mass bin is also available in table 2

#intrinsic scatter
sig_kurc = np.array([0.427,0.273,0.255,0.281,0.220])
sig_err_kurc = np.array([0.011,0.009,0.008,0.017,0.017])

# =============================================================================
# Marmol-Queralto+16
# =============================================================================
# santini data thief 
# z_marm0 = np.array([4.385])
# ssfr_marm0 = np.array([0.7284]) - 9.0
# ssfr_p_marm0 = np.array([1.0286]) - ssfr_marm0 - 9.0
# ssfr_m_marm0 = ssfr_marm0 + 9.0 - np.array([0.433])

z_marm = np.array((1.34, 2.25, 4.55, 4.34))
ssfr_marm = np.log10(np.array((1.7, 4.5, 5.5, 5.1)) * 1e-9)
ssfr_err_marm = np.array((0.2, 0.3, 1.1, 0.6)) * 0.434 / (np.array((1.7, 4.5, 5.5, 5.1)))

# =============================================================================
# Tomczak+16 - CHABRIER IMF
# =============================================================================

# https://www.desmos.com/calculator/2pm5e7auau

# bit more obvious now, but fair amount of effort to get slope and some decisions would need to be made
# taking their star forming parameterisation
# NO ERRORS GIVEN
# santini took specific points for a mass, I'm using STAR FORMING parameterisation at M=9.7

# param based upon Lee 15

z_tomc = np.linspace(0.5, 4.0, 100)

s_sf= 0.448 + 1.220*z_tomc - 0.174*z_tomc**2 
M_sf= 10**(9.458 + 0.865*z_tomc - 0.132*z_tomc**2) # NOT LOG
g_sf= 1.091 

alphaN_tomc = s_sf - np.log10(1+((10**9.7 / M_sf)**-g_sf)) # NOT PLOTTING, INCONSISTENT WITH SLOPE

# linear part (subbing a low mass, M=10^7, based on desmos)
# beta_tomc = g_sf / ((10**1.0 / M_sf)**g_sf + 1) # approaches g_sf as M->0
# beta_tomc = np.full(100, g_sf) # NOT PLOTTING 

ssfr_tomc = alphaN_tomc - 9.7

# can't find anything on scatter

# =============================================================================
# Santini+17 'True' values SALPETER IMF - delayed SFH, SFR from UV slope - errors don't make sense
# =============================================================================
# logSFR = A log(M / M_9.7) + B
# TABLE 2
z_sant = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_sant = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
A_err_sant = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
B_sant = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
B_err_sant = np.array((0.04, 0.03, 0.03, 0.05, 0.13))

# converting normalisation
beta_sant = A_sant
beta_err_sant = A_err_sant
# alphaN_san = alpha_san + (9.7*beta_san) # santini normalised
# alphaN_err_san = (alpha_err_san**2 + (9.7*beta_err_san)**2) ** 0.5 # doesn't make sense
alphaN_sant = A_sant * ((9.7-s2c_M) - 9.7) + B_sant + s2c_SFR
alphaN_err_sant = B_err_sant

# ssfr_san_dt = np.array([0.3139, 0.5188, 0.6712, 0.6665, 1.2953]) - 9.0 # from data thief
ssfr_sant = alphaN_sant - 9.7
ssfr_err_sant = alphaN_err_sant


# fig 4 for scatter BY EYE - not sure if based on true or observed MS - EMMA SAYS TRUE

# Santini scatter can be considered intrinsic convolved with evolution of MS in redshift bin
# mass bins < 8 - 8.8 - 9.6 - 10.3 >
# centres of mass bins 8.4, 9.2, 9.95
# z 1.3 2 3 4 5, centres: 1.65, 2.5, 3.5, 4.5
# blue stars in fig 4: intrinsic scatter 

# M_arr_sant = np.array([8.4, 9.2, 9.95])
z_sig_sant = np.array([1.65, 2.5, 3.5, 4.5])

# BY EYE, can't find numbers in paper
# sig8p4_sant = np.array([0.36,0.27,0.0,0.34])
sig9p2_sant = np.array([0.35,0.1,0.0,0.23])
# sig9p95_sant = np.array([0.0,0.0,0.33,0.04])

'''
for i in range(4):
    plt.plot((M_arr_san[0],M_arr_san[1],M_arr_san[2]), (sig8p4[i],sig9p2[i],sig9p95_san[i]))
    grad = (sig9p2[i]-sig9p95_san[i]) / (M_arr_san[1]-M_arr_san[2])
    c = sig9p95_san[i] - (grad*M_arr_san[2])
    plt.scatter(9.7, grad*9.7 + c)
    print(grad, c, grad*9.7 + c)

plt.show()
sig9p7 = np.array([0.1167,0.0333,0.2200,0.1033])
print(z_arr[:-1])
'''

# =============================================================================
# SANTINI Original values - obtained by eye - delayed SFH, SFR from UV slope - errors don't make sense
# =============================================================================
'''
# logSFR = alpha log(M / M_9.7) + beta
z_san0 = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san0 = np.array((1.05, 1.1, 0.9, 0.75, 0.55))
A_err_san0 = np.array((0.03, 0.03, 0.04, 0.05, 0.18))
B_san0 = np.array((1.0, 1.15, 1.25, 1.2, 1.5))
B_err_san0 = np.array((0.05, 0.03, 0.03, 0.06, 0.12))

# converting normalisation
beta_san0 = A_san0
beta_err_san0 = A_err_san0
alphaN_san0 = B_san0
alphaN_err_san0 = B_err_san0

# from data thief
ssfr_san0 = np.array([0.2662, 0.4282, 0.4949, 0.2424, 0.5615]) - 9.0
ssfr_p_san0 = np.array([0.5855, 0.5855, 0.7379, 0.4235, 1.2144]) - ssfr_san0 - 9.0
ssfr_m_san0 = ssfr_san0 + 9.0 - np.array([-0.0482, 0.2852, 0.2472, 0.0661, -0.0911])
'''
# =============================================================================
# Boogaard+18 - low mass 7 <M<10.5 for 0.11<z<0.91 - CHABRIER IMF
# =============================================================================

# z_boo = np.linspace(0.11, 0.91, 100)
# sig_boo = np.linspace(0.44, 0.44, 100)
# sig_p_boo = np.linspace(0.05, 0.05, 100)
# sig_m_boo = np.linspace(0.04, 0.04, 100)

# their MS eqn in abstract is wrong (I think) and other values in paper aren't clear


# last 2 rows table 1
#log SFR(M yr−1) = a × (logM∗ − 8.5) + b

z_boog = np.array(((0.1 + 0.5)/2, (0.5 + 1.0)/2))

beta_boog = np.array((0.86, 0.84)) # a
beta_err_boog = np.array((0.08, 0.06)) # a


alphaN_boog = (1.2)*beta_boog + np.array((-0.92, -0.73)) # santini normalised
alphaN_err_boog = ((np.array((0.07, 0.06)))**2+(1.2*beta_err_boog)**2)**0.5

ssfr_boog = alphaN_boog - 9.7
ssfr_err_boog = alphaN_err_boog


sig_boog = np.array((0.57, 0.46))
sig_err_boog = np.array((0.06, 0.05))

# =============================================================================
# Pearson+18 - CHABRIER IMF - TABLE 2
# =============================================================================
# S = a [logM - 10.5] + b

z_pear = np.array([0.37,0.66,0.95,1.24,1.59,2.02,2.59,3.23,4.34,5.18]) # medians
z_mid = np.array([(0.2+0.5)/2,(0.5+0.8)/2,(0.8+1.1)/2,(1.1+1.4)/2,(1.4+1.8)/2,(1.8+2.3)/2,(2.3+2.9)/2,(2.9+3.8)/2,(3.8+4.9)/2,(4.9+6.0)/2,])

a = np.array([0.43,0.5,0.46,0.48,0.51,0.74,0.83,0.70,0.93,1.00])
a_err = np.array([0.09,0.10,0.11,0.09,0.09,0.14,0.15,0.09,0.22,0.22])

b = np.array([0.58,0.92,1.10,1.22,1.31,1.39,1.59,1.77,1.87,1.92])
b_err = np.array([0.09,0.08,0.08,0.07,0.08,0.19,0.20,0.08,0.20,0.21])

beta_pear = a
beta_err_pear = a_err

alphaN_pear = a*(9.7 - 10.5) + b
alphaN_err_pear = ((b_err)**2+((9.7 - 10.5)*a_err)**2)**0.5

ssfr_pear = alphaN_pear - 9.7
ssfr_err_pear = alphaN_err_pear

sig_pear = np.array([0.35,0.33,0.34,0.28,0.24,0.29,0.28,0.10,0.15,0.08])
sig_err_pear = np.array([0.06,0.06,0.06,0.06,0.05,0.08,0.07,0.05,0.07,0.05])

# =============================================================================
# Donnari+19 - ILLUSTRIS TNG - TABLE 3
# =============================================================================
# scatter gets complicated, I also don't agree with their Speagle on their fig 7

# EQN 1: logSFR (sf galaxies) = a logM + b


z_donn = np.array([0.0, 0.75, 1.0, 1.75, 2.0])

# original paper
# beta_donn = np.array([0.80, 0.80, 0.81, 0.77, 0.77]) # a
# beta_err_donn = np.array([0.01, 0.01, 0.02, 0.02, 0.03])
# alpha_donn = np.array([-8.15, -7.72, -7.64, -6.97, -6.83]) # b
# alpha_err_donn = np.array([0.11, 0.03, 0.17, 0.24, 0.25])

# amended erratum
beta_donn = np.array([0.80, 0.81, 0.82, 0.78, 0.79]) # a
beta_err_donn = np.array([0.01, 0.01, 0.02, 0.02, 0.03])
alpha_donn = np.array([-8.09, -7.68, -7.63, -7.02, -6.93]) # b
alpha_err_donn = np.array([0.10, 0.03, 0.17, 0.24, 0.25])

alphaN_donn = alpha_donn + (9.7*beta_donn) # santini normalised
alphaN_err_donn = ((alpha_err_donn)**2+(beta_err_donn*9.7)**2)**0.5

ssfr_donn = alphaN_donn - 9.7
ssfr_err_donn = alphaN_err_donn

# =============================================================================
# Leslie+20 TABLE 1 - CHABRIER IMF
# =============================================================================

# adapted param from lee 15
# using SF values from table 1
# EQUATION 6
# log SFR = s0 - a1 t - log (1+ (10^M_t / 10^M))
# M_t = M0 - a2 t,      t in Gyr , M_t is turnover mass

# assumed a low mass slope of unity (ie gamma = 1)

# assuming z 0.3 to 6.0 from figure 4

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z_lesl = np.linspace(0.3, 6.0, 100)
t_lesl = cosmo.age(z_lesl).value # Gyr

# beta_lesl = np.full(100, 1.0)

# STAR FORMING VALUES (as my outlier model should leave us with SF...)
# ERRORS IN PAPER, but complicated
s0 = 2.97
M0 = 11.06
a1 = 0.22
a2 = 0.12

s0_err = 0.08
M0_err = 0.15
a1_err = 0.01
a2_err = 0.02

M_t = M0 - (a2*t_lesl)


alphaN_lesl = s0 - (a1*t_lesl) - np.log10(1+ ((10**M_t)/(10**9.7)))

alphaN_err_lesl = np.zeros(len(alphaN_lesl))

s0_d = np.random.normal(s0, s0_err, int(1e5))
M0_d = np.random.normal(M0, M0_err, int(1e5))
a1_d = np.random.normal(a1, a1_err, int(1e5))
a2_d = np.random.normal(a2, a2_err, int(1e5))

for i in range(len(alphaN_lesl)):

    M_t_d = M0_d - (a2_d*t_lesl[i])
    alphaN_err_lesl[i] = np.std(s0_d - (a1_d*t_lesl[i]) - np.log10(1+ ((10**M_t_d)/(10**9.7))))

ssfr_lesl = alphaN_lesl - 9.7
ssfr_err_lesl = alphaN_err_lesl


# print(ssfr_err_lesl)

'''
# error on a2 t
a = t_lesl
A = a2
sA = a2_err
f = a2*t_lesl
sf = abs(a) * sA

# error on M_t = M0 - a2 t
A = M0
sA = M0_err
B = f
sB = sf
f = A-B
sf = (sA**2 + sB**2)**0.5

# error on 10^M_t
a = 10.0
b = 1.0
A = f
sA = sf
f = a**(b*A)
sf = abs(f) * abs(b*np.log(a)*sA)

# error on (10^M_t / 10^M)
A = f
sA = sf
B = (10**9.7)
sB = 0.0
f = A/B
sf = abs(f) * (((sA/A)**2 + (sB/B)**2)**0.5)

# error on log (1+ (10^M_t / 10^M))
a = 1.0
b = 1.0
A = 1.0 + f
sA = sf
f = a*np.log10(b*A)
sf = abs((a*sA) / (A*np.log(10)))

# error on s0 - a1 t - log (1+ (10^M_t / 10^M))
A = s0
sA = s0_err
B = a1*t_lesl
sB = abs(t_lesl) * a1_err
C = f
sC = sf
f = A - B - C
sf = (sA**2 + sB**2 + sC**2)**0.5


alphaN_err_lesl = sf
print(sf)

# not much mention of scatter
'''

# =============================================================================
# Leja+21 - RIDGE (table 1) - assuming lower mass based on fig 4, z0.2 - 3.0 - CHABRIER IMF
# =============================================================================
# NO ERRORS IN TABLE 

# logSFR = a(M - M_t) + c for M > M_t
# logSFR = b(M - M_t) + c for M < M_t

a0 = np.array([-0.2384,1.204,-0.5929])
b0 = np.array([0.9387,0.005499,-0.02751])
c0 = np.array([0.3257,0.8805,-0.06114])
Mt0 = np.array([10.37,0.06952,0.1252])

z_leja = np.linspace(0.2, 3.0, 100)

a = a0[0] + a0[1]*z_leja + a0[2]*(z_leja**2)
b = b0[0] + b0[1]*z_leja + b0[2]*(z_leja**2)
c = c0[0] + c0[1]*z_leja + c0[2]*(z_leja**2)
Mt = Mt0[0] + Mt0[1]*z_leja + Mt0[2]*(z_leja**2)

alphaN_leja = b*(9.7 - Mt) + c
ssfr_leja = alphaN_leja - 9.7
beta_leja = b

# plt.plot(z_leja, Mt)
# plt.xlabel('z')
# plt.ylabel('mass that changes slope')
# plt.show()

# from fig 6
z_sig_leja = np.array([0.5, 1.5])
sig9p7_leja = np.array([0.38,0.36]) # refers to z 0.5 and 1.5
# sig10p3_leja = np.array([0.40,0.36,0.43])
# sig11p0_leja = np.array([0.34,0.32,0.33])

'''
# print(min(Mt)) # 10.388912, FOR MY EXAMPLE, ALWAYS ABOVE THIS


plt.plot(z, leja_alpha)
plt.show()

z_arr = np.array([0.3, 0.6, 1.0, 1.5, 2.0, 2.7])

for z in z_arr:
    a = a0[0] + a0[1]*z + a0[2]*(z**2)
    b = b0[0] + b0[1]*z + b0[2]*(z**2)
    c = c0[0] + c0[1]*z + c0[2]*(z**2)
    Mt = Mt0[0] + Mt0[1]*z + Mt0[2]*(z**2)

    Mlow = np.linspace(9.0, Mt, 100)
    Mhigh = np.linspace(Mt, 11.5, 100)
    
    plt.plot(Mlow, b*(Mlow - Mt) + c)
    plt.plot(Mhigh, a*(Mhigh - Mt) + c)

    alphaN = b*(9.7 - Mt) + c
    
    plt.scatter(9.7, alphaN)

plt.show()
'''

# =============================================================================
# Lovell+21 - from Lovell via Emma on skype - using high mass
# =============================================================================

z_love = np.array([10.,  9.,  8.,  7.,  6.,  5.])[3:]
x0 = np.array(['9.55', '9.16', '9.19', '9.35', '9.45', '9.60']).astype(float)[3:] # these values are actually x0+9.7 according to paper
# alpha1 = np.array(['1.24', '1.32', '1.33', '1.31', '1.27', '1.23']).astype(float) # mass < x0
alpha = np.array(['0.84', '0.91', '0.80', '0.72', '0.70', '0.62']).astype(float)[3:] # mass > x0
beta = np.array(['1.65', '1.22', '1.23', '1.30', '1.28', '1.30']).astype(float)[3:]

# alphaN_love1 = beta - (alpha1*(x0-9.7)) 
alphaN_love = beta - (alpha*(x0-9.7)) 

# beta_love1 = alpha1
beta_love = alpha

# ssfr_love1 = alphaN_love1 - 9.7
ssfr_love = alphaN_love - 9.7

# =============================================================================
# thorne 21 - CHABRIER IMF
# =============================================================================
#%%
# eqn 5, table D3, https://academic.oup.com/view-large/256108427
# adapted param from lee 15

z_low_thor = np.array([0.02,0.08,0.14,0.20,0.28,0.36,0.45,0.56,0.68,0.82,1.00,1.20,1.45,1.75,2.20,2.60,3.25,3.75,4.25,5.00])[:-2]
z_high_thor = np.array([0.08,0.14,0.20,0.28,0.36,0.45,0.56,0.68,0.82,1.00,1.20,1.45,1.75,2.20,2.60,3.25,3.75,4.25,5.00,9.00])[:-2]
s0 = np.array([0.064,0.139,0.3174,0.7597,1.1795,0.8962,0.5551,0.7267,0.9698,1.1263,1.3363,1.4696,1.5302,1.6857,1.8791,1.8457,2.3238,2.4844,2.6792,2.7513])[:-2]
s0_err = np.array([0.0019,0.00076,0.032,0.02,0.019,0.0087,0.0045,0.0086,0.016,0.028,0.0053,0.013,0.0092,0.0096,0.007,0.007,0.0099,0.0043,0.0057,0.013])[:-2]
M0 = np.array([9.5971,9.452,9.4452,10.1064,10.6332,10.2632,9.4763,9.5132,9.7331,9.8598,10.0925,10.1153,10.1547,10.3276,10.4174,10.1967,10.7865,10.9361,11.1114,10.9998])[:-2]
M0_err = np.array([0.011,0.0047,0.056,0.012,0.04,0.013,0.021,0.02,0.023,0.045,0.019,0.021,0.018,0.021,0.017,0.018,0.013,0.0078,0.0035,0.0098])[:-2]
a0 = np.array([0.9703,1.1515,1.1838,1.0728,0.9618,0.893,1.0229,1.0327,1.0036,0.972,0.9298,0.9463,0.9931,0.9374,0.9589,0.8964,0.9163,0.9013,0.9299,0.9861])[:-2]
a0_err = np.array([0.015,0.0074,0.052,0.013,0.016,0.0065,0.039,0.02,0.027,0.027,0.02,0.035,0.024,0.029,0.016,0.041,0.0098,0.0045,0.0029,0.011])[:-2]
b0 = np.array([0.187,0.1576,0.15,0.15,0.1997,0.2528,0.15,0.15,0.1519,0.1714,0.1902,0.1784,0.1662,0.2101,0.1903,0.2021,0.2117,0.2126,0.2066,0.2019])[:-2]
b0_err = np.array([0.038,0.016,0.0017,0.017,0.039,0.053,0.015,0.0046,0.022,0.02,0.022,0.036,0.02,0.047,0.034,0.042,0.05,0.05,0.042,0.067])[:-2]


z_thor = (z_low_thor+z_high_thor) / 2.0

# eqn 5: log SFR = s0 - log [ (10^M / 10^M0)^-a + (10^M / 10^M0)^-b ]

# beta_thor = a # low mass power law slope (check M0 for turnover, all M0>10 once z>1)
# beta_err_thor = a_err

alphaN_thor = s0 - np.log10((((10**9.7)/(10**M0))**-a0)+(((10**9.7)/(10**M0))**-b0)) # ERRORS HERE IS RIDICULOUS

alphaN_err_thor = np.zeros(len(alphaN_thor))

for i in range(len(s0)):
    s0_d = np.random.normal(s0[i], s0_err[i], int(1e5))
    M0_d = np.random.normal(M0[i], M0_err[i], int(1e5))
    a0_d = np.random.normal(a0[i], a0_err[i], int(1e5))
    b0_d = np.random.normal(b0[i], b0_err[i], int(1e5))

    alphaN_err_thor[i] = np.std(s0_d - np.log10((((10**9.7)/(10**M0_d))**-a0_d)+(((10**9.7)/(10**M0_d))**-b0_d)))

ssfr_thor = alphaN_thor - 9.7
ssfr_err_thor = alphaN_err_thor




'''
alphaN_thor = s0 - np.log10((((10**9.7)/(10**M0))**-a0)+(((10**9.7)/(10**M0))**-b0)) # ERRORS HERE IS RIDICULOUS

# error on 10**M0
a = 10.0
b = 1.0
A = M0
sA = M0_err
f = a**(b*A)
sf1 = abs(f) * abs(b*np.log(a)*sA)
# plt.plot(sf)
# plt.show()

# error on (10**9.7)/(10**M0)
A = (10**9.7)
sA = 0.0
B = f
sB = sf1
f = A/B
sf2 = abs(f) * (((sA/A)**2 + (sB/B)**2)**0.5)
# plt.plot(sf)
# plt.show()

# error on (((10**9.7)/(10**M0))**-a0)
A = f
sA = sf2
B = -a0
sB = a0_err
fa = A**B
sfa = abs(f) * (((B*sA/A)**2 + (sB*np.log(A))**2)**0.5)
# plt.plot(sfa)
# plt.show()

# error on (((10**9.7)/(10**M0))**-b0)
B = -b0
sB = b0_err
fb = A**B
sfb = abs(f) * (((B*sA/A)**2 + (sB*np.log(A))**2)**0.5)
# plt.plot(sfb)
# plt.show()

# error on (((10**9.7)/(10**M0))**-a0) + (((10**9.7)/(10**M0))**-b0)
A = fa
sA = sfa
B = fb
sB = sfb
f = A+B
sf3 = (sA**2 + sB**2)**0.5
# plt.plot(sf)
# plt.show()

# error on np.log10((((10**9.7)/(10**M0))**-a0)+(((10**9.7)/(10**M0))**-b0))
a = 1.0
b = 1.0
A = f
sA = sf3
f = a*np.log10(b*A)
sf4 = abs((a*sA) / (A*np.log(10)))
# plt.plot(sf)
# plt.show()

# error on s0 - np.log10((((10**9.7)/(10**M0))**-a0)+(((10**9.7)/(10**M0))**-b0))
A = s0
sA = s0_err
B = f
sB = sf4
f = A-B
sf5 = (sA**2 + sB**2)**0.5
# plt.plot(sf)
# plt.show()

alphaN_err_thor = sf5

ssfr_thor = alphaN_thor - 9.7
ssfr_err_thor = alphaN_err_thor




for i in range(len(s0)):
    s0_d = np.random.normal(s0[i], s0_err[i], int(1e5))
    M0_d = np.random.normal(M0[i], M0_err[i], int(1e5))
    a0_d = np.random.normal(a0[i], a0_err[i], int(1e5))
    b0_d = np.random.normal(b0[i], b0_err[i], int(1e5))
    
    
    
    d1 = 10**M0_d
    d2 = (10**9.7)/(10**M0_d)
    da = (((10**9.7)/(10**M0_d))**-a0_d)
    db = (((10**9.7)/(10**M0_d))**-b0_d)
    d3 = (((10**9.7)/(10**M0_d))**-a0_d) + (((10**9.7)/(10**M0_d))**-b0_d)
    d4 = np.log10((((10**9.7)/(10**M0_d))**-a0_d)+(((10**9.7)/(10**M0_d))**-b0_d))
    d5 = s0_d - np.log10((((10**9.7)/(10**M0_d))**-a0_d)+(((10**9.7)/(10**M0_d))**-b0_d))
    
        
    print('{:.1g} {:.1g}  {:.5f} {:.5f}  {:.5f} {:.5f}  {:.5f} {:.5f}  {:.5f} {:.5f}  {:.5f} {:.5f}  {:.5f} {:.5f} '.format(np.std(d1), sf1[i], np.std(d2), sf2[i], np.std(da), sfa[i], np.std(db), sfb[i], np.std(d3), sf3[i], np.std(d4), sf4[i], np.std(d5), sf5[i]))

for i in range(len(s0)):
    s0_d = np.random.normal(s0[i], s0_err[i], int(1e5))
    M0_d = np.random.normal(M0[i], M0_err[i], int(1e5))
    a0_d = np.random.normal(a0[i], a0_err[i], int(1e5))
    b0_d = np.random.normal(b0[i], b0_err[i], int(1e5))
    
    
    
    d1 = 10**M0_d
    d2 = (10**9.7)/(10**M0_d)
    da = (((10**9.7)/(10**M0_d))**-a0_d)
    db = (((10**9.7)/(10**M0_d))**-b0_d)
    d3 = (((10**9.7)/(10**M0_d))**-a0_d) + (((10**9.7)/(10**M0_d))**-b0_d)
    d4 = np.log10((((10**9.7)/(10**M0_d))**-a0_d)+(((10**9.7)/(10**M0_d))**-b0_d))
    d5 = s0_d - np.log10((((10**9.7)/(10**M0_d))**-a0_d)+(((10**9.7)/(10**M0_d))**-b0_d))
    
        
    print('{:.2g}   {:.2g}   {:.2g}   {:.2g}   {:.2g}   {:.2g}   {:.2g} '.format(np.std(d1)/ sf1[i], np.std(d2)/ sf2[i], np.std(da)/ sfa[i], np.std(db)/ sfb[i], np.std(d3)/ sf3[i], np.std(d4)/ sf4[i], np.std(d5)/ sf5[i]))   

'''





#%%







#%%




































































# =============================================================================
# plots
# =============================================================================

print(plt.rcParams['axes.prop_cycle'].by_key()['color'])
# 

ms_dot = 10
ms_tri = 6
ms_x = 8

plt.title('intercept')

plt.errorbar(z_donn, alphaN_donn, ls='dashed', marker='.', ms=0, label='Illustris-TNG')
plt.errorbar(z_love, alphaN_love, ls='dashed', marker='.', ms=0, label='FLARES')

plt.errorbar(z_spea, alphaN_spea, ls='solid', marker='.', ms=0, label='Speagle+14')
plt.errorbar(z_schr, alphaN_schr, ls='solid', marker='.', ms=0, label='Schreiber+15')
plt.errorbar(z_leja, alphaN_leja, ls='solid', marker='.', ms=0, label='Leja+21')

plt.errorbar(z_stei, alphaN_stei, yerr=alphaN_err_stei, ls='none', marker='.', ms=ms_dot, label='Steinhardt+14')
plt.errorbar(z_whit, alphaN_whit, yerr=alphaN_err_whit, ls='none', marker='.', ms=ms_dot, label='Whitaker+14')
plt.errorbar(z_salm, alphaN_salm1, yerr=alphaN_err_salm1, ls='none', marker='.', ms=ms_dot, label='Salmon+15')
plt.errorbar(z_kurc, alphaN_kurc, yerr=alphaN_err_kurc, ls='none', marker='.', ms=ms_dot, label='Kurczynski+16')
plt.errorbar(z_sant, alphaN_sant, yerr=alphaN_err_sant, ls='none', marker='.', ms=ms_dot, label='Santini+17')
plt.errorbar(z_boog, alphaN_boog, yerr=alphaN_err_boog, ls='none', marker='.', ms=ms_dot, label='Boogaard+18')
plt.errorbar(z_pear, alphaN_pear, yerr=alphaN_err_pear, ls='none', marker='.', ms=ms_dot, label='Pearson+18')



# plt.fill_between(z_spea, alphaN_spea-alphaN_err_spea, alphaN_spea+alphaN_err_spea, alpha=0.1, color='k')
plt.fill_between(z_schr, alphaN_schr-alphaN_err_schr, alphaN_schr+alphaN_err_schr, alpha=0.1, color='#d62728')

plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.4), ncol=4,fontsize=12)
plt.show()




plt.title('slope')

# plt.errorbar(z_spea, beta_spea, ls='solid', marker='.', ms=0, label='Speagle+14')
# plt.errorbar(z_schr, beta_schr, ls='solid', marker='.', ms=0, label='Schreiber+15')
# plt.errorbar(z_leja, beta_leja, ls='solid', marker='.', ms=0, label='Leja+21')

# plt.errorbar(z_stei, beta_stei, yerr=beta_err_stei, ls='none', marker='.', ms=ms_dot, label='Steinhardt+14')
# plt.errorbar(z_whit, beta_whit, yerr=beta_err_whit, ls='none', marker='.', ms=ms_dot, label='Whitaker+14')
# plt.errorbar(z_salm, beta_salm1, ls='none', marker='x', ms=8, label='Salmon+15')
# plt.errorbar(z_kurc, beta_kurc, yerr=beta_err_kurc, ls='none', marker='.', ms=ms_dot, label='Kurczynski+16')
# plt.errorbar(z_sant, beta_sant, yerr=beta_err_sant, ls='none', marker='.', ms=ms_dot, label='Santini+17')
# plt.errorbar(z_boog, beta_boog, yerr=beta_err_boog, ls='none', marker='.', ms=ms_dot, label='Boogaard+18')
# plt.errorbar(z_pear, beta_pear, yerr=beta_err_pear, ls='none', marker='.', ms=ms_dot, label='Pearson+18')
# plt.errorbar(z_donn, beta_donn, yerr=beta_err_donn, ls='none', marker='.', ms=ms_dot, label='Illustris-TNG')
# plt.errorbar(z_love, beta_love, ls='none', marker='x', ms=8, label='FLARES')


plt.errorbar(z_donn, beta_donn, ls='dashed', marker='.', ms=0, label='Illustris-TNG')
plt.errorbar(z_love, beta_love, ls='dashed', marker='.', ms=0, label='FLARES')

plt.errorbar(z_spea, beta_spea, ls='solid', marker='.', ms=0, label='Speagle+14')
plt.errorbar(z_schr, beta_schr, ls='solid', marker='.', ms=0, label='Schreiber+15')
plt.errorbar(z_leja, beta_leja, ls='solid', marker='.', ms=0, label='Leja+21')

plt.errorbar(z_stei, beta_stei, yerr=beta_err_stei, ls='none', marker='.', ms=ms_dot, label='Steinhardt+14')
plt.errorbar(z_whit, beta_whit, yerr=beta_err_whit, ls='none', marker='.', ms=ms_dot, label='Whitaker+14')
plt.errorbar(z_salm, beta_salm1, yerr=0, ls='none', marker='.', ms=ms_dot, label='Salmon+15')
plt.errorbar(z_kurc, beta_kurc, yerr=beta_err_kurc, ls='none', marker='.', ms=ms_dot, label='Kurczynski+16')
plt.errorbar(z_sant, beta_sant, yerr=beta_err_sant, ls='none', marker='.', ms=ms_dot, label='Santini+17')
plt.errorbar(z_boog, beta_boog, yerr=beta_err_boog, ls='none', marker='.', ms=ms_dot, label='Boogaard+18')
plt.errorbar(z_pear, beta_pear, yerr=beta_err_pear, ls='none', marker='.', ms=ms_dot, label='Pearson+18')

# plt.fill_between(z_spea, beta_spea-beta_err_spea, beta_spea+beta_err_spea, alpha=0.1, color='k')

plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.4), ncol=4,fontsize=12)
plt.show()

# plt.plot(0,1)






plt.title('scatter')

plt.plot(0, 0.3)
plt.plot(0, 0.3)
plt.errorbar(z_spea, sig_spea, ls='solid', marker='.', ms=0, label='Speagle+14')
plt.errorbar(z_sig_schr, sig_schr, ls='solid', marker='.', ms=0, label='Schreiber+15')
plt.errorbar(z_sig_leja, sig9p7_leja, ls='solid', marker='.', ms=0, label='Leja+21')

plt.errorbar(z_stei, sig_stei, ls='none', marker='x', ms=ms_x, label='Steinhardt+14')
plt.plot(0, 0.3)
plt.errorbar(z_salm, sig_salm, ls='none', marker='x', ms=ms_x, label='Salmon+15')

plt.errorbar(z_kurc, sig_kurc, yerr=sig_err_kurc, ls='none', marker='.', ms=ms_dot, label='Kurczynski+16')
plt.errorbar(z_sig_sant, sig9p2_sant, ls='none', marker='x', ms=ms_x, label='Santini+17')
plt.errorbar(z_boog, sig_boog, yerr=sig_err_boog, ls='none', marker='.', ms=ms_dot, label='Boogaard+18')
plt.errorbar(z_pear, sig_pear, yerr=sig_err_pear, ls='none', marker='.', ms=ms_dot, label='Pearson+18')




# plt.fill_between(z_spea, sig_spea-sig_err_spea, sig_spea+sig_err_spea, alpha=0.1, color='k')

plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.4), ncol=4,fontsize=12)
plt.show()
# plt.plot(0, 0.3)






plt.title('ssfr')

plt.errorbar(z_furl, ssfr_furl, ls='dashed', marker='.', ms=0, label='EAGLE')
# plt.errorbar(z_donn, ssfr_donn, yerr=ssfr_err_donn, ls='solid', marker='.', ms=0, label='Illustris-TNG')
plt.errorbar(z_donn, ssfr_donn, ls='dashed', marker='.', ms=0, label='Illustris-TNG')
plt.errorbar(z_love, ssfr_love, ls='dashed', marker='.', ms=0, label='FLARES')

plt.errorbar(z_spea, ssfr_spea, ls='solid', marker='.', ms=0, label='Speagle+14')
plt.errorbar(z_schr, ssfr_schr, ls='solid', marker='.', ms=0, label='Schreiber+15')
plt.errorbar(z_tomc, ssfr_tomc, ls='solid', marker='.', ms=0, label='Tomczak+16')
plt.errorbar(z_lesl, ssfr_lesl, ls='solid', marker='.', ms=0, label='Leslie+20')
plt.errorbar(z_leja, ssfr_leja, ls='solid', marker='.', ms=0, label='Leja+21')

plt.errorbar(z_redd, ssfr_redd, yerr=(ssfr_m_redd, ssfr_p_redd), ls='none', marker='.', ms=ms_dot, label='Reddy+12')
# plt.errorbar(z_barr, ssfr_barr, yerr=(ssfr_m_barr, ssfr_p_barr), ls='none', marker='.', ms=ms_dot, label='de Barros+14')
plt.errorbar(z_gonz, ssfr_gonz, yerr=(ssfr_m_gonz, ssfr_p_gonz), ls='none', marker='.', ms=ms_dot, label='Gonzalez+14')
plt.errorbar(z_smit, ssfr_smit, yerr=([ssfr_m_smit], [ssfr_p_smit]), ls='none', marker='.', ms=ms_dot, label='Smit+14')
plt.errorbar(z_stei, ssfr_stei, yerr=ssfr_err_stei, ls='none', marker='.', ms=ms_dot, label='Steinhardt+14')
plt.errorbar(z_whit, ssfr_whit, yerr=ssfr_err_whit, ls='none', marker='.', ms=ms_dot, label='Whitaker+14')
plt.errorbar(z_lee, ssfr_lee, yerr=ssfr_err_lee, ls='none', marker='.', ms=ms_dot, label='Lee+15')
plt.errorbar(z_salm, ssfr_salm1, yerr=ssfr_err_salm1, ls='none', marker='.', ms=ms_dot, label='Salmon+15')
plt.errorbar(z_kurc, ssfr_kurc, yerr=ssfr_err_kurc, ls='none', marker='.', ms=ms_dot, label='Kurczynski+16')
plt.errorbar(z_marm, ssfr_marm, yerr=ssfr_err_marm, ls='none', marker='.', ms=ms_dot, label='Marmol-Queralto+16')
plt.errorbar(z_sant, ssfr_sant, yerr=ssfr_err_sant, ls='none', marker='.', ms=ms_dot, label='Santini+17')
plt.errorbar(z_boog, ssfr_boog, yerr=ssfr_err_boog, ls='none', marker='x', ms=ms_tri, label='Boogaard+18')
plt.errorbar(z_pear, ssfr_pear, yerr=ssfr_err_pear, ls='none', marker='x', ms=ms_tri, label='Pearson+18')
plt.errorbar(z_thor, ssfr_thor, yerr=ssfr_err_thor, ls='none', marker='x', ms=ms_tri, label='Thorne+21')

# plt.plot(z_spea, ssfr_spea, label='Speagle+14')
# plt.plot(z_furl, ssfr_furl, label='EAGLE')
# plt.plot(z_schr, ssfr_schr, label='Schreiber+15')
# plt.plot(z_tomc, ssfr_tomc, label='Tomczak+16')
# plt.plot(z_donn, ssfr_donn, label='Illustris-TNG')
# plt.plot(z_lesl, ssfr_lesl, label='Leslie+20')
# plt.plot(z_leja, ssfr_leja, label='Leja+21')
# plt.plot(z_love, ssfr_love, label='FLARES')

# plt.fill_between(z_spea, ssfr_spea-ssfr_err_spea, ssfr_spea+ssfr_err_spea, alpha=0.1, color='k')
plt.fill_between(z_schr, ssfr_schr-ssfr_err_schr, ssfr_schr+ssfr_err_schr, alpha=0.1, color='#9467bd')
plt.fill_between(z_lesl, ssfr_lesl-ssfr_err_lesl, ssfr_lesl+ssfr_err_lesl, alpha=0.1, color='#e377c2')

plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.6), ncol=4,fontsize=12)
plt.show()


# plt.errorbar(0, -9, ls='none')



# =============================================================================
# data thief
# =============================================================================


# =============================================================================
# # Salim+07
# =============================================================================
z_s07 = np.array([0.0992])
ssfr_s07 = np.array([-0.8106]) - 9.0
ssfr_p_s07 = np.array([-0.296]) - ssfr_s07 - 9.0
ssfr_m_s07 = ssfr_s07 + 9.0 - np.array([-1.2918])

# =============================================================================
# # Santini+09
# =============================================================================
z_s09 = np.array([0.4429, 0.8013, 1.2547, 2.0007])
ssfr_s09 = np.array([-0.2388, -0.0578, 0.1185, 0.2995]) - 9.0
ssfr_p_s09 = np.array([0.0566, 0.2376, 0.4187, 0.595]) - ssfr_s09 - 9.0
ssfr_m_s09 = ssfr_s09 + 9.0 - np.array([-0.5342, -0.358, -0.1817, 4.1691e-3])

# =============================================================================
# # Stark+13
# =============================================================================
z_s13 = np.array([3.7999, 4.9774, 5.8916, 5.8916, 6.7911, 6.7911])
ssfr_s13 = np.array([0.7522, 0.7332, 0.8332, 0.7808, 1.1191, 0.9524]) - 9.0
ssfr_p_s13 = np.array([1.0715, 1.0476, 1.1048, 1.1048, 1.2859, 1.2859]) - ssfr_s13 - 9.0
ssfr_m_s13 = ssfr_s13 + 9.0 - np.array([0.4663, 0.4378, 0.4997, 0.4997, 0.676, 0.676])

# =============================================================================
# # Bouwens+12
# =============================================================================
z_b12 = np.array([4.0046, 5.0066, 6.0086, 7.208])
ssfr_b12 = np.array([0.7141, 0.6855, 0.4997, 0.7284]) - 9.0
ssfr_p_b12 = np.array([1.0048, 0.9762, 0.7904, 0.8475]) - ssfr_b12 - 9.0
ssfr_m_b12 = ssfr_b12 + 9.0 - np.array([0.4139, 0.3853, 0.1995, 0.6236])

# =============================================================================
# # Menci et al. (2014) SAM
# =============================================================================
z_m14 = np.array([0.1504, 0.3624, 0.5819, 0.8525, 1.0719, 1.4522, 1.8763, 2.403, 2.9149, 3.4487, 4.0632, 4.8969, 5.6429, 6.228, 6.7473])
ssfr_m14 = np.array([-1.0965, -0.9297, -0.8154, -0.6915, -0.5867, -0.458, -0.3389, -0.215, -0.1197, -0.0339, 0.0661, 0.1948, 0.3139, 0.4139, 0.514]) - 9.0



