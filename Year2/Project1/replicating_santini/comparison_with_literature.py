#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 17:07:08 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

norm = 9.7
#norm = 0

# =============================================================================
# Santini+17 'True' values - delayed SFH, SFR from UV slope
# =============================================================================
# logSFR = alpha log(M / M_9.7) + beta
z_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))

# converting normalisation
alpha_san = B_san - 9.7*A_san
alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
beta_san = A_san
beta_err_san = A_err_san

alpha_san_n = alpha_san + (norm*beta_san) # santini normalised
alpha_err_san_n = (alpha_err_san**2 - (norm*beta_err_san)**2) ** 0.5

# =============================================================================
# Santini+17 Original values - obtained by eye - delayed SFH, SFR from UV slope
# =============================================================================
# logSFR = alpha log(M / M_9.7) + beta
z_san0 = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san0 = np.array((1.05, 1.1, 0.9, 0.75, 0.55))
A_err_san0 = np.array((0.03, 0.03, 0.04, 0.05, 0.18))
B_san0 = np.array((1.0, 1.15, 1.25, 1.2, 1.5))
B_err_san0 = np.array((0.05, 0.03, 0.03, 0.06, 0.12))

# converting normalisation
alpha_san0 = B_san0 - 9.7*A_san0
alpha_err_san0 = (B_err_san0**2 + (9.7*A_err_san0)**2) ** 0.5
beta_san0 = A_san0
beta_err_san0 = A_err_san0

alpha_san0_n = alpha_san0 + (norm*beta_san0) # santini normalised
alpha_err_san0_n = (alpha_err_san0**2 - (norm*beta_err_san0)**2) ** 0.5

# =============================================================================
# Speagle+14 - errors calculated dodgily
# =============================================================================
# log SFR(M∗, t) = (0.84 ± 0.02 − 0.026 ± 0.003 × t ) logM∗−(6.51 ± 0.24 − 0.11 ± 0.03 × t ), where t is the age of the universe in Gyr.

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z_speagle = np.linspace(0.5, 6.5, 1000)
t_speagle = cosmo.age(z_speagle).value # Gyr

alpha_speagle = -(6.51 - (0.11*t_speagle) )
alpha_err_speagle = 0.24 + (0.03*t_speagle)
beta_speagle = 0.84 - (0.026*t_speagle)
beta_err_speagle = 0.02 + (0.003*t_speagle)

alpha_speagle_n = alpha_speagle + (norm*beta_speagle) # santini normalised


# =============================================================================
# Schreiber+15 - ignoring high mass
# =============================================================================
#r ≡ log10(1 + z) and m ≡ log10(M∗/10**9 M):
#log10(SFRMS[M/yr]) = m − m0 + a0r − a1 max(0,m − m1 − a2 r)**2, 
#with m0 = 0.5 ± 0.07, a0 = 1.5 ± 0.15, a1 = 0.3 ± 0.08, m1 = 0.36 ± 0.3 and a2 = 2.5 ± 0.6.

z_schreiber = np.linspace(0.5, 4.0, 1000)
r_schreiber = np.log10(1+z_schreiber)

m0_schreiber = 0.5
a0_schreiber = 1.5
#a1 = 0.3
#m1 = 0.36
#a2 = 2.5

# m - m1 - a2r is usually < 0, except high mass, low redshift, IGNORED FOR NOW
#print( np.log10((10**9.9)/(1e9))- m1 - (a2*r_schreiber))

alpha_schreiber = - (9.0 + m0_schreiber - (a0_schreiber*r_schreiber))
beta_schreiber = np.linspace(1.0, 1.0, 1000)

alpha_schreiber_n = alpha_schreiber + (norm*beta_schreiber) # santini normalised


# =============================================================================
# Salmon+15
# =============================================================================
#log(SFR/M yr−1) = a log(M/M) + b
z_salmon = np.array((4.0, 5.0, 6.0))

alpha_salmon = np.array((-5.7, -4.4, -3.9))
alpha_err_salmon = np.array((2.1, 2.6, 1.6))
beta_salmon = np.array((0.7, 0.59, 0.54))
beta_err_salmon = np.array((0.21, 0.26, 0.16))

alpha_salmon_n = alpha_salmon + (norm*beta_salmon) # santini normalised


# =============================================================================
# Steinhardt+14
# =============================================================================
#log SFR(M yr−1) = α × (logM∗/M − 10.5) + β,
z_steinhardt = np.array(((4.0 + 4.8)/2, (4.8 + 6.0)/2))

beta_steinhardt = np.array((0.78, 0.78))
beta_err_steinhardt = np.array((0.02, 0.02))

alpha_steinhardt = np.array((1.976, 2.110)) - (10.5*beta_steinhardt)

alpha_steinhardt_n = alpha_steinhardt + (norm*beta_steinhardt) # santini normalised


# =============================================================================
# Tomczak+16 - Not obvious how to get these, also maybe not just star forming?!
# =============================================================================

# =============================================================================
# Schreiber+16 - single point, hard to find paper
# =============================================================================

# =============================================================================
# Kurczynski+16 - nice table, also mentions purpose of rescaling which I haven't included here
# =============================================================================
#log SFR = a ´ log M* + b + N (0, sIS).
z_kurc = np.array(((0.5 + 1.0)/2, (1.0 + 1.5)/2, (1.5 + 2.0)/2, (2.0 + 2.5)/2, (2.5 + 3.0)/2))

alpha_kurc = np.array((-8.394, -7.474, -7.484, -7.513, -7.729))
alpha_err_kurc = np.array((0.011, 0.010, 0.011, 0.018, 0.015))
beta_kurc = np.array((0.919, 0.825, 0.867, 0.849, 0.899))
beta_err_kurc = np.array((0.017, 0.012, 0.013, 0.021, 0.017))

alpha_kurc_n = alpha_kurc + (norm*beta_kurc) # santini normalised


# =============================================================================
# lester options - santini replication - need to add errors from bootstrapping
# =============================================================================

# H<27.5, relflag==1

# SAME + zLow<redshift_BEAGLE<zHigh
# SAME + mass_BEAGLE_stellar
# Changed limit on BEAGLE mass to account for mass offset
# Adjusting upper limit to account for mass offset
# SAME + sfr_BEAGLE_instant
# all fields

z_lester = np.array((1.65, 2.5, 3.5, 4.5, 5.5))

alpha_lester000000 = np.array((-8.524643477752727, -8.421954970372497, -6.81530836339029, -8.774787326855662, -2.7516185290043182))
beta_lester000000 = np.array((0.9721675837301412, 0.9589392930321784, 0.7864571528369195, 0.9977181422947992, 0.36909767808861654))
alpha_lester000000_n = alpha_lester000000 + (norm*beta_lester000000) # santini normalised
alpha_err_lester000000 = np.array((alpha_lester000000-np.array([-8.965173503834054, -8.62591428306176, -7.3032213480092, -8.929546305130764, -3.759212511161998]),np.array([-7.3499119955525245, -7.303149125924152, -6.337159747807995, -6.2888382320874925, 0.1117672150686556])-alpha_lester000000))
beta_err_lester000000 = np.array((beta_lester000000-np.array([0.8436245497642918, 0.8323755525507877, 0.7349917315432323, 0.7424350984579859, 0.085210921473392]),np.array([1.0210829236927124, 0.9813506148440654, 0.841235584935741, 1.0133281066171462, 0.4805563720021437])-beta_lester000000))

alpha_lester100000 = np.array((-7.938601389499966, -8.648701922459916, -7.325059385052887, -7.238641395772859, -3.175608153493875))
beta_lester100000 = np.array((0.9062358589355383, 0.9776541256160342, 0.8405873243691521, 0.8370590593738213, 0.4145557651335486))
alpha_lester100000_n = alpha_lester100000 + (norm*beta_lester100000) # santini normalised

alpha_lester110000 = np.array((-8.019339196435595, -8.875556878616448, -6.707382637830726, -5.870674944494401, -10.22962113643443))
beta_lester110000 = np.array((0.9491537547631866, 1.044446807379201, 0.8027721806629016, 0.7611033765206862, 1.2796525994363472))
alpha_lester110000_n = alpha_lester110000 + (norm*beta_lester110000) # santini normalised

alpha_lester111000 = np.array((-7.860311124795692, -9.308325656673535, -6.724376607216616, -5.955283195668764, -8.675523184536532))
beta_lester111000 = np.array((0.9323005205934546, 1.0927520726071906, 0.8057271301397584, 0.7710875646785492, 1.1156962100969041))
alpha_lester111000_n = alpha_lester111000 + (norm*beta_lester111000) # santini normalised

alpha_lester111100 = np.array((-7.86634837801764, -9.151894734686179, -6.724376607216616, -5.955283195668764, -8.956691861260962))
beta_lester111100 = np.array((0.9330219556517038, 1.0743744234390307, 0.8057271301397584, 0.7710875646785492, 1.1490799579403068))
alpha_lester111100_n = alpha_lester111100 + (norm*beta_lester111100) # santini normalised

alpha_lester111110 = np.array((-6.098666698061404, -7.418510150450283, -4.930440373595574, -6.198887664149959, -7.858423777063218))
beta_lester111110 = np.array((0.6981389277951645, 0.8669147984936145, 0.6149290921456422, 0.8044709744823518, 1.0204438424091684))
alpha_lester111110_n = alpha_lester111110 + (norm*beta_lester111110) # santini normalised

alpha_lester000010 = np.array((-6.604878482435299, -6.9033940912282326, -6.431762602018932, -8.116330731863638, -1.9152692550467652))
beta_lester000010 = np.array((0.7392697125801366, 0.7988515518864913, 0.7595798660689483, 0.9359176658661896, 0.30443949603609305))
alpha_lester000010_n = alpha_lester000010 + (norm*beta_lester000010) # santini normalised


# =============================================================================
# lester scenarios
# =============================================================================

alpha_scenario_1_1 = np.array([-8.956260029539061, -9.380651514667798, -7.821309930829113, -6.067132815619622, -4.217687552662691])
beta_scenario_1_1 = np.array([1.0267709079513785, 1.087005395242499, 0.9356425100394657, 0.7499678870408129, 0.588139667320449])
alpha_scenario_1_1_n =  alpha_scenario_1_1 + (norm*beta_scenario_1_1) # santini normalised
alpha_err_scenario_1_1 = np.array((alpha_scenario_1_1-np.array([-9.423285040374944, -9.7209970402005, -8.238955370098019, -6.812799365800618, -9.190299261530159]),np.array([-8.258680070373597, -8.83325152641616, -6.748474245042315, -5.216901466538098, -2.0426290937985363])-alpha_scenario_1_1))
beta_err_scenario_1_1 = np.array((beta_scenario_1_1-np.array([0.9496166582341242, 1.0251275259404125, 0.8172957846116141, 0.653879474561798, 0.34195299262532497]),np.array([1.077967133890977, 1.124586775982536, 0.9807128890354087, 0.8309488626359105, 1.1518392959790955])-beta_scenario_1_1))
alpha_err_scenario_1_1_n = abs(((((alpha_err_scenario_1_1[1]+alpha_err_scenario_1_1[0])/2))**2 - (norm*(((beta_err_scenario_1_1[1]+beta_err_scenario_1_1[0])/2)))**2)) ** 0.5

alpha_scenario_5_1 = np.array([-9.19057034844203, -7.785133004659811, -2.7044434817154612, -3.981913516594088, -2.624200119563876])
beta_scenario_5_1 = np.array([1.0382700779886236, 0.9084132284500976, 0.37561526914054105, 0.5264576871099392, 0.404930660880207])
alpha_scenario_5_1_n = alpha_scenario_5_1 + (norm*beta_scenario_5_1) # santini normalised
alpha_err_scenario_5_1 = np.array((alpha_scenario_5_1-np.array([-9.669166643587046, -8.510818960194317, -4.106033004213572, -5.06079122215039, -5.551185437730209]),np.array([-7.739284859809073, -5.678209452172945, -2.14945673207135, -2.6175676306372115, -1.607135544934253])-alpha_scenario_5_1))
beta_err_scenario_5_1 = np.array((beta_scenario_5_1-np.array([0.8847236965393324, 0.6728948549950081, 0.3136947811020121, 0.37149446550606324, 0.2810257798354253]),np.array([1.0914992647885113, 0.988413159886452, 0.5357165301209235, 0.6509161129243018, 0.7637819295034786])-beta_scenario_5_1))
alpha_err_scenario_5_1_n = abs(((((alpha_err_scenario_5_1[1]+alpha_err_scenario_5_1[0])/2))**2 - (norm*(((beta_err_scenario_5_1[1]+beta_err_scenario_5_1[0])/2)))**2)) ** 0.5

alpha_scenario_5_4 = np.array([-8.658012377717828, -9.499261913490278, -6.1476123512757805, -6.90722057031772, -8.755773191518736])
beta_scenario_5_4 = np.array([1.0211198464432367, 1.1334097737716526, 0.7775598178315407, 0.8751543200088223, 1.1219718596670003])
alpha_scenario_5_4_n = alpha_scenario_5_4 + (norm*beta_scenario_5_4) # santini normalised
alpha_err_scenario_5_4 = np.array((alpha_scenario_5_4-np.array([-9.455940092395462, -9.831407500239955, -7.00526503418614, -8.085375760080208, -15.24505535182757]),np.array([-7.608486631047726, -8.486215548677459, -5.723390507991851, -4.831419986219626, -5.6209603628963904])-alpha_scenario_5_4))
beta_err_scenario_5_4 = np.array((beta_scenario_5_4-np.array([0.9025809555046863, 1.019717814648489, 0.7269164443560806, 0.6403181233337241, 0.7498855915249537]),np.array([1.113184457253098, 1.171366270942317, 0.8710337054486679, 1.009959561487698, 1.9106274793060822])-beta_scenario_5_4))
alpha_err_scenario_5_4_n = abs(((((alpha_err_scenario_5_4[1]+alpha_err_scenario_5_4[0])/2))**2 - (norm*(((beta_err_scenario_5_4[1]+beta_err_scenario_5_4[0])/2)))**2)) ** 0.5

alpha_scenario_5_6 = np.array([-8.246235583266522, -8.53548421172312, -2.5242184857457164, -2.9111214934700946, -3.2942397344721535])
beta_scenario_5_6 = np.array([0.9136631943948326, 0.9620181592470548, 0.33733752538536665, 0.3967666491356541, 0.4943958300719963])
alpha_scenario_5_6_n = alpha_scenario_5_6 + (norm*beta_scenario_5_6) # santini normalised
alpha_err_scenario_5_6 = np.array((alpha_scenario_5_6-np.array([-8.631823250034355, -9.479370087735047, -3.4661632111794543, -4.985667041869504, -4.549433358036185]),np.array([-7.00864735472677, -3.8863845191798516, -1.938570401580888, -2.4707739439227296, -1.6590484102113772])-alpha_scenario_5_6))
beta_err_scenario_5_6 = np.array((beta_scenario_5_6-np.array([0.7783640766431933, 0.4449906337512449, 0.2723977299886322, 0.3504463238258511, 0.2983281150364063]),np.array([0.9575148377787933, 1.0625614684774811, 0.44363881908348135, 0.6385330133705875, 0.6422905081531785])-beta_scenario_5_6))
alpha_err_scenario_5_6_n = abs(((((alpha_err_scenario_5_6[1]+alpha_err_scenario_5_6[0])/2))**2 - (norm*(((beta_err_scenario_5_6[1]+beta_err_scenario_5_6[0])/2)))**2)) ** 0.5

alpha_scenario_5_7 = np.array([-5.930074382373524, -7.132979447177554, -5.152714845683838, -6.97357110019294, -7.292248760365191])
beta_scenario_5_7 = np.array([0.6660660396935073, 0.8344657611425867, 0.642689012211567, 0.8860155758731968, 0.9705069028938051])
alpha_scenario_5_7_n = alpha_scenario_5_7 + (norm*beta_scenario_5_7) # santini normalised
alpha_err_scenario_5_7 = np.array((alpha_scenario_5_7-np.array([-8.072994900198138, -9.001270366246127, -5.703576396355865, -7.858338871518982, -10.12272039964168]),np.array([-5.784868520212281, -6.146976997702867, -4.635564939472313, -5.754294678631154, -5.056870837067507])-alpha_scenario_5_7))
beta_err_scenario_5_7 = np.array((beta_scenario_5_7-np.array([0.6534686744940024, 0.7194187266253252, 0.5839456643566201, 0.7491165304052826, 0.7205661119970764]),np.array([0.9266709825219148, 1.0544231938864492, 0.7053717670838697, 0.9872947407919423, 1.3026300735672762])-beta_scenario_5_7))
alpha_err_scenario_5_7_n = abs(((((alpha_err_scenario_5_7[1]+alpha_err_scenario_5_7[0])/2))**2 - (norm*(((beta_err_scenario_5_7[1]+beta_err_scenario_5_7[0])/2)))**2)) ** 0.5

# =============================================================================
# kelly applied to scenario 5
# =============================================================================

alpha_scenario_5_k = np.array([-7.221,-7.967,-6.191,-9.955,-11.618])
beta_scenario_5_k = np.array([0.841,0.930,0.746,1.239,1.464])
alpha_scenario_5_k_n = alpha_scenario_5_k + (norm*beta_scenario_5_k) # santini normalised
alpha_err_scenario_5_k = np.array((alpha_scenario_5_k-np.array([-8.070,-8.746,-7.833,-13.356,-23.747]),np.array([-6.309,-6.942,-5.088,-7.932,-2.776])-alpha_scenario_5_k))
beta_err_scenario_5_k = np.array((beta_scenario_5_k-np.array([0.738,0.812,0.604,0.994,0.472]),np.array([0.938,1.021,0.941,1.642,2.865])-beta_scenario_5_k))
alpha_err_scenario_5_k_n = abs(((((alpha_err_scenario_5_k[1]+alpha_err_scenario_5_k[0])/2))**2 - (norm*(((beta_err_scenario_5_k[1]+beta_err_scenario_5_k[0])/2)))**2)) ** 0.5

alpha_scenario_5_k0 = np.array([-7.005,-7.539,-4.838,-7.460,-15.009])
beta_scenario_5_k0 = np.array([0.814,0.880,0.601,0.957,1.862])
alpha_scenario_5_k0_n = alpha_scenario_5_k0 + (norm*beta_scenario_5_k0) # santini normalised
alpha_err_scenario_5_k0 = np.array((alpha_scenario_5_k0-np.array([-7.899,-8.478,-6.250,-10.087,-31.111]),np.array([-6.067,-6.745,-3.517,-5.224,-3.375])-alpha_scenario_5_k0))
beta_err_scenario_5_k0 = np.array((beta_scenario_5_k0-np.array([0.708,0.790,0.456,0.704,0.521]),np.array([0.918,0.991,0.760,1.255,3.625])-beta_scenario_5_k0))
alpha_err_scenario_5_k0_n = abs(((((alpha_err_scenario_5_k0[1]+alpha_err_scenario_5_k0[0])/2))**2 - (norm*(((beta_err_scenario_5_k0[1]+beta_err_scenario_5_k0[0])/2)))**2)) ** 0.5


# =============================================================================
# scenario 6
# =============================================================================

alpha_scenario_6_1 = np.array([-8.88188305291079, -7.489462097874605, -2.7840924009706445, -3.7403895731346712, -3.682451423164204])
beta_scenario_6_1 = np.array([1.0062692495173964, 0.8748612826591646, 0.38273792660486944, 0.4980572789896197, 0.5225313073066552])
alpha_scenario_6_1_n = alpha_scenario_6_1 + (norm*beta_scenario_6_1) # santini normalised
alpha_err_scenario_6_1 = np.array((alpha_scenario_6_1-np.array([-9.140243258947148, -8.303021484484871, -3.2431165989221102, -4.385140873984328, -4.310607812197381]),np.array([-7.217987671033388, -6.97724224631726, -2.1978074789870727, -3.0243067765805853, -2.3882554437842867])-alpha_scenario_6_1))
beta_err_scenario_6_1 = np.array((beta_scenario_6_1-np.array([0.823122976623502, 0.8141891436940116, 0.3153936498089162, 0.4109505153429864, 0.36270164777384517]),np.array([1.037435303931541, 0.9624116043177267, 0.43548406476809476, 0.5750316004635354, 0.6031075699449371])-beta_scenario_6_1))
alpha_err_scenario_6_1_n = abs(((((alpha_err_scenario_6_1[1]+alpha_err_scenario_6_1[0])/2))**2 - (norm*(((beta_err_scenario_6_1[1]+beta_err_scenario_6_1[0])/2)))**2)) ** 0.5

alpha_scenario_6_4 = np.array([-7.525971990819344, -7.966143975241509, -6.4998776845300785, -5.793122853226448, -4.879173606988697])
beta_scenario_6_4 = np.array([0.8973973280819733, 0.9604700490076792, 0.8190472750726082, 0.7570345593713296, 0.6575691460865578])
alpha_scenario_6_4_n = alpha_scenario_6_4 + (norm*beta_scenario_6_4) # santini normalised
alpha_err_scenario_6_4 = np.array((alpha_scenario_6_4-np.array([-8.251783258798348, -8.318940320774583, -7.056834032236463, -6.261323319379805, -6.325267707101737]),np.array([-6.763751223339136, -7.49179003665181, -5.340492799398845, -4.898013137446472, -3.442679020223564])-alpha_scenario_6_4))
beta_err_scenario_6_4 = np.array((beta_scenario_6_4-np.array([0.8057638194862033, 0.9085025856603692, 0.6846836868751275, 0.6514336935548594, 0.4949203967541863]),np.array([0.9818290333441889, 1.0034914161415693, 0.8827065605086474, 0.8106701299156557, 0.8505761155505599])-beta_scenario_6_4))
alpha_err_scenario_6_4_n = abs(((((alpha_err_scenario_6_4[1]+alpha_err_scenario_6_4[0])/2))**2 - (norm*(((beta_err_scenario_6_4[1]+beta_err_scenario_6_4[0])/2)))**2)) ** 0.5

alpha_scenario_6_6 = np.array([-7.3317759154033, -6.764849821352314, -2.5016269826024464, -2.9655602761772193, -4.822407198810928])
beta_scenario_6_6 = np.array([0.8179685764244055, 0.768665803686183, 0.3311942314595104, 0.4024316609124473, 0.6664592246796924])
alpha_scenario_6_6_n = alpha_scenario_6_6 + (norm*beta_scenario_6_6) # santini normalised
alpha_err_scenario_6_6 = np.array((alpha_scenario_6_6-np.array([-7.894731243357739, -7.5232488213606095, -3.174588892702917, -4.268586155439379, -5.131116470079184]),np.array([-6.2897076928215405, -6.0630196587779075, -1.7494866996814675, -2.577111180035621, -3.0924435448604575])-alpha_scenario_6_6))
beta_err_scenario_6_6 = np.array((beta_scenario_6_6-np.array([0.6977559504667705, 0.6879150341137981, 0.24526029094294188, 0.35825997322419084, 0.4446959295431766]),np.array([0.8817355306921633, 0.8544151517380246, 0.40996650906424087, 0.5572070668471459, 0.6968152309243711])-beta_scenario_6_6))
alpha_err_scenario_6_6_n = abs(((((alpha_err_scenario_6_6[1]+alpha_err_scenario_6_6[0])/2))**2 - (norm*(((beta_err_scenario_6_6[1]+beta_err_scenario_6_6[0])/2)))**2)) ** 0.5

alpha_scenario_6_7 = np.array([-5.131113356957769, -5.655457996946231, -5.030941316094558, -6.325791517032839, -6.7581848744598965])
beta_scenario_6_7 = np.array([0.5834301387908689, 0.666711932410436, 0.6342806091625183, 0.8196649019516414, 0.8813698479494465])
alpha_scenario_6_7_n = alpha_scenario_6_7 + (norm*beta_scenario_6_7) # santini normalised
alpha_err_scenario_6_7 = np.array((alpha_scenario_6_7-np.array([-6.152907847748257, -6.334093396423188, -5.837956354766358, -6.7856054844006755, -8.075188777019154]),np.array([-4.703311224430971, -4.8916316246092695, -4.467279984396902, -5.6262030739003555, -5.914711787609433])-alpha_scenario_6_7))
beta_err_scenario_6_7 = np.array((beta_scenario_6_7-np.array([0.526853605318538, 0.5782400303165434, 0.5701555638296307, 0.7324816516706052, 0.7967395447013543]),np.array([0.7084524656739071, 0.747321942971973, 0.7299006670581536, 0.8758913603661521, 1.052803304617541])-beta_scenario_6_7))
alpha_err_scenario_6_7_n = abs(((((alpha_err_scenario_6_7[1]+alpha_err_scenario_6_7[0])/2))**2 - (norm*(((beta_err_scenario_6_7[1]+beta_err_scenario_6_7[0])/2)))**2)) ** 0.5

# =============================================================================
# kelly applied to scenario 6
# =============================================================================

alpha_scenario_6_k = np.array([-6.570, -6.987, -6.728, -7.290, -12.189])
beta_scenario_6_k = np.array([0.771, 0.829, 0.826, 0.957, 1.556])
alpha_scenario_6_k_n = alpha_scenario_6_k + (norm*beta_scenario_6_k) # santini normalised
alpha_err_scenario_6_k = np.array((alpha_scenario_6_k-np.array([-7.217, -7.595, -7.425, -8.445, -13.971]),np.array([-6.010, -6.346, -5.621, -6.381, -10.152])-alpha_scenario_6_k))
beta_err_scenario_6_k = np.array((beta_scenario_6_k-np.array([0.699, 0.745, 0.688, 0.850, 1.322]),np.array([0.845, 0.899, 0.910, 1.096, 1.778])-beta_scenario_6_k))
alpha_err_scenario_6_k_n = abs(((((alpha_err_scenario_6_k[1]+alpha_err_scenario_6_k[0])/2))**2 - (norm*(((beta_err_scenario_6_k[1]+beta_err_scenario_6_k[0])/2)))**2)) ** 0.5

alpha_scenario_6_k0 = np.array([-6.522, -6.977, -5.586, -6.687, -14.591])
beta_scenario_6_k0 = np.array([0.763, 0.823, 0.690, 0.888, 1.833])
alpha_scenario_6_k0_n = alpha_scenario_6_k0 + (norm*beta_scenario_6_k0) # santini normalised
alpha_err_scenario_6_k0 = np.array((alpha_scenario_6_k0-np.array([-7.216, -7.574, -6.370, -8.163, -17.885]),np.array([-5.920, -6.230, -4.709, -5.415, -10.567])-alpha_scenario_6_k0))
beta_err_scenario_6_k0 = np.array((beta_scenario_6_k0-np.array([0.690, 0.737, 0.590, 0.734, 1.369]),np.array([0.846, 0.896, 0.778, 1.058, 2.230])-beta_scenario_6_k0))
alpha_err_scenario_6_k0_n = abs(((((alpha_err_scenario_6_k0[1]+alpha_err_scenario_6_k0[0])/2))**2 - (norm*(((beta_err_scenario_6_k0[1]+beta_err_scenario_6_k0[0])/2)))**2)) ** 0.5  

# =============================================================================
# kelly applied to scenario 7 delayed uniform logtau
# =============================================================================

#alpha_scenario_7_7_d = np.array([])
#beta_scenario_7_7_d = np.array([])
#alpha_scenario_7_7_d_n = alpha_scenario_7_7_d + (norm*beta_scenario_7_7_d) # santini normalised
#alpha_err_scenario_7_7_d = np.array((alpha_scenario_7_7_d-np.array([]),np.array([])-alpha_scenario_7_7_d))
#beta_err_scenario_7_7_d = np.array((beta_scenario_7_7_d-np.array([]),np.array([])-beta_scenario_7_7_d))
#alpha_err_scenario_7_7_d_n = abs(((((alpha_err_scenario_7_7_d[1]+alpha_err_scenario_7_7_d[0])/2))**2 - (norm*(((beta_err_scenario_7_7_d[1]+beta_err_scenario_7_7_d[0])/2)))**2)) ** 0.5

alpha_scenario_7_k0_d = np.array([-6.602, -7.066, -5.363, -6.037, -14.595])
beta_scenario_7_k0_d = np.array([0.773, 0.834, 0.664, 0.812, 1.83])
alpha_scenario_7_k0_d_n = alpha_scenario_7_k0_d + (norm*beta_scenario_7_k0_d) # santini normalised
alpha_err_scenario_7_k0_d = np.array((alpha_scenario_7_k0_d-np.array([-7.188, -7.704, -6.272, -7.38, -17.152]),np.array([-6.008, -6.426, -4.519, -4.746, -12.044])-alpha_scenario_7_k0_d))
beta_err_scenario_7_k0_d = np.array((beta_scenario_7_k0_d-np.array([0.701, 0.759, 0.57, 0.655, 1.538]),np.array([0.841, 0.911, 0.77, 0.97, 2.139])-beta_scenario_7_k0_d))
alpha_err_scenario_7_k0_d_n = abs(((((alpha_err_scenario_7_k0_d[1]+alpha_err_scenario_7_k0_d[0])/2))**2 - (norm*(((beta_err_scenario_7_k0_d[1]+beta_err_scenario_7_k0_d[0])/2)))**2)) ** 0.5

'''
alpha_scenario_7_7_d = np.array([])
beta_scenario_7_7_d = np.array([])
alpha_scenario_7_7_d_n = alpha_scenario_7_7_d + (norm*beta_scenario_7_7_d) # santini normalised
alpha_err_scenario_7_7_d = np.array((alpha_scenario_7_7_d-np.array([]),np.array([])-alpha_scenario_7_7_d))
beta_err_scenario_7_7_d = np.array((beta_scenario_7_7_d-np.array([]),np.array([])-beta_scenario_7_7_d))
alpha_err_scenario_7_7_d_n = abs(((((alpha_err_scenario_7_7_d[1]+alpha_err_scenario_7_7_d[0])/2))**2 - (norm*(((beta_err_scenario_7_7_d[1]+beta_err_scenario_7_7_d[0])/2)))**2)) ** 0.5

alpha_scenario_7_k0_d = np.array([])
beta_scenario_7_k0_d = np.array([])
alpha_scenario_7_k0_d_n = alpha_scenario_7_k0_d + (norm*beta_scenario_7_k0_d) # santini normalised
alpha_err_scenario_7_k0_d = np.array((alpha_scenario_7_k0_d-np.array([]),np.array([])-alpha_scenario_7_k0_d))
beta_err_scenario_7_k0_d = np.array((beta_scenario_7_k0_d-np.array([]),np.array([])-beta_scenario_7_k0_d))
alpha_err_scenario_7_k0_d_n = abs(((((alpha_err_scenario_7_k0_d[1]+alpha_err_scenario_7_k0_d[0])/2))**2 - (norm*(((beta_err_scenario_7_k0_d[1]+beta_err_scenario_7_k0_d[0])/2)))**2)) ** 0.5  
'''

# =============================================================================
# kelly applied to scenario 7 danger delayed uniform logtau
# =============================================================================

#alpha_scenario_7_7_dd = np.array([])
#beta_scenario_7_7_dd = np.array([])
#alpha_scenario_7_7_dd_n = alpha_scenario_7_7_dd + (norm*beta_scenario_7_7_dd) # santini normalised
#alpha_err_scenario_7_7_dd = np.array((alpha_scenario_7_7_dd-np.array([]),np.array([])-alpha_scenario_7_7_dd))
#beta_err_scenario_7_7_dd = np.array((beta_scenario_7_7_dd-np.array([]),np.array([])-beta_scenario_7_7_dd))
#alpha_err_scenario_7_7_dd_n = abs(((((alpha_err_scenario_7_7_dd[1]+alpha_err_scenario_7_7_dd[0])/2))**2 - (norm*(((beta_err_scenario_7_7_dd[1]+beta_err_scenario_7_7_dd[0])/2)))**2)) ** 0.5

alpha_scenario_7_k0_dd = np.array([-5.862, -7.532, -4.722, -5.029, -13.657])
beta_scenario_7_k0_dd = np.array([0.69, 0.887, 0.594, 0.69, 1.727])
alpha_scenario_7_k0_dd_n = alpha_scenario_7_k0_dd + (norm*beta_scenario_7_k0_dd) # santini normalised
alpha_err_scenario_7_k0_dd = np.array((alpha_scenario_7_k0_dd-np.array([-6.527, -8.093, -5.684, -6.31, -17.332]),np.array([-5.239, -6.92, -3.709, -3.775, -10.197])-alpha_scenario_7_k0_dd))
beta_err_scenario_7_k0_dd = np.array((beta_scenario_7_k0_dd-np.array([0.615, 0.813, 0.48, 0.539, 1.321]),np.array([0.769, 0.95, 0.703, 0.842, 2.162])-beta_scenario_7_k0_dd))
alpha_err_scenario_7_k0_dd_n = abs(((((alpha_err_scenario_7_k0_dd[1]+alpha_err_scenario_7_k0_dd[0])/2))**2 - (norm*(((beta_err_scenario_7_k0_dd[1]+beta_err_scenario_7_k0_dd[0])/2)))**2)) ** 0.5  

# =============================================================================
# kelly applied to scenario 7 danger constant
# =============================================================================

#alpha_scenario_7_7_dc = np.array([])
#beta_scenario_7_7_dc = np.array([])
#alpha_scenario_7_7_dc_n = alpha_scenario_7_7_dc + (norm*beta_scenario_7_7_dc) # santini normalised
#alpha_err_scenario_7_7_dc = np.array((alpha_scenario_7_7_dc-np.array([]),np.array([])-alpha_scenario_7_7_dc))
#beta_err_scenario_7_7_dc = np.array((beta_scenario_7_7_dc-np.array([]),np.array([])-beta_scenario_7_7_dc))
#alpha_err_scenario_7_7_dc_n = abs(((((alpha_err_scenario_7_7_dc[1]+alpha_err_scenario_7_7_dc[0])/2))**2 - (norm*(((beta_err_scenario_7_7_dc[1]+beta_err_scenario_7_7_dc[0])/2)))**2)) ** 0.5

alpha_scenario_7_k0_dc = np.array([-6.779, -7.727, -7.799, -5.205, -15.054])
beta_scenario_7_k0_dc = np.array([0.785, 0.9, 0.932, 0.674, 1.858])
alpha_scenario_7_k0_dc_n = alpha_scenario_7_k0_dc + (norm*beta_scenario_7_k0_dc) # santini normalised
alpha_err_scenario_7_k0_dc = np.array((alpha_scenario_7_k0_dc-np.array([-7.386, -8.241, -8.76, -6.333, -19.045]),np.array([-6.14, -7.173, -6.907, -4.15, -8.832])-alpha_scenario_7_k0_dc))
beta_err_scenario_7_k0_dc = np.array((beta_scenario_7_k0_dc-np.array([0.708, 0.835, 0.83, 0.549, 1.128]),np.array([0.858, 0.961, 1.043, 0.805, 2.334])-beta_scenario_7_k0_dc))
alpha_err_scenario_7_k0_dc_n = abs(((((alpha_err_scenario_7_k0_dc[1]+alpha_err_scenario_7_k0_dc[0])/2))**2 - (norm*(((beta_err_scenario_7_k0_dc[1]+beta_err_scenario_7_k0_dc[0])/2)))**2)) ** 0.5  


# =============================================================================
# kelly applied to scenario 7, delayed, danger delayed and danger constant, with HOGG
# =============================================================================

alpha_scenario_7_k0_d_Hogg = np.array([-7.697, -7.104, -6.291, -6.362, -11.383])
beta_scenario_7_k0_d_Hogg = np.array([0.91, 0.825, 0.762, 0.84, 1.463])
alpha_scenario_7_k0_d_Hogg_n = alpha_scenario_7_k0_d_Hogg + (norm*beta_scenario_7_k0_d_Hogg) # santini normalised
alpha_err_scenario_7_k0_d_Hogg = np.array((alpha_scenario_7_k0_d_Hogg-np.array([-8.596, -7.578, -6.957, -7.826, -16.617]),np.array([-6.806, -6.639, -5.732, -5.13, -8.228])-alpha_scenario_7_k0_d_Hogg))
beta_err_scenario_7_k0_d_Hogg = np.array((beta_scenario_7_k0_d_Hogg-np.array([0.805, 0.77, 0.697, 0.693, 1.07]),np.array([1.021, 0.879, 0.837, 1.016, 2.089])-beta_scenario_7_k0_d_Hogg))
alpha_err_scenario_7_k0_d_Hogg_n = abs(((((alpha_err_scenario_7_k0_d_Hogg[1]+alpha_err_scenario_7_k0_d_Hogg[0])/2))**2 - (norm*(((beta_err_scenario_7_k0_d_Hogg[1]+beta_err_scenario_7_k0_d_Hogg[0])/2)))**2)) ** 0.5

alpha_scenario_7_k0_dd_Hogg = np.array([-5.449, -7.528, -6.021, -5.994, -9.907])
beta_scenario_7_k0_dd_Hogg = np.array([0.615, 0.877, 0.733, 0.77, 1.271])
alpha_scenario_7_k0_dd_Hogg_n = alpha_scenario_7_k0_dd_Hogg + (norm*beta_scenario_7_k0_dd_Hogg) # santini normalised
alpha_err_scenario_7_k0_dd_Hogg = np.array((alpha_scenario_7_k0_dd_Hogg-np.array([-6.089, -8.11, -6.803, -7.492, -15.434]),np.array([-4.902, -7.02, -5.19, -4.508, -6.165])-alpha_scenario_7_k0_dd_Hogg))
beta_err_scenario_7_k0_dd_Hogg = np.array((beta_scenario_7_k0_dd_Hogg-np.array([0.551, 0.817, 0.639, 0.597, 0.84]),np.array([0.693, 0.944, 0.82, 0.947, 1.946])-beta_scenario_7_k0_dd_Hogg))
alpha_err_scenario_7_k0_dd_Hogg_n = abs(((((alpha_err_scenario_7_k0_dd_Hogg[1]+alpha_err_scenario_7_k0_dd_Hogg[0])/2))**2 - (norm*(((beta_err_scenario_7_k0_dd_Hogg[1]+beta_err_scenario_7_k0_dd_Hogg[0])/2)))**2)) ** 0.5  

alpha_scenario_7_k0_dc_Hogg = np.array([-6.204, -7.61, -7.867, -6.582, -8.937])
beta_scenario_7_k0_dc_Hogg = np.array([0.685, 0.867, 0.921, 0.789, 1.139])
alpha_scenario_7_k0_dc_Hogg_n = alpha_scenario_7_k0_dc_Hogg + (norm*beta_scenario_7_k0_dc_Hogg) # santini normalised
alpha_err_scenario_7_k0_dc_Hogg = np.array((alpha_scenario_7_k0_dc_Hogg-np.array([-6.775, -7.964, -8.466, -7.663, -14.083]),np.array([-5.72, -7.33, -7.203, -5.408, -5.028])-alpha_scenario_7_k0_dc_Hogg))
beta_err_scenario_7_k0_dc_Hogg = np.array((beta_scenario_7_k0_dc_Hogg-np.array([0.628, 0.832, 0.848, 0.657, 0.649]),np.array([0.754, 0.908, 0.991, 0.918, 1.758])-beta_scenario_7_k0_dc_Hogg))
alpha_err_scenario_7_k0_dc_Hogg_n = abs(((((alpha_err_scenario_7_k0_dc_Hogg[1]+alpha_err_scenario_7_k0_dc_Hogg[0])/2))**2 - (norm*(((beta_err_scenario_7_k0_dc_Hogg[1]+beta_err_scenario_7_k0_dc_Hogg[0])/2)))**2)) ** 0.5  




















'''
# =============================================================================
# PLOT slope
# =============================================================================
plt.figure(figsize=(6, 6))
plt.title('slope')

plt.scatter(z_san, beta_san, label='Santini+17')
plt.errorbar(z_san, beta_san, yerr=beta_err_san, ls='none')
plt.scatter(z_san0, beta_san0, label='Santini+17 Raw')
plt.errorbar(z_san0, beta_san0, yerr=beta_err_san0, ls='none')
#plt.scatter(z_salmon, beta_salmon, label='Salmon+15')
#plt.errorbar(z_salmon, beta_salmon, yerr=beta_err_salmon, ls='none')
#plt.scatter(z_steinhardt, beta_steinhardt, label='Steinhardt+14')
#plt.errorbar(z_steinhardt, beta_steinhardt, yerr=beta_err_steinhardt, ls='none')
#plt.scatter(z_kurc, beta_kurc, label='Kurczynski+16')
#plt.errorbar(z_kurc, beta_kurc, yerr=beta_err_kurc, ls='none')

plt.scatter(z_lester+0.1, beta_lester000000, label='Lester000000')
#plt.scatter(z_lester, beta_lester100000, label='Lester100000')
#plt.scatter(z_lester, beta_lester110000, label='Lester110000')
#plt.scatter(z_lester, beta_lester111000, label='Lester111000')
#plt.scatter(z_lester, beta_lester111100, label='Lester111100')
#plt.scatter(z_lester, beta_lester111110, label='Lester111110')
#plt.scatter(z_lester, beta_lester000010, label='Lester000010')

#plt.plot(z_lester, beta_lester000000, label='Lester000000')
#plt.plot(z_lester, beta_lester100000, label='Lester100000')
#plt.plot(z_lester, beta_lester110000, label='Lester110000')
#plt.plot(z_lester, beta_lester111000, label='Lester111000')
#plt.plot(z_lester, beta_lester111100, label='Lester111100')
#plt.plot(z_lester, beta_lester111110, label='Lester111110')
#plt.plot(z_lester, beta_lester000010, label='Lester000010')

plt.errorbar(z_lester+0.1, beta_lester000000, yerr=beta_err_lester000000, label='Lester000000', ls='none')
#plt.errorbar(z_lester, beta_lester100000, yerr=beta_err_lester100000, label='Lester100000', ls='none')
#plt.errorbar(z_lester, beta_lester110000, yerr=beta_err_lester110000, label='Lester110000', ls='none')
#plt.errorbar(z_lester, beta_lester111000, yerr=beta_err_lester111000, label='Lester111000', ls='none')
#plt.errorbar(z_lester, beta_lester111100, yerr=beta_err_lester111100, label='Lester111100', ls='none')
#plt.errorbar(z_lester, beta_lester111110, yerr=beta_err_lester111110, label='Lester111110', ls='none')
#plt.errorbar(z_lester, beta_lester000010, yerr=beta_err_lester000010, label='Lester000010', ls='none')

#plt.plot(z_speagle, beta_speagle, label='Speagle+14')
#plt.plot(z_schreiber, beta_schreiber, label='Schreiber+15')
plt.legend()
plt.show()

# =============================================================================
# PLOT intercept
# =============================================================================
plt.figure(figsize=(6, 6))
plt.title('intercept')

plt.scatter(z_san, alpha_san, label='Santini+17')
plt.plot([], [])
plt.scatter(z_san0, alpha_san0, label='Santini+17 Raw')
plt.plot([], [])
#plt.scatter(z_salmon, alpha_salmon, label='Salmon+15')
#plt.scatter(z_steinhardt, alpha_steinhardt, label='Steinhardt+14')
#plt.scatter(z_kurc, alpha_kurc, label='Kurczynski+16')

#plt.scatter(z_lester, alpha_lester000000, label='Lester000000')
#plt.scatter(z_lester, alpha_lester100000, label='Lester100000')
#plt.scatter(z_lester, alpha_lester110000, label='Lester110000')
#plt.scatter(z_lester, alpha_lester111000, label='Lester111000')
#plt.scatter(z_lester, alpha_lester111100, label='Lester111100')
#plt.scatter(z_lester, alpha_lester111110, label='Lester111110')
#plt.scatter(z_lester, alpha_lester000010, label='Lester000010')

plt.plot(z_lester, alpha_lester000000, label='Lester000000')
plt.plot(z_lester, alpha_lester100000, label='Lester100000')
#plt.plot(z_lester, alpha_lester110000, label='Lester110000')
#plt.plot(z_lester, alpha_lester111000, label='Lester111000')
#plt.plot(z_lester, alpha_lester111100, label='Lester111100')
#plt.plot(z_lester, alpha_lester111110, label='Lester111110')
plt.plot(z_lester, alpha_lester000010, label='Lester000010')

#plt.plot(z_speagle, alpha_speagle, label='Speagle+14')
#plt.plot(z_schreiber, alpha_schreiber, label='Schreiber+15')
plt.legend()
plt.show()

# =============================================================================
# PLOT normalised as in Santini - slope unaffected
# =============================================================================
plt.figure(figsize=(6, 6))
plt.title('santini normalised intercept')

plt.scatter(z_san, alpha_san_n, label='Santini+17')
plt.plot([], [])
plt.scatter(z_san0, alpha_san0_n, label='Santini+17 Raw')
plt.plot([], [])
#plt.scatter(z_salmon, alpha_salmon_n, label='Salmon+15')
#plt.scatter(z_steinhardt, alpha_steinhardt_n, label='Steinhardt+14')
#plt.scatter(z_kurc, alpha_kurc_n, label='Kurczynski+16')

#plt.scatter(z_lester, alpha_lester000000_n, label='Lester000000')
#plt.scatter(z_lester, alpha_lester100000_n, label='Lester100000')
#plt.scatter(z_lester, alpha_lester110000_n, label='Lester110000')
#plt.scatter(z_lester, alpha_lester111000_n, label='Lester111000')
#plt.scatter(z_lester, alpha_lester111100_n, label='Lester111100')
#plt.scatter(z_lester, alpha_lester111110_n, label='Lester111110')
#plt.scatter(z_lester, alpha_lester000010_n, label='Lester000010')

plt.plot(z_lester, alpha_lester000000_n, label='Lester000000')
plt.plot(z_lester, alpha_lester100000_n, label='Lester100000')
#plt.plot(z_lester, alpha_lester110000_n, label='Lester110000')
#plt.plot(z_lester, alpha_lester111000_n, label='Lester111000')
#plt.plot(z_lester, alpha_lester111100_n, label='Lester111100')
#plt.plot(z_lester, alpha_lester111110_n, label='Lester111110')
plt.plot(z_lester, alpha_lester000010_n, label='Lester000010')

#plt.plot(z_speagle, alpha_speagle_n, label='Speagle+14')
#plt.plot(z_schreiber, alpha_schreiber_n, label='Schreiber+15')
plt.legend()
plt.show()

'''

# =============================================================================
# UPDATED COMPARISONS SINCE HAVING SANTINI RESULTS
# =============================================================================

inc = 0.03 # increment in offset

# =============================================================================
# PLOT slope
# =============================================================================

plt.figure(figsize=(15, 6))
plt.title('slope')

plt.scatter(z_san, beta_san, label='Santini+17')
plt.errorbar(z_san, beta_san, yerr=beta_err_san, ls='none')
plt.scatter(z_san0, beta_san0, label='Santini+17 Raw')
plt.errorbar(z_san0, beta_san0, yerr=beta_err_san0, ls='none')
#plt.scatter(z_salmon, beta_salmon, label='Salmon+15')
#plt.errorbar(z_salmon, beta_salmon, yerr=beta_err_salmon, ls='none')
#plt.scatter(z_steinhardt, beta_steinhardt, label='Steinhardt+14')
#plt.errorbar(z_steinhardt, beta_steinhardt, yerr=beta_err_steinhardt, ls='none')
#plt.scatter(z_kurc, beta_kurc, label='Kurczynski+16')
#plt.errorbar(z_kurc, beta_kurc, yerr=beta_err_kurc, ls='none')

#plt.scatter(z_lester+ 1*inc, beta_scenario_1_1, label='scenario 1 1 Santini M+SFR', marker='x')
#plt.scatter(z_lester+ 2*inc, beta_scenario_5_1, label='scenario 5 1 Santini M+SFR', marker='x')
#plt.scatter(z_lester+ 3*inc, beta_scenario_5_4, label='scenario 5 4 Santini SFR', marker='x')
#plt.scatter(z_lester+ 4*inc, beta_scenario_5_6, label='scenario 5 6 Santini M', marker='x')
#plt.scatter(z_lester+ 5*inc, beta_scenario_5_7, label='scenario 5 7', marker='x')
#plt.scatter(z_lester+ 6*inc, beta_scenario_5_k0, label='scenario 5 k0', marker='x')
#plt.scatter(z_lester+ 7*inc, beta_scenario_5_k, label='scenario 5 k', marker='x')

#plt.errorbar(z_lester+ 1*inc, beta_scenario_1_1, yerr=beta_err_scenario_1_1, ls='none')
#plt.errorbar(z_lester+ 2*inc, beta_scenario_5_1, yerr=beta_err_scenario_5_1, ls='none')
#plt.errorbar(z_lester+ 3*inc, beta_scenario_5_4, yerr=beta_err_scenario_5_4, ls='none')
#plt.errorbar(z_lester+ 4*inc, beta_scenario_5_6, yerr=beta_err_scenario_5_6, ls='none')
#plt.errorbar(z_lester+ 5*inc, beta_scenario_5_7, yerr=beta_err_scenario_5_7, ls='none')
#plt.errorbar(z_lester+ 6*inc, beta_scenario_5_k0, yerr=beta_err_scenario_5_k0, ls='none')
#plt.errorbar(z_lester+ 7*inc, beta_scenario_5_k, yerr=beta_err_scenario_5_k, ls='none')

#plt.gca().set_prop_cycle(None)
#plt.scatter((),())
#plt.scatter((),())
plt.scatter(z_lester+ 1*inc, beta_scenario_1_1, label='scenario 1 1 Santini M+SFR')
plt.scatter(z_lester+ 2*inc, beta_scenario_6_1, label='scenario 6 1 Santini M+SFR')
plt.scatter(z_lester+ 3*inc, beta_scenario_6_4, label='scenario 6 4 Santini SFR')
plt.scatter(z_lester+ 4*inc, beta_scenario_6_6, label='scenario 6 6 Santini M')
plt.scatter(z_lester+ 5*inc, beta_scenario_6_7, label='scenario 6 7')
plt.scatter(z_lester+ 6*inc, beta_scenario_6_k0, label='scenario 6 k0')
plt.scatter(z_lester+ 7*inc, beta_scenario_6_k, label='scenario 6 k')
plt.scatter(z_lester+ 8*inc, beta_scenario_7_k0_d, label='scenario 7 k0 d')
plt.scatter(z_lester+ 9*inc, beta_scenario_7_k0_dd, label='scenario 7 k0 dd')
plt.scatter(z_lester+ 10*inc, beta_scenario_7_k0_dc, label='scenario 7 k0 dc')
plt.scatter(z_lester+ 11*inc, beta_scenario_7_k0_d_Hogg, label='scenario 7 k0 d Hogg')
plt.scatter(z_lester+ 12*inc, beta_scenario_7_k0_dd_Hogg, label='scenario 7 k0 dd Hogg')
plt.scatter(z_lester+ 13*inc, beta_scenario_7_k0_dc_Hogg, label='scenario 7 k0 dc Hogg')

plt.errorbar(z_lester+ 1*inc, beta_scenario_1_1, yerr=beta_err_scenario_1_1, ls='none')
plt.errorbar(z_lester+ 2*inc, beta_scenario_6_1, yerr=beta_err_scenario_6_1, ls='none')
plt.errorbar(z_lester+ 3*inc, beta_scenario_6_4, yerr=beta_err_scenario_6_4, ls='none')
plt.errorbar(z_lester+ 4*inc, beta_scenario_6_6, yerr=beta_err_scenario_6_6, ls='none')
plt.errorbar(z_lester+ 5*inc, beta_scenario_6_7, yerr=beta_err_scenario_6_7, ls='none')
plt.errorbar(z_lester+ 6*inc, beta_scenario_6_k0, yerr=beta_err_scenario_6_k0, ls='none')
plt.errorbar(z_lester+ 7*inc, beta_scenario_6_k, yerr=beta_err_scenario_6_k, ls='none')
plt.errorbar(z_lester+ 8*inc, beta_scenario_7_k0_d, yerr=beta_err_scenario_7_k0_d, ls='none')
plt.errorbar(z_lester+ 9*inc, beta_scenario_7_k0_dd, yerr=beta_err_scenario_7_k0_dd, ls='none')
plt.errorbar(z_lester+ 10*inc, beta_scenario_7_k0_dc, yerr=beta_err_scenario_7_k0_dc, ls='none')
plt.errorbar(z_lester+ 11*inc, beta_scenario_7_k0_d_Hogg, yerr=beta_err_scenario_7_k0_d_Hogg, ls='none')
plt.errorbar(z_lester+ 12*inc, beta_scenario_7_k0_dd_Hogg, yerr=beta_err_scenario_7_k0_dd_Hogg, ls='none')
plt.errorbar(z_lester+ 13*inc, beta_scenario_7_k0_dc_Hogg, yerr=beta_err_scenario_7_k0_dc_Hogg, ls='none')

#plt.plot(z_speagle, beta_speagle, label='Speagle+14')
#plt.plot(z_schreiber, beta_schreiber, label='Schreiber+15')
plt.ylim(0.0, 1.5)
plt.legend()
plt.show()

# =============================================================================
# PLOT intercept
# =============================================================================
plt.figure(figsize=(15, 6))
plt.title('intercept')

plt.scatter(z_san, alpha_san, label='Santini+17')
plt.errorbar(z_san, alpha_san, yerr=alpha_err_san, ls='none')
#plt.plot([], [])
plt.scatter(z_san0, alpha_san0, label='Santini+17 Raw')
plt.errorbar(z_san0, alpha_san0, yerr=alpha_err_san0, ls='none')
#plt.plot([], [])
#plt.scatter(z_salmon, alpha_salmon, label='Salmon+15')
#plt.scatter(z_steinhardt, alpha_steinhardt, label='Steinhardt+14')
#plt.scatter(z_kurc, alpha_kurc, label='Kurczynski+16')

#plt.scatter(z_lester+ 1*inc, alpha_scenario_1_1, label='scenario 1 1 Santini M+SFR', marker='x')
#plt.scatter(z_lester+ 2*inc, alpha_scenario_5_1, label='scenario 5 1 Santini M+SFR', marker='x')
#plt.scatter(z_lester+ 3*inc, alpha_scenario_5_4, label='scenario 5 4 Santini SFR', marker='x')
#plt.scatter(z_lester+ 4*inc, alpha_scenario_5_6, label='scenario 5 6 Santini M', marker='x')
#plt.scatter(z_lester+ 5*inc, alpha_scenario_5_7, label='scenario 5 7', marker='x')
#plt.scatter(z_lester+ 6*inc, alpha_scenario_5_k0, label='scenario 5 k0', marker='x')
#plt.scatter(z_lester+ 7*inc, alpha_scenario_5_k, label='scenario 5 k', marker='x')

#plt.errorbar(z_lester+ 1*inc, alpha_scenario_1_1, yerr=alpha_err_scenario_1_1, ls='none')
#plt.errorbar(z_lester+ 2*inc, alpha_scenario_5_1, yerr=alpha_err_scenario_5_1, ls='none')
#plt.errorbar(z_lester+ 3*inc, alpha_scenario_5_4, yerr=alpha_err_scenario_5_4, ls='none')
#plt.errorbar(z_lester+ 4*inc, alpha_scenario_5_6, yerr=alpha_err_scenario_5_6, ls='none')
#plt.errorbar(z_lester+ 5*inc, alpha_scenario_5_7, yerr=alpha_err_scenario_5_7, ls='none')
#plt.errorbar(z_lester+ 6*inc, alpha_scenario_5_k0, yerr=alpha_err_scenario_5_k0, ls='none')
#plt.errorbar(z_lester+ 7*inc, alpha_scenario_5_k, yerr=alpha_err_scenario_5_k, ls='none')

#plt.gca().set_prop_cycle(None)
#plt.scatter((),())
#plt.scatter((),())
plt.scatter(z_lester+ 1*inc, alpha_scenario_1_1, label='scenario 1 1 Santini M+SFR')
plt.scatter(z_lester+ 2*inc, alpha_scenario_6_1, label='scenario 6 1 Santini M+SFR')
plt.scatter(z_lester+ 3*inc, alpha_scenario_6_4, label='scenario 6 4 Santini SFR')
plt.scatter(z_lester+ 4*inc, alpha_scenario_6_6, label='scenario 6 6 Santini M')
plt.scatter(z_lester+ 5*inc, alpha_scenario_6_7, label='scenario 6 7')
plt.scatter(z_lester+ 6*inc, alpha_scenario_6_k0, label='scenario 6 k0')
plt.scatter(z_lester+ 7*inc, alpha_scenario_6_k, label='scenario 6 k')
plt.scatter(z_lester+ 8*inc, alpha_scenario_7_k0_d, label='scenario 7 k0 d')
plt.scatter(z_lester+ 9*inc, alpha_scenario_7_k0_dd, label='scenario 7 k0 dd')
plt.scatter(z_lester+ 10*inc, alpha_scenario_7_k0_dc, label='scenario 7 k0 dc')
plt.scatter(z_lester+ 11*inc, alpha_scenario_7_k0_d_Hogg, label='scenario 7 k0 d Hogg')
plt.scatter(z_lester+ 12*inc, alpha_scenario_7_k0_dd_Hogg, label='scenario 7 k0 dd Hogg')
plt.scatter(z_lester+ 13*inc, alpha_scenario_7_k0_dc_Hogg, label='scenario 7 k0 dc Hogg')

plt.errorbar(z_lester+ 1*inc, alpha_scenario_1_1, yerr=alpha_err_scenario_1_1, ls='none')
plt.errorbar(z_lester+ 2*inc, alpha_scenario_6_1, yerr=alpha_err_scenario_6_1, ls='none')
plt.errorbar(z_lester+ 3*inc, alpha_scenario_6_4, yerr=alpha_err_scenario_6_4, ls='none')
plt.errorbar(z_lester+ 4*inc, alpha_scenario_6_6, yerr=alpha_err_scenario_6_6, ls='none')
plt.errorbar(z_lester+ 5*inc, alpha_scenario_6_7, yerr=alpha_err_scenario_6_7, ls='none')
plt.errorbar(z_lester+ 6*inc, alpha_scenario_6_k0, yerr=alpha_err_scenario_6_k0, ls='none')
plt.errorbar(z_lester+ 7*inc, alpha_scenario_6_k, yerr=alpha_err_scenario_6_k, ls='none')
plt.errorbar(z_lester+ 8*inc, alpha_scenario_7_k0_d, yerr=alpha_err_scenario_7_k0_d, ls='none')
plt.errorbar(z_lester+ 9*inc, alpha_scenario_7_k0_dd, yerr=alpha_err_scenario_7_k0_dd, ls='none')
plt.errorbar(z_lester+ 10*inc, alpha_scenario_7_k0_dc, yerr=alpha_err_scenario_7_k0_dc, ls='none')
plt.errorbar(z_lester+ 11*inc, alpha_scenario_7_k0_d_Hogg, yerr=alpha_err_scenario_7_k0_d_Hogg, ls='none')
plt.errorbar(z_lester+ 12*inc, alpha_scenario_7_k0_dd_Hogg, yerr=alpha_err_scenario_7_k0_dd_Hogg, ls='none')
plt.errorbar(z_lester+ 13*inc, alpha_scenario_7_k0_dc_Hogg, yerr=alpha_err_scenario_7_k0_dc_Hogg, ls='none')

#plt.plot(z_speagle, alpha_speagle, label='Speagle+14')
#plt.plot(z_schreiber, alpha_schreiber, label='Schreiber+15')
plt.ylim(-11, -2)
plt.legend()
plt.show()

# =============================================================================
# PLOT normalised as in Santini - slope unaffected
# =============================================================================
plt.figure(figsize=(15, 6))
plt.title('santini normalised intercept')

plt.scatter(z_san, alpha_san_n, label='Santini+17')
plt.errorbar(z_san, alpha_san_n, yerr=alpha_err_san_n, ls='none')
#plt.plot([], [])
plt.scatter(z_san0, alpha_san0_n, label='Santini+17 Raw')
plt.errorbar(z_san0, alpha_san0_n, yerr=alpha_err_san0_n, ls='none')
#plt.plot([], [])
#plt.scatter(z_salmon, alpha_salmon_n, label='Salmon+15')
#plt.scatter(z_steinhardt, alpha_steinhardt_n, label='Steinhardt+14')
#plt.scatter(z_kurc, alpha_kurc_n, label='Kurczynski+16')

#plt.scatter(z_lester+ 1*inc, alpha_scenario_1_1_n, label='scenario 1 1 Santini M+SFR', marker='x')
#plt.scatter(z_lester+ 2*inc, alpha_scenario_5_1_n, label='scenario 5 1 Santini M+SFR', marker='x')
#plt.scatter(z_lester+ 3*inc, alpha_scenario_5_4_n, label='scenario 5 4 Santini SFR', marker='x')
#plt.scatter(z_lester+ 4*inc, alpha_scenario_5_6_n, label='scenario 5 6 Santini M', marker='x')
#plt.scatter(z_lester+ 5*inc, alpha_scenario_5_7_n, label='scenario 5 7', marker='x')
#plt.scatter(z_lester+ 6*inc, alpha_scenario_5_k0_n, label='scenario 5 k0', marker='x')
#plt.scatter(z_lester+ 7*inc, alpha_scenario_5_k_n, label='scenario 5 k', marker='x')

#plt.errorbar(z_lester+ 1*inc, alpha_scenario_1_1_n, yerr=alpha_err_scenario_1_1_n, ls='none')
#plt.errorbar(z_lester+ 2*inc, alpha_scenario_5_1_n, yerr=alpha_err_scenario_5_1_n, ls='none')
#plt.errorbar(z_lester+ 3*inc, alpha_scenario_5_4_n, yerr=alpha_err_scenario_5_4_n, ls='none')
#plt.errorbar(z_lester+ 4*inc, alpha_scenario_5_6_n, yerr=alpha_err_scenario_5_6_n, ls='none')
#plt.errorbar(z_lester+ 5*inc, alpha_scenario_5_7_n, yerr=alpha_err_scenario_5_7_n, ls='none')
#plt.errorbar(z_lester+ 6*inc, alpha_scenario_5_k0_n, yerr=alpha_err_scenario_5_k0_n, ls='none')
#plt.errorbar(z_lester+ 7*inc, alpha_scenario_5_k_n, yerr=alpha_err_scenario_5_k_n, ls='none')

#plt.gca().set_prop_cycle(None)
#plt.scatter((),())
#plt.scatter((),())
plt.scatter(z_lester+ 1*inc, alpha_scenario_1_1_n, label='scenario 1 1 Santini M+SFR')
plt.scatter(z_lester+ 2*inc, alpha_scenario_6_1_n, label='scenario 6 1 Santini M+SFR')
plt.scatter(z_lester+ 3*inc, alpha_scenario_6_4_n, label='scenario 6 4 Santini SFR')
plt.scatter(z_lester+ 4*inc, alpha_scenario_6_6_n, label='scenario 6 6 Santini M')
plt.scatter(z_lester+ 5*inc, alpha_scenario_6_7_n, label='scenario 6 7')
plt.scatter(z_lester+ 6*inc, alpha_scenario_6_k0_n, label='scenario 6 k0')
plt.scatter(z_lester+ 7*inc, alpha_scenario_6_k_n, label='scenario 6 k')
plt.scatter(z_lester+ 8*inc, alpha_scenario_7_k0_d_n, label='scenario 7 k0 d')
plt.scatter(z_lester+ 9*inc, alpha_scenario_7_k0_dd_n, label='scenario 7 k0 dd')
plt.scatter(z_lester+ 10*inc, alpha_scenario_7_k0_dc_n, label='scenario 7 k0 dc')
plt.scatter(z_lester+ 11*inc, alpha_scenario_7_k0_d_Hogg_n, label='scenario 7 k0 d Hogg')
plt.scatter(z_lester+ 12*inc, alpha_scenario_7_k0_dd_Hogg_n, label='scenario 7 k0 dd Hogg')
plt.scatter(z_lester+ 13*inc, alpha_scenario_7_k0_dc_Hogg_n, label='scenario 7 k0 dc Hogg')

plt.errorbar(z_lester+ 1*inc, alpha_scenario_1_1_n, yerr=alpha_err_scenario_1_1_n, ls='none')
plt.errorbar(z_lester+ 2*inc, alpha_scenario_6_1_n, yerr=alpha_err_scenario_6_1_n, ls='none')
plt.errorbar(z_lester+ 3*inc, alpha_scenario_6_4_n, yerr=alpha_err_scenario_6_4_n, ls='none')
plt.errorbar(z_lester+ 4*inc, alpha_scenario_6_6_n, yerr=alpha_err_scenario_6_6_n, ls='none')
plt.errorbar(z_lester+ 5*inc, alpha_scenario_6_7_n, yerr=alpha_err_scenario_6_7_n, ls='none')
plt.errorbar(z_lester+ 6*inc, alpha_scenario_6_k0_n, yerr=alpha_err_scenario_6_k0_n, ls='none')
plt.errorbar(z_lester+ 7*inc, alpha_scenario_6_k_n, yerr=alpha_err_scenario_6_k_n, ls='none')
plt.errorbar(z_lester+ 8*inc, alpha_scenario_7_k0_d_n, yerr=alpha_err_scenario_7_k0_d_n, ls='none')
plt.errorbar(z_lester+ 9*inc, alpha_scenario_7_k0_dd_n, yerr=alpha_err_scenario_7_k0_dd_n, ls='none')
plt.errorbar(z_lester+ 10*inc, alpha_scenario_7_k0_dc_n, yerr=alpha_err_scenario_7_k0_dc_n, ls='none')
plt.errorbar(z_lester+ 11*inc, alpha_scenario_7_k0_d_Hogg_n, yerr=alpha_err_scenario_7_k0_d_Hogg_n, ls='none')
plt.errorbar(z_lester+ 12*inc, alpha_scenario_7_k0_dd_Hogg_n, yerr=alpha_err_scenario_7_k0_dd_Hogg_n, ls='none')
plt.errorbar(z_lester+ 13*inc, alpha_scenario_7_k0_dc_Hogg_n, yerr=alpha_err_scenario_7_k0_dc_Hogg_n, ls='none')

#plt.plot(z_speagle, alpha_speagle_n, label='Speagle+14')
#plt.plot(z_schreiber, alpha_schreiber_n, label='Schreiber+15')
plt.ylim(-2, 3)
#plt.ylim(-11, -2)
plt.legend()
plt.show()




print(alpha_err_scenario_1_1)
print(alpha_err_scenario_1_1[0])
print(alpha_err_scenario_1_1[1])
print((alpha_err_scenario_1_1[1]-alpha_err_scenario_1_1[0])/2)







