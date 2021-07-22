#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: lester
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm,truncnorm,multivariate_normal

def calc_sfr_surface(X, beta_a, beta_b, alphaN_a, alphaN_b):
    x,z = X
    return ((beta_a + beta_b*z)*(x-alphaNorm)) + (alphaN_a + alphaN_b*z) # error if alphaN_a <0

def calc_sigsqr(xi, sig0, k=1.0, xi_min=8.5, xi_max=10.0):
    sigsqr = ( sig0 * ( ((1.0-k)*(xi-xi_max)/(xi_max-xi_min)) + 1.0 ) ) ** 2.0
    return sigsqr

xi = 8.7
zeta = 1.5
beta_a = 1.0
beta_b = 0.0
alphaN_a = 1.6
alphaN_b = 0.0
sig0=0.3
alphaNorm = 9.7

outlier_mean = 0.69
outlier_sigma = 2.0
pbad = 0.2

#update each object in turn
eta_curr = 0.69+0.365
eta_prop = 0.69-0.365

# P(eta|xi,zeta,theta) - this is only really true for the proportionality, I'm not calculating the correct normalization
p_eta_xi_zeta_theta_curr = norm.pdf(eta_curr, scale=np.sqrt(calc_sigsqr(xi, sig0)), loc=calc_sfr_surface((xi, zeta), beta_a, beta_b, alphaN_a, alphaN_b))
p_eta_xi_zeta_theta_prop = norm.pdf(eta_prop, scale=np.sqrt(calc_sigsqr(xi, sig0)), loc=calc_sfr_surface((xi, zeta), beta_a, beta_b, alphaN_a, alphaN_b))

print(np.sqrt(calc_sigsqr(xi, sig0)))
print(calc_sfr_surface((xi, zeta), beta_a, beta_b, alphaN_a, alphaN_b))
print(p_eta_xi_zeta_theta_curr)
print(p_eta_xi_zeta_theta_prop)

# outlier distribution
p_eta_om_os_curr = norm.pdf(eta_curr, loc=outlier_mean, scale=outlier_sigma)
p_eta_om_os_prop = norm.pdf(eta_prop, loc=outlier_mean, scale=outlier_sigma)

print(p_eta_om_os_curr)
print(p_eta_om_os_prop)

# P(eta|xi,zeta,x,y)
# mean, cov = get_3d_mean_cov(i)
mean_curr = np.array([xi, eta_curr, zeta])
mean_prop = np.array([xi, eta_prop, zeta])
cov = np.array([[1,0,0],[0,1,0],[0,0,1]])

log_p_eta_xi_zeta_x_y_curr = multivariate_normal.logpdf([xi,eta_curr,zeta],mean_curr,cov)
log_p_eta_xi_zeta_x_y_prop = multivariate_normal.logpdf([xi,eta_prop,zeta],mean_prop,cov)

print(log_p_eta_xi_zeta_x_y_curr)
print(log_p_eta_xi_zeta_x_y_prop)

log_target_curr = log_p_eta_xi_zeta_x_y_curr + np.log((((1.0-pbad)*p_eta_xi_zeta_theta_curr) + (pbad*p_eta_om_os_curr)))
log_target_prop = log_p_eta_xi_zeta_x_y_prop + np.log((((1.0-pbad)*p_eta_xi_zeta_theta_prop) + (pbad*p_eta_om_os_prop)))

print(log_target_curr)
print(log_target_prop)

ratio1 = np.exp(log_target_prop - log_target_curr)

#working out separately

curr = np.exp(log_p_eta_xi_zeta_x_y_curr) * ((1.0-pbad)*(p_eta_xi_zeta_theta_curr) + (pbad)*(p_eta_om_os_curr))

prop = np.exp(log_p_eta_xi_zeta_x_y_prop) * ((1.0-pbad)*(p_eta_xi_zeta_theta_prop) + (pbad)*(p_eta_om_os_prop))

ratio2 = prop / curr

print(ratio1)
print(ratio2)


#%%






acceptanceProb = min(1,np.exp(log_target_prop - log_target_curr)) # 1 if prop more likely than curr

if acceptanceProb == 0:
    print('eta acceptanceProb 0', ichain, i, log_target_prop, log_target_curr, np.exp(log_target_prop - log_target_curr))

u = rng.uniform() # random between 0 and 1
if u <= acceptanceProb: # accept proposal
    eta[i] = eta_prop[i]
    accept_eta[i] = accept_eta[i] + 1
else:
    reject_eta[i] = reject_eta[i] + 1

test = accept_eta[i]+reject_eta[i]
if ichain >= nBurn and test > 0 and test % testIter == 0:
    if float(accept_eta[i])/float(accept_eta[i]+reject_eta[i]) > 0.5:
        proposalscale_eta[i] = proposalscale_eta[i]*1.1
    elif float(accept_eta[i])/float(accept_eta[i]+reject_eta[i]) < 0.4:
        proposalscale_eta[i] = proposalscale_eta[i]*0.9

    accept_eta[i] = 0
    reject_eta[i] = 0
            
            
            


# =============================================================================
# mass
# =============================================================================


'''
xi_curr = 9.7                                                             # current xi values
xi_prop = self.rng.normal(size=len(xi_curr),scale=self.proposalscale_xi)+xi_curr    # new proposal for each xi
muArr = []                                                                      # mu are means of xi gaussians
tausqrArr = []                                                                  # tausqr are sigsqr of xi gaussians
for i in range(len(self.xi)):
    tempIdx = np.where(self.G[i] == 1)[0]                                       # [0], [1] or [2] decides which xi gaussian
    muArr.append(self.mu[tempIdx][0])                                           # chooses mu from [mu1, mu2, mu3]
    tausqrArr.append(self.tausqr[tempIdx][0])                                   # chooses tau from [ta1, ta2, ta3]
# mass distribution
log_p_xi_psi_curr = norm.logpdf(xi_curr,loc=muArr,scale=np.sqrt(tausqrArr))     # gives log( gaussian value at xi_curr ) per object
log_p_xi_psi_prop = norm.logpdf(xi_prop,loc=muArr,scale=np.sqrt(tausqrArr))     # gives log( gaussian value at xi_prop ) per object
# P(xi|eta,zeta,theta) - this is only really true for the proportionality, I'm not calculating the correct normalization
# MS relation
p_xi_eta_zeta_theta_curr = norm.pdf(self.eta, scale=np.sqrt(self.calc_sigsqr(xi_curr)), loc=self.calc_sfr_surface((xi_curr,self.zeta), self.beta_a, self.beta_b, self.alphaN_a, self.alphaN_b))
p_xi_eta_zeta_theta_prop = norm.pdf(self.eta, scale=np.sqrt(self.calc_sigsqr(xi_prop)), loc=self.calc_sfr_surface((xi_prop,self.zeta), self.beta_a, self.beta_b, self.alphaN_a, self.alphaN_b))
# outlier distribution
p_xi_om_os_curr = norm.pdf(self.eta, loc=self.outlier_mean, scale=self.outlier_sigma)
p_xi_om_os_prop = norm.pdf(self.eta, loc=self.outlier_mean, scale=self.outlier_sigma)   
for i in range(len(self.xi)):
#                P(xi|eta,zeta,x,y)
    mean, cov = self.get_3d_mean_cov(i)
    # GMM measurement uncertainties
    log_p_xi_eta_zeta_x_y_curr = multivariate_normal.logpdf([xi_curr[i],self.eta[i],self.zeta[i]],mean,cov)
    log_p_xi_eta_zeta_x_y_prop = multivariate_normal.logpdf([xi_prop[i],self.eta[i],self.zeta[i]],mean,cov)
    log_target_curr = log_p_xi_eta_zeta_x_y_curr + log_p_xi_psi_curr[i] + np.log((((1.0-self.pbad)*p_xi_eta_zeta_theta_curr[i]) + (self.pbad*p_xi_om_os_curr[i])))
    log_target_prop = log_p_xi_eta_zeta_x_y_prop + log_p_xi_psi_prop[i] + np.log((((1.0-self.pbad)*p_xi_eta_zeta_theta_prop[i]) + (self.pbad*p_xi_om_os_prop[i])))
    acceptanceProb = min(1,np.exp(log_target_prop - log_target_curr)) # 1 if prop more likely than curr
    if acceptanceProb == 0:
        print('xi acceptanceProb 0', self.ichain, i, log_target_prop, log_target_curr, np.exp(log_target_prop - log_target_curr))
    u = self.rng.uniform() # random between 0 and 1
    if u <= acceptanceProb: # accept proposal
        self.xi[i] = xi_prop[i]
        self.accept_xi[i] = self.accept_xi[i] + 1
    else:
        self.reject_xi[i] = self.reject_xi[i] + 1
    test = self.accept_xi[i]+self.reject_xi[i]
    if self.ichain >= self.nBurn and test > 0 and test % self.testIter == 0:
        if float(self.accept_xi[i])/float(self.accept_xi[i]+self.reject_xi[i]) > 0.5:
            self.proposalscale_xi[i] = self.proposalscale_xi[i]*1.1
        elif float(self.accept_xi[i])/float(self.accept_xi[i]+self.reject_xi[i]) < 0.4:
            self.proposalscale_xi[i] = self.proposalscale_xi[i]*0.9
        self.accept_xi[i] = 0
        self.reject_xi[i] = 0

'''







'''
# =============================================================================
# the model (HOGG + redshift dependent alpha, NO mass dependent scatter (k=1))
# values taken from full scenario 29 fit (ssfr alpha z1.25 - 6.0)
# =============================================================================

n = 100000
mass = np.random.normal(8.5, 0.5, n)
z = np.random.uniform(1.25, 2.0, n)

alpha_a = 0.9
alpha_b = 0.0

beta_a = 1.0
beta_b = 0.0

alpha = alpha_a + alpha_b*z
beta = beta_a + beta_b*z

sig0 = 0.25

pbad = 0.2
nBad = np.random.binomial(n,pbad)
outlier_mean = -2.0
outlier_sigma = 1.2

tempIdx = np.random.choice(np.fromiter(range(n),np.int),size=nBad,replace=False)
good = np.ones(n,np.bool)
good[tempIdx] = False

sfr = np.zeros_like(mass)
sfr[good] = alpha[good] + (mass[good]-9.7)*beta[good] + np.random.normal(0, sig0, sum(good))
sfr[~good] = np.random.normal(outlier_mean, outlier_sigma, nBad)


plt.hist(mass)
plt.show()

plt.hist(sfr)
plt.show()

plt.hist(z)
plt.show()

# =============================================================================
# plot
# =============================================================================
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(mass, z, sfr, c=z)
plt.show()


plt.scatter(mass, sfr, c=z)
plt.colorbar()
plt.show()

# =============================================================================
# kelly input
# =============================================================================
x_GMM_3d = np.array([mass]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))
y_GMM_3d = np.array([sfr]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))
z_GMM_3d = np.array([z]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))

xsig_GMM_3d = np.random.normal(0.15, 0.05, (n,3))
xsig_GMM_3d = np.where(xsig_GMM_3d<0.01, 0.01, xsig_GMM_3d)
ysig_GMM_3d = np.random.normal(0.2, 0.05, (n,3))
ysig_GMM_3d = np.where(ysig_GMM_3d<0.01, 0.01, ysig_GMM_3d)
zsig_GMM_3d = np.random.normal(0.07, 0.05, (n,3))
zsig_GMM_3d = np.where(zsig_GMM_3d<0.01, 0.01, zsig_GMM_3d)

xycov_GMM_3d = np.random.normal(-0.3, 0.4, (n,3))
xycov_GMM_3d = np.where(np.abs(xycov_GMM_3d)>1.0, -0.3, xycov_GMM_3d) * xsig_GMM_3d * ysig_GMM_3d
xzcov_GMM_3d = np.random.normal(0.3, 0.4, (n,3))
xzcov_GMM_3d = np.where(np.abs(xzcov_GMM_3d)>1.0, 0.3, xzcov_GMM_3d) * xsig_GMM_3d * zsig_GMM_3d
yzcov_GMM_3d = np.random.normal(0.3, 0.4, (n,3))
yzcov_GMM_3d = np.where(np.abs(yzcov_GMM_3d)>1.0, 0.3, yzcov_GMM_3d) * ysig_GMM_3d * zsig_GMM_3d

amp_GMM_3d = np.random.uniform(0, 1, (n,3))
amp_GMM_3d_norm = np.sum(amp_GMM_3d, axis=1)
amp_GMM_3d = amp_GMM_3d/amp_GMM_3d_norm[:, np.newaxis]


# =============================================================================
# sorting out the eigenvalue issue
# =============================================================================
count = 0
for i in range(len(x_GMM_3d)):
    for j in np.arange(3):
        
        cov = [[xsig_GMM_3d[i][j]**2,xycov_GMM_3d[i][j],xzcov_GMM_3d[i][j]],\
                  [xycov_GMM_3d[i][j],ysig_GMM_3d[i][j]**2,yzcov_GMM_3d[i][j]],\
                  [xzcov_GMM_3d[i][j],yzcov_GMM_3d[i][j],zsig_GMM_3d[i][j]**2]]
            
        while not np.all(np.linalg.eigvals(cov) > 0):
            count+=1
            print(count)
            
            xycov_GMM_3d[i][j] = np.random.normal(-0.3, 0.4)
            xycov_GMM_3d[i][j] = np.where(np.abs(xycov_GMM_3d[i][j])>1.0, -0.3, xycov_GMM_3d[i][j]) * xsig_GMM_3d[i][j] * ysig_GMM_3d[i][j]
            xzcov_GMM_3d[i][j] = np.random.normal(0.3, 0.4)
            xzcov_GMM_3d[i][j] = np.where(np.abs(xzcov_GMM_3d[i][j])>1.0, 0.3, xzcov_GMM_3d[i][j]) * xsig_GMM_3d[i][j] * zsig_GMM_3d[i][j]
            yzcov_GMM_3d[i][j] = np.random.normal(0.3, 0.4)
            yzcov_GMM_3d[i][j] = np.where(np.abs(yzcov_GMM_3d[i][j])>1.0, 0.3, yzcov_GMM_3d[i][j]) * ysig_GMM_3d[i][j] * zsig_GMM_3d[i][j]                
            
            cov = [[xsig_GMM_3d[i][j]**2,xycov_GMM_3d[i][j],xzcov_GMM_3d[i][j]],\
                  [xycov_GMM_3d[i][j],ysig_GMM_3d[i][j]**2,yzcov_GMM_3d[i][j]],\
                  [xzcov_GMM_3d[i][j],yzcov_GMM_3d[i][j],zsig_GMM_3d[i][j]**2]]

# =============================================================================
# perturbing means based on covariance matrices
# =============================================================================
for i in range(len(x_GMM_3d)):
    for j in np.arange(3):
        
        cov = [[xsig_GMM_3d[i][j]**2,xycov_GMM_3d[i][j],xzcov_GMM_3d[i][j]],\
                  [xycov_GMM_3d[i][j],ysig_GMM_3d[i][j]**2,yzcov_GMM_3d[i][j]],\
                  [xzcov_GMM_3d[i][j],yzcov_GMM_3d[i][j],zsig_GMM_3d[i][j]**2]]
        mean = [x_GMM_3d[i,j],y_GMM_3d[i,j],z_GMM_3d[i,j]]
        
        xyz = np.random.multivariate_normal(mean, cov)
        
        x_GMM_3d[i][j] = xyz[0]
        y_GMM_3d[i][j] = xyz[1]
        z_GMM_3d[i][j] = xyz[2]


plt.scatter(x_GMM_3d.flatten(), y_GMM_3d.flatten(), c=z_GMM_3d.flatten())
plt.colorbar()
plt.show()


data = {'x_GMM_3d':x_GMM_3d, 'y_GMM_3d':y_GMM_3d, 'z_GMM_3d':z_GMM_3d, 'xsig_GMM_3d':xsig_GMM_3d, 'ysig_GMM_3d':ysig_GMM_3d, 'zsig_GMM_3d':zsig_GMM_3d, 'xycov_GMM_3d':xycov_GMM_3d, 'xzcov_GMM_3d':xzcov_GMM_3d, 'yzcov_GMM_3d':yzcov_GMM_3d, 'amp_GMM_3d':amp_GMM_3d}

#pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/mock_z1_001.p','w'))

#WINDOWS
# pickle.dump(data, open('/Users/LSand/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis_windows/kelly_input/mock_z1_0{}.p'.format(str(m+11)),'wb'))


print(x_GMM_3d[0])
print(y_GMM_3d[0])
print(z_GMM_3d[0])
print(amp_GMM_3d[0])
print(nBad)

#test = np.polyfit(mass, sfr, 1)
#print(test)


'''

