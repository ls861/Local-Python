#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 10:04:18 2021

@author: lester
"""



from scipy.stats import norm,truncnorm,multivariate_normal
import numpy as np
import matplotlib.pyplot as plt
'''
xi = np.linspace(0, 20, 1000)
eta = 70.0
zeta = 0.0

m_x = 10.0
m_y = 100.0
m_z = 0.0

xsig = 1.0
ysig = 10.0
zsig = 10.0

xvar = xsig**2
yvar = ysig**2
zvar = zsig**2

xycov = 0.8*xsig*ysig
xzcov = 0.0*xsig*zsig
yzcov = 0.0*ysig*zsig

mean = np.array([m_x, m_y, m_z])
cov = np.matrix([[xvar, xycov, xzcov],[xycov, yvar, yzcov],[xzcov, yzcov, zvar]])

Y = np.matrix([eta, zeta]).T
m_X = np.matrix([m_x]).T
m_Y = np.matrix([m_y, m_z]).T
E_XX = np.matrix([xvar])
E_XY = np.matrix([xycov, xzcov])
E_YX = E_XY.T
E_YY = np.matrix([[yvar, yzcov],[yzcov, zvar]])

m_X_Y = m_X + E_XY*np.linalg.inv(E_YY)*(Y-m_Y)
E_XX_Y = E_XX - E_XY*np.linalg.inv(E_YY)*E_YX



# =============================================================================
# normalisation
# =============================================================================
xi_test = 7.0
pdf1d = norm.pdf(xi_test, m_X_Y[0,0], np.sqrt(E_XX_Y[0,0]))
pdf3d = multivariate_normal.pdf(np.array([xi_test, eta, zeta]), mean, cov)
print(pdf1d, pdf3d)

# =============================================================================
# 1d vaues from 3d
# =============================================================================
pdf = []
for i in range(len(xi)):
    pdf.append(multivariate_normal.pdf(np.array([xi[i], eta, zeta]), mean, cov))
pdf = np.array(pdf)
plt.plot(xi, pdf*(pdf1d/pdf3d), label='3d')


# =============================================================================
# 1d values from 1d
# =============================================================================
plt.plot(xi, norm.pdf(xi, m_X_Y[0,0], np.sqrt(E_XX_Y[0,0])), linestyle='dashed', label='1d')
plt.legend()
plt.show()


print(pdf3d/pdf1d)


# =============================================================================
# sampling from the 1d gaussian
# =============================================================================
'''
#%%
# =============================================================================
# should make version for 3 clouds to test...
# =============================================================================

def GMM_normalisation((eta, zeta), (m_x, m_y, m_z), (xsig, ysig, zsig), (xycov, xzcov, yzcov)):
    
    xvar = xsig**2
    yvar = ysig**2
    zvar = zsig**2
    
    xycov = xycov*xsig*ysig
    xzcov = xzcov*xsig*zsig
    yzcov = yzcov*ysig*zsig
    
    mean = np.array([m_x, m_y, m_z])
    cov = np.matrix([[xvar, xycov, xzcov],[xycov, yvar, yzcov],[xzcov, yzcov, zvar]])
    
    print(mean)
    print(cov)
    
    Y = np.matrix([eta, zeta]).T
    m_X = np.matrix([m_x]).T
    m_Y = np.matrix([m_y, m_z]).T
    E_XX = np.matrix([xvar])
    E_XY = np.matrix([xycov, xzcov])
    E_YX = E_XY.T
    E_YY = np.matrix([[yvar, yzcov],[yzcov, zvar]])
    
    m_X_Y = m_X + E_XY*np.linalg.inv(E_YY)*(Y-m_Y)
    E_XX_Y = E_XX - E_XY*np.linalg.inv(E_YY)*E_YX
    
    pdf1d = norm.pdf(m_X_Y[0,0], m_X_Y[0,0], np.sqrt(E_XX_Y[0,0]))
    pdf3d = multivariate_normal.pdf(np.array([m_X_Y[0,0], eta, zeta]), mean, cov)
    
    return pdf3d/pdf1d, mean, cov, m_X_Y[0,0], E_XX_Y[0,0]


#%%
    

# =============================================================================
# get normalisation, mean, cov and conditionals 
# =============================================================================
    
eta = 15.0
zeta = 0.0

n1, mean1, cov1, m_X_Y1, E_XX_Y1 = GMM_normalisation((eta, zeta), (10.0, 20.0, 0.0), (3.0, 3.0, 3.0), (0.9, 0.0, 0.0))
n2, mean2, cov2, m_X_Y2, E_XX_Y2 = GMM_normalisation((eta, zeta), (30.0, 20.0, 0.0), (3.0, 3.0, 3.0), (-0.9, 0.0, 0.0))
n3, mean3, cov3, m_X_Y3, E_XX_Y3 = GMM_normalisation((eta, zeta), (20.0, 10.0, 0.0), (3.0, 3.0, 3.0), (0.0, 0.0, 0.0))
N1 = n1 / np.sum((n1, n2, n3))
N2 = n2 / np.sum((n1, n2, n3))
N3 = n3 / np.sum((n1, n2, n3))

# =============================================================================
# 2d histogram of arbitrary 3x 3D gaussians
# =============================================================================

t1 = np.random.multivariate_normal(mean1, cov1, size=1000000)
t2 = np.random.multivariate_normal(mean2, cov2, size=1000000)
t3 = np.random.multivariate_normal(mean3, cov3, size=1000000)
t = np.concatenate((t1, t2, t3))

plt.hist2d(t[:,0], t[:,1], bins=50)
plt.plot((min(t[:,0]), max(t[:,0])), (5, 5), color='r')
plt.plot((min(t[:,0]), max(t[:,0])), (10, 10), color='r')
plt.plot((min(t[:,0]), max(t[:,0])), (15, 15), color='r')
plt.plot((min(t[:,0]), max(t[:,0])), (20, 20), color='r')
plt.plot((min(t[:,0]), max(t[:,0])), (25, 25), color='r')
plt.show()

# =============================================================================
# TRUE sample from 3x 3D gaussians within conditional bounds for histogram
# =============================================================================

idx = (abs(t[:,1] - eta) < 0.5) & (abs(t[:,2] - zeta) < 0.5)
plt.hist(t[:,0][idx], bins=50, density=True, label='true samples') # samples from GMMs

# =============================================================================
# CALCULATED sample from 1D conditionals for histogram
# =============================================================================

h1 = np.random.normal(loc=m_X_Y1, scale=np.sqrt(E_XX_Y1), size=int(30000*N1))
h2 = np.random.normal(loc=m_X_Y2, scale=np.sqrt(E_XX_Y2), size=int(30000*N2))
h3 = np.random.normal(loc=m_X_Y3, scale=np.sqrt(E_XX_Y3), size=int(30000*N3))
h = np.concatenate((h1, h2, h3))
plt.hist(h, bins=50, density=True, histtype='step', linewidth=2, label='conditional samples') # samples from adjusted 1D conditionals

# =============================================================================
# CALCULATED gaussians from 1D conditionals
# =============================================================================

x_coord1 = np.linspace(min(h), max(h), 1000)

y_coord1 = norm.pdf(x_coord1, loc=m_X_Y1, scale=np.sqrt(E_XX_Y1))*N1
y_coord2 = norm.pdf(x_coord1, loc=m_X_Y2, scale=np.sqrt(E_XX_Y2))*N2
y_coord3 = norm.pdf(x_coord1, loc=m_X_Y3, scale=np.sqrt(E_XX_Y3))*N3

x_coord = np.concatenate((x_coord1, x_coord1, x_coord1))
y_coord = np.concatenate((y_coord1, y_coord2, y_coord3))

plt.plot(x_coord, y_coord, color='k', label='single gaussians')
plt.plot(x_coord1, y_coord1+y_coord2+y_coord3, color='r', linewidth=2, label='conditional distribution')
plt.legend()
plt.show()










#%%
# =============================================================================
# 2d sanity check as something isn't right - turned out I was using E instead of sqrt E when sampling!!!
# =============================================================================
'''
def GMM_normalisation_2d((xi, eta), (m_x, m_y), (xsig, ysig), (xycov)):
    
    xvar = xsig**2
    yvar = ysig**2
#    zvar = zsig**2
    
    xycov = xycov*xsig*ysig
#    xzcov = xzcov*xsig*zsig
#    yzcov = yzcov*ysig*zsig
    
    mean = np.array([m_x, m_y])
    cov = np.matrix([[xvar, xycov],[xycov, yvar]])
    
    print(mean)
    print(cov)
    
    Y = np.matrix([eta]).T
    m_X = np.matrix([m_x]).T
    m_Y = np.matrix([m_y]).T
    E_XX = np.matrix([xvar])
    E_XY = np.matrix([xycov])
    E_YX = E_XY.T
    E_YY = np.matrix([yvar])
    
    m_X_Y = m_X + E_XY*np.linalg.inv(E_YY)*(Y-m_Y)
    E_XX_Y = E_XX - E_XY*np.linalg.inv(E_YY)*E_YX
    
    
    print(E_XX)
    print(E_XY)
    print(np.linalg.inv(E_YY))
    print(E_YX)
    print(E_XY*np.linalg.inv(E_YY)*E_YX)
    
#    print(Y)
#    print(m_X)
#    print(m_Y)
#    print(E_XX)
#    print(E_XY)
#    print(E_YX)
#    print(E_YY)
#    print(m_X_Y)
#    print(E_XX_Y)
    
    pdf1d = norm.pdf(xi, m_X_Y[0,0], np.sqrt(E_XX_Y[0,0]))
    pdf3d = multivariate_normal.pdf(np.array([xi, eta]), mean, cov)
    
    return pdf3d/pdf1d, mean, cov, m_X_Y[0,0], E_XX_Y[0,0]

#n = GMM_normalisation((0.0, 70.0, 0.0), (10.0, 100.0, 0.0), (1.0, 10.0, 10.0), (0.8, 0.0, 0.0))
#print(n)
n4, mean4, cov4, m_X_Y4, E_XX_Y4 = GMM_normalisation_2d((20.0, 15.0), (20.0, 10.0), (3.0, 3.0), (0.0))


t4 = np.random.multivariate_normal(mean4, cov4, size=10000000)


plt.hist2d(t4[:,0], t4[:,1], bins=50)
plt.show()

print(mean4, cov4, m_X_Y4, E_XX_Y4)

plt.hist(t4[:,0], bins=50, histtype='step')
plt.hist(t4[:,0][abs(t4[:,1]-15.0)<0.5], bins=50, histtype='step')
plt.show()

plt.hist(t4[:,0], bins=50, density=True, histtype='step')
plt.hist(t4[:,0][abs(t4[:,1]-15.0)<0.5], bins=50, density=True, histtype='step')
plt.show()


print(np.var(t4[:,0]))
print(np.var(t4[:,0][abs(t4[:,1]-15.0)<0.5]))


'''








