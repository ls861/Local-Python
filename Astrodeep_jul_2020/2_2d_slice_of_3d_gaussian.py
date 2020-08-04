#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 19:44:13 2020

@author: lester
"""


from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

'''
This finds a 2d slice of a 3d multivariate gaussian, given fixed z
'''


xmean = 2.0
ymean = 30.0
zmean = 4000.0

xvar = 1.0
yvar = 1.0
zvar = 100.0

xycov = 0.
xzcov = 0.
yzcov = -0.9

mean_3d = np.array([xmean, ymean, zmean])
cov_3d = np.array([[xvar, xycov, xzcov],[xycov, yvar, yzcov],[xzcov, yzcov, zvar]])

mvn_3d = np.random.multivariate_normal(mean_3d, cov_3d, 5000)



ax = plt.axes(projection='3d')
ax.scatter(mvn_3d[:,0], mvn_3d[:,1], mvn_3d[:,2])
plt.show()



a = np.array([3990]) # fixed z

mean_2d = np.array([xmean,ymean]).T + ( np.array([xzcov,yzcov]).T * np.array([1/zvar]) * (a - np.array([zmean])) )

cov_2d = np.array([[xvar,xycov],[xycov,yvar]]) - ( np.array([1/zvar]) * np.outer(np.array([xzcov,yzcov]).T, np.array([xzcov,yzcov]) ) )

mvn_2d = np.random.multivariate_normal(mean_2d, cov_2d, 500)


a_full = np.full(len(mvn_2d[:,0]), a)

plt.scatter(mvn_2d[:,0], mvn_2d[:,1])
plt.show()

ax = plt.axes(projection='3d')
ax.scatter(mvn_2d[:,0], mvn_2d[:,1], a_full)
plt.show()



# =============================================================================
# again
# =============================================================================

a1s = np.array([3850, 3900, 3950, 4000, 4050, 4100, 4150]) # fixed z
ax = plt.axes(projection='3d')

for a1 in a1s:
    print(a1)
    
    

    mean_2d1 = np.array([xmean,ymean]).T + ( np.array([xzcov,yzcov]).T * np.array([1/zvar]) * (a1 - np.array([zmean])) )

    cov_2d1 = np.array([[xvar,xycov],[xycov,yvar]]) - ( np.array([1/zvar]) * np.outer(np.array([xzcov,yzcov]).T, np.array([xzcov,yzcov]) ) )
    
    mvn_2d1 = np.random.multivariate_normal(mean_2d1, cov_2d1, 500)

    
    a_full1 = np.full(len(mvn_2d1[:,0]), a1)
    
    ax.scatter(mvn_2d1[:,0], mvn_2d1[:,1], a_full1)
plt.show()







