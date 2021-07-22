#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 16:28:49 2021

@author: lester
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm,truncnorm,multivariate_normal
import time
from astropy.io import fits
import copy


#%%
# =============================================================================
# real stuff
# =============================================================================

def pdf(s, idx, G, n, dintegration = 1E-3 ):

    # =============================================================================
    # gmm pdfs
    # =============================================================================
    n_gmm = 1.0
    
    #The GMM parameters are separate from the grid
    pi_gmm = s['amp_GMM_3d'][idx,G]
    mean = np.array([s['x_GMM_3d'][idx,G],s['y_GMM_3d'][idx,G],s['z_GMM_3d'][idx,G]])
    cov = np.array([[np.power(s['xsig_GMM_3d'][idx,G],2), s['xycov_GMM_3d'][idx,G], s['xzcov_GMM_3d'][idx,G]],[s['xycov_GMM_3d'][idx,G], np.power(s['ysig_GMM_3d'][idx,G],2), s['yzcov_GMM_3d'][idx,G]],[s['xzcov_GMM_3d'][idx,G], s['yzcov_GMM_3d'][idx,G], np.power(s['zsig_GMM_3d'][idx,G],2)]])

    # =============================================================================
    # setting up the optimal coordinate grid for this Gaussian
    #
    # Doing this a simple way for now - start with n x n x n grid points
    # =============================================================================
    #Why not call them mass, sfr and redshift?
    #Because we need to do the integral, I'm going to actually evaluate the GMM
    #at the mid-point of the initial grid points
    x = np.linspace(0, 20, n)
    dx = (x[1]-x[0])/2.
    xMid = x[:-1]+dx
    
    y = np.linspace(-10, 10, n)
    dy = (y[1]-y[0])/2.
    yMid = y[:-1]+dy
    
    z = np.linspace(0, 10, n)
    dz = (z[1]-z[0])/2.
    zMid = z[:-1]+dz
    
    midPointsMesh = np.meshgrid(xMid,yMid,zMid)
    nPts = len(xMid)*len(yMid)*len(zMid)
    midPoints = np.array([midPointsMesh[0].reshape(nPts),\
                          midPointsMesh[1].reshape(nPts),\
                          midPointsMesh[2].reshape(nPts)]).transpose()
                          
    #The next one is to keep track of the volume
    dGridMesh = np.meshgrid(x[1:]-x[:-1],y[1:]-y[:-1],z[1:]-z[:-1])
    dGrid = np.array([dGridMesh[0].reshape(nPts),\
                      dGridMesh[1].reshape(nPts),\
                      dGridMesh[2].reshape(nPts)]).transpose()
                      
    #And let's keep track of the grid points which will become irregular!
    gridPointsLowerMesh = np.meshgrid(x[:-1],y[:-1],z[:-1])
    gridPointsLower = np.array([gridPointsLowerMesh[0].reshape(nPts),\
                                gridPointsLowerMesh[1].reshape(nPts),\
                                gridPointsLowerMesh[2].reshape(nPts)]).transpose()
    gridPointsUpperMesh = np.meshgrid(x[1:],y[1:],z[1:])
    gridPointsUpper = np.array([gridPointsUpperMesh[0].reshape(nPts),\
                                gridPointsUpperMesh[1].reshape(nPts),\
                                gridPointsUpperMesh[2].reshape(nPts)]).transpose()
    
    #And let's track which points to sub-divide or not; we start off with
    #dividing all the points!
    divideGrid = np.ones_like(midPoints[:,0],np.bool)
    gridMarker = np.fromiter((x for x in range(len(divideGrid))), np.int)
    
    adjustGrid = True
    counter = 1
    while (adjustGrid == True) and (counter < 1000):
    
      #Implementing in an initially simple way to figure out how to do it!
      newGridPointsLower = []
      newGridPointsUpper = []
      newdGrid = []
      newGridMarker = []
      for i in range(len(divideGrid)):
        if divideGrid[i] == True:
          #We're essentially sub-dividing a cube, which gives 8 (2^3) points for the original point
          newSubX = [gridPointsLower[i,0],gridPointsLower[i,0]+dGrid[i,0]/2.]
          newSubY = [gridPointsLower[i,1],gridPointsLower[i,1]+dGrid[i,1]/2.]
          newSubZ = [gridPointsLower[i,2],gridPointsLower[i,2]+dGrid[i,2]/2.]
          newLowerMesh = np.meshgrid(newSubX,newSubY,newSubZ)
          newUpperMesh = np.meshgrid(newSubX+dGrid[i,0]/2.,newSubY+dGrid[i,1]/2,\
                                     newSubZ+dGrid[i,2]/2.)
          newLower = np.array([newLowerMesh[0].reshape(8), newLowerMesh[1].reshape(8), \
                      newLowerMesh[2].reshape(8)]).transpose()
          newUpper = np.array([newUpperMesh[0].reshape(8), newUpperMesh[1].reshape(8), \
                      newUpperMesh[2].reshape(8)]).transpose()
          #OK, I know this bit's horrible, sorry!
          for j in range(8):
            newGridPointsLower.append(newLower[j,:])
            newGridPointsUpper.append(newUpper[j,:])
            newdGrid.append(dGrid[i,:]/2.)
            newGridMarker.append(gridMarker[i])
        else:
          newGridPointsLower.append(gridPointsLower[i,:])
          newGridPointsUpper.append(gridPointsUpper[i,:])
          newdGrid.append(dGrid[i,:])
          newGridMarker.append(gridMarker[i])
          
      newGridPointsLower = np.array(newGridPointsLower)
      newGridPointsUpper = np.array(newGridPointsUpper)
      newdGrid = np.array(newdGrid)
      newGridMarker = np.array(newGridMarker)
      newMidPoints = newGridPointsLower + newdGrid/2.
          
                             
      pdf_gmm_orig = multivariate_normal.pdf(midPoints, mean, cov)
      integration_orig = pdf_gmm_orig
#      print('pdf_gmm_orig: ', pdf_gmm_orig)
      for i in range(3):
        integration_orig = integration_orig*dGrid[:,i]
      
      pdf_gmm_new = multivariate_normal.pdf(newMidPoints, mean, cov)
      integration_new = pdf_gmm_new
      for i in range(3):
        integration_new = integration_new*newdGrid[:,i]
        
      #Then you have to determine which integrations of the individual grid points actually changed....
      maxDiff = -1.E5
      keepIdx_new = np.ones_like(newGridMarker,np.bool)
      keepIdx_orig = np.zeros_like(gridMarker,np.bool)
      for i in range(np.max(gridMarker)+1):
        tempIdx1 = np.where(gridMarker == gridMarker[i])[0]
        int1 = np.sum(integration_orig[tempIdx1])
        tempIdx2 = np.where(newGridMarker == gridMarker[i])[0]
        int2 = np.sum(integration_new[tempIdx2])
#        print(int1, int2, np.abs(int1-int2), dintegration)
        if np.abs(int1-int2) < dintegration: #don't accept the new sub-divided grid
          keepIdx_new[tempIdx2] = False
          keepIdx_orig[tempIdx1] = True
        else:
#          print(np.abs(int1-int2))
          if np.abs(int1-int2) > maxDiff:
            maxDiff = np.abs(int1-int2)
#      print('maxDiff: ', maxDiff)
#      print(keepIdx_orig)
        
      #Form the new grid from this information, again, this is super horrible
      adjustedGridLower = []
      adjustedGridUpper = []
      adjusteddGrid = []
      adjustedDivide = []
      for i in range(len(keepIdx_orig)):
        if keepIdx_orig[i] == True:
          adjustedGridLower.append(gridPointsLower[i,:])
          adjustedGridUpper.append(gridPointsUpper[i,:])
          adjusteddGrid.append(dGrid[i])
          adjustedDivide.append(False)
        else:
          tempIdx = np.where(newGridMarker == gridMarker[i])[0]
          for j in range(len(tempIdx)):
            adjustedGridLower.append(newGridPointsLower[tempIdx[j],:])
            adjustedGridUpper.append(newGridPointsUpper[tempIdx[j],:])
            adjusteddGrid.append(newdGrid[tempIdx[j]])
            adjustedDivide.append(True) #Assume we haven't converged here
          
      gridPointsLower = np.array(adjustedGridLower)
      gridPointsUpper = np.array(adjustedGridUpper)
      divideGrid = np.array(adjustedDivide)
      gridMarker = np.fromiter((x for x in range(len(divideGrid))), np.int)
      dGrid = np.array(adjusteddGrid)
      midPoints = np.array(gridPointsLower+dGrid/2.)
      
      tempIdx = np.where(adjustedDivide)[0]
      if len(tempIdx) == 0:
        adjustGrid = False
        
      counter = counter+1.
      
    pdf_gmm = multivariate_normal.pdf(midPoints, mean, cov)
    return(pdf,gridPointsLower,gridPointsUpper,midPoints,dGrid)
    
    
#    integration = pdf_gmm
#    for i in range(3):
#      integration = integration*dcoords[:,i]
#
#    #sub-divide each grid point identified to refine
#    #First define a logical array that identifies which grid points to subdivide
#    divide = np.ones([len(xMid)*len(yMid)*len(zMid),3], np.boolean)
#    print(divide.shape, dcoords.shape)
#    print(divide[:10,0])
#
#    #Setting this up the simple (probably inefficient way)
#    xNew = []
#    yNew = []
#    zNew = []
#    for i
##    print('pdf_gmm', sum(pdf_gmm))
#
#    # =============================================================================
#    # ms pdfs
#    # =============================================================================
#    n_ms = 1.0
#    alpha = np.full(np.shape(coords[:,0]), 0.0)
#    beta = np.full(np.shape(coords[:,0]), 0.0)
#    sigma = np.full(np.shape(coords[:,0]), 0.0)
#
#    # ms per z bin
#    alpha = np.where((coords[:,2]>-1.0)&(coords[:,2]<2.0), 0.98, alpha)
#    beta = np.where((coords[:,2]>-1.0)&(coords[:,2]<2.0), 0.86, beta)
#    sigma = np.where((coords[:,2]>-1.0)&(coords[:,2]<2.0), 0.31, sigma)
#
#    alpha = np.where((coords[:,2]>2.0)&(coords[:,2]<3.0), 1.02, alpha)
#    beta = np.where((coords[:,2]>2.0)&(coords[:,2]<3.0), 0.89, beta)
#    sigma = np.where((coords[:,2]>2.0)&(coords[:,2]<3.0), 0.21, sigma)
#
#    alpha = np.where((coords[:,2]>3.0)&(coords[:,2]<4.0), 1.10, alpha)
#    beta = np.where((coords[:,2]>3.0)&(coords[:,2]<4.0), 0.72, beta)
#    sigma = np.where((coords[:,2]>3.0)&(coords[:,2]<4.0), 0.22, sigma)
#
#    alpha = np.where((coords[:,2]>3.0)&(coords[:,2]<5.0), 1.68, alpha)
#    beta = np.where((coords[:,2]>3.0)&(coords[:,2]<5.0), 0.93, beta)
#    sigma = np.where((coords[:,2]>3.0)&(coords[:,2]<5.0), 0.33, sigma)
#
#    alpha = np.where((coords[:,2]>4.0)&(coords[:,2]<99.0), 2.67, alpha)
#    beta = np.where((coords[:,2]>4.0)&(coords[:,2]<99.0), 1.54, beta)
#    sigma = np.where((coords[:,2]>4.0)&(coords[:,2]<99.0), 0.46, sigma)
#
#    pdf_ms = norm.pdf(coords[:,1], beta*(coords[:,0]-9.7) + alpha, sigma)
##    print('pdf_ms', sum(pdf_ms))
#
#    # =============================================================================
#    # mass pdfs
#    # =============================================================================
#    n_m = 1.0
#    mass_mu = np.full((len(alpha),3), 0.0)
#    mass_sigma = np.full((len(alpha),3), 0.0)
#    mass_pi = np.full((len(alpha),3), 0.0)
#
#    # mass per z bin
#    mass_mu[:,0] = np.where(coords[:,2]<99, 8.0, mass_mu[:,0])
#    mass_mu[:,1] = np.where(coords[:,2]<99, 8.0, mass_mu[:,1])
#    mass_mu[:,2] = np.where(coords[:,2]<99, 8.0, mass_mu[:,2])
#    mass_sigma[:,0] = np.where(coords[:,2]<99, 0.5, mass_sigma[:,0])
#    mass_sigma[:,1] = np.where(coords[:,2]<99, 0.5, mass_sigma[:,1])
#    mass_sigma[:,2] = np.where(coords[:,2]<99, 0.5, mass_sigma[:,2])
#    mass_pi[:,0] = np.where(coords[:,2]<99, 0.3, mass_pi[:,0])
#    mass_pi[:,1] = np.where(coords[:,2]<99, 0.3, mass_pi[:,1])
#    mass_pi[:,2] = np.where(coords[:,2]<99, 0.3, mass_pi[:,2])
#
#    pdf_m = mass_pi[:,0] * norm.pdf(coords[:,0], mass_mu[:,0], mass_sigma[:,0]) + \
#            mass_pi[:,1] * norm.pdf(coords[:,0], mass_mu[:,1], mass_sigma[:,1]) + \
#            mass_pi[:,2] * norm.pdf(coords[:,0], mass_mu[:,2], mass_sigma[:,2])
#
#    pdf_m = 1.0
#
#    # =============================================================================
#    # Total
#    # =============================================================================
#    pdf = (n_gmm * n_ms * n_m * pi_gmm * pdf_gmm * pdf_ms * pdf_m) # volume of cube **3 not needed, can be factored into nomalisations if necessary
#    return sum(pdf)


fileName = 'scenario_29_clusters_z1p25-2p0.fits'
s = fits.open(fileName)[1].data

distances = []

#for i in range(len(s)):
for i in [264]:
    for G in range(1):
        #Make sure you set a high enough n so that it catches
        #your gaussian in the grid somewhere!
        pdf,lowerGrid,upperGrid,midGrid,dGrid = pdf(s, i, G, 20)
    print(' ')
