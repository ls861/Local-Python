#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 16:06:01 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import lambertw



### ### ### ### ### ### ###
### Delayed Exponential ###
### ### ### ### ### ### ### 

### ### mass vs tau ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5
msa = 5E8                   # 10**8.70
t0 = t - msa

tau_arr = 10**np.arange(7.5, 10.1, 0.1)
mass_arr = A*(-tau_arr * ( ((np.exp(-(t-t0)/tau_arr))*(-t0+tau_arr+t)) - tau_arr) )

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(tau_arr, mass_arr / max(mass_arr))
ax.plot((msa, msa), (0, 1.1), color='k', linestyle=':')   
ax.set_xlim(1E10-t0, 0-t0)
ax.set_ylim(0, 1.1)
ax.set_xlabel('TAU')
ax.set_ylabel('MASS / A')
plt.show()
    


### ### ### ### ###  ###
### Double Power Law ###
### ### ### ### ###  ###  

### ### mass vs tau ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5 / 10**9.049218
#alpha = 10.**np.arange(-1, 4, 1) 
alpha = [0.5, 50]
beta = 1

tau_arr = 10**np.arange(7.5, 10.1, 0.0005)
mass_arr = np.zeros(len(tau_arr))

fig, ax = plt.subplots(figsize=(8, 4))
    
for j in range(len(alpha)):

    for i in range(len(tau_arr)):
        
        integrand = lambda T: 1 / (((T/tau_arr[i])**alpha[j])+((T/tau_arr[i])**-beta))
        integral  = quad(integrand, 0, t)
    
    
        mass_arr[i] = A * integral[0] 
    
    # https://www.wolframalpha.com/input/?i=%28%28t%2FC%29%5EB+%28B+-+A+%28t%2FC%29%5E%28A+%2B+B%29%29%29%2F%28t+%281+%2B+%28t%2FC%29%5E%28A+%2B+B%29%29%5E2%29&assumption=%22ClashPrefs%22+-%3E+%7B%22Math%22%7D
    
    #t_peak = np.exp( (   (alpha*np.log(tau)) - np.log(alpha) + (beta*np.log(tau)) + np.log(beta))    /  (alpha + beta))
    
    
    tau_peak = np.exp(    ( (alpha[j]*np.log(t)) + np.log(alpha[j]) + (beta*np.log(t)) - np.log(beta))  /  (alpha[j] + beta) )
    

    ax.plot(tau_arr, mass_arr / max(mass_arr), color='C%s' % (j), label=r'$\alpha$ = {%.2g}, $\beta$ = {%.2g}' % (alpha[j], beta))
    ax.plot((t, t), (0, 1.1), color='k', linestyle=':')   
    ax.plot((tau_peak, tau_peak), (0, 1.1), color='C%s' % (j), linestyle=':')   
    ax.set_xlim(5E9, 0)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel('TAU')
    ax.set_ylabel('MASS / A')

plt.legend()
plt.show()

### ### mass vs alpha ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5 / 10**9.049218
alpha_arr = 10.**np.arange(-1, 4, 0.1)
beta = 50.

tau = [10**8.5, 10**9]
mass_arr = np.zeros(len(alpha_arr))

fig, ax = plt.subplots(figsize=(8, 4))
    
for j in range(len(tau)):
    
    for i in range(len(alpha_arr)):
        
        integrand = lambda T: 1 / (((T/tau[j])**alpha_arr[i])+((T/tau[j])**-beta))
        integral  = quad(integrand, 0, t)
    
        mass_arr[i] = A * integral[0] 
    
    # https://www.wolframalpha.com/input/?i=%28%28t%2FC%29%5EB+%28B+-+A+%28t%2FC%29%5E%28A+%2B+B%29%29%29%2F%28t+%281+%2B+%28t%2FC%29%5E%28A+%2B+B%29%29%5E2%29&assumption=%22ClashPrefs%22+-%3E+%7B%22Math%22%7D
    #t_peak = np.exp( (   (alpha*np.log(tau)) - np.log(alpha) + (beta*np.log(tau)) + np.log(beta))    /  (alpha + beta))
    #tau_peak = np.exp(    ( (alpha[j]*np.log(t)) + np.log(alpha[j]) + (beta*np.log(t)) - np.log(beta))  /  (alpha[j] + beta) )
    #Solve[((t/C)^B (B - A (t/C)^(A + B)))/(t (1 + (t/C)^(A + B))^2) == 0, A]
    
    alpha_peak = np.real(lambertw(  beta * ((t/tau[j])**(-beta)) * np.log(t/tau[j])   ))  / np.log(t/tau[j])  
    
    ax.plot(alpha_arr, mass_arr / max(mass_arr), color='C%s' % (j), label=r'$\tau$ = {%.2g}, $\beta$ = {%.2g}' % (tau[j], beta))
    ax.plot((alpha_peak, alpha_peak), (0, 1.1), color='C%s' % (j), linestyle=':')   
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel('ALPHA')
    ax.set_ylabel('MASS / A')
    ax.set_xscale('log')

plt.legend()
plt.show()

### ### mass vs beta ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5 / 10**9.049218
beta_arr = 10.**np.arange(-1, 4, 0.1)
alpha = 50.

tau = [10**8.5, 10**9]
mass_arr = np.zeros(len(beta_arr))

fig, ax = plt.subplots(figsize=(8, 4))
    
for j in range(len(tau)):
    
    for i in range(len(beta_arr)):
        
        integrand = lambda T: 1 / (((T/tau[j])**alpha)+((T/tau[j])**-beta_arr[i]))
        integral  = quad(integrand, 0, t)
    
        mass_arr[i] = A * integral[0] 
    
    # https://www.wolframalpha.com/input/?i=%28%28t%2FC%29%5EB+%28B+-+A+%28t%2FC%29%5E%28A+%2B+B%29%29%29%2F%28t+%281+%2B+%28t%2FC%29%5E%28A+%2B+B%29%29%5E2%29&assumption=%22ClashPrefs%22+-%3E+%7B%22Math%22%7D
    #t_peak = np.exp( (   (alpha*np.log(tau)) - np.log(alpha) + (beta*np.log(tau)) + np.log(beta))    /  (alpha + beta))
    #tau_peak = np.exp(    ( (alpha[j]*np.log(t)) + np.log(alpha[j]) + (beta*np.log(t)) - np.log(beta))  /  (alpha[j] + beta) )
    #Solve[((t/C)^B (B - A (t/C)^(A + B)))/(t (1 + (t/C)^(A + B))^2) == 0, A]
    #alpha_peak = np.real(lambertw(  beta * ((t/tau[j])**(-beta)) * np.log(t/tau[j])   ))  / np.log(t/tau[j])  
    
    beta_peak = -np.real(lambertw(-alpha * ((t/tau[j])**alpha) * np.log(t/tau[j])))  / np.log(t/tau[j])
    
    ax.plot(beta_arr, mass_arr / max(mass_arr), color='C%s' % (j), label=r'$\tau$ = {%.2g}, $\alpha$ = {%.2g}' % (tau[j], alpha))
    ax.plot((beta_peak, beta_peak), (0, 1.1), color='C%s' % (j), linestyle=':')   
    ax.set_xlim(1000, 0.1)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel('BETA')
    ax.set_ylabel('MASS / A')
    ax.set_xscale('log')
    print(beta_peak)

plt.legend()
plt.show()



### ### ### ### ### ###  ###
### Linear / Exp ###
### ### ### ### ### ###  ###

def lin_exp(t_plot, t0_value, tp_value, tau_value):
    sfr_plot = np.zeros(len(t_plot))
    for i in range(len(t_plot)):
        if t_plot[i] < t0_value:
            sfr_plot[i] = 0
        elif t_plot[i] < tp_value:
            sfr_plot[i] = (t_plot[i] - t0_value) / (tp_value - t0_value)
        else:
            sfr_plot[i] = np.exp(-t_plot[i]/tau_value) * np.exp(tp_value/tau_value)
    sfr_plot = sfr_plot / max(sfr_plot)
    return sfr_plot

### ### mass vs tau rising ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5
msa = 5E8                   # 10**8.70
t0 = t - msa
tp = t + msa/2

tau_arr = np.linspace(1E7, 1E10, 1000)
mass_arr = np.zeros(len(tau_arr))

for i in range(len(mass_arr)):
    if t <= tp:
        mass_arr[i] = ((t-t0)**2) / (2*(tp-t0))
    elif t > tp:
        mass_arr[i] = ((tp - t0) / 2) + tau_arr[i] * (1 - np.exp(-((t-tp)/tau_arr[i])))

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(tau_arr, mass_arr / max(mass_arr))
#ax.plot((t0, t0), (0, 1.1), color='k', linestyle=':')   
#ax.plot((tp, tp), (0, 1.1), color='k', linestyle=':')   
ax.set_xlim(0, 1E10)
ax.set_ylim(0, 1.1)
ax.set_xlabel('TAU')
ax.set_ylabel('MASS / A')
plt.show()

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, t0, tp, min(tau_arr)), label=r'$\tau$ = {%.2g}' % (min(tau_arr)))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, t0, tp, max(tau_arr)))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###


### ### mass vs tau falling ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5
msa = 5E8                   # 10**8.70
t0 = t - msa
tp = t - msa/2

tau_arr = np.linspace(1E7, 1E10, 1000)
mass_arr = np.zeros(len(tau_arr))

for i in range(len(mass_arr)):
    if t <= tp:
        mass_arr[i] = ((t-t0)**2) / (2*(tp-t0))
    elif t > tp:
        mass_arr[i] = ((tp - t0) / 2) + tau_arr[i] * (1 - np.exp(-((t-tp)/tau_arr[i])))

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(tau_arr, mass_arr / max(mass_arr))
#ax.plot((t0, t0), (0, 1.1), color='k', linestyle=':')   
#ax.plot((tp, tp), (0, 1.1), color='k', linestyle=':')   
ax.set_xlim(0, 1E10)
ax.set_ylim(0, 1.1)
ax.set_xlabel('TAU')
ax.set_ylabel('MASS / A')
plt.show()

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, t0, tp, min(tau_arr)), label=r'$\tau$ = {%.2g}' % (min(tau_arr)))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, t0, tp, max(tau_arr)))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###


### ### mass vs t0 rising ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5


tp = t + 5E8
t0_arr = np.linspace(0.1, tp, 1000)

tau = 1E9
mass_arr = np.zeros(len(t0_arr))

for i in range(len(mass_arr)):
    if t0_arr[i] >= t:
        mass_arr[i] = 0
    else:
        if t <= tp:
            mass_arr[i] = ((t-t0_arr[i])**2) / (2*(tp-t0_arr[i]))
        elif t > tp:
            mass_arr[i] = ((tp - t0_arr[i]) / 2) + tau * (1 - np.exp(-((t-tp)/tau)))

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t0_arr, mass_arr / max(mass_arr))
#ax.plot((t0, t0), (0, 1.1), color='k', linestyle=':')   
#ax.plot((tp, tp), (0, 1.1), color='k', linestyle=':')   
ax.set_xlim(0, tp)
ax.set_ylim(0, 1.1)
ax.set_xlabel('T0')
ax.set_ylabel('MASS / A')
plt.show()

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, min(t0_arr), tp, tau), label=r'$\tau$ = {%.2g}' % (min(tau_arr)))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, max(t0_arr), tp, tau))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###


### ### mass vs t0 falling ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5
tp = t - 5E8

t0_arr = np.linspace(0.1, tp, 1000)



tau = 1E9
mass_arr = np.zeros(len(t0_arr))

for i in range(len(mass_arr)):
    if t <= tp:
        print('error')
        mass_arr[i] = ((t-t0_arr[i])**2) / (2*(tp-t0_arr[i]))
    elif t > tp:

        mass_arr[i] = ((tp - t0_arr[i]) / 2) + tau * (1 - np.exp(-((t-tp)/tau)))

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t0_arr, mass_arr / max(mass_arr))
#ax.plot((t0, t0), (0, 1.1), color='k', linestyle=':')   
#ax.plot((tp, tp), (0, 1.1), color='k', linestyle=':')   
ax.set_xlim(0, tp)
ax.set_ylim(0, 1.1)
ax.set_xlabel('T0')
ax.set_ylabel('MASS / A')
plt.show()

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, min(t0_arr), tp, tau), label=r'$\tau$ = {%.2g}' % (min(tau_arr)))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, max(t0_arr), tp, tau))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###


### ### mass vs tp rising ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5
tau = 1E9
t0 = t - 5E8

tp_arr = np.linspace(t, 1E10, 1000)


mass_arr = np.zeros(len(tp_arr))

for i in range(len(mass_arr)):
    if t <= tp_arr[i]:
        mass_arr[i] = ((t-t0)**2) / (2*(tp_arr[i]-t0))
    elif t > tp:
        print('error')
        mass_arr[i] = ((tp - t0_arr[i]) / 2) + tau * (1 - np.exp(-((t-tp)/tau)))

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(tp_arr, mass_arr / max(mass_arr))
#ax.plot((t0, t0), (0, 1.1), color='k', linestyle=':')   
#ax.plot((tp, tp), (0, 1.1), color='k', linestyle=':')   
ax.set_xlim(t, 1E10)
ax.set_ylim(0, 1.1)
ax.set_xlabel('TP')
ax.set_ylabel('MASS / A')
plt.show()

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, t0, min(tp_arr), tau), label=r'$\tau$ = {%.2g}' % (min(tau_arr)))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, t0, max(tp_arr), tau))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###

### ### mass vs tp falling ### ###

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5
tau = 1E9
t0 = t - 5E8

tp_arr = np.linspace(t0, t, 1000)


mass_arr = np.zeros(len(tp_arr))

for i in range(len(mass_arr)):
    if t < tp_arr[i]:
        print('error')
        mass_arr[i] = ((t-t0)**2) / (2*(tp_arr[i]-t0))
    elif t >= tp_arr[i]:

        mass_arr[i] = ((tp_arr[i] - t0) / 2) + tau * (1 - np.exp(-((t-tp_arr[i])/tau)))

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(tp_arr, mass_arr / max(mass_arr))
#ax.plot((t0, t0), (0, 1.1), color='k', linestyle=':')   
#ax.plot((tp, tp), (0, 1.1), color='k', linestyle=':')   
ax.set_xlim(t0, t)
ax.set_ylim(0, 1.1)
ax.set_xlabel('TP')
ax.set_ylabel('MASS / A')
plt.show()

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, t0, min(tp_arr), tau), label=r'$\tau$ = {%.2g}' % (min(tau_arr)))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###

### PLOT ###
t_plot = np.linspace(0., 1E10, 1000)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_plot, lin_exp(t_plot, t0, max(tp_arr), tau))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')    
ax.set_xlim(0, 3E9)
ax.set_ylim(0, 1.1)
ax.set_xlabel('t')
ax.set_ylabel('SFR')
plt.show()
### ###  ###






















