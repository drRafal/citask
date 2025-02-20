#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 19:52:47 2022

@author: rafalperczynski
"""
'''
Pierwszy przykład z doktoratu, równanie ODE u'(t)=tu(t)+cos(om*t)*u(t); u(0)=1;
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA

import matplotlib.pyplot as plot

log = lambda x: np.log(x)/np.log(10)



T = 1; 


Nvec = 2.0**np.arange(1,11)



data1=0

#prawa strona równania
om = 1
#jawne rozwiązanie równania
sol = np.exp(np.sin(om*T)/om+T**2/2)
f = lambda y,t: t*y+np.cos(om*t)*y
for n in Nvec:
    h = 1/n; t0=0; y0=1; nplots = round(T/(h))
    for i in np.arange(1, nplots+1):
        k1 = h*f(y0,t0)
        k2 = h*f(y0+k1/2,t0+h/2)
        k3 = h*f(y0+k2/2,t0+h/2)
        k4 = h*f(y0+k3,t0+h)
        y1 = y0+(k1+2*k2+2*k3+k4)/6
        t0 = t0+h
        y0=y1
    er = LA.norm(abs(sol-y0))
    data1 = np.append(data1,er)
    
data1 = np.delete(data1,0)


data2=0

#prawa strona równania
om = 1000
#jawne rozwiązanie równania
sol = np.exp(np.sin(om*T)/om+T**2/2)
f = lambda y,t: t*y+np.cos(om*t)*y
for n in Nvec:
    h = 1/n; t0=0; y0=1; nplots = round(T/(h))
    for i in np.arange(1, nplots+1):
        k1 = h*f(y0,t0)
        k2 = h*f(y0+k1/2,t0+h/2)
        k3 = h*f(y0+k2/2,t0+h/2)
        k4 = h*f(y0+k3,t0+h)
        y1 = y0+(k1+2*k2+2*k3+k4)/6
        t0 = t0+h
        y0=y1
    er = LA.norm(abs(sol-y0))
    data2 = np.append(data2,er)
    
data2 = np.delete(data2,0)

data3=0

#prawa strona równania
om = 10000
#jawne rozwiązanie równania
sol = np.exp(np.sin(om*T)/om+T**2/2)
f = lambda y,t: t*y+np.cos(om*t)*y
for n in Nvec:
    h = 1/n; t0=0; y0=1; nplots = round(T/(h))
    for i in np.arange(1, nplots+1):
        k1 = h*f(y0,t0)
        k2 = h*f(y0+k1/2,t0+h/2)
        k3 = h*f(y0+k2/2,t0+h/2)
        k4 = h*f(y0+k3,t0+h)
        y1 = y0+(k1+2*k2+2*k3+k4)/6
        t0 = t0+h
        y0=y1
    er = LA.norm(abs(sol-y0))
    data3 = np.append(data3,er)
    
data3 = np.delete(data3,0)

data4=0

#prawa strona równania
om = 100
#jawne rozwiązanie równania
sol = np.exp(np.sin(om*T)/om+T**2/2)
f = lambda y,t: t*y+np.cos(om*t)*y
for n in Nvec:
    h = 1/n; t0=0; y0=1; nplots = round(T/(h))
    for i in np.arange(1, nplots+1):
        k1 = h*f(y0,t0)
        k2 = h*f(y0+k1/2,t0+h/2)
        k3 = h*f(y0+k2/2,t0+h/2)
        k4 = h*f(y0+k3,t0+h)
        y1 = y0+(k1+2*k2+2*k3+k4)/6
        t0 = t0+h
        y0=y1
    er = LA.norm(abs(sol-y0))
    data4 = np.append(data4,er)
    
data4 = np.delete(data4,0)




plt.plot(log(T/Nvec),log(data1), '*-', label='$\omega$=1')

plt.plot(log(T/Nvec),log(data4), '<-',label='$\omega=10^2$')

plt.plot(log(T/Nvec),log(data2), '.-',label='$\omega=10^3$')

plt.plot(log(T/Nvec),log(data3), '-',label='$\omega=10^4$')

plt.plot(log(T/Nvec), log(Nvec**(-4)), '--' ,label='$h^{-4}$')

plt.xlabel('time step $h$')
plt.ylabel('error')

plt.legend()












