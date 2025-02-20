#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 09:31:48 2020

@author: rafalperczynski
"""

from matplotlib import cm
import numpy as np
from scipy.linalg import expm
from D2toeplitz import D2toeplitz
import matplotlib.pyplot as plt




N=128; h=2* np.pi/N; 
[x,D2] = D2toeplitz(N)

#przeskalowanie 
L = 10; x = L*(x-np.pi)/np.pi; D2 = (np.pi/L)**2*D2


#war początkowe
t= 0; dt = h/4; ep = 1
tmax = 1; tplot = 0.01;
p0 =np.zeros(N); q0 = x**2*np.exp(-x**2/2);
pq00 = np.append(p0,q0)

#macierze pomocnicze + macierz antydiagonalna z układu
I = np.identity(N); Z = np.zeros((N,N))
A = np.block([[Z,I],[D2-I*x**2,Z]])




c1 = -np.sqrt(3)/3; c2 = np.sqrt(3)/3
AA = expm(A*dt); AA0 = expm(A*dt/2)
AA1 = expm(-dt/2*c1*A); AA2 = expm(-dt/2*c2*A)
fx = -ep**2/2*x**4*np.exp(-x**2/2)
V0f = lambda c: np.append(np.zeros(N),fx*np.cos(c))

   


dt = tplot
nplots = round(tmax/tplot)

data = np.vstack([pq00, np.zeros((nplots,2*N))]); tdata = t;
for i in np.arange(1, nplots+1):
    t = t+dt 
#    cos1 = np.cos((i-1)*dt+dt/2*(1-c1))
#    cos2 = np.cos((i-1)*dt+dt/2*(1-c2))

    V1 = V0f(dt/2*c1+t-dt/2)
    V2 = V0f(dt/2*c2+t-dt/2)
    EC1 = AA1.dot(V1); EC2 = AA2.dot(V2)
    QW = dt/2*AA0.dot(EC1+EC2)
    pqnew = AA.dot(pq00)+QW; pq00 = pqnew
    data[i,:] = pq00; tdata = np.append(tdata,t)
  

omega = 1000
x,t = np.meshgrid(x,tdata)    
'''
p00 =  np.exp(-0.5*x**2)*np.cos(t)
p21 = ep*x**2/2*np.exp(-x**2/2)*np.cos(t)
U_as = p00+1/(omega**2)*(data[:,0:N]+p21*2*np.cos(omega*t))
'''




#3Dplot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(x,t, data[:,0:N],rstride=1, cstride=1,
             cmap=cm.coolwarm , edgecolor='none')   
ax.set_xlim(L, -L); ax.set_ylim(0,tmax); ax.set_zlim(-0.3,0.5);
ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('u')
ax.view_init(70, -50)
plt.title('r = 2, EXMP, dt ='+str(dt))