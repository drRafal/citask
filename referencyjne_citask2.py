#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 14:36:09 2021

@author: rafalperczynski
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 11:39:21 2021

@author: rafalperczynski
"""

'''
CITASK
Rozwiązanie referencyjne dla równania 
u_t = i*u_xx-i*V(x,t)*u t\in [0,5],  x\in[-10,10]
V(x,t) = V_D(x)*(1+exp(i*omega*t))
V_D(x) = (x-1/4)^2
Rozwinięcie Magnusa, dwie pierwsze omegi (czwartego rzędu)
exponent macierzowy liczymy za pomocą funkcji EXPM
 dt = 1/sto_tysięcy; N = 300
x_0 = 0; delta = 0.2
'''

from scipy.linalg import expm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz
import time
import csv
def D2toeplitz(N):
    
   
    h=2* np.pi/N; x= h*np.arange(0,N)
    k= np.arange(1,N)
    csc = lambda t: 1/np.sin(t)
    column = 0.5*(-1)**(k+1) *csc(k*h/2)**2
    column = np.append(-np.pi**2/(3*h**2)-1/6, column)
    D2 = toeplitz(column, column)
    return x, D2


def Magnus_ref(N,dt, t_max,L,omega):

    s = time.time()
    [x,D2] = D2toeplitz(N)
    
    #przeskalowanie -
    x = L*(x-np.pi)/np.pi; D2 = (np.pi/L)**2*D2
    
    #VD = x**4-1/2*x**2
    VD = (x-1/4)**2
    #war początkowy
    #x_0 = -np.sqrt(2)/4; 
    delta = 0.2
    x_0 = 0
    u_0 = (delta*np.pi)**(-1/4)*np.exp((-(x-x_0)**2)/(2*delta))
    
    I = np.identity(N); #Z = np.zeros((N,N))
    
    nplots = 100
    t = 0
    
    #zapisujemy dane w macierzy data wymiaru (nplots  x N)
    data = np.vstack([u_0, np.zeros((nplots,N),dtype=complex)]); tdata = t
    plotgap = t_max/(dt*nplots)
    
    #komutator i stała które możemy wyznaczyc poza pętlą
    com = D2.dot(np.diag(VD))-np.diag(VD).dot(D2)
    const = 1/(2*omega**2)*com*(2-2*np.exp(1j*omega*dt)+1j*omega*dt*(1+np.exp(1j*omega*dt)))
       
    for i in np.arange(1,nplots+1):
        for j in np.arange(1,plotgap+1):
            OM1 = dt*1j*D2-I*VD*(1j*dt+np.exp(1j*omega*t)/omega*(np.exp(1j*omega*dt)-1))
            OM2 = const*np.exp(1j*omega*t)
            #OM1 =  dt*1j*D2-I*VD*1j*dt
            
            
            u_1 = expm(OM1+OM2).dot(u_0)
            #u_1 = expm(OM1).dot(u_0)
            u_0 = u_1
            t = t+dt
        data[i,:] = u_0; tdata = np.append(tdata,t)
        
            # siatka meshgrid potrzebna do rysowania wykresu       
    xx,tt = np.meshgrid(x,tdata)
    data = abs(data)    
    el = time.time() - s  
    el = float(el)
    return data, xx, tt, el


 
#N = 300; dt = 0.00001; t_max = 5; L =10; omega = 50; 

#N = 100; dt = 0.001; t_max = 0.5; L =10; omega = 30; 

#[data,xx,tt,el] = Magnus_ref(N,dt, t_max,L,omega)


'''
  #3Dplot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(xx,tt, data,rstride=1, cstride=1,
             cmap='coolwarm' , edgecolor='none')   

ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('u')
'''


#plt.plot(xx[-1,:], data[-1,:])


#zapisanie danych w zewnętrznym pliku

'''
df = pd.DataFrame(data)
df1 = pd.DataFrame(xx)
df2 = pd.DataFrame(tt)


df.to_csv(r'DATA_om'+str(omega)+'.csv',index=False, header=False)
df1.to_csv(r'SPACE_om'+str(omega)+'.csv',index=False, header=False)
df2.to_csv(r'TIME_om'+str(omega)+'.csv',index=False, header=False)

'''


'''
with open("test.csv", "a+", newline='') as file:
    writer = csv.writer(file, delimiter=';')
    writer.writerow([el])

'''




df1 = pd.read_csv (r'/Users/rafalperczynski/Desktop/doktorat/CITASK/Dane_ t_max=5, L=10, dt=0.00001/SPACE_om25.csv')
df2 = pd.read_csv (r'/Users/rafalperczynski/Desktop/doktorat/CITASK/Dane_ t_max=5, L=10, dt=0.00001/TIME_om25.csv')
df3 = pd.read_csv (r'/Users/rafalperczynski/Desktop/doktorat/CITASK/Dane_ t_max=5, L=10, dt=0.00001/DATA_om25.csv')

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(df1,df2,df3,rstride=1, cstride=1,
             cmap='coolwarm' , edgecolor='none')   

ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('u')

df=pd.DataFrame(df1.iloc[-1,:].values)


