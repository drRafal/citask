
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 16:08:16 2020

@author: rafalperczynski
"""

'''
Rozwiązanie dla 
g(x,t) = x^2+epsilon*cos(omega*t)*x^2
r = 0,1,2,3
Praca paper wave equations, rozwinięcie asymptotyczne, metoda drugiego rzędu
Strang splitting

'''
from scipy.linalg import expm
import pandas as pd
import numpy as np
from D2toeplitz import D2toeplitz
import matplotlib.pyplot as plt
from SSWE import SSWE
from mpl_toolkits.mplot3d import Axes3D


def wave_asimptotic(N,L,t,tmax,dt,epsilon,omega, plotgap):



#funkcja macierzy drugiej pochodnej
    [x,D2] = D2toeplitz(N)

#przeskalowanie -
    x = L*(x-np.pi)/np.pi; D2 = (np.pi/L)**2*D2




    #warunek początkowy dla równania r = 2
    p0 = -epsilon*x**2*np.exp(-x**2/2); q0 = np.zeros(N);
    pq00 = np.append(p0,q0)

    I = np.identity(N); Z = np.zeros((N,N))


    nplots = round(tmax/(plotgap*dt))


    #współcznniki kwadratury
    c1 = -np.sqrt(3)/3; c2 = np.sqrt(3)/3
    #eksponenty z kwadratury wyznaczane przy pomocy Strang Splitting (zewnętrna funkcja SSWE)
    #które liczymy poza pętlą
    
    B = lambda h: np.block([[ np.diag(np.cos(h/2*x)),np.diag(h/2*np.sinc(h/2*x/np.pi)) ],[-np.diag(x*np.sin(h/2*x)),np.diag(np.cos(h/2*x))]  ])
    C = lambda h: np.block([[I,Z ],[h*D2,I]])
    OM = lambda h: B(h).dot(C(h)).dot(B(h))
   
    AAA = np.block([[Z,I],[D2-I*x**2,Z]])
    EX = expm(dt*AAA)
    AA =  OM(dt); AA0 = OM(dt/2)
    AA1 = OM(-dt/2*c1); AA2 = OM(-dt/2*c2)
    fx = -epsilon**2/2*x**4*np.exp(-x**2/2)
    #lokalna funkcja pomocnicza
    V0f = lambda c: np.append(np.zeros(N),fx*np.cos(c))
    
    er = abs(EX-OM(dt)).max()


    # algorytm wyznaczający współczynnik p20
    # ( dane w macierzy data zapisujemy co plotgap)
    data = np.vstack([pq00, np.zeros((nplots,2*N))]); tdata = t;
    for i in np.arange(1, nplots+1):
        for j in np.arange(1,plotgap+1):
    
            t = t+dt 
            V1 = V0f(dt/2*c1+t-dt/2)
            V2 = V0f(dt/2*c2+t-dt/2)
            EC1 = AA1.dot(V1); EC2 = AA2.dot(V2)
            QW = dt/2*AA0.dot(EC1+EC2)
            pqnew = AA.dot(pq00)+QW; pq00 = pqnew
        data[i,:] = pq00; tdata = np.append(tdata,t)
  



#wspólczynniki rozwinięcia asymptotycznego do r = 3


    xx,tt = np.meshgrid(x,tdata)    
    p00 =  np.exp(-xx**2/2)*np.cos(tt)
    As0 = p00;
    p20 = data[:,0:N]
    p21 = epsilon*xx**2/2*np.exp(-xx**2/2)*np.cos(tt)
    As2 = 1/(omega**2)*(p20+p21*2*np.cos(omega*tt))
    p31 = epsilon*xx**2*np.exp(-xx**2/2)*np.sin(tt)
    As3 = 1/(omega**3)*2*np.sin(omega*tt)*p31
    U_as = As0 + As2 + As3

    return U_as, xx, tt,er
    #return  As0, As2, As3, xx, tt




N = 100; L =10 ; t = 0; tmax = 1; dt = 0.001; epsilon = 1; omega = 1000; plotgap =10


[U_as, xx,tt,er] = wave_asimptotic(N,L,t,tmax,dt,epsilon,omega, plotgap)

plt.plot(xx[-1,:], U_as[-1,:])


#plt.title(r'Paper Wave, $g(x,t) = (1+ \epsilon(\cos(\omega t))$, $\epsilon$ = '+str(epsilon)+ ', $\omega=$'+str(omega))
#max = U_as[-1,:].max()

#plt.savefig('wave.png')

'''
#3Dplot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(xx,tt, U_as,rstride=1, cstride=1,
             cmap='coolwarm' , edgecolor='none')   

ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('u')
'''



#plt.title(' $U_{asym}, \ r \leq 3, \   g(x,t) = x^2+\epsilon cos(\omega t)x^2$, $\epsilon=$'+str(epsilon)+', $\omega=$'+str(omega))     


#plt.savefig('porownanie.png')

#zapisanie danych w zewnętrznym pliku
'''
df = pd.DataFrame(U_as)
df1 = pd.DataFrame(xx)
df2 = pd.DataFrame(tt)
df.to_csv(r'/Users/rafalperczynski/Desktop/CITASK/DATA.csv', index=False, header=True)
df1.to_csv(r'/Users/rafalperczynski/Desktop/CITASK/SPACE.csv', index=False, header=True)
df2.to_csv(r'/Users/rafalperczynski/Desktop/CITASK/TIME.csv', index=False, header=True)
'''
'''
df.rename(columns={0: 'DATA: EPSILON = '+str(ep)+',  OMEGA = '+str(omega)+',  dt = '+str(dt)+',  T = '+str(tmax)+',  N = '+str(N)+'   |||'}, inplace=True)
df1 = pd.DataFrame(t)
df1.rename(columns={0:'CZAS t: dt = '+str(dt)+',  T = '+str(tmax)+', N = '+str(N)+'   |||'}, inplace = True)
df2 = pd.DataFrame(x)
df2.rename(columns={0:'PRZESTRZEN x: dt = '+str(dt)+',  T = '+str(tmax)+', N = '+str(N)+'   |||'}, inplace=True)

#data
df.to_csv(r'/Users/rafalperczynski/Desktop/CITASK/CITASK.csv', index=False, header=True)
#df1 - czas t

#df1.to_csv(r'/Users/rafalperczynski/Desktop/data/T=10;dt=0.001;N=128/dt=0.001.csv',index=False, header=True)
# df2 - przestrzeń x
#df2.to_csv(r'/Users/rafalperczynski/Desktop/data/T=10;dt=0.001;N=128/x.csv',index=False, header=True)

#ax.set_xlim(-L, +L); ax.set_ylim(0,tmax); 
'''

'''
DATA = pd.read_csv (r'/Users/rafalperczynski/Desktop/CITASK/CITASK.csv')
TIME = pd.read_csv (r'/Users/rafalperczynski/Desktop/CITASK/SPACE.csv')
SPACE = pd.read_csv (r'/Users/rafalperczynski/Desktop/CITASK/TIME.csv')

'''
'''
 #3Dplot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(SPACE,TIME, DATA,rstride=1, cstride=1,
             cmap='coolwarm' , edgecolor='none')   
ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('u')

'''


