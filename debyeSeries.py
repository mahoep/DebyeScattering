# -*- coding: utf-8 -*-
"""
Created on Sun May 26 15:38:41 2024

@author: Matt

This script computes the Debye series for a homogoneous sphere. Alogrithim is 
derived from the following publication
https://doi.org/10.1364/AO.45.001260
"""

import numpy as np

def angularFunctions(theta, N):
    """
    Calculates the Angular functions using recurrence formulation
    Reference:
    https://opg.optica.org/ao/fulltext.cfm?uri=ao-49-13-2422&id=198319
    """

    if np.ndim(theta) == 0:
        pi = np.zeros((N+1, ), dtype=np.longdouble)
        tau = np.zeros((N+1, ), dtype=np.longdouble)
        
        pi[0] = 0.
        pi[1] = 1.
        
        tau[0] = 0.
        tau[1] = 1*np.cos(theta) * pi[1]- (1+1)*pi[1-1]
        
        if N == 1:
            return pi[1], tau[1]
        
        else:
            for n in range(2,N+1):
                pi[n] = (2*n-1) / (n-1) * np.cos(theta) * pi[n-1] - n/(n-1) * pi[n-2]
                tau[n] = n * np.cos(theta) * pi[n] - (n+1) * pi[n-1]
        
        return pi[-1], tau[-1]

    elif np.ndim(theta) == 1:
        pi = np.zeros((N+1, np.size(theta)), dtype=np.longdouble)
        tau = np.zeros((N+1, np.size(theta)), dtype=np.longdouble)
        
        pi[0,:] = 0.
        pi[1,:] = 1.
        
        tau[0,:] = 0.
        tau[1,:] = 1.*np.cos(theta) * pi[1,:] - (1+1)*pi[1-1,:]
        
        if N == 1:
            return pi[1,:], tau[1,:]
        
        else:
            for n in range(2,N+1):
                pi[n,:] = (2*n-1) / (n-1) * np.cos(theta) * pi[n-1,:] - n/(n-1) * pi[n-2,:]
                tau[n,:] = n * np.cos(theta) * pi[n,:] - (n+1) * pi[n-1,:]
        
        return pi[-1,:], tau[-1,:]
    
    else:
        TypeError("Theta must be 0 or 1 dimensional")
    
    
def A3f(z, N):
    """
    Computes the Riccati Bessel Functions using recurrence
    Reference: https://opg.optica.org/viewmedia.cfm?r=1&rwjcode=ao&uri=ao-45-6-1260&html=true
    """
    
    A3 = np.zeros( (N+1, ), dtype=np.complex128)
    
    A3[0] = 1j
    
    for n in range(1,N+1):
        A3[n] = 1 / (n/z - A3[n-1]) - n/z
        
    return A3


def A4f(z, N):
    """
    Computes the Riccati Bessel Functions
    Reference: https://opg.optica.org/viewmedia.cfm?r=1&rwjcode=ao&uri=ao-45-6-1260&html=true
    """
    
    A4 = np.zeros((N+1, ), dtype=np.complex128)
    
    A4[0] = -1j
    
    for n in range(1,N+1):
        A4[n] = 1 / (n/z - A4[n-1]) - n/z        
    
    return A4


def Bf(z, N):
    from scipy import special as spc
    
    B = np.zeros( (N+1, ), dtype=np.complex128)
    
    for n in range(1, N+1):
        B[n] = (spc.riccati_jn(n, z)[0][-1] + 1j*spc.riccati_yn(n, z)[0][-1]) / (spc.riccati_jn(n, z)[0][-1] - 1j*spc.riccati_yn(n, z)[0][-1])

    return B
          
    
if __name__ == "__main__":
    
    m1 = 1.333+0.00j
    m2 = np.complex128(1.0)
    a = 10e-6
    wavelength = 532e-9
    k = (2*np.pi / wavelength)

    x = 100
    y = m1 / m2 * x
    # N = 250
    N = round(x + 4.05*x**(1/3) + 2)
    P = 15
    Pc = 1
    
    theta = np.linspace(0,180,1801) * np.pi/180

    
    S1 = np.zeros(np.shape(theta), dtype=np.complex128)
    S2 = np.zeros(np.shape(theta), dtype=np.complex128)
    S1p = np.zeros(np.shape(theta), dtype=np.complex128)
    S2p = np.zeros(np.shape(theta), dtype=np.complex128)
    
    R212 = np.array([0,0], dtype=np.complex128)
    R121 = np.array([0,0], dtype=np.complex128)
    T1   = np.array([0,0], dtype=np.complex128)
    ray = np.zeros((2, P+1), dtype=np.complex128)
    
    A3x = A3f(x, N)
    A3y = A3f(y, N)
    A4x = A4f(x, N)
    A4y = A4f(y, N)
    Bx = Bf(x, N)
    By = Bf(y, N)
    
    
    for n in range(1,N):
        for i in range(0,2):
            if i == 0:
                alpha = m1/m2
                beta = 1
            elif i == 1:
                alpha = 1
                beta = m1/m2
                
        
            R212[i] = 1/Bx[n] * (alpha*A4x[n] - beta*A4y[n]) / (-alpha*A3x[n] + beta*A4y[n])
                
                
            R121[i] = By[n] * (alpha*A3x[n] - beta*A3y[n]) / (-alpha*A3x[n] + beta*A4y[n])
                    
                    
            T1[i] = m1/m2 * By[n]/Bx[n] * (A4y[n] - A3y[n]) / (-alpha*A3x[n] + beta*A4y[n]) * (A4x[n] - A3x[n]) / (-alpha*A3x[n] + beta*A4y[n])
        
            for p in range(1,P+1):
                ray[i,p] = T1[i] * (R121[i])**(p-1)
        raysum = np.sum(ray, axis=1)
            
        a1 = 0.5*(1 - R212[0] - raysum[0])
        b1 = 0.5*(1 - R212[1] - raysum[1])
        
        
        pi, tau = angularFunctions(theta, n)

        S1 += (2*n+1)/(n*(n+1)) * (a1*pi + b1*tau)
        S2 += (2*n+1)/(n*(n+1)) * (a1*tau + b1*pi)

        
        if P > 1:
            a1p = -0.5 * ray[0,Pc]
            b1p = -0.5 * ray[1,Pc]
            S1p += (2*n+1)/(n*(n+1)) * (a1p*pi + b1p*tau)
            S2p += (2*n+1)/(n*(n+1)) * (a1p*tau + b1p*pi)
        
        #end loop
            
    S1a = np.abs(S1)**2
    S2a = np.abs(S2)**2
    
    S1p = np.abs(S1p)**2
    S2p = np.abs(S2p)**2


    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    mpl.rcParams['figure.dpi'] = 300
    plt.semilogy(theta/np.pi*180, S1a, linewidth=1)
    # plt.semilogy(theta/np.pi*180, S1p, linewidth=1)
    plt.ylim([1e0, 1e8])
    plt.xlim([0,180])

                
        
            
                       
            