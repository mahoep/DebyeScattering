# -*- coding: utf-8 -*-
"""
Created on Wed May 29 13:24:04 2024

@author: Matt

This script computes the Debye series for a homogoneous sphere. Alogrithim is 
derived from the following publication
https://doi.org/10.1364/AO.49.002422
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
        
        
        
def computeBessel(x, y, N):
    """
    This function pre computes all the permuations of spherical bessel functions 
    outlined in (https://doi.org/10.1364/AO.49.002422). All the computations are
    performed via recursion and never call scipy methods for bessel functions.
    This formulation ensures stability for an arbitrary function order, n.

    """
    nu = N + 0.5

    ajx = np.zeros((N+1,), dtype=np.complex128)
    ajy = np.zeros((N+1,), dtype=np.complex128)
    
    ajx[1] = 2 * nu/x
    ajy[1] = 2 * nu/y
    
    for j in range(2,N+1):
        ajx[j] = (-1.0) ** (j-1) * 2.0*(nu + np.real(j-1)) / x
        ajy[j] = (-1.0) ** (j-1) * 2.0*(nu + np.real(j-1)) / y
        
    D1x0 = ajx[-1]
    D1y0 = ajy[-1]
    
    for j in range(1, N):
        D1x0 = ajx[N-j] + 1/D1x0
        D1y0 = ajy[N-j] + 1/D1y0
        
    D1x = np.zeros( (N+1,), dtype=np.complex128)
    D1y = np.zeros( (N+1,), dtype=np.complex128)
    
    D1x[-1] = -N/x + D1x0
    D1y[-1] = -N/y + D1y0
    
    for n in range(1, N+1):
        en = N - n + 1
        D1x[N-n] = en/x - 1/(D1x[N-n+1] + en/x)
        D1y[N-n] = en/y - 1/(D1y[N-n+1] + en/y)
    
    D3x = np.zeros( (N+1,), dtype=np.complex128)
    D3y = np.zeros( (N+1,), dtype=np.complex128)
    
    D3x[0] = 1j
    D3y[0] = 1j
    
    w31 = np.zeros( (N+1,), dtype=np.complex128)
    w31[0] = -1/np.tan(y) + 1j
    
    for n in range(1,N+1):
        D3x[n] = -n/x + 1 / (n/x - D3x[n-1])
        w31[n] = w31[n-1] / ( (D3y[n-1] - n/y) * (D1y[n-1] - n/y) )
        D3y[n] = D1y[n] + w31[n]
        
        
    A13x = np.zeros( (N+1,), dtype=np.complex128)
    lnA13y = np.zeros( (N+1,), dtype=np.complex128)
    
    A13x[0] = 2/(1-1j/np.tan(x))
    lnA13y[0] = 2 * np.imag(y) + np.log( np.exp(-2*np.imag(y)) - np.cos(2*np.real(y)) + 1j*np.sin(2*np.real(y)))
    for n in range(1, N+1):
        A13x[n] = A13x[n-1] * (D1x[n-1] - n/x) / (D3x[n-1] - n/x)
        lnA13y[n] = lnA13y[n-1] + np.log( (D1y[n-1] - n/y) / (D3y[n-1] - n/y) )
    
    return D1x, D1y, D3x, D3y, A13x, lnA13y



    
if __name__ == "__main__":
    
    m1 = 1.33+0.00j # index of refraction, supports complex arguments
    m2 = np.complex128(1.0) # index of refraction of surrounding medium
    a = 10e-6 # radius of particle
    wavelength = 532e-9 # wavelength of light
    k = (2*np.pi / wavelength) # wavenumber of light

    x = 59.05249 # size parameter
    y = m1 / m2 * x
    N = (int(x + 4.05 * x**(1/3)) + 2) # number of terms to sum over
    N = 250
    P = 500 # number of Debye modes to consider
    Pc = 0 # Debye mode to consider by itself, must be 1 or greater
    
    theta = np.linspace(0,180,1801) * np.pi/180

    
    S1 = np.zeros(np.shape(theta), dtype=np.complex128)
    S2 = np.zeros(np.shape(theta), dtype=np.complex128)
    S1p = np.zeros(np.shape(theta), dtype=np.complex128)
    S2p = np.zeros(np.shape(theta), dtype=np.complex128)
    
    R212 = np.array([0,0], dtype=np.complex128)
    R121 = np.array([0,0], dtype=np.complex128)
    T1   = np.array([0,0], dtype=np.complex128)
    ray = np.zeros((2, P+1), dtype=np.complex128)
    
    
    D1x, D1y, D3x, D3y, A13x, lnA13y = computeBessel(x, y, N)
        
    
    for n in range(1,N+1):
        for i in range(0,2):
            if i == 0:
                alpha = m1/m2
                beta = 1
            elif i == 1:
                alpha = 1
                beta = m1/m2
            
            u0 = (D1x[n] - D3x[n]) * (D1y[n] - D3y[n])
            u11 = alpha*D1x[n] - beta*D1y[n]
            u31 = alpha*D3x[n] - beta*D1y[n]
            u13 = alpha*D1x[n] - beta*D3y[n]
            u33 = alpha*D3x[n] - beta*D3y[n]
            
            if (np.real(lnA13y[n]) < 0.0) or (np.imag(m1) == 0.0):
                A13y = np.exp(lnA13y[n])
                T1[i] = m1 * A13x[n] * A13y * u0 / (u33 - A13y*u31)**2
                
                umR212 = A13x[n] * (u13 - A13y*u11) / (u33 - A13y*u31)
                umR121 = -A13y * u31 / (u33 - A13y*u31)
                
                R212[i] = 1 - A13x[n] * (u13 - A13y*u11) / (u33 - A13y*u31)
                R121[i] = 1 + A13y * u31 / (u33 - A13y*u31)
            else:
                A13y = np.exp(-lnA13y[n])
                T1[i] = m1 * A13x[n] * A13y * u0 / (u33 - A13y*u31)**2
                
                umR212 = A13x[n] * (u13 - A13y*u11) / (u33 - A13y*u31)
                umR121 = -A13y * u31 / (u33 - A13y*u31)
                
                R212[i] = 1 - A13x[n] * (u13 - A13y*u11) / (u33 - A13y*u31)
                R121[i] = 1 + A13y * u31 / (u33 - A13y*u31)
                
            for p in range(1,P+1):
                ray[i,p] = T1[i] * (R121[i])**(p-1)
        raysum = np.sum(ray, axis=1)
            
        a1 = 0.5*(1 - R212[0] - raysum[0])
        b1 = 0.5*(1 - R212[1] - raysum[1])
        
        
        pi, tau = angularFunctions(theta, n)

        S1 += (2*n+1)/(n*(n+1)) * (a1*pi + b1*tau)
        S2 += (2*n+1)/(n*(n+1)) * (a1*tau + b1*pi)

            
        if Pc == 0:
            a1p = -0.5 + 0.5*R212[0]
            b1p = -0.5 + 0.5*R212[1]
            S1p += (2*n+1)/(n*(n+1)) * (a1p*pi + b1p*tau)
            S2p += (2*n+1)/(n*(n+1)) * (a1p*tau + b1p*pi)
        else:
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
    plt.semilogy(theta/np.pi*180, S1p, linewidth=1)
    plt.ylim([1e0, 1e8])
    plt.xlim([0,180])

                
        
            
                       
            