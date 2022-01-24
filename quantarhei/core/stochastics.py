# -*- coding: utf-8 -*-
"""
Entirely experimental - if I delete it, I would forget about. That's why it is here

"""
import numpy as np
import time
import matplotlib.pyplot as plt


# correlation functions
def cor(lam, g, a, b, t):
    #return lam*(a-1j*b)*np.exp(-g*np.abs(t))
    #om = a
    #T = b
    #re = lam*(1.0/np.tanh(om/(2.0*T)))*np.cos(om*t)
    #im = -lam*np.sin(om*t)
    #return re + 1j*im
    gam = a
    T = b
    re = 2.0*lam*T*np.exp(-gam*np.abs(t))
    im = -lam*gam*np.exp(-gam*np.abs(t))

    return re + 1j*im


def run():
    g11 = 1.0/200
    l11 = 0.2
    l12 = 0.2
    dt = 0.2
    
    a = 1.0/100.0
    b = a/3.0
    
    # time
    iTmax = 1000
    tt = np.zeros(iTmax,dtype=np.float64)
    
    # mean and covariance matrix
    mean = np.zeros(2*iTmax,dtype=np.float64)
    cov = np.zeros((2*iTmax,2*iTmax),dtype=np.float64)
    
    # values of correlation function
    val = np.zeros(2*iTmax,dtype=np.complex128)
    for tau in range(2*iTmax):
        val[tau] = cor(l11,g11,a,b,tau*dt-iTmax*dt)
    
    t1 = time.time()
    
    for it1 in range(iTmax):
        tt1 = it1*dt
        tt[it1] = tt1
        for it2 in range(iTmax):
            itau = it1+iTmax-it2
            cov[it1,it2] = 3.0*np.real(val[itau])
            cov[iTmax+it1,it2] = np.imag(val[itau])
            cov[it1,iTmax+it2] = np.imag(val[itau])
            cov[iTmax+it1,iTmax+it2] = np.real(val[itau])   
    cov = cov/2.0
    
    t2 = time.time()
      
    # get realizations of noise    
    Nreal = 10000
    xi = np.random.multivariate_normal(mean, cov, Nreal).T
    
    # construct complex noise
    zi = np.zeros((iTmax,Nreal), dtype=np.complex128)
    for i in range(iTmax):
        zi[i,:] = xi[i,:]+1j*xi[iTmax+i,:]
        
    t3 = time.time()
    print(zi.shape)
    
    # correlation function for verification
    ct = np.zeros(iTmax,dtype=np.complex128)
    ctC = np.zeros(iTmax,dtype=np.complex128)
    for it in range(iTmax):
        ct[it] = np.sum(zi[0,:]*zi[it,:])/xi.shape[1]
        ctC[it] = np.sum(zi[0,:]*np.conj(zi[it,:]))/xi.shape[1]
        
    print(ct.shape)
    print(t2-t1)
    print(t3-t1)

    # propagation
    rhoav = np.zeros(iTmax, dtype=np.complex128)
    rhoeg = np.zeros(iTmax, dtype=np.complex128)
    #iTmax = 500
    for k in range(Nreal):
        rhoeg[0] = 1.0
        for it in range(1,iTmax):
            rhoeg[it] = np.exp(-1j*np.real(zi[it-1,k])*dt)*rhoeg[it-1]
        rhoav += rhoeg
    rhoav = rhoav/(Nreal)
    
    plt.plot(tt,np.real(ct))
    plt.plot(tt,np.real(cor(l11,g11,a,b,tt)))
    plt.plot(tt,np.imag(ct))
    plt.plot(tt,np.imag(cor(l11,g11,a,b,tt)))
    #plt.plot(tt,np.real(ctC))
    #plt.plot(tt,np.imag(ctC))
    
    plt.show()

    return zi, tt, rhoav

#zi, tt, rhoeg = run()



##import pyximport
##pyximport.install()
##
##import cy_stochastics as cs
##
##cs.run()
#print(zi.shape)
#print(tt.shape)
#plt.plot(tt, np.real(zi[:,0]))
#plt.plot(tt, np.imag(zi[:,0]))
#plt.show()
##plt.plot(tt, np.real(zi[:,1]))
##plt.plot(tt, np.imag(zi[:,1]))
##plt.show()
#plt.plot(tt, np.real(rhoeg))
#plt.plot(tt, np.imag(rhoeg))
#plt.show()