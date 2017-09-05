import numpy as np
import time
import matplotlib.pyplot as plt


# correlation functions
cdef double complex cor(float lam, float g, float a,float b,t):
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
    cdef int iTmax
    cdef int it1, it2, itau
    
    g11 = 1.0/200
    l11 = 1.0
    l12 = 0.2
    dt = 1.0
    
    a = 1.0/100.0
    b = a/3.0
    
    iTmax = 800
    
    mean = np.zeros(2*iTmax,dtype=np.float64)
    cov = np.zeros((2*iTmax,2*iTmax),dtype=np.float64)
    tt = np.zeros(iTmax,dtype=np.float64)
    
    val = np.zeros(2*iTmax,dtype=np.complex128)
    for itau in range(2*iTmax):
        val[itau] = cor(l11,g11,a,b,itau*dt-iTmax)
    tt = [it1*dt for it1 in range(iTmax)]
    
    t1 = time.time()
    for it1 in range(iTmax):
        for it2 in range(iTmax):
            itau = it1+iTmax-it2
            cov[it1,it2] = 3.0*np.real(val[itau])
            cov[iTmax+it1,it2] = np.imag(val[itau])
            cov[it1,iTmax+it2] = np.imag(val[itau])
            cov[iTmax+it1,iTmax+it2] = np.real(val[itau])
     
    cov = cov/2.0
    t2 = time.time()
          
    xi = np.random.multivariate_normal(mean, cov, 1000).T
    
    zi = np.zeros(xi.shape, dtype=np.complex128)
    for i in range(iTmax):
        zi[i,:] = xi[i,:]+1j*xi[iTmax+i,:]
        
    t3 = time.time()
    print(zi.shape)
    
    ct = np.zeros(iTmax,dtype=np.complex128)
    ctC = np.zeros(iTmax,dtype=np.complex128)
    for it in range(iTmax):
        ct[it] = np.sum(zi[0,:]*zi[it,:])/xi.shape[1]
        ctC[it] = np.sum(zi[0,:]*np.conj(zi[it,:]))/xi.shape[1]
        
    print(ct.shape)
    print(t2-t1)
    print(t3-t1)
    
    plt.plot(tt,np.real(ct))
    plt.plot(tt,np.real(cor(l11,g11,a,b,tt)))
    plt.plot(tt,np.imag(ct))
    plt.plot(tt,np.imag(cor(l11,g11,a,b,tt)))
    #plt.plot(tt,np.real(ctC))
    #plt.plot(tt,np.imag(ctC))
    
    plt.show()
