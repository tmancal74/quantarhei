# -*- coding: utf-8 -*-
_show_plots_ = False

print("""
***********************************************************      
*
*
*          Integrodifferential Propagator Demo
*
*
***********************************************************
""")

import time
import matplotlib.pyplot as plt
import quantarhei as qr
import numpy

##########################################################
#
#   DEMONSTRATION OF THE FFT METHOD
#
##########################################################

Om = 2.0*numpy.pi/30.0
g1 = 1.0/100.0
a0 = 1.0

#
# This is how to problem is solved
#
tt = qr.TimeAxis(0.0, 300, 1.0)
gg = 5.0/tt.data[tt.length-1]
#
om = (2.0*numpy.pi)*numpy.fft.fftfreq(tt.length, tt.step)
aom = numpy.zeros(tt.length, dtype=qr.COMPLEX)
for ii in range(tt.length):
    aom[ii] = a0/(-1j*om[ii] + 1j*Om + g1 + gg)
    

at = numpy.fft.fft(aom)*numpy.exp(gg*tt.data)/tt.data[tt.length-1]

#
#  Expected solution
#
atsol = a0*numpy.exp(-1j*Om*tt.data)*numpy.exp(-g1*tt.data)
at[0] = a0

#
# Plot the results
#
if _show_plots_:
    tshow = tt.data[:tt.length//2]
    atshow = at[:tt.length//2]
    atsolshow = atsol[:tt.length//2]
    
    plt.plot(tshow, numpy.real(atshow))
    plt.plot(tshow, numpy.real(atsolshow))
    plt.show()


###############################################################
#  
#   SOLVING LIOUVILLE-VON NEUMANN EQUATION BY FFT METHOD
#
###############################################################
def test_kernel(timeaxis, ham, operators, rates, ctime):
    """Returns a simple kernel for tests
    
    """
    
    dim = ham.dim
    Nt = timeaxis.length
    
    MM = numpy.zeros((Nt, dim, dim, dim, dim), dtype=qr.COMPLEX)
    gamma = 1.0/ctime

    if dim == 2:
        
        sys_ops = operators
        sbi = qr.qm.SystemBathInteraction(sys_operators=sys_ops, rates=rates)
        lbf = qr.qm.LindbladForm(ham, sbi, as_operators=False)
        #return lbf.data
        
        for ti in range(Nt):
            tm = timeaxis.data[ti]
            MM[ti,:,:,:,:] = -lbf.data*numpy.exp(-gamma*tm)
            
        return MM

from quantarhei.qm.liouvillespace.integrodiff.integrodiff \
     import IntegrodiffPropagator

timea = qr.TimeAxis(0.0, 200, 0.5)
Nt = timea.length
ham = qr.Hamiltonian(data=[[0.0, 0.1], 
                           [0.1, 0.01]])

#
# Propagation without relaxation
#

ip1 = IntegrodiffPropagator(timea, ham, 
                            timefac=3, decay_fraction=2.0)
np = qr.ReducedDensityMatrixPropagator(timea, ham)

rhoi = qr.ReducedDensityMatrix(data=[[0.0, 0.0],[0.0, 1.0]])

t1 = time.time()
rhot_i = ip1.propagate(rhoi)
t2 = time.time()
print("Propagated in frequency domain in:", t2-t1)

t1 = time.time()
rhot_n = np.propagate(rhoi)
t2 = time.time()
print("Propagated in time domain in:", t2-t1)

if _show_plots_:
    plt.plot(timea.data, numpy.real(rhot_n.data[:,0,0]),"-b")
    plt.plot(timea.data, numpy.real(rhot_n.data[:,1,1]),"-r")
    plt.plot(timea.data, numpy.real(rhot_i.data[:,0,0]),"--g")
    plt.plot(timea.data, numpy.real(rhot_i.data[:,1,1]),"--k")
    #plt.axis([0,10,0,1])
    plt.show()
    #rhot_i.plot(coherences=False, show=True)
    #rhot_n.plot(coherences=False, show=True)


#
#  Propagation with relaxation kernel
#
    
K01 = qr.qm.ProjectionOperator(0,1,ham.dim)
K10 = qr.qm.ProjectionOperator(1,0,ham.dim)

sys_ops = [K01, K10]
rates = [1.0/30.0, 1.0/20.0]

ker = test_kernel(timea, ham, sys_ops, rates, ctime=20.0)

# time domain propagator
ip2 = IntegrodiffPropagator(timea, ham, kernel=ker, 
                            fft=False, cutoff_time=80)

# frequency domain propagator
ip3 = IntegrodiffPropagator(timea, ham, kernel=ker, 
                            fft=True, timefac=3, decay_fraction=2.0)

t1 = time.time()
rhot_k = ip2.propagate(rhoi)
t2 = time.time()
print("Propagated in time domain in:", t2-t1)

t1 = time.time()
rhot_k3 = ip3.propagate(rhoi)
t2 = time.time()
print("Propagated in frequency domain in:", t2-t1)

#plotit = True
if _show_plots_:
    plt.plot(timea.data, numpy.real(rhot_k.data[:,0,0]),"-b")
    plt.plot(timea.data, numpy.real(rhot_k.data[:,1,1]),"-r")
    plt.plot(timea.data, numpy.real(rhot_k3.data[:,0,0]),"--b")
    plt.plot(timea.data, numpy.real(rhot_k3.data[:,1,1]),"--r")

