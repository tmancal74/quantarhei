# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import quantarhei as qr
import numpy

plotit = False

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
if plotit:
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

from quantarhei.qm.liouvillespace.integrodiff.integrodiff \
     import IntegrodiffPropagator

timea = qr.TimeAxis(0.0, 400, 0.1)
Nt = timea.length
ham = qr.Hamiltonian(data=[[0.0, 0.1], [0.1, 0.01]])

ip = IntegrodiffPropagator(timea, ham, timefac=3, decay_fraction=2.0)

ker = ip._test_kernel()

ip2 = IntegrodiffPropagator(timea, ham, kernel=ker, 
                            fft=False)
ip3 = IntegrodiffPropagator(timea, ham, kernel=ker, 
                            fft=True, timefac=3, decay_fraction=2.0)

np = qr.ReducedDensityMatrixPropagator(timea, ham)

rhoi = qr.ReducedDensityMatrix(data=[[0.0, 0.0],[0.0, 1.0]])

rhot_i = ip.propagate(rhoi)
rhot_n = np.propagate(rhoi)

if plotit:
    plt.plot(timea.data, numpy.real(rhot_n.data[:,0,0]),"-b")
    plt.plot(timea.data, numpy.real(rhot_n.data[:,1,1]),"-r")
    plt.plot(timea.data, numpy.real(rhot_i.data[:,0,0]),"--g")
    plt.plot(timea.data, numpy.real(rhot_i.data[:,1,1]),"--k")
    #plt.axis([0,10,0,1])
    plt.show()
    #rhot_i.plot(coherences=False, show=True)
    #rhot_n.plot(coherences=False, show=True)



rhot_k = ip2.propagate(rhoi)

rhot_k3 = ip3.propagate(rhoi)

#plotit = True
if plotit:
    plt.plot(timea.data, numpy.real(rhot_k.data[:,0,0]),"-b")
    plt.plot(timea.data, numpy.real(rhot_k.data[:,1,1]),"-r")
    plt.plot(timea.data, numpy.real(rhot_k3.data[:,0,0]),"--b")
    plt.plot(timea.data, numpy.real(rhot_k3.data[:,1,1]),"--r")
