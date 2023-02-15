# -*- coding: utf-8 -*-
_show_plots_ = False

import time

import numpy

import quantarhei as qr
from quantarhei.qm.liouvillespace.integrodiff.integrodiff \
     import IntegrodiffPropagator

print("")
print("***********************************************************")
print("*                                                         *")
print("*          Quantarhei's HEOM implementation demo          *")
print("*                                                         *")
print("***********************************************************")
###############################################################################
#
#   Model system definition
#
###############################################################################

#   Three molecules
with qr.energy_units("1/cm"):
    m1 = qr.Molecule([0.0, 10100.0])
    m2 = qr.Molecule([0.0, 10300.0])
    m3 = qr.Molecule([0.0, 10000.0])

#   Aggregate is built from the molecules    
agg = qr.Aggregate([m1, m2, m3])

#   Couplings between them are set
with qr.energy_units("1/cm"):
    agg.set_resonance_coupling(0,1,80.0)
    agg.set_resonance_coupling(0,2,100.0)

#   Interaction with the bath is set through bath correlation functions
timea = qr.TimeAxis(0.0, 500, 1.0)
cpar1 = dict(ftype="OverdampedBrownian-HighTemperature", reorg=50,
            cortime=50, T=300)
cpar2 = dict(ftype="OverdampedBrownian-HighTemperature", reorg=50,
            cortime=50, T=300)

with qr.energy_units("1/cm"):
    cfce1 = qr.CorrelationFunction(timea, cpar1)
    cfce2 = qr.CorrelationFunction(timea, cpar2)
    
m1.set_transition_environment((0, 1), cfce1)
m2.set_transition_environment((0, 1), cfce1)
m3.set_transition_environment((0, 1), cfce2)
    
#    Aggregate is built
agg.build()

###############################################################################
#
#    Definition of the hierarchy
#
###############################################################################

#    Hamiltonian and the system-bath interaction operator is needed to 
#    define the Kubo-Tanimura hierarchy
ham = agg.get_Hamiltonian()
sbi = agg.get_SystemBathInteraction() 

#    We define the hierarchy
#Hy3 = qr.KTHierarchy(ham, sbi, 3)
#Hy4 = qr.KTHierarchy(ham, sbi, 4)
#Hy5 = qr.KTHierarchy(ham, sbi, 5)
Hy6 = qr.KTHierarchy(ham, sbi, 3)
print("Size of hierarchy of depth",Hy6.depth,"is",Hy6.hsize)
Hy7 = qr.KTHierarchy(ham, sbi, 4)
print("Size of hierarchy of depth",Hy7.depth,"is",Hy7.hsize)
# testing generation of hierarchy indices
#print(Hy.generate_indices(4, level=4))
#
#raise Exception()


###############################################################################
#
#    Propagation of the HEOM
#
###############################################################################

#   Initial density matrix
rhoi = qr.ReducedDensityMatrix(dim=ham.dim)

with qr.eigenbasis_of(ham):
    rhoi.data[2,2] = 0.8
    rhoi.data[1,1] = 0.1
    rhoi.data[3,3] = 0.1

#print(rhoi)
    
#   Definition of the HEOM propagator
#kprop3 = qr.KTHierarchyPropagator(timea, Hy3)
#kprop4 = qr.KTHierarchyPropagator(timea, Hy4)
#kprop5 = qr.KTHierarchyPropagator(timea, Hy5)
kprop6 = qr.KTHierarchyPropagator(timea, Hy6)
kprop7 = qr.KTHierarchyPropagator(timea, Hy7)

#   Propagation of the hierarchy and saving the density operator
t1 = time.time()
#rhot3 = kprop3.propagate(rhoi, report_hierarchy=False, free_hierarchy=False)
#rhot4 = kprop4.propagate(rhoi, report_hierarchy=False, free_hierarchy=False)
#rhot5 = kprop5.propagate(rhoi, report_hierarchy=False, free_hierarchy=False)
rhot6 = kprop6.propagate(rhoi, report_hierarchy=False, free_hierarchy=False)
t2 = time.time()
print("Propagated in", t2-t1,"s")
t1 = time.time()
rhot7 = kprop7.propagate(rhoi, report_hierarchy=False, free_hierarchy=False)
t2 = time.time()
print("Propagated in", t2-t1,"s")

###############################################################################
#
#    Graphical output of the results
#
###############################################################################

if _show_plots_:
    
    import matplotlib.pyplot as plt
    N = timea.length
    with qr.eigenbasis_of(ham):
    #    plt.plot(timea.data[0:N], rhot3.data[0:N,1,1],"-b")
    #    plt.plot(timea.data[0:N], rhot3.data[0:N,2,2],"-r")
    #    plt.plot(timea.data[0:N], rhot3.data[0:N,3,3],"-k")
    #    plt.plot(timea.data[0:N], rhot4.data[0:N,2,2],"-r")
    #    plt.plot(timea.data[0:N], rhot4.data[0:N,1,1],"-b")
    #    plt.plot(timea.data[0:N], rhot4.data[0:N,3,3],"-k")
    #    plt.plot(timea.data[0:N], rhot5.data[0:N,1,1],"-b")
    #    plt.plot(timea.data[0:N], rhot5.data[0:N,2,2],"-r")
    #    plt.plot(timea.data[0:N], rhot5.data[0:N,3,3],"-k")
        plt.plot(timea.data[0:N], rhot6.data[0:N,0,0])
        plt.plot(timea.data[0:N], rhot6.data[0:N,1,3],"-b")
        plt.plot(timea.data[0:N], rhot6.data[0:N,2,3],"-r")
        plt.plot(timea.data[0:N], rhot6.data[0:N,1,2],"-k")    
        plt.plot(timea.data[0:N], rhot7.data[0:N,1,3],"--b")
        plt.plot(timea.data[0:N], rhot7.data[0:N,2,3],"--r")
        plt.plot(timea.data[0:N], rhot7.data[0:N,1,2],"--k")    
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,1], "-k")
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,2], "-k")
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,3], "-b")
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,4], "-b")
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,5], "-b")
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,6], "-r")
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,7], "-r")
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,8], "-r")
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,9], "-r")
        #plt.plot(timea.data[0:N], Hy.hpop[0:N,10], "-g")
    plt.show()

print("Kernel generation")

ker = Hy6.get_kernel(timea)

ip8 = IntegrodiffPropagator(timea, ham, kernel=ker,
                            fft=True, timefac=3, decay_fraction=2.0)
   #fft=False) #, cutoff_time=100)

rhot8 = ip8.propagate(rhoi)

trc = numpy.zeros(timea.length, dtype=qr.REAL)
for ti in range(timea.length):
    trc[ti] = numpy.real(numpy.trace(rhot8.data[ti,:,:]))

if _show_plots_:
    N = timea.length
    with qr.eigenbasis_of(ham):
        #plt.plot(timea.data[0:N], rhot8.data[0:N,0,0])
        #plt.plot(timea.data[0:N], trc[0:N],"-m")
        plt.plot(timea.data[0:N], ker[0:N,1,1,1,1],"-m")
        plt.plot(timea.data[0:N], ker[0:N,1,2,1,2],"-m")
        plt.plot(timea.data[0:N], ker[0:N,2,2,2,2],"-m")
        plt.show()
        plt.plot(timea.data[0:N], rhot8.data[0:N,1,1],"-b")
        plt.plot(timea.data[0:N], rhot8.data[0:N,2,2],"-r")
        plt.plot(timea.data[0:N], rhot8.data[0:N,1,2],"-k") 
        plt.plot(timea.data[0:N], rhot6.data[0:N,1,1],"--b")
        plt.plot(timea.data[0:N], rhot6.data[0:N,2,2],"--r")
        plt.plot(timea.data[0:N], rhot6.data[0:N,1,2],"--k")
    plt.show()

print("")
print("***********************************************************")
print("*                                                         *")
print("*               Demo finished successfully                 *")
print("*                                                         *")
print("***********************************************************")

