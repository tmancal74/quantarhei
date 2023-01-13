# -*- coding: utf-8 -*-
_show_plots_ = False

import time
import quantarhei as qr
from quantarhei.qm.liouvillespace.integrodiff.integrodiff \
     import IntegrodiffPropagator

print("")
print("***************************************************************")
print("*                                                             *")
print("*  Quantarhei's HEOM kernel and rate extraction demo          *")
print("*                                                             *")
print("***************************************************************")
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
timea = qr.TimeAxis(0.0, 1000, 1.0)
cpar1 = dict(ftype="OverdampedBrownian-HighTemperature", reorg=70,
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
