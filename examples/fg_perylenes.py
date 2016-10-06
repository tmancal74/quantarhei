# -*- coding: utf-8 -*-

import numpy
import matplotlib.pylab as plt
import time as tm

from quantarhei import TimeAxis
from quantarhei.qm import ReducedDensityMatrix
from quantarhei.qm import RedfieldRelaxationTensor, RedfieldRateMatrix
##from cu.oqs.liouvillespace import TDRedfieldRelaxationTensor                                 
from quantarhei.qm import RDMPropagator

from quantarhei import CorrelationFunction
from quantarhei.qm.corfunctions import CorrelationFunctionMatrix
from quantarhei import Aggregate
from quantarhei import Molecule
from quantarhei import AbsSpect
from quantarhei import energy_units

from quantarhei.utils import normalize2 
from quantarhei.core.units import cm2int, kB_intK

""" Create monomers """
# state energies
en1 = [0.0,12100] #*cm2int]
en2 = [0.0,12000] #*cm2int]

# transition dipole moment
""" d0 = [0.0,0.0,0.0] """
d1 = normalize2([0.0,1.0,0.0],1.0)

with energy_units("1/cm"):
    m1 = Molecule("M1",en1) #,[d0,d1])
    m1.set_dipole(0,1,d1)
    m2 = Molecule("M2",en1) #,[d0,d1])
    m2.set_dipole(0,1,d1)
    m3 = Molecule("M3",en2) #,[d0,d1])
    m3.set_dipole(0,1,d1)
    m4 = Molecule("M4",en1) #,[d0,d1])
    m4.set_dipole(0,1,d1)
    m5 = Molecule("M5",en1) #,[d0,d1])
    m5.set_dipole(0,1,d1)
    m6 = Molecule("M6",en1) #,[d0,d1])
    m6.set_dipole(0,1,d1)
    m7 = Molecule("M7",en1) #,[d0,d1])
    m7.set_dipole(0,1,d1)
    m8 = Molecule("M8",en1) #,[d0,d1])
    m8.set_dipole(0,1,d1)
    m9 = Molecule("M9",en1) #,[d0,d1])
    m9.set_dipole(0,1,d1)
    m10 = Molecule("M10",en1) #,[d0,d1])
    m10.set_dipole(0,1,d1)
    m11 = Molecule("M11",en1) #,[d0,d1])
    m11.set_dipole(0,1,d1)
    m12 = Molecule("M12",en1) #,[d0,d1])
    m12.set_dipole(0,1,d1)
    
#possitions in space
a = 0.1   # C-C distance
b = 2 * a * numpy.sin(numpy.pi/3.)
dist1 = 9 * a
dist2 = 4 * b
r1 = [0,0,0]
r2 = [0,dist1, 0]
r3 = [0, 2 * dist1,0]
r4 = [0,3 * dist1,0]
r5 = [0,4 * dist1,0]
r6 = [0,5 * dist1, 0]
r7 = [dist2, 0, 0]
r8 = [dist2, dist1, 0]
r9 = [dist2, 2 * dist1, 0]
r10 = [dist2, 3 * dist1, 0]
r11 = [dist2, 4 * dist1, 0]
r12 = [dist2, 5 * dist1, 0]

m1.position = r1
m2.position = r2
m3.position = r3
m4.position = r4
m5.position = r5
m6.position = r6
m7.position = r7
m8.position = r8
m9.position = r9
m10.position = r10
m11.position = r11
m12.position = r12

""" Correlation functions """

temperature = 300.0 # in Kelvins
# Time axis on which everything is calculated
time = TimeAxis(0.0, 2000, 1.0) # in fs

cfce_params1 = dict(ftype="OverdampedBrownian",
                   reorg=30.0,
                   cortime=60.0,
                   T=temperature)
cfce_params2 = dict(ftype="OverdampedBrownian",
                   reorg=30.0,
                   cortime=60.0,
                   T=temperature)

with energy_units("1/cm"):
    cfce1 = CorrelationFunction(time,cfce_params1)
    cfce2 = CorrelationFunction(time,cfce_params2)

explicit_mapping = False

if explicit_mapping:
    
    cm = CorrelationFunctionMatrix(time,12)
    i = cm.set_correlation_function(cfce1, [(1,1),(2,2)]) #,
                            #iof=1)
    cm.set_correlation_function(cfce1, [(4,4),(6,6),(7,7),(11,11)]) 
    cm.set_correlation_function(cfce2, [(0,0),(3,3),(5,5),(8,8),(9,9),(10,10)]) #,
                            #iof=2)
    # Mapping of the correlation functions on the transitions in monomers 
    m1.set_egcf_mapping((0,1),cm,0)
    m2.set_egcf_mapping((0,1),cm,1)
    m3.set_egcf_mapping((0,1),cm,2)
    m4.set_egcf_mapping((0,1),cm,3)
    m5.set_egcf_mapping((0,1),cm,4)
    m6.set_egcf_mapping((0,1),cm,5)
    m7.set_egcf_mapping((0,1),cm,6)
    m8.set_egcf_mapping((0,1),cm,7)
    m9.set_egcf_mapping((0,1),cm,8)
    m10.set_egcf_mapping((0,1),cm,9)
    m11.set_egcf_mapping((0,1),cm,10)
    m12.set_egcf_mapping((0,1),cm,11)
    
else:
    
    m1.set_transition_environment((0,1),cfce2)
    m2.set_transition_environment((0,1),cfce1)
    m3.set_transition_environment((0,1),cfce1)
    m4.set_transition_environment((0,1),cfce2)
    m5.set_transition_environment((0,1),cfce1)
    m6.set_transition_environment((0,1),cfce2)
    m7.set_transition_environment((0,1),cfce1)
    m8.set_transition_environment((0,1),cfce1)
    m9.set_transition_environment((0,1),cfce2)
    m10.set_transition_environment((0,1),cfce2)
    m11.set_transition_environment((0,1),cfce2)
    m12.set_transition_environment((0,1),cfce1)    
    

# create an aggregate
AG = Aggregate("TestAggregate")

if explicit_mapping:
    AG.set_egcf_matrix(cm)

# fill the cluster with monomers
AG.add_Molecule(m1)
AG.add_Molecule(m2)
AG.add_Molecule(m3)
AG.add_Molecule(m4)
AG.add_Molecule(m5)
AG.add_Molecule(m6)
AG.add_Molecule(m7)
AG.add_Molecule(m8)
AG.add_Molecule(m9)
AG.add_Molecule(m10)
AG.add_Molecule(m11)
AG.add_Molecule(m12)


print(AG._has_egcf_matrix)
print(AG._has_system_bath_interaction)


# setting coupling by dipole-dipole formula
AG.set_coupling_by_dipole_dipole(prefac=0.0147520827152)

print (AG.resonance_coupling[1,2]/cm2int)

# building aggregate
print("Building aggregate... ")
start = tm.time()
AG.build()
print("...done in ",tm.time()-start)

with energy_units("1/cm"):
    anth = Molecule("Anth",[0.0,11380]) #*cm2int],[d0,d1])
anth.set_dipole(0,1,d1)
anth.set_egcf((0,1),cfce1)

# FIXME: use with construct
aAnth = AbsSpect(time,anth)
aAnth.calculate(12000*cm2int)
aAnth.normalize2()

a1 = AbsSpect(time,AG)
a1.calculate(12000*cm2int)
a1.normalize2()
a2 = AbsSpect(time,m1)
a2.calculate(12000*cm2int)
a2.normalize2()

plt.plot(a1.frequency/cm2int,a1.data,'-r')
plt.plot(a2.frequency/cm2int,a2.data,'-g')
plt.plot(aAnth.axis.data/cm2int,aAnth.data,'-c')

p = plt.gca()
p.set_autoscale_on(False)
p.axis([10000,15000,0,1.1])
en1st = str(en1[1]/cm2int)
plt.title('L12per_pbl --- en1 = ' + en1st)

#plt.title('aggregate: 12 perylenes, linear dense zigzag ')
#plt.savefig('L12per_pbl_spec.pdf')

plt.show()

with energy_units("1/cm"):
    aAnth.plot(axis=[10000,15000,0.0,1.1])


# Redfield Tensor 
print("Calculating relaxation tensor ...")
start = tm.time()
# System bath interaction
#sbi = AG.get_system_bath_coupling()
sbi = AG.get_SystemBathInteraction()
# Hamiltonian 
HH = AG.get_Hamiltonian()

# Relaxation tensor
RT = RedfieldRelaxationTensor(HH,sbi)
#RT = TDRedfieldRelaxationTensor(HH,sbi)
RT.secularize()
print("... done in ",tm.time()-start)

#for n in range(4):
#    print(numpy.real(RT.data[n,n,5,5]))
#
#RR = RedfieldRateMatrix(HH,sbi)
#for n in range(4):
#    print(numpy.real(RR.data[n,5]))    

"""

    IMPLEMENT AG.get_initial_RDM()

"""
# Sets an initial population of sites according to the transition dipole
#AG.set_impulsive_population()
rho = ReducedDensityMatrix(dim=HH.dim)
#rho.data = AG.rho0
#rho.data[0,0] = 1.0
rho.data[HH.dim-1,HH.dim-1] = 1.0
#rho.data[1,0] = 0.5
#rho.data[0,1] = 0.5
rho.normalize2()


print("Propagating density matrix ...")
start = tm.time()
# Convert to exciton basis in which RT was calculated 
ss = HH.diagonalize()  # Relaxation tensors are defined in the eigenstate basis
prop = RDMPropagator(time,HH,RTensor=RT)
#prop.setDtRefinement(10)

pr2 = prop.propagate(rho,method="short-exp")
print("... done in ", tm.time()-start)
#pr2.transform(ss.T)

""" 

Plotting the time evolution

"""
print("\nEvolution of populations")
pr2.plot(coherences=False,how='-')

"""

    IMPLEMENT AG.get_cannonical_equilibrium_RDM()
              AG.get_cannonical_populations()
              
"""
""" 

Calculate termal cannonical distribution for comparison

"""
pop = numpy.zeros((time.length,13),dtype=numpy.float64)
tpop = 0.0
for i in range(1,HH.dim):
    pop[:,i] = numpy.exp(-HH.data[i,i]/(temperature*kB_intK))
    tpop += pop[1,i]

prtot = 0.0
for i in range(1,HH.dim):
    prtot += pr2.data[time.length-1,i,i]
    
pop[:,:] = numpy.real((pop[:,:]/tpop)*prtot)

# plot the termal distrubution
plt.title('L12per_pbl')
plt.plot(time.data,pop[:,1],'--r')
plt.plot(time.data,pop[:,2],'--b')
plt.plot(time.data,pop[:,3],'--g')
plt.plot(time.data,pop[:,4],'--m')
plt.plot(time.data,pop[:,5],'--y')
plt.plot(time.data,pop[:,6],'--c')
plt.plot(time.data,pop[:,7],'--k')
plt.plot(time.data,pop[:,8],'--',color='deeppink')
plt.plot(time.data,pop[:,9],'--',color='darkviolet')
plt.plot(time.data,pop[:,10],'--',color='dimgrey')
plt.plot(time.data,pop[:,11],'--',color='firebrick')
plt.plot(time.data,pop[:,12],'--',color='gold')
plt.ylim((0.0,1.0)) 

#plt.savefig('L12per_pbl.pdf')
