# -*- coding: utf-8 -*-
import numpy

from quantarhei import *

e_units = energy_units("1/cm")

ta = TimeAxis(0.0, 1000, 2.0)

cfce_params1 = dict(ftype="OverdampedBrownian",
                   reorg=20.0,
                   cortime=60.0,
                   T=300,
                   matsubara=20)

en = 12000.0

with e_units:
    cfce1 = CorrelationFunction(ta,cfce_params1)
    cfce2 = CorrelationFunction(ta,cfce_params2)
    m1 = Molecule("M1",[0.0, 12100])
    m1.set_dipole(0,1,[0.0,3.0,0.0])
    m1.set_transition_environment((0,1), cfce1)
    m2 = Molecule("M2",[0.0, 11800])
    m2.set_dipole(0,1,[0.0,1.0,2.0])
    m2.set_transition_environment((0,1), cfce1)
    m3 = Molecule("M3",[0.0, 12500])
    m3.set_dipole(0,1,[0.0,1.0,1.0])    
    m3.set_transition_environment((0,1), cfce1)    


m1.position = [0.0,0.0,0.0]
m2.position = [7.0,0.0,0.0]
m3.position = [0.0,7.0,0.0]

# create an aggregate
AG = Aggregate("TestAggregate")


# fill the cluster with monomers
AG.add_Molecule(m1)
AG.add_Molecule(m2)
AG.add_Molecule(m3)

# setting coupling by dipole-dipole formula
#AG.set_coupling_by_dipole_dipole(epsr=1.0)
with energy_units("1/cm"):
    AG.set_coupling_matrix([[0.0, 80, 50],
                            [80, 0.0, 80],
                            [50, 80,  0.0]])

AG.build()

HH = AG.get_Hamiltonian()
with e_units:
    print(HH)
    
(RRf,hamf) = AG.get_RelaxationTensor(ta,
                                   relaxation_theory="standard_Foerster",
                                   time_dependent=True)
(RRr,hamr) = AG.get_RelaxationTensor(ta,
                                   relaxation_theory="standard_Redfield",
                                   time_dependent=True)

with e_units:
    (RRc,hamc) = AG.get_RelaxationTensor(ta,
                                   relaxation_theory="combined_RedfieldFoerster",
                                   time_dependent=True,
                                   coupling_cutoff=30.0)
    

a1 = AbsSpect(ta, AG, relaxation_tensor=RRf, effective_hamiltonian=hamf)
a2 = AbsSpect(ta, AG, relaxation_tensor=RRr, effective_hamiltonian=hamr)
a3 = AbsSpect(ta, AG, relaxation_tensor=RRc, effective_hamiltonian=hamc)

with e_units:
    a2.calculate(rwa=12100)
    a3.calculate(rwa=12100)
    a1.calculate(rwa=12100)
    a1.plot(show=False)
    a3.plot(show=False)
    a2.plot(axis=[11000,13000,0,numpy.max(a2.data)*1.1])

