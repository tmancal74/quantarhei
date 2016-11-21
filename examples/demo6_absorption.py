# -*- coding: utf-8 -*-

import quantarhei as qr

import numpy

ta = qr.TimeAxis(0.0, 1000, 1.0)

"""

    Absorption of a monomeric two-level molecule


"""
cfce_params1 = dict(ftype="OverdampedBrownian",
                   reorg=20.0,
                   cortime=100.0,
                   T=300,matsubara=20)

en = 12000.0

e_units = qr.energy_units("1/cm")

with e_units:
    m = qr.Molecule("Molecule",[0.0,en])
    with qr.energy_units("1/cm"):
        cfce1 = qr.CorrelationFunction(ta,cfce_params1)
    
m.set_egcf((0,1),cfce1)   
m.set_dipole(0,1,[0.0, 1.0, 0.0])

a1 = qr.AbsSpect(ta,m) 

with qr.energy_units("1/cm"):    
    a1.calculate(rwa=en)

HH = m.get_Hamiltonian()
with qr.frequency_units("1/cm"):
    print(HH)
    a1.plot(axis=[11500,12500,0,numpy.max(a1.data)*1.1])
    
save_load = False
if save_load:

    filename = "abs_1mol_20cm_100fs_300K_m20"
    with qr.frequency_units("1/cm"):
        a1.save(filename,ext="dat")
    
    with qr.frequency_units("1/cm"):
    
        f = qr.DFunction()
    
        f.load(filename,ext="dat",axis="frequency")
    
        f.plot()
    
"""

    Absorption of a simple dimer aggregate of two-level molecules


"""

ta = qr.TimeAxis(0.0, 2000, 5.0)

cfce_params1 = dict(ftype="OverdampedBrownian",
                   reorg=30.0,
                   cortime=60.0,
                   T=300,
                   matsubara=10)
cfce_params2 = dict(ftype="OverdampedBrownian",
                   reorg=30.0,
                   cortime=60.0,
                   T=300,
                   matsubara=10)

with qr.energy_units("1/cm"):
    cfce1 = qr.CorrelationFunction(ta,cfce_params1)
    cfce2 = qr.CorrelationFunction(ta,cfce_params2)
    m1 = qr.Molecule("M1",[0.0, 12100])
    m1.set_dipole(0,1,[0.0,3.0,0.0])
    m1.set_transition_environment((0,1),cfce1)
    m1.position = [0.0,0.0,0.0]
    m2 = qr.Molecule("M1",[0.0, 12000])
    m2.set_dipole(0,1,[0.0,1.0,1.0])
    m2.position = [5.0,0.0,0.0]
    m2.set_transition_environment((0,1),cfce2)    

# create an aggregate
AG = qr.Aggregate("TestAggregate")
#AG.set_egcf_matrix(cm)

# fill the cluster with monomers
AG.add_Molecule(m1)
AG.add_Molecule(m2)

# setting coupling by dipole-dipole formula
AG.set_coupling_by_dipole_dipole()

AG.build()

HH = AG.get_Hamiltonian()

a2 = qr.AbsSpect(ta, AG)

with e_units:
    print(HH)
    a2.calculate(rwa=12000)
    a2.plot(axis=[11500,12500,0,numpy.max(a2.data)*1.1])

save_load = False    
if save_load:
    with e_units:
        a2.save("abs_2mol_30cm_60fs_300K_m10",ext="dat")
    
        f = qr.DFunction()
        f.load("abs_2mol_30cm_60fs_300K_m10",ext="dat",axis="frequency")
        f.plot()
  

"""

    Absorption of a simple trimeric aggregate of two-level molecules


""" 

cfce_params1 = dict(ftype="OverdampedBrownian",
                   reorg=20.0,
                   cortime=60.0,
                   T=300,
                   matsubara=20)
cfce_params2 = dict(ftype="OverdampedBrownian",
                   reorg=20.0,
                   cortime=60.0,
                   T=300,
                   matsubara=20)

with qr.energy_units("1/cm"):
    cfce1 = qr.CorrelationFunction(ta,cfce_params1)
    cfce2 = qr.CorrelationFunction(ta,cfce_params2)
    m1 = qr.Molecule("M1",[0.0, 12100])
    m1.set_dipole(0,1,[0.0,3.0,0.0])
    m1.set_transition_environment((0,1), cfce2)
    m2 = qr.Molecule("M1",[0.0, 12000])
    m2.set_dipole(0,1,[0.0,1.0,2.0])
    m2.set_transition_environment((0,1), cfce1)
    m3 = qr.Molecule("M1",[0.0, 12000])
    m3.set_dipole(0,1,[0.0,1.0,1.0])    
    m3.set_transition_environment((0,1), cfce2)    


m1.position = [0.0,0.0,0.0]
m2.position = [5.0,0.0,0.0]
m3.position = [0.0,5.0,0.0]

# create an aggregate
AG = qr.Aggregate("TestAggregate")


# fill the cluster with monomers
AG.add_Molecule(m1)
AG.add_Molecule(m2)
AG.add_Molecule(m3)

# setting coupling by dipole-dipole formula
AG.set_coupling_by_dipole_dipole()

AG.build()

HH = AG.get_Hamiltonian()
a2 = qr.AbsSpect(ta,AG)

with e_units:
    print(HH)
    a2.calculate(rwa=12000)
    a2.plot(axis=[11500,12500,0,numpy.max(a2.data)*1.1])

save_load = False    
if save_load:
    with e_units:
        a2.save("abs_3mol_20cm_60fs_300K_m20",ext="dat")
    
        f = qr.DFunction()
        f.load("abs_3mol_20cm_60fs_300K_m20",ext="dat",axis="frequency")
        f.plot()
        
