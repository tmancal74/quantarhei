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
                   T=100,matsubara=20)

en = 12000.0

e_units = qr.energy_units("1/cm")

with e_units:
    m = qr.Molecule("Molecule",[0.0,en])
    with qr.energy_units("1/cm"):
        cfce1 = qr.CorrelationFunction(ta,cfce_params1)
    
m.set_egcf((0,1),cfce1)   
m.set_dipole(0,1,[0.0, 1.0, 0.0])

ac = qr.AbsSpectrumCalculator(ta,m) 

with qr.energy_units("1/cm"): 
    ac.bootstrap(rwa=en)
    a1 = ac.calculate()

HH = m.get_Hamiltonian()
with qr.frequency_units("1/cm"):
    print(HH)
    a1.plot(axis=[11500,12500,0,numpy.max(a1.data)*1.1])
    
save_load = False
if save_load:

    filename = "abs_1mol_20cm_100fs_100K_m20"
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
                   reorg=10.0,
                   cortime=60.0,
                   T=100,
                   matsubara=20)
cfce_params2 = dict(ftype="OverdampedBrownian",
                   reorg=10.0,
                   cortime=60.0,
                   T=100,
                   matsubara=20)

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
(RR,ham) = AG.get_RelaxationTensor(ta,
                                   relaxation_theory="standard_Redfield")

ac2 = qr.AbsSpectrumCalculator(ta, AG, relaxation_tensor=RR)
ac3 = qr.AbsSpectrumCalculator(ta, AG)

with e_units:
    print(HH)
    ac2.bootstrap(rwa=12000)
    ac3.bootstrap(rwa=12000)
    a2 = ac2.calculate()
    a3 = ac3.calculate()

if True:   
#with e_units:
    a3.save("spectrum_a3")  
    a4 = qr.AbsSpectrum()
    a4.load("spectrum_a3.npz")

with e_units:
    #a3.plot(show=False)
    a4.plot(show=False)
    a2.plot(axis=[11500,12500,0,numpy.max(a3.data)*1.1])

save_load = True    
if save_load:
    with e_units:
        a3.save_data("abs_2mol_10cm_60fs_100K_m20",ext="dat")
    
        f = qr.DFunction()
        f.load_data("abs_2mol_10cm_60fs_100K_m20",ext="dat",axis="frequency")
        f.plot(axis=[11500,12500,0,numpy.max(a3.data)*1.1])

"""

    Absorption of a simple trimeric aggregate of two-level molecules


""" 

cfce_params1 = dict(ftype="OverdampedBrownian",
                   reorg=20.0,
                   cortime=60.0,
                   T=100,
                   matsubara=20)
cfce_params2 = dict(ftype="OverdampedBrownian",
                   reorg=20.0,
                   cortime=60.0,
                   T=100,
                   matsubara=20)

with qr.energy_units("1/cm"):
    cfce1 = qr.CorrelationFunction(ta,cfce_params1)
    cfce2 = qr.CorrelationFunction(ta,cfce_params2)
    m1 = qr.Molecule("M1",[0.0, 12100])
    m1.set_dipole(0,1,[0.0,3.0,0.0])
    m1.set_transition_environment((0,1), cfce2)
    m2 = qr.Molecule("M2",[0.0, 11800])
    m2.set_dipole(0,1,[0.0,1.0,2.0])
    m2.set_transition_environment((0,1), cfce1)
    m3 = qr.Molecule("M3",[0.0, 12500])
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
AG.set_coupling_by_dipole_dipole(epsr=3.0)

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

with qr.energy_units("1/cm"):
    (RRc,hamc) = AG.get_RelaxationTensor(ta,
                                   relaxation_theory="combined_RedfieldFoerster",
                                   time_dependent=True,
                                   coupling_cutoff=50.0)
    

ac1 = qr.AbsSpectrumCalculator(ta, AG, relaxation_tensor=RRf,
                               effective_hamiltonian=hamf)

ac2 = qr.AbsSpectrumCalculator(ta, AG, relaxation_tensor=RRr,
                               effective_hamiltonian=hamr)
ac3 = qr.AbsSpectrumCalculator(ta, AG, relaxation_tensor=RRc,
                               effective_hamiltonian=hamc)

with e_units:
    ac2.bootstrap(rwa=12000)
    ac3.bootstrap(rwa=12000)
    ac1.bootstrap(rwa=12000)
    
    
a2 = ac2.calculate()
fa = a2.axis
ACont = qr.AbsSpectrumContainer(fa)
ACont.set_spectrum(a2,tag=1)

a3 = ac3.calculate()
ACont.set_spectrum(a3,tag=2)

a1 = ac1.calculate()
ACont.set_spectrum(a1,tag=0)

with e_units:
    a1 = ACont.get_spectrum(tag=0)
    a1.plot(show=False)
    a3 = ACont.get_spectrum(tag=2)
    a3.plot(show=False)
    a1 = ACont.get_spectrum(tag=1)
    a2.plot(axis=[11000,13000,0,numpy.max(a2.data)*1.1])

    

save_load = True   
if save_load:
    with e_units:
        a2.save_data("abs_3mol_20cm_60fs_100K_m20",ext="dat")
    
        f = qr.DFunction()
        f.load_data("abs_3mol_20cm_60fs_100K_m20",ext="dat",axis="frequency")
        f.plot(axis=[11000,13000,0,numpy.max(a2.data)*1.1])
        
        ACont.save("container.hdf5")
        
        Cont2 = qr.AbsSpectrumContainer()
        Cont2.load("container.hdf5")
        
        a6 = Cont2.get_spectrum(1)
        a6.plot(axis=[11000,13000,0,numpy.max(a2.data)*1.1])
        
#        a5 = qr.AbsSpectrum()
#        a5.load("container_2.npz")
#        a5.plot(axis=[11000,13000,0,numpy.max(a2.data)*1.1])