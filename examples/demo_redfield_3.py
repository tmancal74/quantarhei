# -*- coding: utf-8 -*-
"""

    Propagation with Redfield theory using Aggregate object




"""

import os.path
import pickle
import numpy

from quantarhei import *

print("""
*******************************************************************************
*                                                                             *
*                         Redfield Theory Demo                                *
*                                                                             *                  
*******************************************************************************
""")

time = TimeAxis(0.0, 5000, 1.0)

if not os.path.exists("aggregate.hdf5"):
    
    with energy_units("1/cm"):
    
        m1 = Molecule(name="Mol 1", elenergies=[0.0, 10100.0])
        m2 = Molecule(name="Mol 2", elenergies=[0.0, 10050.0])
        m3 = Molecule(name="Mol 3", elenergies=[0.0, 10000.0])
        
        m1.position = [0.0, 0.0, 0.0]
        m2.position = [10.0, 0.0, 0.0]
        m3.position = [5.0, 5.0, 0.0]
        m1.set_dipole(0,1,[5.8, 0.0, 0.0])
        m2.set_dipole(0,1,[5.8, 0.0, 0.0])
        m3.set_dipole(0,1,[numpy.sqrt(12.0), 0.0, 0.0])
    
        agg = Aggregate(name="Trimer")
        agg.add_Molecule(m1)
        agg.add_Molecule(m2)
        agg.add_Molecule(m3)
        
        agg.set_coupling_by_dipole_dipole(epsr=1.92921)
    
        params = dict(ftype="OverdampedBrownian", reorg=20, cortime=300, T=100)
        cf = CorrelationFunction(time, params)
    
    m1.set_transition_environment((0,1), cf)
    m2.set_transition_environment((0,1), cf)
    m3.set_transition_environment((0,1), cf)
    
    agg.save("unbuilt.hdf5")
    fid = open("unbuilt.pkl","wb")
    pickle.dump(agg, fid)
    fid.close()  
    
    agg.build()
    
    print("Number: ", isinstance(agg._built,bool))
    print("bool  : ", isinstance(agg._built,bool))
    print(agg._built)
    

else:
    
    agg = Aggregate()
    agg.load("aggregate.hdf5")
    


#
# Aggregate object can return a propagator
#
prop = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="standard_Redfield",
                           time_dependent=False,
                           secular_relaxation=True)

prop_nonsec = agg.get_ReducedDensityMatrixPropagator(time,
                           time_dependent=False,
                           relaxation_theory="standard_Redfield")
         
ham = agg.get_Hamiltonian()                  
with eigenbasis_of(ham):
#if True:

    rho_i1 = ReducedDensityMatrix(dim=4, name="Initial DM")
    rho_i1.data[3,3] = 1.0   
    
    rho_t1 = prop.propagate(rho_i1,
             name="Redfield evolution from aggregate, non-secular")
    rho_t1.plot(coherences=False)
    
    rho_i2 = ReducedDensityMatrix(dim=4, name="Initial DM")
    rho_i2.data[3,3] = 1.0   
    
#    rho_t2 = prop_nonsec.propagate(rho_i2,
#             name="Redfield evolution from aggregate, non-secular")
#    rho_t2.plot(coherences=False)    

    rho_i1.save_data("/Users/tomas/rho_i1.dat")
    rho_t1.save_data("/Users/tomas/rho_t1.dat")

    rho = ReducedDensityMatrixEvolution(time)
    rho.load_data("/Users/tomas/rho_t1.dat")


    rho.plot(coherences=False, axis=[0,3000,0.0,1.0])
  
agg.save("aggregate.hdf5")

fid = open("aggregate.pkl","wb")
pickle.dump(agg, fid)
fid.close()

