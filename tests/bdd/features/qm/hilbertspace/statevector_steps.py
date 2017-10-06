from aloe import step
from aloe import world

import numpy

from quantarhei import StateVector, Molecule, Aggregate
from quantarhei import energy_units, eigenbasis_of

@step(r'state vector with parameters')
def state_vector(self):
    
    for row in self.hashes:
        ncomp = int(row['nocomp'])
        oneat = int(row['oneat'])
        
    sva = numpy.zeros(ncomp, dtype=numpy.complex128)
    sva[oneat] = 1.0
    sv = StateVector(data=sva)
    world.sv = sv

@step(r'homodimer hamiltonian is created')
def hamiltonian_is_created(self):
    
    #
    # Homodimer
    #
    with energy_units(world.h_units):
        en = world.senergy
        m1 = Molecule([0.0, en], "mol1")
        m2 = Molecule([0.0, en], "mol2")
        
    agg = Aggregate(name="Homodimer")
    
    agg.add_Molecule(m1)
    agg.add_Molecule(m2)
        
    with energy_units(world.r_units):
        agg.set_resonance_coupling(0,1,world.r_coupl)
        
    agg.build()
    
    world.HH = agg.get_Hamiltonian()
    
    
@step(r'state vector is transformed in eigenbasis of the hamiltonian')
def state_vector_transformed_in_eig_of_hamiltonian(self):
    
    with eigenbasis_of(world.HH):
        vec = world.sv.data
    world.sv_transformed = vec
    

@step(r'I get correctly transformed homodimer state')
def assert_correct_homodimer_state(self):
    
    
    for row in self.hashes:
        ncomp = int(row['nocomp'])
        oneat = int(row['oneat'])
        
    S = numpy.zeros((3,3), dtype=numpy.complex128)
    S[0,0] = 1.0
    S[1,1] = -1.0
    S[2,2] = 1.0
    S[1,2] = 1.0
    S[2,1] = 1.0
    S = S/numpy.sqrt(2.0)
    
    v = numpy.zeros(ncomp, dtype=numpy.complex128)
    v[oneat] = 1.0
    
    vtr = numpy.dot(S,v)
    
    numpy.testing.assert_allclose(world.sv_transformed, vtr)
    
    