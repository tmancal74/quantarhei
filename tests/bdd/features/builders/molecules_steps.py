# -*- coding: utf-8 -*-

from aloe import world
from aloe import step
from nose.tools import assert_equal 

import numpy

from quantarhei import Molecule
from quantarhei import energy_units

""" ********************************************************************* """

""" Feature: Molecule creation """


"""Given"""
@step(r'I define a two-level molecule')
def trans_energies_for_a_molecule(self): #,en0,en1,units):
    # Initialize from the sentence
    # Regular expressions give strings!
    for row in self.hashes: 
        e0 = float(row['ground_state_energy'])
        e1 = float(row['excited_state_energy'])
        units = row['units']

    # save it for other steps
    world.e0 = e0 
    world.e1 = e1
    world.units = units

    print("\nunits ",units)
    
"""When"""
@step(r'molecule is created')
def when_Molecule_is_created(self):
    with energy_units(world.units):
        world.m1 = Molecule([world.e0,world.e1], "Mol1") 

"""Then"""
@step(r"molecule has a correct Hamiltonian with values (\d+(?:\.\d+)?) and (\d+(?:\.\d+)?) in internal units")
def molecule_has_correct_Hamiltonian(self,val0,val1):
    e0 = float(val0)
    e1 = float(val1)
    # get the Hamiltonian from the Molecule
    h = world.m1.get_Hamiltonian()

    # This is how the Hamiltonian data shoud look like
    data = numpy.zeros((2,2),dtype=numpy.float32)
    data[0,0] = e0
    data[1,1] = e1 

    numpy.testing.assert_allclose(data,h.data)


"""Given"""
@step(r'I define two-level molecules')
def i_define_a_two_level_molecule(self):
    world.e0 = []
    world.e1 = []
    world.units = []
    for row in self.hashes:
        world.e0.append(float(row['ground_state_energy']))
        world.e1.append(float(row['excited_state_energy']))
        world.units.append(row['units'])

"""And"""
@step(r'I create transition dipole moment vectors')
def i_create_tdm_vectors(self):
    world.dvec = []
    for row in self.hashes:
        dx = float(row['dx'])
        dy = float(row['dy'])
        dz = float(row['dz'])
        vec = [dx, dy, dz]
        world.dvec.append(vec)

"""When"""
@step(r'molecules are created')
def when_molecules_are_created(self):
    N = len(world.e0)
    molecules = []
    for i in range(N):
        e0 = world.e0[i]
        e1 = world.e1[i]
        units = world.units[i]
        m1 = Molecule("Mol",[e0,e1])
        molecules.append(m1)

    world.molecules = molecules

"""And"""
@step(r'transition dipole moments are set to molecules')
def tdm_are_set_to_molecules(self):
    N = len(world.molecules)
    for i in range(N):
        m = world.molecules[i]
        vec = world.dvec[i]
        m.set_dipole(0,1,vec)

"""Then"""
@step(r'molecules return transition dipole moment vectors')
def molecules_return_tdm_vectors(self):
    pass


 
""" ********************************************************************** """



