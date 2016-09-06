# -*- coding: utf-8 -*-

from aloe import world
from aloe import step
from nose.tools import assert_equal 

import numpy

from quantarhei import Molecule


""" ********************************************************************* """

@step(r'When Molecule is created')
def when_Molecule_is_created(self):
    world.m1 = Molecule("Mol1",[world.e0,world.e1]) 

@step(r'Given ground state energy is (\d+(?:\.\d+)?) and excited state energy is (\d+(?:\.\d+)?) "([^"]*)"')
def trans_energies_for_a_molecule(self,en0,en1,units):
    # Initialize from the sentence
    # Regular expressions give strings!
    e0 = float(en0)
    e1 = float(en1)

    # save it for other steps
    world.e0 = e0 
    world.e1 = e1
    world.units = units


@step(r"Then it has a correct Hamiltonian with values (\d+(?:\.\d+)?) and (\d+(?:\.\d+)?) in internal units")
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
 
""" ********************************************************************** """



