# -*- coding: utf-8 -*-

import aloe

@aloe.step(r"When Molecule is created")
def when_Molecule_is_created(self):
    pass

@aloe.step(r"Given ground state energy is (\d+(?:\.\d+)? and excited state energy is (\d+(?:\.\d+)? 1/cm")
def trans_energies_for_a_molecule(self,en0,en1):
    print(en0, en1)

@aloe.step(r"Then it has a correct Hamiltonian")
def molecule_has_correct_Hamiltonian(self):
    pass

