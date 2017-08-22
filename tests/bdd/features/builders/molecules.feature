Feature: Molecule creation

Scenario: A user creates a Molecule
          Given transition energy and transition dipole moment
          When Molecule is created
          Then it has a correct Hamiltonian and Dipole Moment Operator

