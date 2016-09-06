Feature: Molecule creation

Scenario: A user creates a Molecule
    Given ground state energy is 0.0 and excited state energy is 12000 1/cm
    When Molecule is created
    Then it has a correct Hamiltonian

