Feature: Molecule creation

    As a user I want to create molecules using different units of energy
    and asign transition dipole moments to them
    and asign position to them

Scenario Outline: A user creates a Molecule using different units 
    Given ground state energy is <ground_state_energy> and excited state energy is <excited_state_energy> "<units>"
    When Molecule is created
    Then it has a correct Hamiltonian with values <gse_internal> and <ese_internal> in internal units

    Examples:
        |units | ground_state_energy | excited_state_energy | gse_internal | ese_internal |
        | 1/cm |                 0.0 |              12000.0 |          0.0 |      12000.0 |
        | Thz  |                 0.0 |                300.0 |          0.0 |        300.0 |


