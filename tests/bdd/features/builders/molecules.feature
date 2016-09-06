Feature: Molecule creation

    As a user I want to create molecules using different units of energy
    and asign transition dipole moments to them
    and asign position to them

Scenario Outline: A user creates a Molecule using different units 
    Given I define a two-level molecule:
        |  units  |  ground_state_energy  |  excited_state_energy  |
        | <units> | <ground_state_energy> | <excited_state_energy> |
    When molecule is created
    Then molecule has a correct Hamiltonian with values <gse_internal> and <ese_internal> in internal units

    Examples:
        |units | ground_state_energy | excited_state_energy | gse_internal | ese_internal |
        | 1/cm |                 0.0 |              12000.0 |          0.0 |      12000.0 |
        | Thz  |                 0.0 |                300.0 |          0.0 |        300.0 |

@in_development
Scenario Outline: A user adds a transition dipole moment to a molecule
    Given I define two-level molecules:
        |ground_state_energy | excited_state_energy | units |
        |                0.0 |              12000.0 |  1/cm |
        |                0.0 |                300.0 |  THz  |
    And I create transition dipole moment vectors:
        | dx  |  dy |  dz | units |
        | 0.0 | 0.0 | 1.0 | Debye |
        | 1.0 | 0.0 | 0.0 | a.u.  |
    When molecules are created
    And transition dipole moments are set to molecules
    Then molecules return transition dipole moment vectors:
        | dx  |  dy |  dz | units |
        | 0.0 | 0.0 | 1.0 | Debye |
        | 1.0 | 0.0 | 0.0 | a.u.  |



