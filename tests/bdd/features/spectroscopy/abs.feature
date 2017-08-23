Feature: Molecular absorption spectrum

@absorption
Scenario Outline: A user calculates absorption spectrum of a two-level molecule
    Given reorganization energy <reorg> "<e_units>"
    And correlation time <ctime> "<t_units>" 
    And temperature <temp> "<T_units>"
    And number of Matsubara frequencies <mats>
    And upper-half TimeAxis with parameters:
        | start | number_of_steps | step | units |
        | 0.0   |    1000         |  1.0 | fs    |
    When I calculate the <ctype> correlation function
    And I calculate absorption spectrum of a molecule
    Then I get absorption spectrum from the file <file> in 1/cm

    Examples:
        | ctype              | reorg | e_units | ctime | t_units | temp   | T_units | mats | file                             |
        | OverdampedBrownian | 20.0  | 1/cm    | 100   |   fs    | 300    | K       | 20   | abs_1mol_20cm_100fs_300K_m20.dat |
        | OverdampedBrownian | 20.0  | 1/cm    | 100   |   fs    | 100    | K       | 20   | abs_1mol_20cm_100fs_100K_m20.dat |

@absorption
Scenario Outline: A user calculates absorption spectrum of an aggregate of two two-level molecule
    Given reorganization energy <reorg> "<e_units>"
    And correlation time <ctime> "<t_units>" 
    And temperature <temp> "<T_units>"
    And number of Matsubara frequencies <mats>
    And upper-half TimeAxis with parameters:
        | start | number_of_steps | step | units |
        | 0.0   |    2000         |  5.0 | fs    |
    When I calculate the <ctype> correlation function
    And I calculate absorption spectrum of a dimer aggregate
    Then I get absorption spectrum from the file <file> in 1/cm

    Examples:
        | ctype              | reorg | e_units | ctime | t_units | temp   | T_units | mats | file                            |
        | OverdampedBrownian | 30.0  | 1/cm    | 60   |   fs     | 300    | K       | 10   | abs_2mol_30cm_60fs_300K_m10.dat |
        | OverdampedBrownian | 10.0  | 1/cm    | 60   |   fs     | 100    | K       | 20   | abs_2mol_10cm_60fs_100K_m20.dat |

@absorption
Scenario Outline: A user calculates absorption spectrum of an aggregate of three two-level molecule
    Given reorganization energy <reorg> "<e_units>"
    And correlation time <ctime> "<t_units>" 
    And temperature <temp> "<T_units>"
    And number of Matsubara frequencies <mats>
    And upper-half TimeAxis with parameters:
        | start | number_of_steps | step | units |
        | 0.0   |    2000         |  5.0 | fs    |
    When I calculate the <ctype> correlation function
    And I calculate absorption spectrum of a trimer aggregate
    Then I get absorption spectrum from the file <file> in 1/cm

    Examples:
        | ctype              | reorg | e_units | ctime | t_units | temp   | T_units | mats | file                            |
        | OverdampedBrownian | 20.0  | 1/cm    | 60   |   fs     | 300    | K       | 20   | abs_3mol_20cm_60fs_300K_m20.dat |
        | OverdampedBrownian | 20.0  | 1/cm    | 60   |   fs     | 100    | K       | 20   | abs_3mol_20cm_60fs_100K_m20.dat |
