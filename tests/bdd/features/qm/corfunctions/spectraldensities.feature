Feature: Spectral density definition

@spectraldensity
Scenario Outline: A user creates spectral density with various parameters
    Given reorganization energy <reorg> "<e_units>"
    And correlation time <ctime> "<t_units>" 
    And temperature <temp> "<T_units>"
    And number of Matsubara frequencies <mats>
    And upper-half TimeAxis with parameters:
        | start | number_of_steps | step  | units |
        | 0.0   |    1000         |  1.0  | fs    |
    When I calculate the <ctype> spectral density 
    Then spectral density corresponds to file <file> in internal units

    Examples:
        | ctype              | reorg | e_units | ctime | t_units | temp   | T_units | mats | file                          |
        | OverdampedBrownian | 20.0  | 1/cm    | 100   |   fs    | 300    | K       | 20   | ob_sd_20cm_100fs_300K_m20.dat |
        | OverdampedBrownian | 30.0  | 1/cm    | 200   |   fs    | 300    | K       | 20   | ob_sd_30cm_200fs_300K_m20.dat |


