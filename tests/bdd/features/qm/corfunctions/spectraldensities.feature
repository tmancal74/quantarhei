Feature: Spectral density and correlation function relations

    Spectral density and correlation function are related by Fourier transform,
    and the objects representing them should be created from one another

@spectraldensity
Scenario Outline: Spectral density creation from correlation function
    Given correlation function parameters:
        | cf_type      | temp     | T_units | reorg      | e_units | ctime      | t_units | mats |
        | <cf_type>    | <temp>   | K       | <reorg>    | 1/cm    | <ctime>    | fs      |   20 |
    And upper-half TimeAxis with parameters:
        | start | number_of_steps | step  | units |
        | 0.0   |    1000         |  1.0  | fs    |
    When I calculate the <cf_type> correlation function
    And spectral density is created from correlation function
    Then spectral density corresponds to file <file> in internal units

    Examples:
    | cf_type            | temp    | reorg | ctime | file                      |
    | OverdampedBrownian |  300.   |  30.  | 200.  | ob_sd_30cm_200fs_300K.dat |
#    | Lorenzian          |  300.   |  30.  | 300.  | lr_sd_30cm_300fs_300K.dat |

       
@spectraldensity
Scenario Outline: Spectral density from correlation function compared to analytical
    Given correlation function parameters:
        | cf_type      | temp     | T_units | reorg      | e_units | ctime      | t_units | mats |
        | <cf_type>    | <temp>   | K       | <reorg>    | 1/cm    | <ctime>    | fs      |   20 |
    And upper-half TimeAxis with parameters:
        | start | number_of_steps | step  | units |
        | 0.0   |    1000         |  1.0  | fs    |
    When I calculate the <cf_type> correlation function
    And spectral density is created from correlation function
    Then spectral density corresponds to analytical result for <cf_type> in internal units

    Examples:
    | cf_type            | temp    | reorg | ctime | 
    | OverdampedBrownian |  300.   |  30.  | 200.  | 
    | OverdampedBrownian |  100.   |  100. | 100.  |

@in_development
Scenario Outline: Correlation function creation from spectral density
    Given spectral density parameters:
        | cf_type   | temp   | T_units | reorg    | e_units | ctime   | t_units | mats |
        | <cf_type> | <temp> | K       | <reorg>  | 1/cm    | <ctime> | fs      |   20 |
    And TimeAxis:
        | units | start | number_of_steps | step |
        | fs    | 0.0   | 1000            | 1.0  | 
    When I calculate the <cf_type> spectral density
    And correlation function is created from spectral density
    Then correlation function corresponds to file <file>

    Examples:
    | cf_type            | temp   | reorg | ctime | file                       |
    | OverdampedBrownian |  300   |  20   | 100.  | ob_20cm_100fs_300K_m20.dat |
#    | Lorenzian          |  300   |  30   | 300.  | lr_300K_30cm_300fs.dat |       

@in_development
Scenario Outline: Odd FT of correlation function equals spectral density
    Given correlation function parameters:
        | cf_type   | temp   | T_units | reorg    | e_units | ctime   | t_units | mats |
        | <cf_type> | <temp> | K       | <reorg>  | 1/cm    | <ctime> | fs      |   20 |
    And TimeAxis:
        | units | start | number_of_steps | step |
        | fs    | 0.0   | 1000            | 1.0  | 
    When I calculate the <cf_type> correlation function
    And I calculate the <cf_type> spectral density
    And calculate odd FT of the correlation function
    Then odd FT correlation function corresponds to spectral density

    Examples:
    | cf_type            | temp   | reorg | ctime | file                       |
    | OverdampedBrownian |  300   |  20   | 100.  | ob_20cm_100fs_300K_m20.dat |
#    | Lorenzian          |  300   |  30   | 300.  | lr_300K_30cm_300fs.dat |