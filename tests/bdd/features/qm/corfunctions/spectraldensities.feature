Feature: Spectral density and correlation function relations

    Spectral density and correlation function are related by Fourier transform,
    and the objects representing them should be created from one another

@in_development
Scenario Outline: Spectral density creation from correlation function
    Given correlation function parameters:
        | cf_type      | temp     | T_units | reorg      | e_units | ctime      | t_units | mats |
        | <cf_type> | <temp> | K          | <reorg>  | 1/cm    | <ctime> | fs         |     20 |
    And TimeAxis:
        | units | start | number_of_steps | step |
        | fs      | 0.0    |   1000                   | 1.0    | 
    When correlation function is created
    And spectral density is created from correlation function
    Then spectral density corresponds to file <file>

    Examples:
    | cf_type                           | temp   | reorg | ctime | file                                               |
    | OverdampedBrownian |  300.   |  30.    | 200.   | ob_sd_300K_30cm_200fs.dat |
    | Lorenzian                       |  300.   |.  30.   | 300.  |  lr_sd_300K_30cm_300fs.dat   |         

@in_development
Scenario Outline: Correlation function creation from spectral density
    Given spectral density parameters:
        | cf_type      | temp     | T_units | reorg      | e_units | ctime      | t_units | mats |
        | <cf_type> | <temp> | K          | <reorg>  | 1/cm    | <ctime> | fs         |     20 |
    And TimeAxis:
        | units | start | number_of_steps | step |
        | fs      | 0.0    | 1000                     | 1.0    | 
    When spectral density is created
    And correlation function is created from spectral density
    Then correlation function corresponds to file <file>

    Examples:
    | cf_type                           | temp   | reorg | ctime | file                                         |
    | OverdampedBrownian |  300.   |  30.    | 200.   | ob_300K_30cm_200fs.dat |
    | Lorenzian                       |  300.   |.  30.   | 300.  |  lr_300K_30cm_300fs.dat   |       

