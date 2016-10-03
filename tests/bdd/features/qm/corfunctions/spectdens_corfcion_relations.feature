#
# Tests of the mutual relation of the correlation functions and spectral densities
#
# Correlation functions and spectral densities can be specified by a corresponding set of parameters,
# and, usually, only one of the two can be given by a simple formula. The other of the two
# has to be calculated numerically. 
#
# In these tests we test the mutual consistency of the classes representing correlation functions
# and spectral density. 
#
# The following two relations have to be tested:
#
# 1) When spectral density is specified for positive frequencies (the function is odd), one can obtain
# a Fourier transform of the correlation function by multiplying the spectral density by a 
# factor containing temperature. By inverse Fourier transform we can obtain correlation function
#
# 2) When correlation function is speciefied, the Fourier transform of its imagonary part gives
# the spectral density. 
#
# The two relations can be achieved analytically, when formulae are known, or numerically in the case
# that they are not known. Consistency of the numerics can and should be tested on the analytical
# results


Feature: Spectral density and correlation function relations

    Spectral density and correlation function are related by Fourier transform,
    and it should be possible to create the objects representing them from one another

#
# Creation of spectral density using correlation function and test against the 
# expected result specified in a file
#
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

#
# Creation of a spectral density from correlation function and a check
# against the expected analytical result.
#       
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

#
# Correlation function and spectral density created from the same parameters
# and the spectral density obtained by FFT is compared to the spectral density 
#
@spectraldensity
Scenario Outline: Odd FT of correlation function equals spectral density
    Given correlation function parameters:
        | cf_type   | temp   | T_units | reorg    | e_units | ctime   | t_units | mats |
        | <cf_type> | <temp> | K       | <reorg>  | 1/cm    | <ctime> | fs      |   20 |
    And upper-half TimeAxis with parameters:
        | units | start | number_of_steps | step |
        | fs    | 0.0   | 1000            | 1.0  | 
    When I calculate the <cf_type> correlation function
    And I calculate the <cf_type> spectral density
    And I calculate odd FT of the correlation function
    Then odd FT correlation function corresponds to spectral density

    Examples:
    | cf_type            | temp   | reorg | ctime |
    | OverdampedBrownian |  300   |  20   | 100.  |
#    | Lorenzian          |  300   |  30   | 300.  |


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
Scenario Outline: Correlation function creation from spectral density and back
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
Scenario Outline: Spectral density creation from correlation function and back
    Given correlation function parameters:
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
