#
#  Tests of the standard and time-dependent rate theories
#
#
#
#
#


Feature: Calculation of Redfield relaxation rates in a homodimer

    As a user I want to calculate relaxation rates in a homodimer by Redfield theory

#
#  Comparison against results saved in a file
#
@in_development
Scenario Outline: Refield rates for a homodimer against saved results
    Given a homodimer with site transition energies <tr_en> "<e_units>"
    And resonance coupling <coupl_en> "<e_units>"
    And correlation function parameters:
        | cf_type            | reorg   | e_units| ctime    | t_units| temp   | T_units | mats |
        | OverdampedBrownian | <reorg> | 1/cm   | <ctime>  |  fs    | <temp> | K       | 20   |
    And upper-half TimeAxis with parameters: 
        | start | number_of_steps | step | units |
        | 0.0   |   1000          | 1.0  | 1/cm  |
    When I calculate Redfield relaxation rates
    Then I get Redfield relaxation rates from file <file>

    Examples:
        | tr_en | coupl_en | e_units | reorg | ctime | temp | file                                          |
        | 12000 |  300     | 1/cm    | 30    | 200   | 300  | rdf_homo2_e12000cm_J300cm_30cm_200fs_300K.dat |
 
#
#  Comparison against analytical results
#
@in_development
Scenario Outline: Refield rates for a homodimer against analytical result
    Given a homodimer with site transition energies <tr_en> "<e_units>"
    And resonance coupling <coupl_en> "<e_units>"
    And correlation function parameters:
        | cf_type            | reorg   | e_units| ctime    | t_units| temp   | T_units | mats |
        | OverdampedBrownian | <reorg> | 1/cm   | <ctime>  |  fs    | <temp> | K       | 20   |
    And upper-half TimeAxis with parameters: 
        | start | number_of_steps | step | units |
        | 0.0   |   1000          | 1.0  | 1/cm  |
    When I calculate Redfield relaxation rates
    Then I get Redfield relaxation rates corresponding to analytical results for a homodimer

    Examples:
        | tr_en | coupl_en | e_units | reorg | ctime  | temp | 
        | 12000 |  300     | 1/cm    | 300    | 100   | 300  |    