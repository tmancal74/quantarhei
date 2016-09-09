Feature: Calculation of Redfield relaxation rates in a homodimer

    As a user I want to calculate relaxation rates in a homodimer by Redfield theory

Scenario Outline: Refield rates for a homodimer 
    Given a homodimer with site transition energies <tr_en> "<e_units>"
    And resonance coupling <coupl_en> "<e_units>"
    And correlation function with parameters:
        | cf_type            | reorg_en | e_units| cor_time | t_units| temp   | T_units | matsubara |
        | OverdampedBrownian | <reorg>  | 1/cm   | <ctime>  |  fs    | <temp> | K       | 20        |
    And TimeAxis 
        | start | number_of_steps | step | units |
        | 0.0   |   1000          | 1.0  | 1/cm  |
    When I calculate Redfield relaxation rates
    Then I get rates from file <file>

    Examples:
        | tr_en | coupl_en | e_units | reorg | ctime | temp |
        | 12000 |  300     | 1/cm    | 30    | 200   | 300  |
 
