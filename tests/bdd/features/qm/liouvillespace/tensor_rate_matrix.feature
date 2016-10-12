Feature: Redfield tensor and Redfield rate matrix contain the same rates

    The two implementations of the Redfield theory have to give exactly the
    same results

@samerates
Scenario Outline: Comparing relaxation tensor with the rate matrix
    Given a homodimer with site transition energies <tr_en> "<e_units>"
    And resonance coupling <coupl_en> "<e_units>"
    And correlation function parameters:
        | cf_type            | reorg   | e_units| ctime    | t_units| temp   | T_units | mats |
        | OverdampedBrownian | <reorg> | 1/cm   | <ctime>  |  fs    | <temp> | K       | 20   |
    And upper-half TimeAxis with parameters: 
        | start | number_of_steps | step | units |
        | 0.0   |   1000          | 1.0  | 1/cm  |
    When I calculate Redfield relaxation tensor
    And I calculate Redfield relaxation rates
    Then Redfield relaxation tensor has the same rates as Redfield rate matrix with rtol 0.005

    Examples:
        | tr_en | coupl_en | e_units | reorg | ctime  | temp | 
        | 12000 |  100     | 1/cm    | 300   | 100    | 300  |  
        | 12000 |  30      | 1/cm    | 20    | 100    | 300  |  
        | 12000 |  100     | 1/cm    | 20    | 50     | 300  |
        | 12000 |  30      | 1/cm    | 100   | 50     | 300  |