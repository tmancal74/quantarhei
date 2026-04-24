Feature: Bose-Einstein distribution

@in_development
Scenario Outline: Creation of Bose-Einstein distribution
    Given temperature <temp> "<T_units>"
    And FrequencyAxis:
       | start   | number_of_steps | step   | units     |
       | <start> | <nrst>          | <step> | <f_units> |
    When I create Bose-Einstein distribution
    Then it satisfies 1 + n(om) = - n(-om)

    Examples:
       | temp | T_units | f_units | start  | nrst | step  |
       | 300  | K       | 2pi/fs  | -2.0   | 4000 | 0.001 |
       | 100  | K       | 1/cm    | -20000 | 1000 | 40    |

