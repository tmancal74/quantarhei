Feature: Correlation function definition

Scenario Outline: A user creates correlation function with various parameters
    Given reorganization energy <reorg> "<e_units>"
    And correlation time <ctime> "<t_units>" 
    And temperature <temp> "<T_units>"
    And time interval:
        | start | number_of_steps | step | units |
        | 0.0   |    100          |  10  | fs    |
    When I calculate the <ctype> correlation function
    Then I get data from the file <file>

    Examples:
        | ctype              | reorg | e_units | ctime | t_units | temp   | T_units | file                   |
        | OverdampedBrownian | 20.0  | 1/cm    | 100   |   fs    | 300    | K       | ob_20cm_100fs_300K.dat |
        | OverdampedBrownian | 20.0  | 1/cm    | 100   |   fs    | 100    | K       | ob_20cm_100fs_100K.dat |

