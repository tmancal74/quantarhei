Feature: Time axis can be converted into frequency axis

@time
Scenario Outline: Frequency axis creation from a half time axis
    Given upper-half TimeAxis with parameters:
        | start     | number_of_steps   | step     | units   | 
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When FrequencyAxis is created from TimeAxis
    Then FrequencyAxis has correct properties

Examples:
        | start | number_of_steps | step | units |
        | 0.0   |    100          |  10  | fs    |
        | 0.0   |    101          |  1.0 | fs    |


@time
Scenario Outline: Half time axis self-creation through frequency axis
    Given upper-half TimeAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When FrequencyAxis is created from TimeAxis
    Then TimeAxis can be recreated from FrequencyAxis

Examples:
        | start | number_of_steps | step | units |
        | 0.0   |    100          |  10  | fs    |
        | 0.0   |    101          |  1.0 | fs    |

@time
Scenario Outline: Frequency axis creation from a complete time axis
    Given complete TimeAxis with parameters:
        | start     | number_of_steps   | step     | units   | 
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When FrequencyAxis is created from TimeAxis
    Then FrequencyAxis has correct properties

Examples:
        | start | number_of_steps | step | units |
        | 0.0   |    100          |  10  | fs    |
        | 0.0   |    101          |  1.0 | fs    |

@time
Scenario Outline: Comlete time axis self-creation through frequency axis
    Given complete TimeAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When FrequencyAxis is created from TimeAxis
    Then TimeAxis can be recreated from FrequencyAxis

Examples:
        | start | number_of_steps | step | units |
        | 0.0   |    100          |  10  | fs    |
        | 0.0   |    101          |  1.0 | fs    |