Feature: Frequency axis can be converted into time axis

@frequency
Scenario Outline: Time axis creation from a frequency axis
    Given complete FrequencyAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When TimeAxis is created from FrequencyAxis
    Then TimeAxis has correct properties

Examples:
        | start | number_of_steps | step | units  |
        | 0.0   |    100          |  10  | 2pi/fs |
        | 0.0   |    101          |  1.0 | 2pi/fs |

@frequency
Scenario Outline: Frequency axis recreation through time axis
    Given complete FrequencyAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When TimeAxis is created from FrequencyAxis
    Then FrequencyAxis can be recreated from TimeAxis

Examples:
        | start | number_of_steps | step | units  |
        | 0.0   |    100          |  10  | 2pi/fs |
        | 0.0   |    101          |  1.0 | 2pi/fs |


@frequency
Scenario Outline: Time axis creation from a frequency axis
    Given upper-half FrequencyAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When TimeAxis is created from FrequencyAxis
    Then TimeAxis has correct properties

Examples:
        | start | number_of_steps | step | units  |
        | 0.0   |    100          |  10  | 2pi/fs |
        | 0.0   |    100          |  1.0 | 2pi/fs |

@frequency
Scenario Outline: Frequency axis recreation through time axis
    Given upper-half FrequencyAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When TimeAxis is created from FrequencyAxis
    Then FrequencyAxis can be recreated from TimeAxis

Examples:
        | start | number_of_steps | step | units  |
        | 0.0   |    100          |  10  | 2pi/fs |
        | 0.0   |    100          |  1.0 | 2pi/fs |