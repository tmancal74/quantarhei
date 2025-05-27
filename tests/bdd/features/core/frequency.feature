Feature: Frequency axis can be converted into time axis

@frequency
Scenario Outline: Complete time axis creation from a frequency axis
    Given complete FrequencyAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When TimeAxis is created from FrequencyAxis
    Then TimeAxis has correct properties

Examples:
        | start | number_of_steps | step | units  |
        | 0.0   |    100          |  10  | 2pi/fs |
        | 0.0   |    101          |  1.0 | 2pi/fs |
        | 10.0  |    200          |  0.1 | 2pi/fs |

@frequency
Scenario Outline: Comlete frequency axis recreation through time axis
    Given complete FrequencyAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When TimeAxis is created from FrequencyAxis
    Then FrequencyAxis can be recreated from TimeAxis

Examples:
        | start | number_of_steps | step | units  |
        | 0.0   |    100          |  10  | 2pi/fs |
        | 0.0   |    101          |  1.0 | 2pi/fs |
        | 10.0  |    201          |  0.1 | 2pi/fs |


@frequency
Scenario Outline: Upper-half time axis creation from a frequency axis
    Given upper-half FrequencyAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When TimeAxis is created from FrequencyAxis
    Then TimeAxis has correct properties

Examples:
        | start | number_of_steps | step | units  |
        | 0.0   |    100          |  10  | 2pi/fs |
        | 0.0   |    100          |  1.0 | 2pi/fs |
        | 20.0  |    290          | 0.1  | 2pi/fs |

@frequency
Scenario Outline: Upper-half frequency axis recreation through time axis
    Given upper-half FrequencyAxis with parameters:
        | start     | number_of_steps   | step     | units   |
        | <start>   | <number_of_steps> |  <step>  | <units> |
    When TimeAxis is created from FrequencyAxis
    Then FrequencyAxis can be recreated from TimeAxis

Examples:
        | start | number_of_steps | step | units  |
        | 0.0   |    100          |  10  | 2pi/fs |
        | 0.0   |    100          |  1.0 | 2pi/fs |
        | 10.0  |    200          | 0.1  | 2pi/fs |