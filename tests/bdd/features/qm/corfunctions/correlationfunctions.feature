Feature: Correlation function definition

Scenario Outline: A user creates OverdampedBrownian correlation function
    Given reorganization energy <reorg> and correlation time <ctime>
    When I calculate the <ctype> correlation function
    Then I get data from the file <file>

    Examples:
        | ctype             | reorg | ctime | file          |
        | OverdampedBrownian| 20.0  | 100   | ob_20_100.dat |
