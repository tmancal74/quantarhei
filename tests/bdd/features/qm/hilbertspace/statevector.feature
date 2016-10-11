Feature: Transformation of a state vector

    A user creates a state vector in one basis
    and wants to transform it into a different one.

@statevector
Scenario Outline: State vector is transformed 
    Given a homodimer with site transition energies <tr_en> "<e_units>"
    And resonance coupling <coupl_en> "<e_units>"
    And a state vector with parameters:
        |  nocomp  | oneat   |
        | <nocomp> | <oneat> |
    When homodimer hamiltonian is created
    And state vector is transformed in eigenbasis of the hamiltonian
    Then I get correctly transformed homodimer state:
        |  nocomp  | oneat   |
        | <nocomp> | <oneat> |

    Examples:
        | tr_en | coupl_en | e_units | nocomp | oneat  |  
        | 12000 |  300     | 1/cm    | 3      |  2     |   
        | 12000 |  30      | 1/cm    | 3      |  1     |  
