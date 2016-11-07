
#
#  Some types of relaxation theories result in canonical equilibrum
#
#
#

Feature: Foerster relaxation theory result in cannonical equilibrium 

    Theories such as Foerster relaxation tensor etc. should result in
    cannonical equilibrium

@foerster
#@in_development
Scenario Outline: Foerster rates for a small chain of molecules


    Given time interval from zero with time step <time_step> fs and <nsteps> steps
    And a small homochain of N molecules with parameters
       | tr_energy | neighbor_coupling | corfce | reorg | ctime | matsubara | temperature | e_units |
       | 12000     |  100              | OB     | 40    | 100   |  <matsu>  |  <temp>     | 1/cm    | 
       | 12200     |  50               | OB     | 40    | 80    |  <matsu>  |  <temp>     | 1/cm    |
       | 12100     |  200              | OB     | 50    | 100   |  <matsu>  |  <temp>     | 1/cm    |
    When I calculate aggregate Foerster relaxation tensor
    Then long time Foerster populations will correspond to canonical equilibrium with temperature <temp> K with atol <atol>
    Examples:
       | temp | time_step | nsteps | matsu  | atol   |
       | 300  | 0.1       | 10000  | 100    |  0.01  |
       | 200  | 0.1       | 10000  | 100    |  0.02  |
       | 100  | 0.1       | 10000  | 100    |  0.02  |
       | 50   | 0.1       | 10000  | 100    |  0.02  |


