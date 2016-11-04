
#
#  Secular Redfield relaxation theory results in canonical equilibrum
#
#
#

Feature: Secular Redfield thgeory results in cannonical equilibrium 

    Theories such as Redfield relaxation tensor etc. should result in
    cannonical equilibrium

#@redfield
@in_development
Scenario Outline: Redfield rates for a small chain of molecules


    Given time interval from zero with time step <time_step> fs and <nsteps> steps
    And a small homochain of N molecules with parameters
       | tr_energy | neighbor_coupling | corfce | reorg | ctime | matsubara | temperature | e_units |
       | 12000     |  100              | OB     | 20    | 100   |  <matsu>  |  <temp>     | 1/cm    | 
       | 12200     |  50               | OB     | 20    | 80    |  <matsu>  |  <temp>     | 1/cm    |
       | 12100     |  200              | OB     | 20    | 100   |  <matsu>  |  <temp>     | 1/cm    |
    When I calculate aggregate Redfield relaxation tensor
    Then long time Redfield populations will correspond to canonical equilibrium with temperature <temp> K with atol <atol>
    Examples:
       | temp | time_step | nsteps | matsu  | atol   |
       | 300  | 1.0       | 10000  | 100    |  0.01  |
       | 200  | 1.0       | 10000  | 100    |  0.02  |
       | 100  | 1.0       | 10000  | 100    |  0.02  |
       | 50   | 1.0       | 10000  | 100    |  0.02  |


