
#
#  Combined Redfield and Foerster relaxation theories 
#
#
#

Feature: Combined Redfield-Foerster relaxation theory 

    Combination of Redfield and Foerster relaxation tensors 

#@redfield-foerster
@in_development
Scenario Outline: Redfield-Foerster rates for a small chain of molecules


    Given time interval from zero with time step <time_step> fs and <nsteps> steps
    And a small homochain of N molecules with parameters
       | tr_energy | neighbor_coupling | corfce | reorg | ctime | matsubara | temperature | e_units |
       | 12000     |  100              | OB     | 40    | 100   |  <matsu>  |  <temp>     | 1/cm    |
       | 12000     |  80               | OB     | 40    | 100   |  <matsu>  |  <temp>     | 1/cm    | 
       | 12200     |  30               | OB     | 40    | 80    |  <matsu>  |  <temp>     | 1/cm    |
       | 12100     |  200              | OB     | 50    | 100   |  <matsu>  |  <temp>     | 1/cm    |
       | 12100     |  100              | OB     | 50    | 100   |  <matsu>  |  <temp>     | 1/cm    |
    When I calculate aggregate Redfield-Foerster relaxation tensor with cutoff coupling of <cutoff_coup> 1/cm
    Then long time Redfield-Foerster populations will correspond to canonical equilibrium with temperature <temp> K with atol <atol>
    Examples:
       | temp | time_step | nsteps | matsu  | cutoff_coup |atol   |
       | 300  | 0.1       | 10000  | 100    |   50        | 0.01  |
       | 200  | 0.1       | 10000  | 100    |   50        | 0.02  |
       | 100  | 0.1       | 10000  | 100    |   50        | 0.02  |
       | 50   | 0.1       | 10000  | 100    |   50        | 0.02  |

