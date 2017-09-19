
#
#  Population progator `propagate` method and the propagation matrix 
#  from `get_PropagationMatrix` method must aggree
#
#
#

Feature: Result of propagate method and the propagation matrix agree 

    We check that the propagate method gives the same populations as the
    propagation matrix

@poppropagator
Scenario Outline: Propagator for a simple relaxation matrix


    Given a propagation time interval from zero with time step <time_step> fs and <nsteps> steps
    And a subset time interval from zero with time step <time_step_2> fs and <nsteps_2> steps 
    And a relaxation matrix with uphill time <t_up> fs and downhill time <t_down> fs
    When I calculate density matrix propagation and propagation matrix
    Then density matrix propagation and propagation matrix aggree
    Examples:
       | time_step | nsteps  | time_step_2 | nsteps_2 | t_up  | t_down |
       | 0.1       | 10000   | 10.0        |   50     | 100.0 | 30.0   |
       | 0.1       | 10000   | 10.0        |   50     | 200.0 | 100.0  |
       | 0.1       | 10000   | 10.0        |   50     | 100.0 | 100.0  |
       | 0.1       | 10000   | 10.0        |   50     | 60.0  | 30.0   |

