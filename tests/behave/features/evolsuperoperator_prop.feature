#
#  EvolutionSuperOperator testing 
#

Feature: Evolution super operator for large step propagation
    Evolution super operators enables to propagate evolution of density
    matrix with a large time step. 

    Scenario Outline: Evolution superoperator can be calculated with a time step which is integer times larger than required for numerics
        Given I have a Hamiltonian H, Lidblad form L and initial density matrix R 
         When I calculate evolution superoperator U with time step <N_dense> times large than needed for numerics
          And I calculate dynamics of R using H and L to get R1
          And I apply the evolution superoperator to R to get R2 at times <t_prop>
         Then R1 equals R2 at times <t_prop>         

       Examples: Propagation times
           | N_dense | t_prop  |
           | 10      | 0.0     |
           | 10      | 50.0    |
           | 20      | 100.0   |
           | 20      | 200.0   |
           | 15      | 300.0   |


    Scenario Outline: Evolution superoperator returns a superoperator whose application is equivalent to propagation by evolution super operator
        Given I have a Hamiltonian H, Lidblad form L and initial density matrix R 
         When I calculate evolution superoperator U with time step <N_dense> times large than needed for numerics
          And I apply the evolution superoperator to R to get R2 at times <t_prop>
          And I get the evolution super operator S at a single time <t_prop>
          And I apply the evolution superoperator S to R to get R3 at times <t_prop>
         Then R2 equals R3 at times <t_prop>         

       Examples: Propagation times
           | N_dense | t_prop  |
           | 10      | 0.0     |
           | 10      | 50.0    |
           | 20      | 100.0   |
           | 20      | 200.0   |
           | 15      | 300.0   |

