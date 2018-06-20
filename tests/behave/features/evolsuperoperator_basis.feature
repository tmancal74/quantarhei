#
#  EvolutionSuperOperator testing 
#

Feature: Evolution super operator basis transformations
    Evolution super operators need to transform correctly under basis
    transformation.

    Scenario Outline: Application of an evolution superoperator in different basis leads to the same result
        Given I have a general evolution superoperator U, initial density matrix R and a Hamiltonian H
         When I apply U to R in in time <t_prop> to get density matrix Rt
          And I transform U and R to the eigenbasis of H to get U_ and R_, respectively
          And I apply U_ to R_ to get operator Rt_trans
          And I transform Rt into eigenbasis of H to get Rt_
         Then Rt_ equals Rt_trans

       Examples: Propagation times
           | t_prop  |
           | 0.0     |
           | 50.0    |
           | 100.0   |
           | 200.0   |
           | 300.0   |
           
    Scenario Outline: Evolution superoperator is calculated correctly in site basis
        Given I have a Hamiltonian H, Lidblad form L and initial density matrix R
        When I calculate evolution superoperator using H and L in sitebasis
        And I calculate dynamics of R using H and L to get R1
        And I apply the evolution superoperator to R to get R2 at times <t_prop>
        Then R1 equals R2 at times <t_prop>
        
       Examples: Propagation times
           | t_prop  |
           | 0.0     |
           | 10.0    |
           | 20.0    |
           | 50.0    |
           | 30.0    |
                   

    Scenario Outline: Evolution superoperator is calculated correctly in exciton basis
        Given I have a Hamiltonian H, Lidblad form L and initial density matrix R
        When I calculate evolution superoperator using H and L in exciton basis
        And I calculate dynamics of R using H and L to get R1
        And I apply the evolution superoperator to R to get R2 at times <t_prop>
        Then R1 equals R2 at times <t_prop>
        
       Examples: Propagation times
           | t_prop  |
           | 0.0     |
           | 10.0    |
           | 20.0    |
           | 50.0    |
           | 30.0    |

           
    Scenario Outline: Evolution superoperator is calculated correctly in combination of bases
        Given I have a Hamiltonian H, Lidblad form L and initial density matrix R
        When I calculate evolution superoperator in site basis using H and L in exciton basis
        And I calculate dynamics of R using H and L to get R1
        And I apply the evolution superoperator to R to get R2 at times <t_prop>
        Then R1 equals R2 at times <t_prop>
        
       Examples: Propagation times
           | t_prop  |
           | 0.0     |
           | 10.0    |
           | 20.0    |
           | 50.0    |
           | 30.0    | 