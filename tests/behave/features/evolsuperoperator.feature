#
#  EvolutionSuperOperator testing 
#

Feature: Evolution super operator basis transformations
    Evolution super operators need to transform correctly under basis
    transformation.

    Scenario Outline: Application of an evolution superoperator different basis leads to the same result
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