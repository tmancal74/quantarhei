#
#
#  Test of Lorentzian and Gaussian dephasing in density matrix evolution
#
#

Feature: EvolutionSuperOperator can be calculated with pure dephasing
    EvolutionSuperOperator accepts additional argument PDeph which allows
    one to describe pure dephasing in the eigen state basis of the system
    
    Scenario Outline: EvolutionSuperOperator can be calculated with pure dephasing only
        Given I have PureDephasing object D with dephasing time <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator using only PureDephasing D with <time_step> and <N_dense>
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I multiply each coherence element by corresponding exponential decay with dephasing time <dtime> to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  10.0   |    10.0     |    10     |    50.0  |
        |  10.0   |    10.0     |    10     |   100.0  |
        |  10.0   |    10.0     |    10     |   200.0  |
        
         
    
         
        