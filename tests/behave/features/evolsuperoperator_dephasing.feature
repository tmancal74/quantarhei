#
#
#  Test of Lorentzian and Gaussian dephasing in density matrix evolution
#
#

Feature: EvolutionSuperOperator can be calculated with pure dephasing
    EvolutionSuperOperator accepts additional argument PDeph which allows
    one to describe pure dephasing in the eigen state basis of the system
    
    Scenario Outline: EvolutionSuperOperator can be calculated with Lorenzian pure dephasing only
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
        

#
#
#   Gaussian dephasing
#
#         
    
    Scenario Outline: EvolutionSuperOperator can be calculated with Gaussian pure dephasing only
        Given I have PureDephasing object D with dephasing constant <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator using only PureDephasing D with <time_step> and <N_dense>
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I multiply each coherence element by corresponding Gaussian decay with dephasing time <dtime> to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  50.0   |    10.0     |    10     |    50.0  |
        |  50.0   |    10.0     |    10     |   100.0  |
        |  50.0   |    10.0     |    10     |   200.0  |
        |  10.0   |    10.0     |    10     |    50.0  |
        |  10.0   |    10.0     |    10     |   100.0  |
        |  10.0   |    10.0     |    10     |   200.0  |         
        |  10.0   |    10.0     |    10     |   500.0  |
        |  10.0   |    10.0     |    10     |   800.0  |
        |  10.0   |    10.0     |    10     |   990.0  |         
        
        
    Scenario Outline: EvolutionSuperOperator can be calculated with Gaussian pure dephasing step by step
        Given I have PureDephasing object D with dephasing constant <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator step by step using only PureDephasing D with <time_step> and <N_dense> 
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I multiply each coherence element by corresponding Gaussian decay with dephasing time <dtime> to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  50.0   |    10.0     |    10     |    50.0  |
        |  50.0   |    10.0     |    10     |   100.0  |
        |  50.0   |    10.0     |    10     |   200.0  |
        
        
    Scenario Outline: EvolutionSuperOperator can be calculated with Lindbladform and Gaussian pure dephasing step by step
        Given I have a Hamiltonian H, Lidblad form L, PureDephasing object D with dephasing constant <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator step by step using PureDephasing D and Lindblad form L with <time_step> and <N_dense> 
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I propagate with Lindblad form L and PureDephasing D to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |     0.0  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  50.0   |    10.0     |    10     |    50.0  |
        |  50.0   |    10.0     |    10     |   100.0  |
        |  50.0   |    10.0     |    10     |   200.0  |


    Scenario Outline: EvolutionSuperOperator can be calculated with Gaussian pure dephasing in one shot
        Given I have PureDephasing object D with dephasing constant <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator in one shot using only PureDephasing D with <time_step> and <N_dense> 
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I multiply each coherence element by corresponding Gaussian decay with dephasing time <dtime> to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  50.0   |    10.0     |    10     |    50.0  |
        |  50.0   |    10.0     |    10     |   100.0  |
        |  50.0   |    10.0     |    10     |   200.0  |


    Scenario Outline: EvolutionSuperOperator can be calculated with Lindbladform and Gaussian pure dephasing in one shot
        Given I have a Hamiltonian H, Lidblad form L, PureDephasing object D with dephasing constant <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator in one shot using PureDephasing D and Lindblad form L with <time_step> and <N_dense> 
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I propagate with Lindblad form L and PureDephasing D to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |     0.0  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  50.0   |    10.0     |    10     |    50.0  |
        |  50.0   |    10.0     |    10     |   100.0  |
        |  50.0   |    10.0     |    10     |   200.0  |


#
#
# Lorentzian pure dephasing
#
#


    Scenario Outline: EvolutionSuperOperator can be calculated with Lorentzian pure dephasing step by step
        Given I have PureDephasing object D with dephasing time <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator step by step using only PureDephasing D with <time_step> and <N_dense> 
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I multiply each coherence element by corresponding exponential decay with dephasing time <dtime> to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  50.0   |    10.0     |    10     |    50.0  |
        |  50.0   |    10.0     |    10     |   100.0  |
        |  50.0   |    10.0     |    10     |   200.0  |
        
        
    Scenario Outline: EvolutionSuperOperator can be calculated with Lindbladform and Lorentzian pure dephasing step by step
        Given I have a Hamiltonian H, Lidblad form L, PureDephasing object D with dephasing time <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator step by step using PureDephasing D and Lindblad form L with <time_step> and <N_dense> 
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I propagate with Lindblad form L and PureDephasing D to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |     0.0  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  50.0   |    10.0     |    10     |    50.0  |
        |  50.0   |    10.0     |    10     |   100.0  |
        |  50.0   |    10.0     |    10     |   200.0  |


    Scenario Outline: EvolutionSuperOperator can be calculated with Lorentzian pure dephasing in one shot
        Given I have PureDephasing object D with dephasing time <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator in one shot using only PureDephasing D with <time_step> and <N_dense> 
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I multiply each coherence element by corresponding exponential decay with dephasing time <dtime> to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  50.0   |    10.0     |    10     |    50.0  |
        |  50.0   |    10.0     |    10     |   100.0  |
        |  50.0   |    10.0     |    10     |   200.0  |


    Scenario Outline: EvolutionSuperOperator can be calculated with Lindbladform and Lorentzian pure dephasing in one shot
        Given I have a Hamiltonian H, Lidblad form L, PureDephasing object D with dephasing time <dtime> and initial density matrix R
         When I calculate EvolutionSuperOperator in one shot using PureDephasing D and Lindblad form L with <time_step> and <N_dense> 
          And I apply the EvolutionSuperOperator to R to get RD at time <t_prop>
          And I propagate with Lindblad form L and PureDephasing D to get RE at time <t_prop>
         Then RD equals RE at times <t_prop>
         
    Examples: Parameters 
        |  dtime  |  time_step  |  N_dense  |  t_prop  |
        |  100.0  |    10.0     |    10     |     0.0  |
        |  100.0  |    10.0     |    10     |    50.0  |
        |  100.0  |    10.0     |    10     |   100.0  |
        |  100.0  |    10.0     |    10     |   200.0  |
        |  50.0   |    10.0     |    10     |    50.0  |
        |  50.0   |    10.0     |    10     |   100.0  |
        |  50.0   |    10.0     |    10     |   200.0  |
