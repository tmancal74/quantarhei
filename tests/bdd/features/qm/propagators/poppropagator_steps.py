# -*- coding: utf-8 -*-

"""

    Test of the population propagator
    


"""
from aloe import step
from aloe import world

import numpy

from quantarhei.testing.feature import FeatureFileGenerator
from quantarhei.testing.feature import match_number 

from quantarhei import TimeAxis
from quantarhei import PopulationPropagator


feature = """
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

"""

example_1 = """
    Examples:
       | time_step | nsteps  | time_step_2 | nsteps_2 | t_up  | t_down |
       | 0.1       | 10000   | 10.0        |   50     | 100.0 | 30.0   |
       | 0.1       | 10000   | 10.0        |   50     | 200.0 | 100.0  |
       | 0.1       | 10000   | 10.0        |   50     | 100.0 | 100.0  |
       | 0.1       | 10000   | 10.0        |   50     | 60.0  | 30.0   |

"""


@step(r'a propagation time interval from zero with time step '+match_number+r' fs and '
      +match_number+r' steps')
def time_interval_prop(self, time_step, nsteps):
    """Matching the following
    
    #begin_feature    
    Given a propagation time interval from zero with time step <time_step> fs and <nsteps> steps
    #end_feature
    
    """
    world.time = TimeAxis(0.0, int(nsteps), float(time_step))
    print("Setting time")

@step(r'a subset time with time step '+match_number+r' fs and '
      +match_number+r' steps')
def time_interval_sub(self, time_step, nsteps):
    """Matching the following
    
    #begin_feature    
    And a subset time with time step <time_step_2> fs and <nsteps_2> steps
    #end_feature
    
    """
    world.subtime = TimeAxis(0.0, int(nsteps), float(time_step))
    print("Setting subtime")


@step(r'a relaxation matrix with uphill time '+match_number+r' fs and downhill time '
      +match_number+r' fs')
def relaxation_matrix(self, uphill, downhill):
    """Matching the following
    
    #begin_feature    
    And a relaxation matrix with uphill time <t_up> fs and downhill time <t_down> fs
    #end_feature
    
    """
    world.KK = numpy.zeros((2,2), dtype=numpy.float64)
    Kup = 1.0/float(uphill)
    world.KK[0,0] = -Kup
    world.KK[1,0] = Kup
    Kdn = 1.0/float(downhill)
    world.KK[0,1] = Kdn
    world.KK[1,1] = -Kdn
    
    
@step(r' I calculate density matrix propagation and propagation matrix')
def calculation_of_propagation(self):
    """Matching the following
    
    #begin_feature    
    When I calculate density matrix propagation and propagation matrix
    #end_feature
    
    """ 
    
    prop = PopulationPropagator(world.time, rate_matrix=world.KK)
    
    pop_ini = numpy.array([1.0, 0.0])
    
    pop_t = prop.propagate(pop_ini)
    
    sta = world.subtime
    
    U = prop.get_PropagationMatrix(sta)
    
    pop_sub = numpy.zeros((2,sta.length))
    
    for i in range(sta.length):
        pop_sub[:,i] = numpy.dot(U[:,:,i],pop_ini)  
        
    world.pop_t = pop_t
    world.pop_sub = pop_sub
    
@step(r' density matrix propagation and propagation matrix aggree')
def test_of_agreement(self):
    """Matching the following
    
    #begin_feature    
    Then density matrix propagation and propagation matrix aggree
    #end_feature
    
    """ 
            
    pop_t = world.pop_t 
    pop_sub = world.pop_sub   
    
    
    Ns = world.subtime.length
    
    dt = world.time.step
    ds = world.subtime.step
    N = round(ds/dt)
    for i in range(Ns):
        numpy.testing.assert_allclose(pop_sub[:,i],pop_t[i*N,:])
        

if __name__ == '__main__':
    """
    
    Feature file: poppropagator_1.feature
    
    """
    gen = FeatureFileGenerator(feature+example_1)
    
    gen.add_Given(time_interval_prop)
    gen.add_Given(time_interval_sub)
    gen.add_Given(relaxation_matrix)
    gen.add_When(calculation_of_propagation)
    gen.add_Then(test_of_agreement)
    
    gen.generate_feature_file("poppropagator_1.feature")
