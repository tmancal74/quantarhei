# -*- coding: utf-8 -*-

from aloe import step
from aloe import world

from ..stepslib import upper_half_time_axis
from ..stepslib import complete_time_axis

import numpy

@step('FrequencyAxis is created from TimeAxis')
def fa_from_ta(self):
    
    ta = world.ta
    world.wa = ta.get_FrequencyAxis()
    

@step('FrequencyAxis has correct properties')
def fa_correct(self):
    wa = world.wa
        
    if wa.atype == "complete":
        # dw = 2*pi/(n*dt)
        numpy.testing.assert_allclose(wa.data[1]-wa.data[0],
                                  2.0*numpy.pi/(world.ta.length*world.ta.dt))

        # minimum value
        # -2*pi*n/(2*n*dt) if n is even
        if world.ta.length%2 == 0:
            numpy.testing.assert_allclose(wa.omin,
            -2.0*numpy.pi/(world.ta.dt*2))
        # when odd 2*pi*(n-1)/(2*n*dt)
        else:
            numpy.testing.assert_allclose(wa.omin,
            -2.0*numpy.pi*(world.ta.length-1)/(world.ta.dt*world.ta.length*2))
            
        # array length is the same as in time axis
        numpy.testing.assert_equal(wa.length,world.ta.length)        
        
    if wa.atype == "upper-half":
        
        # dw = 2*pi/(2*n*dt)
        numpy.testing.assert_allclose(wa.data[1]-wa.data[0],
                                  2.0*numpy.pi/(2*world.ta.length*world.ta.dt)) 
                                  
 
        # minimum value
        numpy.testing.assert_allclose(wa.omin,-2.0*numpy.pi/(2.0*world.ta.dt))
        # maximum value
        numpy.testing.assert_allclose(wa.omax,
            2.0*numpy.pi*(world.ta.length-1)/(2.0*world.ta.dt*world.ta.length))
       
            
        # array length is twice the one of the time axis
        numpy.testing.assert_equal(wa.length,2*world.ta.length)                                 
        
    numpy.testing.assert_equal(wa.omin,wa.data[0])
   


    
@step('TimeAxis can be recreated from FrequencyAxis')
def ta_correct_thought_fa(self):
        
    ta = world.wa.get_TimeAxis()
    
    # dw = 2*pi/(n*dt)
    numpy.testing.assert_allclose(world.ta.dt,
                                        ta.dt)

    # minimum value
    numpy.testing.assert_allclose(world.ta.tmin,ta.tmin,atol=1.0e-7)
    
    #maximum
    numpy.testing.assert_allclose(world.ta.tmax,ta.tmax)
   

    # array length is the same as in time axis
    numpy.testing.assert_allclose(world.ta.length,ta.length)  