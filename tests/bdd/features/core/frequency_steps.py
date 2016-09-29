# -*- coding: utf-8 -*-

from aloe import step
from aloe import world

from ..stepslib import upper_half_frequency_axis
from ..stepslib import complete_frequency_axis

import numpy

@step('TimeAxis is created from FrequencyAxis')
def ta_from_fa(self):
    
    wa = world.wa
    world.ta = wa.get_TimeAxis()
    

@step('TimeAxis has correct properties')
def ta_correct(self):
    ta = world.ta
        
    # dt = 2*pi/(n*domega)
    numpy.testing.assert_allclose(ta.data[1]-ta.data[0],
                    (2.0*numpy.pi)/(world.wa.length*world.wa.step))

    numpy.testing.assert_allclose(ta.max,ta.data[ta.length-1])

    
@step('FrequencyAxis can be recreated from TimeAxis')
def ta_correct_through_fa(self):
        
    wa = world.ta.get_FrequencyAxis()
    #print(wa.omin,world.wa.omin)    
    
    # dw = 2*pi/(n*dt)
    numpy.testing.assert_allclose(world.wa.step,
                                        wa.step)

    # minimum value
    numpy.testing.assert_allclose(world.wa.min,wa.min,atol=1.0e-7)
    
    #maximum
    numpy.testing.assert_allclose(world.wa.max,wa.max)
   

    # array length is the same as in time axis
    numpy.testing.assert_allclose(world.wa.length,wa.length)  