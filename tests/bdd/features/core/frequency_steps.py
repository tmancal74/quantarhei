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
                    (2.0*numpy.pi)/(world.wa.length*world.wa.domega))

    numpy.testing.assert_allclose(ta.tmax,ta.data[ta.length-1])

#    # minimum value
#    # 2*pi*n/2 if n is even
#    if world.ta.length%2 == 0:
#        numpy.testing.assert_allclose(wa.omin,
#            -2.0*numpy.pi/(world.ta.dt*2))
#    # when odd
#    else:
#        numpy.testing.assert_allclose(wa.omin,
#            -2.0*numpy.pi*(world.ta.length-1)/(world.ta.dt*world.ta.length*2))
#    numpy.testing.assert_equal(wa.omin,wa.data[0])
#   
#
#    # array length is the same as in time axis
#    numpy.testing.assert_equal(wa.length,world.ta.length)
    
@step('FrequencyAxis can be recreated from TimeAxis')
def ta_correct_through_fa(self):
        
    wa = world.ta.get_FrequencyAxis()
    #print(wa.omin,world.wa.omin)    
    
    # dw = 2*pi/(n*dt)
    numpy.testing.assert_allclose(world.wa.domega,
                                        wa.domega)

    # minimum value
    numpy.testing.assert_allclose(world.wa.omin,wa.omin,atol=1.0e-7)
    
    #maximum
    numpy.testing.assert_allclose(world.wa.omax,wa.omax)
   

    # array length is the same as in time axis
    numpy.testing.assert_allclose(world.wa.length,wa.length)  