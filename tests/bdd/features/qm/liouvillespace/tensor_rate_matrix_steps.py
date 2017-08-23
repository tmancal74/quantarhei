# -*- coding: utf-8 -*-
from aloe import step
from aloe import world

import numpy


@step(r'Redfield relaxation tensor has the same rates as Redfield rate matrix with rtol (\d+(?:\.\d+)?)')
def compare_rates_and_tensor(self, rtols):
    
    rtol = float(rtols)

    rates_T = world.rates_T
    rates_M = world.rates_M
    
    print(rates_T)
    print(rates_M)
            
    numpy.testing.assert_allclose(rates_M, rates_T, rtol=rtol)    