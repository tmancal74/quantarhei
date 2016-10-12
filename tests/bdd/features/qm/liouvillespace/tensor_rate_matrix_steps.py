# -*- coding: utf-8 -*-
from aloe import step
from aloe import world

import numpy

@step(r'Redfield relaxation tensor has the same rates as Redfield rate matrix with rtol (\d+(?:\.\d+)?)')
def compare_rates_and_tensor(self, rtols):
    
    rtol = float(rtols)
    
    RRT = world.RRT
    RRM = world.RRM
    
    dim = world.RRM.data.shape[0]
    rates_T = numpy.zeros(dim*dim)
    rates_M = numpy.zeros(dim*dim)
    
    k = 0
    for i in range(dim):
        for j in range(dim):
            rates_T[k] = numpy.real(RRT.data[i,i,j,j])
            rates_M[k] = RRM.data[i,j]
            k += 1
            
    print(rates_T)
    print(rates_M)
            
    numpy.testing.assert_allclose(rates_M, rates_T, rtol=rtol)    