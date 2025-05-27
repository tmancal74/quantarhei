# -*- coding: utf-8 -*-

import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.propagators.poppropagator module


*******************************************************************************
"""


from quantarhei import PopulationPropagator
from quantarhei import TimeAxis

class TestPopulationPropagator(unittest.TestCase):
    """Tests population propagator module
    
    
    """
    
    def test_of_population_evolution_1(self):
        """Testing population evolution matrix 2x2 starting from t = 0"""
        

        KK = numpy.array([[-1.0/100.0,  1.0/100.0],
                          [ 1.0/100.0, -1.0/100.0]])
        
        U0 = numpy.eye(2)
        Ntd = 10
        
        t = TimeAxis(0.0, 1000, 1.0)
        prop = PopulationPropagator(t, rate_matrix=KK)
        
        
        td = TimeAxis(0.0, Ntd, 10.0)
        U = prop.get_PropagationMatrix(td)
        
        # analytical result
        
        Ucheck = numpy.zeros((2,2,Ntd))
        Ucheck[0,0,:] = 0.5*(1.0+numpy.exp(2.0*KK[0,0]*td.data))
        Ucheck[1,1,:] = Ucheck[0,0,:]
        Ucheck[1,0,:] = 0.5*(1.0-numpy.exp(2.0*KK[0,0]*td.data))
        Ucheck[0,1,:] = Ucheck[1,0,:]
        
        for n in range(Ntd):
            numpy.testing.assert_allclose(U[:,:,n],Ucheck[:,:,n])


    def test_of_population_evolution_2(self):
        """Testing population evolution matrix 2x2 starting from t > 0"""
        

        KK = numpy.array([[-1.0/100.0,  1.0/100.0],
                          [ 1.0/100.0, -1.0/100.0]])
        
        U0 = numpy.eye(2)
        Ntd = 10
        
        t = TimeAxis(0.0, 1000, 1.0)
        prop = PopulationPropagator(t, rate_matrix=KK)


        td = TimeAxis(2.0, Ntd, 10.0)
        U = prop.get_PropagationMatrix(td)
        
        # analytical result
        
        Ucheck = numpy.zeros((2,2,Ntd))
        Ucheck[0,0,:] = 0.5*(1.0+numpy.exp(2.0*KK[0,0]*td.data))
        Ucheck[1,1,:] = Ucheck[0,0,:]
        Ucheck[1,0,:] = 0.5*(1.0-numpy.exp(2.0*KK[0,0]*td.data))
        Ucheck[0,1,:] = Ucheck[1,0,:]
        
        for n in range(Ntd):
            numpy.testing.assert_allclose(U[:,:,n],Ucheck[:,:,n])
        
