# -*- coding: utf-8 -*-

import unittest
import numpy

import quantarhei as qr

"""
*******************************************************************************


    Tests of the quantarhei.qm.propagators.poppropagator module


*******************************************************************************
"""

from quantarhei import ReducedDensityMatrixPropagator
from quantarhei import TimeAxis
from quantarhei.qm import Liouvillian

class TestRDMPropagator(unittest.TestCase):
    """Tests reduced density matrix propagator module
    
    
    """
    
    def test_rdm_evolution_time_independent_tensor(self):
        """Testing evolution of reduced density matrix with time-independent relaxation tensor
        
        """
        HH = qr.Hamiltonian(data=[[0.0, 1.0],[1.0, 0.2]])
        P01 = qr.qm.ProjectionOperator(0, 1, dim=2)
        P10 = qr.qm.ProjectionOperator(1, 0, dim=2)
        rates = [1.0/100.0, 1.0/600.0]
        
        sbi = qr.qm.SystemBathInteraction(sys_operators=[P01,P10],
                                          rates=rates)
        LL = qr.qm.LindbladForm(HH, sbi)
        LL.convert_2_tensor()
        Nt = 500
        dt = 0.01
        time = qr.TimeAxis(0.0, Nt, dt)
        
        #
        # Propagation by propagator
        #
        prop = qr.ReducedDensityMatrixPropagator(time, Ham=HH, RTensor=LL)
        rho_ini = qr.ReducedDensityMatrix(data=[[0.0, 0.0],[0.0, 1.0]])
        
        
        #
        # WE TEST THIS
        #
        rhot_1 = prop.propagate(rho_ini)
        
        
        #
        # Propagation by matrix exponentiation
        #
        Li = Liouvillian(HH)
        # exponentiate matrix with dt
        Leff = (-1j*Li.data + LL.data)*dt
        
        Leff.shape = [4,4]
        
        Ld, SS = numpy.linalg.eig(Leff)
        Ld1 = numpy.diag(Ld)
        
        Lexpd = numpy.diag(numpy.exp(Ld))
        
        S1 = numpy.linalg.inv(SS)
        
        # test
        Ld2 = numpy.dot(S1, numpy.dot(Leff, SS))
        
        numpy.testing.assert_allclose(Ld1, Ld2, rtol=1.0e-7, atol=1.0e-7)
        
        Ut = numpy.dot(SS, numpy.dot(Lexpd, S1))
        Ut.shape=[2,2,2,2]
        
        rhot_2 = qr.qm.ReducedDensityMatrixEvolution(time, rhoi=rho_ini)
        
        for i_t in range(1, time.length):
            rhot_2.data[i_t,:,:] = numpy.tensordot(Ut, rhot_2.data[i_t-1,:,:])
        
        
        #
        # Compare the two
        #
        numpy.testing.assert_allclose(rhot_1.data, rhot_2.data, rtol=1.0e-6,
                                      atol=1.0e-9)
        
        
        


    def test_rdm_evolution_Saveable(self):
        pass
