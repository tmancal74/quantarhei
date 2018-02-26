# -*- coding: utf-8 -*-
import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.EvolutionSuperOperator class


*******************************************************************************
"""


class TestEvSupOp(unittest.TestCase):
    """Tests for the EvolutionSuperOperator class
    
    
    """
    
    def setUp(self,verbose=False):
        """Initializes the calculation
        
        """
        
        self.verbose = verbose
                  
        
        
    def test_redfield_dynamics_comp(self):
        """Compares Redfield dynamics calculated from propagator and superoperator
        
        
        
        """
        # Aggregate
        import quantarhei as qr
        import quantarhei.models.modelgenerator as mgen
        
        time = qr.TimeAxis(0.0, 1000, 1.0)
        
        # create a model
        mg = mgen.ModelGenerator()
        
        agg = mg.get_Aggregate_with_environment(name="trimer-1_env", 
                                                timeaxis=time )
        
        agg.build()
        
        sbi = agg.get_SystemBathInteraction()
        ham = agg.get_Hamiltonian()

        # calculate relaxation tensor
        ham.protect_basis()
        with qr.eigenbasis_of(ham):
    
            RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi)
            RRT.secularize()
            
        ham.unprotect_basis()
        with qr.eigenbasis_of(ham):
 
            #
            # Evolution of reduced density matrix
            #
            prop = qr.ReducedDensityMatrixPropagator(time, ham, RRT)        
            
            #
            # Evolution by superoperator
            #
            
            eSO = qr.qm.EvolutionSuperOperator(time, ham, RRT)
            eSO.set_dense_dt(5)
            eSO.calculate()
            
            # compare the two propagations
            pairs = [(1,3), (3,2), (3,3)]
            for p in pairs:
                
            
                rho_i1 = qr.ReducedDensityMatrix(dim=ham.dim)
                rho_i1.data[p[0],p[1]] = 1.0
            
                rho_t1 = prop.propagate(rho_i1)

            
                exp_rho_t2 = eSO.data[:,:,:,p[0],p[1]]


                #import matplotlib.pyplot as plt
        
                #plt.plot(rho_t1.TimeAxis.data, numpy.real(rho_t1.data[:,p[0],p[1]]))
                #plt.plot(rho_t1.TimeAxis.data, numpy.real(exp_rho_t2[:,p[0],p[1]]))
                #plt.show()

                #numpy.testing.assert_allclose(RRT.data, rtd)
                numpy.testing.assert_allclose(numpy.real(rho_t1.data[:,:,:]),
                                          numpy.real(exp_rho_t2[:,:,:]), 
                                          rtol=1.0e-7,
                                          atol=1.0e-6)
                
                
