# -*- coding: utf-8 -*-
import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.EvolutionSuperOperator class


*******************************************************************************
"""
from quantarhei import TimeAxis
from quantarhei import Molecule, Mode
from quantarhei import CorrelationFunction
from quantarhei import energy_units, eigenbasis_of
from quantarhei import EvolutionSuperOperator
from quantarhei import ReducedDensityMatrix

from quantarhei import COMPLEX


def _test_as_op(as_op, obj):
    """Test calculation with teme independent Redfield tensor
    
    """
    
    time = obj.time
    mol1 = obj.mol1
   
    #######################################################################
    #
    #  Direct calculation - no RWA
    #
    #######################################################################
    
    # Relaxation tensor and Hamiltonian to calculate the evolution
    (RTi, HH) = mol1.get_RelaxationTensor(time, relaxation_theory="stR",
                                   time_dependent=False, 
                                   as_operators=as_op)

    RTi_noRWA = RTi.deepcopy()

    # a very simple initial condition
    rhoi_t = ReducedDensityMatrix(dim=HH.dim)
    rhoi_a = ReducedDensityMatrix(dim=HH.dim)
    with eigenbasis_of(HH):
        rhoi_t.data[1,2] = 1.0
        rhoi_a.data[0,2] = 1.0
    
    # store initial condition independently
    rhoi_c = rhoi_t.deepcopy()
    
    # want to make sure that Molecule returns the right Hamiltonian
    H1 = mol1.get_Hamiltonian()
    obj.assertIs(H1, HH)
    

    #
    # Control propagation of the above initial condition
    #
    prop_i = mol1.get_ReducedDensityMatrixPropagator(time,
                                    relaxation_theory="stR",
                                    time_dependent=False, 
                                    as_operators=as_op)
    
    rhot_i = prop_i.propagate(rhoi_t)
    rhot_a = prop_i.propagate(rhoi_a)

    #
    # calculation of evolution superoperator
    #
    time_eop = TimeAxis(0.0, 100, 10.0)
    UUi = EvolutionSuperOperator(time=time_eop, ham=HH, relt=RTi)
    UUi.set_dense_dt(10)

    UUi.calculate()     

    #
    #
    #  Verification that the results are the same
    #
    #
    
    # verify that the initial condition was not changed
    numpy.testing.assert_allclose(rhoi_t.data, rhoi_c.data)
    numpy.testing.assert_allclose(rhot_i.data[0,:,:], rhoi_t.data[:,:])
    with eigenbasis_of(HH):
        numpy.testing.assert_allclose(rhot_i.data[0,:,:], 
                                      rhoi_t.data[:,:])
        
    # verifying some time points
    attime = UUi.apply(0.0, rhoi_t)
    attime_a = UUi.apply(0.0, rhoi_a)
    numpy.testing.assert_allclose(rhot_i.data[0,:,:], attime.data)
    numpy.testing.assert_allclose(rhot_a.data[0,:,:], attime_a.data)
    attime = UUi.apply(10.0, rhoi_t)
    attime_a = UUi.apply(10.0, rhoi_a)
    numpy.testing.assert_allclose(rhot_i.data[10,:,:], attime.data)
    numpy.testing.assert_allclose(rhot_a.data[10,:,:], attime_a.data)
    attime = UUi.apply(100.0, rhoi_t)
    attime_a = UUi.apply(100.0, rhoi_a)
    numpy.testing.assert_allclose(rhot_i.data[100,:,:], attime.data)
    numpy.testing.assert_allclose(rhot_a.data[100,:,:], attime_a.data)
    attime = UUi.apply(900.0, rhoi_t)
    attime_a = UUi.apply(900.0, rhoi_a)
    numpy.testing.assert_allclose(rhot_i.data[900,:,:], attime.data)
    numpy.testing.assert_allclose(rhot_a.data[900,:,:], attime_a.data)        
    
    _plot_ = False
    n1 = 1
    n2 = 2
    if _plot_:
        import matplotlib.pyplot as plt
        
        with eigenbasis_of(HH):

            #plt.plot(time_eop.data, numpy.real(dat_t[:, n1, n2, n1, n2]))
            plt.plot(time_eop.data, 
                     numpy.real(UUi.data[:, n1, n2, n1, n2]),"-")
            plt.plot(time.data, numpy.real(rhot_i.data[:, n1, n2]),"--")
            # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,2,2,2]))
            # plt.plot(time_eop.data, numpy.real(UU.data[:, 3,3,3,3]))
            # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,3,2,3]))
            # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,0,2,0]))
        
        plt.show()


    # saving results for comparison
    saved_data = numpy.zeros((2, rhot_i.data.shape[0]), dtype=COMPLEX)
    with eigenbasis_of(HH):
        saved_data[0,:] = rhot_i.data[:, n1, n2]
        saved_data[1,:] = rhot_a.data[:, 0, n2]

    #######################################################################
    #
    #  RWA calculation 
    #
    #######################################################################
    
    # enforcing new Hamiltonian obejct and setting RWA
    HH = mol1.get_Hamiltonian(recalculate=True)

    HH.set_rwa([0,1])
    
    # Relaxation tensor and Hamiltonian to calculate the evolution
    (RTi, HH) = mol1.get_RelaxationTensor(time, relaxation_theory="stR",
                                   time_dependent=False, 
                                   as_operators=as_op)
 


    #
    # Control propagation of the above initial condition
    #
    prop_i = mol1.get_ReducedDensityMatrixPropagator(time,
                                    relaxation_theory="stR",
                                    time_dependent=False, 
                                    as_operators=as_op)

    
    rhot_i = prop_i.propagate(rhoi_t)
    rhot_a = prop_i.propagate(rhoi_a)
    obj.assertTrue(rhot_i.is_in_rwa)
    obj.assertTrue(rhot_a.is_in_rwa)
    rhot_i.convert_from_RWA(HH)
    rhot_a.convert_from_RWA(HH)

    #
    # calculation of evolution superoperator
    #
    time_eop = TimeAxis(0.0, 100, 10.0)
    UUi = EvolutionSuperOperator(time=time_eop, ham=HH, relt=RTi)
    UUi.set_dense_dt(10)
    obj.assertFalse(UUi.is_in_rwa)

    UUi.calculate()
    obj.assertTrue(UUi.is_in_rwa)
    UUi.convert_from_RWA()
    obj.assertFalse(UUi.is_in_rwa)

    #
    #
    #  Verification that the results are the same
    #
    #

    # make sure that RWA and non-RWA relaxation tensors are the same
    if RTi.as_operators:
        RTi.convert_2_tensor()
    if RTi_noRWA.as_operators:
        RTi_noRWA.convert_2_tensor()
    numpy.testing.assert_allclose(RTi.data, RTi_noRWA.data)

    
    # verify that the initial condition was not changed
    numpy.testing.assert_allclose(rhoi_t.data, rhoi_c.data)
    numpy.testing.assert_allclose(rhot_i.data[0,:,:], rhoi_t.data[:,:])
    with eigenbasis_of(HH):
        numpy.testing.assert_allclose(rhot_i.data[0,:,:], rhoi_t.data[:,:])
        
    #verifying some time points
    attime = UUi.apply(0.0, rhoi_t)
    attime_a = UUi.apply(0.0, rhoi_a)
    numpy.testing.assert_allclose(rhot_i.data[0,:,:], attime.data)
    numpy.testing.assert_allclose(rhot_a.data[0,:,:], attime_a.data)
    attime = UUi.apply(10.0, rhoi_t)
    attime_a = UUi.apply(10.0, rhoi_a)
    numpy.testing.assert_allclose(rhot_i.data[10,:,:], attime.data)
    numpy.testing.assert_allclose(rhot_a.data[10,:,:], attime_a.data)
    attime = UUi.apply(100.0, rhoi_t)
    attime_a = UUi.apply(100.0, rhoi_a)
    numpy.testing.assert_allclose(rhot_i.data[100,:,:], attime.data)
    numpy.testing.assert_allclose(rhot_a.data[100,:,:], attime_a.data)
    attime = UUi.apply(900.0, rhoi_t)
    attime_a = UUi.apply(900.0, rhoi_a)
    numpy.testing.assert_allclose(rhot_i.data[900,:,:], attime.data)
    numpy.testing.assert_allclose(rhot_a.data[900,:,:], attime_a.data)


    #
    #
    # Verification against data without RWA
    #
    #
    with eigenbasis_of(HH):
        numpy.testing.assert_allclose(numpy.real(rhot_i.data[:,n1,n2]), 
                                      numpy.real(saved_data[0,:]),
                                      atol=2.0e-4)
        numpy.testing.assert_allclose(numpy.real(rhot_a.data[:,0,n2]), 
                                      numpy.real(saved_data[1,:]),
                                      atol=2.0e-4)        

    _plot_ = False
    if _plot_:
        import matplotlib.pyplot as plt
        
        with eigenbasis_of(HH):
            n1 = 1
            n2 = 2
            #plt.plot(time_eop.data, numpy.real(dat_t[:, n1, n2, n1, n2]))
            plt.plot(time_eop.data, 
                     numpy.real(UUi.data[:, 0, n2, 0, n2]),"-b")
            plt.plot(time.data, numpy.real(rhot_a.data[:, 0, n2]),"--")
            plt.plot(time.data, numpy.real(saved_data[1, :]),"-k")
            # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,2,2,2]))
            # plt.plot(time_eop.data, numpy.real(UU.data[:, 3,3,3,3]))
            # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,3,2,3]))
            # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,0,2,0]))
        
        plt.show()





class TestEvSupOp(unittest.TestCase):
    """Tests for the EvolutionSuperOperator class
    
    
    """
    
    def setUp(self,verbose=False):
        """Initializes the calculation
        
        """
        
        self.verbose = verbose

        
        with energy_units("1/cm"):
            
            #
            # The molecule has a small energy gap so that we can compare
            # RWA and non-RWA calculations
            #
            
            mol1 = Molecule([0.0, 1200.0])
            mod1 = Mode(300.0)
            mol1.add_Mode(mod1)
            mod1.set_nmax(0, 1)
            mod1.set_nmax(1, 4)
            mod1.set_HR(1,0.1)
            
            mol1.set_dipole((0,1), [1.0, 0.0, 0.0])
            
            time = TimeAxis(0.0, 1000, 1.0)
            params = dict(ftype="OverdampedBrownian",
                          cortime=30.0, reorg=20, T=300,
                          matsubara=30)
            
            cfce = CorrelationFunction(time, params=params)
            
            params = dict(ftype="OverdampedBrownian",
                          cortime=30.0, reorg=50, T=300,
                          matsubara=30)
            
            cfc1 = CorrelationFunction(time, params=params)
            
        mol1.set_mode_environment(0, 0, corfunc=cfce)
        mol1.set_mode_environment(0, 1, corfunc=cfce)
        mol1.set_transition_environment((0,1), cfc1)                  
        
        self.time = time
        self.mol1 = mol1
        
        
        
    def test_Redfield_calculation_time_independent_op(self):
        """(EvolutionSuperOperator) Calculation with time-independent Redfield tensor in operator form
        
        """

        as_op = True
        
        _test_as_op(as_op, self)
        

    def test_Redfield_calculation_time_independent_ten(self):
        """(EvolutionSuperOperator) Calculation with time-independent Redfield tensor in tensor form
        
        """

        as_op = False
        
        _test_as_op(as_op, self)




    def test_Redfield_calculation_time_dep(self):
        """(EvolutionSuperOperator) Test calculation with time-dependent Redfield tensor
        
        """
        #
        # 1) HOW ARE DIFFERENT TIME AXES HANDLED?
        # 2) Time dependent Redfield gives different results with different axes
        # 3) Is RWA handled correctly in EvolutionSuperOperator ? 
        # 4) We have to add is_in_rwa property and convert_to and from RWA
        #

        as_op = False
        tdep = True

        time = self.time
        mol1 = self.mol1
       
        #######################################################################
        #
        #  Direct calculation - no RWA
        #
        #######################################################################
        
        # Relaxation tensor and Hamiltonian to calculate the evolution
        (RTi, HH) = mol1.get_RelaxationTensor(time, relaxation_theory="stR",
                                       time_dependent=tdep, 
                                       as_operators=as_op)

        RTi_noRWA = RTi.deepcopy()

        # a very simple initial condition
        rhoi_t = ReducedDensityMatrix(dim=HH.dim)
        rhoi_a = ReducedDensityMatrix(dim=HH.dim)
        with eigenbasis_of(HH):
            rhoi_t.data[1,2] = 1.0
            rhoi_a.data[0,2] = 1.0
        
        # store initial condition independently
        rhoi_c = rhoi_t.deepcopy()
        
        # want to make sure that Molecule returns the right Hamiltonian
        H1 = mol1.get_Hamiltonian()
        self.assertIs(H1, HH)

        #HH.set_rwa([0,1])

        #
        # Control propagation of the above initial condition
        #
        prop_i = mol1.get_ReducedDensityMatrixPropagator(time,
                                        relaxation_theory="stR",
                                        time_dependent=tdep, 
                                        as_operators=as_op)
        
        rhot_i = prop_i.propagate(rhoi_t)
        rhot_a = prop_i.propagate(rhoi_a)

        #
        # calculation of evolution superoperator
        #
        time_eop = TimeAxis(0.0, 1000, 1.0)
        UUi = EvolutionSuperOperator(time=time_eop, ham=HH, relt=RTi)
        UUi.set_dense_dt(1)

        UUi.calculate()     

        #
        #
        #  Verification that the results are the same
        #
        #
        
        # verify that the initial condition was not changed
        numpy.testing.assert_allclose(rhoi_t.data, rhoi_c.data)
        numpy.testing.assert_allclose(rhot_i.data[0,:,:], rhoi_t.data[:,:])
        with eigenbasis_of(HH):
            numpy.testing.assert_allclose(rhot_i.data[0,:,:], 
                                          rhoi_t.data[:,:])
            
        # verifying some time points
        attime = UUi.apply(0.0, rhoi_t)
        attime_a = UUi.apply(0.0, rhoi_a)
        numpy.testing.assert_allclose(rhot_i.data[0,:,:], attime.data)
        numpy.testing.assert_allclose(rhot_a.data[0,:,:], attime_a.data)
        attime = UUi.apply(10.0, rhoi_t)
        attime_a = UUi.apply(10.0, rhoi_a)
        numpy.testing.assert_allclose(rhot_i.data[10,:,:], attime.data)
        numpy.testing.assert_allclose(rhot_a.data[10,:,:], attime_a.data)
        attime = UUi.apply(100.0, rhoi_t)
        attime_a = UUi.apply(100.0, rhoi_a)
        numpy.testing.assert_allclose(rhot_i.data[100,:,:], attime.data)
        numpy.testing.assert_allclose(rhot_a.data[100,:,:], attime_a.data)
        attime = UUi.apply(900.0, rhoi_t)
        attime_a = UUi.apply(900.0, rhoi_a)
        numpy.testing.assert_allclose(rhot_i.data[900,:,:], attime.data)
        numpy.testing.assert_allclose(rhot_a.data[900,:,:], attime_a.data)        
        
        _plot_ = False
        n1 = 0
        n2 = 2
        if _plot_:
            import matplotlib.pyplot as plt
            
            with eigenbasis_of(HH):

                #plt.plot(time_eop.data, numpy.real(dat_t[:, n1, n2, n1, n2]))
                plt.plot(time_eop.data, 
                         numpy.real(UUi.data[:, n1, n2, n1, n2]),"-")
                plt.plot(time.data, numpy.real(rhot_a.data[:, n1, n2]),"--")
                # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,2,2,2]))
                # plt.plot(time_eop.data, numpy.real(UU.data[:, 3,3,3,3]))
                # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,3,2,3]))
                # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,0,2,0]))
            
            plt.show()


        # saving results for comparison
        saved_data = numpy.zeros((2, rhot_i.data.shape[0]), dtype=COMPLEX)
        with eigenbasis_of(HH):
            saved_data[0,:] = rhot_i.data[:, n1, n2]
            saved_data[1,:] = rhot_a.data[:, 0, n2]

        #######################################################################
        #
        #  RWA calculation 
        #
        #######################################################################
        
        # enforcing new Hamiltonian obejct and setting RWA
        HH = mol1.get_Hamiltonian(recalculate=tdep)

        HH.set_rwa([0,1])
        
        # Relaxation tensor and Hamiltonian to calculate the evolution
        (RTi, HH) = mol1.get_RelaxationTensor(time, relaxation_theory="stR",
                                       time_dependent=tdep, 
                                       as_operators=as_op)
     


        #
        # Control propagation of the above initial condition
        #
        prop_i = mol1.get_ReducedDensityMatrixPropagator(time,
                                        relaxation_theory="stR",
                                        time_dependent=tdep, 
                                        as_operators=as_op)

        
        rhot_i = prop_i.propagate(rhoi_t)
        rhot_a = prop_i.propagate(rhoi_a)
        self.assertTrue(rhot_i.is_in_rwa)
        self.assertTrue(rhot_a.is_in_rwa)
        rhot_i.convert_from_RWA(HH)
        rhot_a.convert_from_RWA(HH)

        #
        # calculation of evolution superoperator
        #
        time_eop = TimeAxis(0.0, 100, 10.0)
        UUi = EvolutionSuperOperator(time=time_eop, ham=HH, relt=RTi)
        UUi.set_dense_dt(10)
        self.assertFalse(UUi.is_in_rwa)

        UUi.calculate()
        self.assertTrue(UUi.is_in_rwa)
        UUi.convert_from_RWA()
        self.assertFalse(UUi.is_in_rwa)

        #
        #
        #  Verification that the results are the same
        #
        #

        # make sure that RWA and non-RWA relaxation tensors are the same
        if RTi.as_operators:
            RTi.convert_2_tensor()
        if RTi_noRWA.as_operators:
            RTi_noRWA.convert_2_tensor()
        numpy.testing.assert_allclose(RTi.data, RTi_noRWA.data)

        
        # verify that the initial condition was not changed
        numpy.testing.assert_allclose(rhoi_t.data, rhoi_c.data)
        numpy.testing.assert_allclose(rhot_i.data[0,:,:], rhoi_t.data[:,:])
        with eigenbasis_of(HH):
            numpy.testing.assert_allclose(rhot_i.data[0,:,:], rhoi_t.data[:,:])
            
        #verifying some time points
        attime = UUi.apply(0.0, rhoi_t)
        attime_a = UUi.apply(0.0, rhoi_a)
        numpy.testing.assert_allclose(rhot_i.data[0,:,:], attime.data)
        numpy.testing.assert_allclose(rhot_a.data[0,:,:], attime_a.data)
        attime = UUi.apply(10.0, rhoi_t)
        attime_a = UUi.apply(10.0, rhoi_a)
        numpy.testing.assert_allclose(rhot_i.data[10,:,:], attime.data)
        numpy.testing.assert_allclose(rhot_a.data[10,:,:], attime_a.data)
        attime = UUi.apply(100.0, rhoi_t)
        attime_a = UUi.apply(100.0, rhoi_a)
        numpy.testing.assert_allclose(rhot_i.data[100,:,:], attime.data)
        numpy.testing.assert_allclose(rhot_a.data[100,:,:], attime_a.data)
        attime = UUi.apply(900.0, rhoi_t)
        attime_a = UUi.apply(900.0, rhoi_a)
        numpy.testing.assert_allclose(rhot_i.data[900,:,:], attime.data)
        numpy.testing.assert_allclose(rhot_a.data[900,:,:], attime_a.data)


        #
        #
        # Verification against data without RWA
        #
        #
        with eigenbasis_of(HH):
            numpy.testing.assert_allclose(numpy.real(rhot_i.data[:,n1,n2]), 
                                          numpy.real(saved_data[0,:]),
                                          atol=3.0e-4)
            numpy.testing.assert_allclose(numpy.real(rhot_a.data[:,0,n2]), 
                                          numpy.real(saved_data[1,:]),
                                          atol=3.0e-4)        

        _plot_ = False
        if _plot_:
            import matplotlib.pyplot as plt
            
            with eigenbasis_of(HH):
                n1 = 1
                n2 = 2
                #plt.plot(time_eop.data, numpy.real(dat_t[:, n1, n2, n1, n2]))
                plt.plot(time_eop.data, 
                         numpy.real(UUi.data[:, 0, n2, 0, n2]),"-b")
                plt.plot(time.data, numpy.real(rhot_a.data[:, 0, n2]),"--")
                #plt.plot(time.data, numpy.real(saved_data[1, :]),"-k")
                # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,2,2,2]))
                # plt.plot(time_eop.data, numpy.real(UU.data[:, 3,3,3,3]))
                # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,3,2,3]))
                # plt.plot(time_eop.data, numpy.real(UU.data[:, 2,0,2,0]))
            
            plt.show()
            

    def test_block_calculation(self):
        """(EvolutionSuperOperator) Test block calculation with Redfield tensor
        
        """
        #
        # 1) HOW ARE DIFFERENT TIME AXES HANDLED?
        # 2) Time dependent Redfield gives different results with different axes
        # 3) Is RWA handled correctly in EvolutionSuperOperator ? 
        # 4) We have to add is_in_rwa property and convert_to and from RWA
        #

        # U1 = EvolutionSuperOperator(time=time_eop, ham=HH, 
        #                             relt=RTt, block=(0,1))
        # U1.set_dense_dt(10)
        
        # U1.calculate()

        pass


        
        
    def test_redfield_dynamics_comp(self):
        """(EvolutionSuperOperator) Compares Redfield dynamics calculated from propagator and superoperator
        
        
        
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
                
                
    def test_Lindblad_dynamics_comp(self):
        """(EvolutionSuperOperator) Compares Lindblad dynamics calculated from propagator and superoperator
        
        
        
        """
        # Aggregate
        import quantarhei as qr
        import quantarhei.models.modelgenerator as mgen
        
        time = qr.TimeAxis(0.0, 1000, 1.0)
        
        # create a model
        mg = mgen.ModelGenerator()
        
        agg = mg.get_Aggregate(name="trimer-1")
        
        agg.build()
        
        ham = agg.get_Hamiltonian()

        # calculate relaxation tensor

        with qr.eigenbasis_of(ham):
        
            #
            # Operator describing relaxation
            #
            from quantarhei.qm import Operator

            K = Operator(dim=ham.dim,real=True)
            K.data[1,2] = 1.0        

            #
            # System bath interaction with prescribed rate
            #
            from quantarhei.qm import SystemBathInteraction

            sbi = SystemBathInteraction(sys_operators=[K], rates=(1.0/100.0,))
            agg.set_SystemBathInteraction(sbi)    

            #
            # Corresponding Lindblad form
            #
            from quantarhei.qm import LindbladForm

            LF = LindbladForm(ham, sbi, as_operators=False)       
 
            #
            # Evolution of reduced density matrix
            #
            prop = qr.ReducedDensityMatrixPropagator(time, ham, LF)        
            
            #
            # Evolution by superoperator
            #
            
            eSO = qr.qm.EvolutionSuperOperator(time, ham, LF)
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
                
                
    def test_LindbladWithVibrations_dynamics_comp(self):
        """(EvolutionSuperOperator) Compares Lindblad dynamics of a system with vibrations calculated from propagator and superoperator
        
        
        
        """
        # Aggregate
        import quantarhei as qr
        
        time = qr.TimeAxis(0.0, 1000, 1.0)
        
        # create a model
        with qr.energy_units("1/cm"):
            
            me1 = qr.Molecule([0.0, 12100.0])
            me2 = qr.Molecule([0.0, 12000.0])
            me3 = qr.Molecule([0.0, 12900.0])
            
            agg_el = qr.Aggregate([me1, me2, me3])
            
            agg_el.set_resonance_coupling(0, 1, qr.convert(150, "1/cm", to="int"))
            agg_el.set_resonance_coupling(1, 2, qr.convert(50, "1/cm", to="int"))
 
            m1 = qr.Molecule([0.0, 12100.0])
            m2 = qr.Molecule([0.0, 12000.0])
            m3 = qr.Molecule([0.0, 12900.0])
           
            mod1 = qr.Mode(frequency=qr.convert(100, "1/cm", "int"))
            m1.add_Mode(mod1)
            mod1.set_HR(1, 0.01)
            
            agg = qr.Aggregate([m1, m2, m3])
            
            agg.set_resonance_coupling(0, 1, qr.convert(150, "1/cm", to="int"))
            agg.set_resonance_coupling(1, 2, qr.convert(50, "1/cm", to="int"))
            
        
        agg_el.build()
        agg.build()
        
        hame = agg_el.get_Hamiltonian()
        ham = agg.get_Hamiltonian()

        # calculate relaxation tensor

        with qr.eigenbasis_of(hame):
        
            #
            # Operator describing relaxation
            #

            K = qr.qm.ProjectionOperator(1, 2, dim=hame.dim)      

            #
            # System bath interaction with prescribed rate
            #
            from quantarhei.qm import SystemBathInteraction

            sbi = SystemBathInteraction(sys_operators=[K], rates=(1.0/100.0,))
            sbi.set_system(agg) #agg.set_SystemBathInteraction(sbi) 
            
            
        with qr.eigenbasis_of(ham):

            #
            # Corresponding Lindblad form
            #
            from quantarhei.qm import ElectronicLindbladForm

            LF = ElectronicLindbladForm(ham, sbi, as_operators=True)       
            
            
            #
            # Evolution of reduced density matrix
            #
            prop = qr.ReducedDensityMatrixPropagator(time, ham, LF)        
            
            
            
            #
            # Evolution by superoperator
            #
            
            eSO = qr.qm.EvolutionSuperOperator(time, ham, LF)
            eSO.set_dense_dt(5)
            eSO.calculate()
            
            # compare the two propagations
            pairs = [(5,4), (5,5), (6,5), (7,5)]
            for p in pairs:
                
            
                rho_i1 = qr.ReducedDensityMatrix(dim=ham.dim)
                rho_i1.data[p[0],p[1]] = 1.0
            
                rho_t1 = prop.propagate(rho_i1)

            
                exp_rho_t2 = eSO.data[:,:,:,p[0],p[1]]


                #import matplotlib.pyplot as plt
        
                #plt.plot(rho_t1.TimeAxis.data, numpy.real(rho_t1.data[:,p[0],p[1]]), "-r")
                #plt.plot(rho_t1.TimeAxis.data, numpy.real(exp_rho_t2[:,p[0],p[1]]), "--g")
                #plt.show()
                
                #for kk in range(rho_t1.TimeAxis.length):
                #    print(kk, numpy.real(rho_t1.data[kk,p[0],p[1]]),numpy.real(exp_rho_t2[kk,p[0],p[1]]))

                #numpy.testing.assert_allclose(RRT.data, rtd)
                numpy.testing.assert_allclose(numpy.real(rho_t1.data[:,:,:]),
                                          numpy.real(exp_rho_t2[:,:,:]), 
                                          rtol=5.0e-2,
                                          atol=1.0e-3)
                


if __name__ == '__main__':
    unittest.main()
    
    
    