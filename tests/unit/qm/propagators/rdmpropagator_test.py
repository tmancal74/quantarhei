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

from quantarhei import Molecule, Mode
from quantarhei import CorrelationFunction
from quantarhei import energy_units, eigenbasis_of

from quantarhei import LabSetup, convert
from quantarhei import ReducedDensityMatrix

class TestRDMPropagator(unittest.TestCase):
    """Tests reduced density matrix propagator module
    
    
    """
    
    def test_rdm_evolution_time_independent_tensor(self):
        """(RDMPropagator) Testing evolution of reduced density matrix with Lindblad relaxation tensor
        
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
        
        
    def test_tdp(self):
        """(RDMPropagator) Testing with time dependent Redfield tensor
        
        """
        
        expstring = """
            0 4 4 -3.96678714134e-19 0.0
            25 4 4 0.0123395314303 1.10216297674e-16
            50 4 4 0.0521821949151 4.41386645353e-16
            75 4 4 0.0375170438965 3.03228186587e-16
            0 3 3 -2.15151658263e-18 0.0
            25 3 3 0.0542304209663 3.27264685377e-16
            50 3 3 0.133087497663 5.56816047005e-16
            75 3 3 0.114486463162 5.08756535778e-16
            0 0 5 -0.855304990279 0.0
            25 0 5 0.338783123139 0.205218959392
            50 0 5 -0.0404678523217 -0.0413353312194
            75 0 5 -0.0747762812157 0.00493941699829
        """
        
        with energy_units("1/cm"):
            mol1 = Molecule([0.0, 12000.0])
            mod1 = Mode(300.0)
            mol1.add_Mode(mod1)
            mod1.set_nmax(0, 1)
            mod1.set_nmax(1, 8)
            mod1.set_HR(1,0.2)
            
            mol1.set_dipole((0,1), [1.0, 0.0, 0.0])
            
            time = TimeAxis(0.0, 100, 1.0)
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
        
        rhoi = mol1.get_excited_density_matrix(condition="delta")
        
        prop = mol1.get_ReducedDensityMatrixPropagator(time,
                                                relaxation_theory="stR",
                                                time_dependent=True)
        HH = mol1.get_Hamiltonian()
        #rhoi.data[:,:] = 0.0
        with eigenbasis_of(HH):
            rhoi.data[0,4] = 1.0
            rhoi.data[0,3] = 1.0
            
        rhot1 = prop.propagate(rhoi, Nref=1)

        HH.set_rwa([0,1])    
            
        rhot2 = prop.propagate(rhoi, Nref=1)
        
        rhot2.convert_from_RWA(HH)
        
        #
        # checking that rotating wave approximation is the same as standard 
        #
        #numpy.testing.assert_allclose(rhot1.data, rhot2.data, rtol=8.0e-2)
        
        
        _show_plot_ = False
        if _show_plot_:
            import matplotlib.pyplot as plt
            with eigenbasis_of(HH):
                for ii in range(HH.dim):
                    plt.plot(time.data, numpy.real(rhot1.data[:,ii,ii]))
                    plt.plot(time.data, numpy.real(rhot2.data[:,ii,ii]),"--")
            
            plt.show()
        
        _create_data_ = False
        if _create_data_:
            
            elements = [(4,4), (3,3), (0,5)]
            for el in elements:
                for tt in range(0, 100, 25):
                    print(tt, el[0], el[1],
                              numpy.real(rhot2.data[tt,el[0],el[1]]),
                              numpy.imag(rhot2.data[tt,el[0],el[1]]))
                    
        #  
        # compare to precalculated data
        #
        expected = {}
        il = 0
        for line in expstring.splitlines():
            rel = line.strip()
            if len(rel) > 0:
                nmbrs = rel.split(" ")
                tmi = int(nmbrs[0])
                i1 = int(nmbrs[1])
                i2 = int(nmbrs[2])
                re = float(nmbrs[3])
                im = float(nmbrs[4])
                expected[il] = (tmi, i1, i2, re, im)
                il += 1

        _perform_test_ = True
        if _perform_test_:
            for ks in expected:
                dats = expected[ks]
                tt = dats[0]
                i1 = dats[1]
                i2 = dats[2]
                numpy.testing.assert_allclose(dats[3], 
                                       numpy.real(rhot2.data[tt,i1,i2]),
                                       atol=1.0e-6)
                numpy.testing.assert_allclose(dats[4], 
                                       numpy.imag(rhot2.data[tt,i1,i2]),
                                       atol=1.0e-6)



    #def test_rdm_evolution_Saveable(self):
    #    pass


    def test_rdm_field(self):
        """(RDMPropagator) Propagation with external field
        
        """
        
        with energy_units("1/cm"):
            mol1 = Molecule([0.0, 12000.0])
            mol1.set_dipole((0,1), [1.0, 0.0, 0.0])
            time = TimeAxis(0.0, 100, 1.0, atype="complete")
            params = dict(ftype="OverdampedBrownian",
                          cortime=30.0, reorg=20, T=300,
                          matsubara=30)
            
            cfce = CorrelationFunction(time, params=params)
            mol1.set_transition_environment((0,1), cfce)

        HH = mol1.get_Hamiltonian()
        DD = mol1.get_TransitionDipoleMoment()
        
        dd = DD.data[1,:,:]
           
        with self.assertRaises(Exception) as context:
            prop = ReducedDensityMatrixPropagator(time, Ham=HH, Trdip=dd)
        
        self.assertTrue("Operator or None expected here."
                        in str(context.exception))
            
        lab = LabSetup(nopulses=3)
        om = convert(12000.0, "1/cm", "int")
        omegas = [om, om, om]
        lab.set_pulse_frequencies(omegas)
        lab.set_pulse_arrival_times([100.0, 100.0, 100.0])
        lab.set_pulse_polarizations([[1.0, 0.0, 0.0],[1.0, 0.0, 0.0],
                                     [1.0, 0.0, 0.0]])
        prms = dict(ptype="Gaussian", FWHM=50.0, amplitude=1.0)
        lab.set_pulse_shapes(time, params=[prms, prms, prms])
        
        f1 = lab.get_labfield(0)
        
        
        #
        # No relaxation, field as array
        #
        print("\n+++ 14")
        
        prop = ReducedDensityMatrixPropagator(time, Ham=HH,
                                              Trdip=DD, Efield=f1.field)
        
        rhoi = ReducedDensityMatrix(dim=HH.dim)
        rhoi.data[0,0] = 1.0
        
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi) #, method="short-exp-4")
        self.assertTrue("NOT IMPLEMENTED"
                        in str(context.exception))
        
        
        
        #
        # No relaxation, field as object
        #
        
        print("\n+++ 15")
        prop = ReducedDensityMatrixPropagator(time, Ham=HH,
                                              Trdip=DD, Efield=f1)
        
        rhoi = ReducedDensityMatrix(dim=HH.dim)
        rhoi.data[0,0] = 1.0
        
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi) #, method="short-exp-4")
        self.assertTrue("NOT IMPLEMENTED"
                        in str(context.exception))        
        
        
        
        #
        # No field, no relaxation
        #
        
        print("\n+++ 1")
        HH.set_rwa([0, 1])
        
        prop = ReducedDensityMatrixPropagator(time, Ham=HH)
        rhoi = ReducedDensityMatrix(dim=HH.dim)
        rhoi.data[:,:] = 0.5

        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="runge-kutta")
        self.assertTrue("Unknown propagation method"
                        in str(context.exception))            

        rhot = prop.propagate(rhoi, method="short-exp-4")  
        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        rhot.convert_from_RWA(HH)
        

        # import matplotlib.pyplot as plt
        
        # plt.plot(time.data, rhot.data[:,0,1])
        # plt.show()
        


if __name__ == '__main__':
    unittest.main()
    
    