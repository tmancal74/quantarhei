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
from quantarhei import ReducedDensityMatrix
from quantarhei import TimeAxis
from quantarhei.qm import Liouvillian, ProjectionOperator
from quantarhei import LabSetup

from quantarhei import Molecule, Mode, Aggregate
from quantarhei import CorrelationFunction
from quantarhei import energy_units, eigenbasis_of


class TestRDMPropagatorUsage(unittest.TestCase):
    """Tests reduced density matrix propagator module
    
    
    """
    
    def setUp(self):
        
        with energy_units("1/cm"):
            mol1 = Molecule([0.0, 12000.0])
            mol1.set_dipole((0,1), [1.0, 0.0, 0.0])
            mol2 = Molecule([0.0, 12100.0])
            mol2.set_dipole((0,1), [1.0, 0.0, 0.0])
            
            agg = Aggregate(molecules=[mol1, mol2])
            
            agg.set_resonance_coupling(0, 1, 50.0)
            
        agg.build()
        
        
        HH = agg.get_Hamiltonian()
        DD = agg.get_TransitionDipoleMoment()
        P12 = ProjectionOperator(1, 2, dim=HH.dim)
        P21 = ProjectionOperator(2, 1, dim=HH.dim)
        rates = [1.0/100.0, 1.0/600.0]
        
        sbi = qr.qm.SystemBathInteraction(sys_operators=[P12,P21],
                                          rates=rates)
        LL = qr.qm.LindbladForm(HH, sbi)
        LL.convert_2_tensor()
        
        L1 = qr.qm.LindbladForm(HH, sbi)
        
        Nt = 500
        dt = 1.0
        time = TimeAxis(0.0, Nt, dt,atype="complete")
        
        lab = LabSetup(nopulses=3)
        ppar = dict(ptype="Gaussian", FWHM=20, amplitude=0.1)
        params = (ppar, ppar, ppar)
        lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
        lab.set_pulse_shapes(time, params)
        with eigenbasis_of(HH):
            om = HH.data[1,1] - HH.data[0,0]
        lab.set_pulse_frequencies([om, om, om])
        lab.set_pulse_polarizations([[1.0, 0.0, 0.0],
                                     [1.0, 0.0, 0.0],
                                     [1.0, 0.0, 0.0]])
        
        
        Efield = lab.get_labfield(0)
        
        self.HH = HH
        self.DD = DD
        self.LL = LL
        self.L1 = L1
        self.time = time
        self.lab = lab
        self.Efield = Efield
        self.agg = agg
        
        
    def test_prop_init_const(self):
        """(RDMPropagator) Testing initialization of a propagator (const. R)
        
        """
        
        time = self.time
        HH = self.HH
        LL = self.LL
        L1 = self.L1
        Efield = self.Efield
        DD = self.DD
        agg = self.agg
        
        with self.assertRaises(Exception) as context:
            prop = ReducedDensityMatrixPropagator()
        self.assertTrue("TimeAxis and Hamiltonian are required" 
                        in str(context.exception))
        
        with self.assertRaises(Exception) as context:
            prop = ReducedDensityMatrixPropagator([1.0, 2.0])
        self.assertTrue("TimeAxis expected here." 
                        in str(context.exception))        
        
        with self.assertRaises(Exception) as context:
            prop = ReducedDensityMatrixPropagator(timeaxis=[1.0, 2.0])
        self.assertTrue("TimeAxis expected here." 
                        in str(context.exception))        
            
        with self.assertRaises(Exception) as context:
            prop = ReducedDensityMatrixPropagator(timeaxis=time)
        self.assertTrue("Hamiltonian is required." 
                        in str(context.exception))

        with self.assertRaises(Exception) as context:
            prop = ReducedDensityMatrixPropagator(timeaxis=time,
                                                  Ham=[1.0, 2.0])
        self.assertTrue("Hamiltonian represented by a wrong type." 
                        in str(context.exception))
        
        with self.assertRaises(Exception) as context:
            prop = ReducedDensityMatrixPropagator(timeaxis=time,
                                                  Ham=HH, RTensor=[1.0, 2.0])
        self.assertTrue("RelaxationTensor or None expected here."
                        in str(context.exception))

        with self.assertRaises(Exception) as context:
            prop = ReducedDensityMatrixPropagator(timeaxis=time,
                                                  Ham=HH, Iterm=LL)
        self.assertTrue("RelaxationTensor has to be set first."
                        in str(context.exception))
        
        #
        # Call all predefined methods with EField object
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=LL, Efield=Efield,
                                              Trdip=DD)
        
        rhoi = agg.get_thermal_ReducedDensityMatrix()
        
        print("\n+++ 11")
        
        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))    
        
 
        #
        # Call all predefined methods with field as an array
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=LL, Efield=Efield.field,
                                              Trdip=DD)
        
        rhoi = agg.get_thermal_ReducedDensityMatrix()

        print("\n+++ 10")
        
        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))
        


        #
        # Call all predefined methods with EField object
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=L1, Efield=Efield,
                                              Trdip=DD)
        
        rhoi = agg.get_thermal_ReducedDensityMatrix()
        
        print("\n+++ 13")
        
        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))    
        
 
        #
        # Call all predefined methods with field as an array
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=L1, Efield=Efield.field,
                                              Trdip=DD)
        
        rhoi = agg.get_thermal_ReducedDensityMatrix()

        print("\n+++ 12")
        
        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))

        
        #
        # Call relaxation without field; all methods
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=LL, Efield=None,
                                              Trdip=DD)        
        
        rhoi = ReducedDensityMatrix(dim=HH.dim)
        rhoi.data[2,2] = 1.0

        print("\n+++ 2")

        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))        
        

        #
        # Call relaxation without field; all methods
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=L1, Efield=None,
                                              Trdip=DD)        
        
        rhoi = ReducedDensityMatrix(dim=HH.dim)
        rhoi.data[2,2] = 1.0

        print("\n+++ 3")

        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))        


        
        
    def test_prop_init_tdep(self):
        """(RDMPropagator) Testing initialization of a propagator (time-dep. R)
        
        """
        
        time = self.time
        HH = self.HH
        LL = self.LL
        Efield = self.Efield
        DD = self.DD
        agg = self.agg    
        
        params = dict(ftype="OverdampedBrownian", reorg=20.0, cortime=30.0,
                      T=300, matsubara=30)
        with energy_units("1/cm"):
            cf = CorrelationFunction(time, params)
            
        for mol in agg.monomers:
            
            mol.set_transition_environment((0,1), cf)
            
        agg.rebuild()
        
        RR, ham = agg.get_RelaxationTensor(time, relaxation_theory="stR",
                                        time_dependent=True)
        
        
        
        #
        # Call time-dependent relaxation without field; all methods
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=RR, Efield=None,
                                              Trdip=DD)        
        
        rhoi = ReducedDensityMatrix(dim=HH.dim)
        rhoi.data[2,2] = 1.0

        print("\n+++ 4")

        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))


        RR, ham = agg.get_RelaxationTensor(time, relaxation_theory="stR",
                                           time_dependent=True, 
                                           as_operators=True)
        
        
        
        #
        # Call time-dependent operator form relaxation without field
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=RR, Efield=None,
                                              Trdip=DD)        

        print("\n+++ 5")
        
        rhoi = ReducedDensityMatrix(dim=HH.dim)
        rhoi.data[2,2] = 1.0
        rhot = prop.propagate(rhoi, method="short-exp-4")
        
        #
        # Now the same with cutoff-time
        # 
        RR.set_cutoff_time(300.0)
        #
        # Call time-dependent operator form relaxation without field
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=RR, Efield=None,
                                              Trdip=DD)        
        
        rhoi = ReducedDensityMatrix(dim=HH.dim)
        rhoi.data[2,2] = 1.0
        rhot = prop.propagate(rhoi, method="short-exp-4")        
        


        #
        # Call time-dependent relaxation with field; all methods
        #
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=RR, Efield=Efield,
                                              Trdip=DD)        
        
        rhoi = agg.get_thermal_ReducedDensityMatrix()

        print("\n+++ 9")

        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))        


        #
        # Call time-dependent relaxation with field ; as an array
        # all methods
        #
        
        RR, ham = agg.get_RelaxationTensor(time, relaxation_theory="stR",
                                           time_dependent=True, 
                                           as_operators=False)        
        
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=RR, 
                                              Efield=Efield,
                                              Trdip=DD)        
        
        rhoi = agg.get_thermal_ReducedDensityMatrix()

        print("\n+++ 8")

        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))         
        


        #
        # Call time-dependent relaxation with field; all methods
        #
        RR, ham = agg.get_RelaxationTensor(time, relaxation_theory="stR",
                                           time_dependent=True, 
                                           as_operators=True)                
        
        
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=RR, 
                                              Efield=Efield.field,
                                              Trdip=DD)        
        
        rhoi = agg.get_thermal_ReducedDensityMatrix()

        print("\n+++ 7")

        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))        


        #
        # Call time-dependent relaxation with field ; as an array
        # all methods
        #
        
        RR, ham = agg.get_RelaxationTensor(time, relaxation_theory="stR",
                                           time_dependent=True, 
                                           as_operators=False)        
        
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              RTensor=RR, 
                                              Efield=Efield.field,
                                              Trdip=DD)        
        
        rhoi = agg.get_thermal_ReducedDensityMatrix()

        print("\n+++ 6")

        rhot = prop.propagate(rhoi, method="short-exp-2")
        rhot = prop.propagate(rhoi, method="short-exp-4")
        rhot = prop.propagate(rhoi, method="short-exp-6")
        
        with self.assertRaises(Exception) as context:
            rhot = prop.propagate(rhoi, method="short-exp-8")
        self.assertTrue("Unknown propagation method:"
                        in str(context.exception))         

    
    
        
        _show_plot_ = False
        
        if _show_plot_:
            
            import matplotlib.pyplot as plt
                
            plt.plot(time.data, numpy.real(rhot.data[:, 1, 1]))
            plt.show()
            
            
        
        
if __name__ == '__main__':
    unittest.main()