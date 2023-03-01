# -*- coding: utf-8 -*-

import unittest

import numpy
import matplotlib.pyplot as plt

import quantarhei as qr
from quantarhei.utils.vectors import X
from quantarhei.utils.vectors import Y
from quantarhei.utils.vectors import Z

from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import energy_units
from quantarhei import ReducedDensityMatrixPropagator


class TestLabSetup(unittest.TestCase):
    """Test of the laboratory setup 
    
    """
    
    def setUp(self):
        
        self._plot_ = False
        
        ################################################################################
        #
        # Set up laboratory/experiment configuration
        #
        ################################################################################
        
        # three pulses will be used
        lab = qr.LabSetup(nopulses=3)
        
        # on a time axis starting at a specified time, with a certain number of steps
        # and a step size
        time = qr.TimeAxis(-500.0, 1500, 1.0, atype="complete")

        # pulse shapes are specified below
        pulse2 = dict(ptype="Gaussian", FWHM=150, amplitude=1.0)
        params = (pulse2, pulse2, pulse2)
        
        with self.assertRaises(Exception) as context:
        
            lab.set_pulse_shapes(time, params)
            
        self.assertTrue("Pulse arrival times have to specified" 
                        in str(context.exception))
        
        
        # each pulse has a defined frequency
        lab.set_pulse_polarizations(pulse_polarizations=(X,Y,Z), detection_polarization=X)
        
        # time of arrival
        lab.set_pulse_arrival_times([0.0, 0.0, 100.0])
        
        # and polarization
        lab.set_pulse_frequencies([1.0, 1.0, 1.0])
        
        # additional phases can be also controlled
        lab.set_pulse_phases([0.0, 1.0, 0.0])    
        
        lab.set_pulse_shapes(time, params)

        self.lab = lab
        self.time = time

        with energy_units("1/cm"):
            m1 = Molecule([0.0, 10000.0])
            m1.set_dipole((0,1), [1.0, 0.0, 0.0])
            m2 = Molecule([0.0, 10000.0])
            m2.set_dipole((0,1), [1.0, 0.0, 0.0])
            
            agg = Aggregate(molecules=[m1, m2])
        
            agg.set_resonance_coupling(0,1, 200.0)
            
        agg.build()
        
        
        self.agg = agg


    
    def test_lab_pulse_setters(self):
        """(LabSetup) Testing LabSetup pulse properties setters
        
        """
        lab = self.lab
        time = self.time

        fields = lab.get_labfields()
        
        fld = lab.get_labfield(1)
        
        self.assertTrue(fld.om == fields[1].om)
        self.assertTrue(fld.phi == fields[1].phi)
        
        
        _plot_ = self._plot_
        
        if _plot_:
            
            fields[2].set_rwa(0.9)
            fields[2].phi = numpy.pi

            fld = fields[2].get_field()
            plt.plot(time.data, numpy.real(fld))
            plt.show()
            
            fields[2].tc = 300.0

            fld = fields[2].get_field()
            plt.plot(time.data, numpy.real(fld))
            plt.show()
            

    def test_dm_propagation_with_fields(self):
        """(LabSetup) Time evolution with explicit electric field
        
        """
        from quantarhei.qm import LindbladForm, SystemBathInteraction
        from quantarhei.qm import Operator
        
        lab = self.lab
        time = self.time
        agg = self.agg


        
        HH = agg.get_Hamiltonian()
        DD = agg.get_TransitionDipoleMoment()
        
        ops = []
        KK = Operator(dim=HH.dim)
        KK.data[1,2] = 1.0
        ops.append(KK)
        rates = []
        rates.append(1.0/100.0)
        
        SBI = SystemBathInteraction(sys_operators=ops, rates=rates, system=agg)
        
        LT = LindbladForm(HH, SBI, as_operators=False)
        print(LT.data[1,1,2,2], 1.0/100.0)
        print(LT.data[1,2,1,2], 1.0/200.0)
        
        ef = lab.get_labfield(0)
        
        rhoi = agg.get_thermal_ReducedDensityMatrix()
        
        prop = ReducedDensityMatrixPropagator(timeaxis=time, Ham=HH,
                                              Efield=ef, Trdip=DD, 
                                              RTensor=LT)
        
        
        #
        # propagation has to be reimplemented with LabFields
        #
        #rhot = prop.propagate(rhoi)
        
        
        #raise Exception("STOP HERE")
        



if __name__ == '__main__':
    unittest.main()