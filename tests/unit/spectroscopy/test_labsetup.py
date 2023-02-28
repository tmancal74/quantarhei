# -*- coding: utf-8 -*-

import unittest

import numpy
import matplotlib.pyplot as plt

import quantarhei as qr
from quantarhei.utils.vectors import X
from quantarhei.utils.vectors import Y
from quantarhei.utils.vectors import Z


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
        time = qr.TimeAxis(-500.0, 1000, 1.0, atype="complete")

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

    
    def test_lab_pulse_setters(self):
        """(LabSetup) Testing LabSetup pulse properties setters
        
        """
        lab = self.lab
        time = self.time


        

    def test_get_field(self):
        """ (LabSetup) Testing get_field with rwa 
        """

        _plot_ = self._plot_
        lab = self.lab
        time = self.time
        
        numpy.testing.assert_array_equal(lab.e[0,:], X)
        numpy.testing.assert_array_equal(lab.e[1,:], Y)
        
        Ef = lab.get_field(2, rwa=0.9)
        ed = Ef.field()
        
        if _plot_:

            plt.plot(time.data, ed)
            plt.show()
        



if __name__ == '__main__':
    unittest.main()