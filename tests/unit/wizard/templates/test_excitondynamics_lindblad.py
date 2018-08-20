# -*- coding: utf-8 -*-

import pkg_resources

import unittest
import numpy

from quantarhei import load_parcel as load

class TestLindbladDynamics(unittest.TestCase):
    """Tests for the units package
    
    
    """

    def test_that_ExcitonDynamics_Lindblad_template_works(self):
        """Test that ExcitonDynamics_Lindblad template works
        
        
        """
        
        from quantarhei.wizard.templates.excitondynamics_lindblad \
             import mySimulation
        
        ms = mySimulation()
        ms.run()

        # get file for comparison
        filename = pkg_resources.resource_filename(__package__,   
                                    "excitondynamics_lindblad_rhot.qrp_test")
        # load object for comparison
        rhot = load(filename)

        # test that data aggree
        numpy.testing.assert_array_equal(rhot.data, ms._rhot.data)

        