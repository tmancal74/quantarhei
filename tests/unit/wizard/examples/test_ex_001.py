# -*- coding: utf-8 -*-
import pkg_resources
import unittest



class TestLindbladDynamics(unittest.TestCase):
    """Tests for the units package
    
    
    """

    def test_example_001(self):
        """Testing example 001
        
        """
        

        # get file for comparison
        efile = "ex_001_Molecule_Hamiltonian.py"
        resource_path = '/'.join(('wizard', 'examples', efile))
        filename = pkg_resources.resource_filename("quantarhei", resource_path)
        exec(open(filename).read())
        
