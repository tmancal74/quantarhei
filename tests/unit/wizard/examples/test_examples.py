# -*- coding: utf-8 -*-
import pkg_resources
import unittest



class TestExamples(unittest.TestCase):
    """Tests for the units package
    
    
    """

    def test_examples_001_005(self):
        """Testing examples 001 - 005
        """

        # example files
        efiles = ["ex_001_Molecule_Hamiltonian.py",
                  "ex_002_Molecule_Aggregate.py",
                 ]
        
        do_the_work(efiles)

    def test_examples_006_009(self):
        """Testing examples 006 - 009
        """

        # example files
        efiles = [
                  "ex_006_Absorption_1.py",
                 ]
        
        do_the_work(efiles)
        

    def test_examples_010_019(self):
        """Testing examples 010 - 019
        """

        # example files
        efiles = [
                  "ex_010_RedfieldTheory_1.py",
                  "ex_011_LindbladForm_1.py"
                 ]
        
        do_the_work(efiles)
        
        
def do_the_work(efiles):

    for efile in efiles:
        resource_path = '/'.join(('wizard', 'examples', efile))
        filename = pkg_resources.resource_filename("quantarhei", resource_path)
        try:
            exec(open(filename).read())
        except:
            raise Exception("Example "+efile+" failed")
        
        
        
