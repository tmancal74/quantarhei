# -*- coding: utf-8 -*-
import pkg_resources
import unittest



class TestExamples(unittest.TestCase):
    """Tests for the units package
    
    
    """

    def test_examples(self):
        """Testing examples
        
        """
        

        # example files
        efiles = ["ex_001_Molecule_Hamiltonian.py",
                  "ex_002_Molecule_Aggregate.py"
                 ]
        
        for efile in efiles:
            do_the_work(efile)
        
        
        
def do_the_work(efile):
    
    resource_path = '/'.join(('wizard', 'examples', efile))
    filename = pkg_resources.resource_filename("quantarhei", resource_path)
    try:
        exec(open(filename).read())
    except:
        raise Exception("Example "+efile+" failed")
        
        
        
