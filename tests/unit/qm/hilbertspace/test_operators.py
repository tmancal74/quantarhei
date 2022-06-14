# -*- coding: utf-8 -*-

import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.Molecule class


*******************************************************************************
"""

from quantarhei import Molecule
from quantarhei.qm import ReducedDensityMatrix

class TestReducedDensityMatrix(unittest.TestCase):
    """Tests for the Manager class
    
    
    """
    
    def setUp(self):
        self.en = [0.0, 1.0, 2.0]
        self.m = Molecule(name="Molecule",elenergies=self.en)
        self.m.set_dipole(0,1,[1.0,0.0,0.0])
        self.m.set_dipole(0,2,[0.5,0.0,0.0])
    
        self.rho_eq = self.m.get_thermal_ReducedDensityMatrix()
    
    def test_excitation_by_delta(self):
        """Testing reduced density matrix excitation by delta-pulse
        
        
        """

        
        
        rho_exp = numpy.zeros((3,3),dtype=numpy.float)
        rho_exp[0,0] = 1.0
        
        self.assertTrue(numpy.allclose(self.rho_eq._data,rho_exp))
        
        dd = self.m.get_TransitionDipoleMoment()
        epol = [1.0, 0.0, 0.0]        
        
        rho_ex = self.rho_eq.excite_delta(dmoment=dd,epolarization=epol)
            
        rho_exp = numpy.array([[ 0.00+0.j,  0.00+0.j,  0.00+0.j],
                               [ 0.00+0.j,  1.00+0.j,  0.50+0.j],
                               [ 0.00+0.j,  0.50+0.j,  0.25+0.j]])
                               
        self.assertTrue(numpy.allclose(rho_exp,rho_ex._data))
        
        
        
        rho_ex2 = self.rho_eq.excite_delta(dmoment=dd,
                                      epolarization=[1.0/numpy.sqrt(2.0),
                                                     1.0/numpy.sqrt(2.0),0.0])
                                                     
        rho_exp2 = numpy.array([[ 0.000+0.j,  0.000+0.j,  0.000+0.j],
                                [ 0.000+0.j,  0.500+0.j,  0.250+0.j],
                                [ 0.000+0.j,  0.250+0.j,  0.125+0.j]])
                                
        self.assertTrue(numpy.allclose(rho_ex2._data,rho_exp2))
        
        
    def test_conversion_2_populations(self):
        """Testing conversion of reduced density matrix to population vector
        
        
        """
        
        rdm = ReducedDensityMatrix(data=[[0.5, 0.0, 0.1],
                                         [0.0, 0.3, 0.0],
                                         [0.1, 0.0, 0.2]])
                                         
        pop = rdm.get_populations()
        
        self.assertTrue(numpy.allclose(pop,[0.5,0.3,0.2]))
        
        
    def test_saveable(self):
        """Testing reduced density matrix as Saveable
        
        
        """
        
        rdm = ReducedDensityMatrix(data=[[0.5, 0.0, 0.1],
                                         [0.0, 0.3, 0.0],
                                         [0.1, 0.0, 0.2]])
                         
        #import h5py
              
        #with h5py.File("test_file_operators",driver="core", 
        #                   backing_store=False) as fid:
        import tempfile
        with tempfile.TemporaryFile() as fid:
            
            rdm.save(fid) #, test=True)
            fid.seek(0)
            
            rdm2 = ReducedDensityMatrix()
            
            rdm2 = rdm2.load(fid) #, test=True)
        
        self.assertTrue(numpy.allclose(rdm2.data,rdm.data))
                                    