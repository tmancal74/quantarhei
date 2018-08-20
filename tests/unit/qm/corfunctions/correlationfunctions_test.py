# -*- coding: utf-8 -*-


import unittest
import numpy
#import h5py


"""
*******************************************************************************


    Tests of the quantarhei.qm.corfunctions.correlationfunctions module


*******************************************************************************
"""


from quantarhei import CorrelationFunction
from quantarhei import TimeAxis
from quantarhei import energy_units
from quantarhei import Manager

class TestCorrelationFunction(unittest.TestCase):
    """Tests corelation functions module
    
    
    """
    
    def test_of_correlation_function_addition(self):
        """(CorrelationFunction) Testing addition of CorrelationFunction objects """
        
        t = TimeAxis(0.0, 1000, 1.0)
        params1 = dict(ftype="OverdampedBrownian",
                       reorg = 30.0,
                       cortime = 100.0,
                       T = 300.0)
        params2 = dict(ftype="OverdampedBrownian",
                       reorg = 40.0,
                       cortime = 100.0,
                       T = 300.0)
        
        with energy_units("1/cm"):
            f1 = CorrelationFunction(t, params1)
            f2 = CorrelationFunction(t, params2)
        
        #
        # normal addition
        #
        f = f1 + f2
        
        sum_data = f1.data + f2.data
        sum_lamb = f1.lamb + f2.lamb
        sum_cutoff = max(f1.cutoff_time, f2.cutoff_time)
        sum_temp = f1.temperature
        
        self.assertEqual(f.lamb, sum_lamb)
        numpy.testing.assert_allclose(f.data, sum_data)
        self.assertEqual(f.cutoff_time, sum_cutoff)
        self.assertEqual(f.temperature, sum_temp)
        
        #self.assertFalse(f.is_analytical())
        self.assertTrue(f._is_composed)
        self.assertFalse(f._is_empty)
        
        
        #
        # inplace addition by function
        #
        f1.add_to_data(f2)
        
        self.assertEqual(f1.lamb, sum_lamb)
        numpy.testing.assert_allclose(f1.data, sum_data)       
        self.assertEqual(f1.cutoff_time, sum_cutoff)  
        self.assertEqual(f1.temperature, sum_temp)

        #self.assertFalse(f1.is_analytical())
        self.assertTrue(f1._is_composed)
        self.assertFalse(f1._is_empty)
        
        #
        # inplace addition by += operator
        #
        with energy_units("1/cm"):        
            f1 = CorrelationFunction(t, params1) # new object  
        fs = f1
        f1 += f2

        self.assertEqual(f1.lamb, sum_lamb)
        numpy.testing.assert_allclose(f1.data, sum_data)
        self.assertEqual(f1.cutoff_time, sum_cutoff)
        self.assertEqual(f1.temperature, sum_temp)  
        
        #self.assertFalse(f1.is_analytical())
        self.assertTrue(f1._is_composed)
        self.assertFalse(f1._is_empty)
        
        # test if inplace addition really happend
        self.assertEqual(fs, f1)
        
        
    def test_of_correlation_function_addition_T(self):
        """(CorrelationFunction) Testing that addition with different temperatures raises Exception
        
        """
        t = TimeAxis(0.0, 1000, 1.0)
        params1 = dict(ftype="OverdampedBrownian",
                       reorg = 30.0,
                       cortime = 100.0,
                       T = 300.0)
        params2 = dict(ftype="OverdampedBrownian",
                       reorg = 40.0,
                       cortime = 100.0,
                       T = 100.0)    
        

        with energy_units("1/cm"):        
            f1 = CorrelationFunction(t, params1)
            f2 = CorrelationFunction(t, params2)   
        
        with self.assertRaises(Exception):
            f = f1 + f2
        
        
    def test_of_multiple_addition(self):
        """(CorrelationFunction) Testing multiple addition of objects """
        
        t = TimeAxis(0.0, 1000, 1.0)
        
        params1 = dict(ftype="OverdampedBrownian",
                       reorg = 30.0,
                       cortime = 100.0,
                       T = 300.0)
        params2 = dict(ftype="OverdampedBrownian",
                       reorg = 40.0,
                       cortime = 100.0,
                       T = 300.0)
        params3 = dict(ftype="OverdampedBrownian",
                       reorg = 15.0,
                       cortime = 200.0,
                       T = 300.0)
        params4 = dict(ftype="OverdampedBrownian",
                       reorg = 10.0,
                       cortime = 50.0,
                       T = 300.0)     
        
        with energy_units("1/cm"):
            f1 = CorrelationFunction(t, params1)
            f2 = CorrelationFunction(t, params2)
            f3 = CorrelationFunction(t, params3)
            f4 = CorrelationFunction(t, params4)
        
        f = f1 + f2 + f3 + f4
        
        sum_data = f1.data + f2.data + f3.data + f4.data
        sum_lamb = f1.lamb + f2.lamb + f3.lamb + f4.lamb
        sum_cutoff = max(f1.cutoff_time, f2.cutoff_time,
                         f3.cutoff_time, f4.cutoff_time)
        sum_temp = f1.temperature        

        self.assertEqual(f.lamb, sum_lamb)
        numpy.testing.assert_allclose(f.data, sum_data)
        self.assertEqual(f.cutoff_time, sum_cutoff)
        self.assertEqual(f.temperature, sum_temp)
        
        #self.assertFalse(f.is_analytical())
        self.assertTrue(f._is_composed)
        self.assertFalse(f._is_empty)    
        
        #
        # Inplace
        #  
        with energy_units("1/cm"):
            f1 = CorrelationFunction(t, params1)
            f2 = CorrelationFunction(t, params2)
            f3 = CorrelationFunction(t, params3)
            f4 = CorrelationFunction(t, params4)
        fs = f1
        
        sum_data = f1.data + f2.data + f3.data + f4.data
        sum_lamb = f1.lamb + f2.lamb + f3.lamb + f4.lamb
        sum_cutoff = max(f1.cutoff_time, f2.cutoff_time,
                         f3.cutoff_time, f4.cutoff_time)
        sum_temp = f3.temperature 
        
        f1 += f2 + f3
        f1 += f4
        
        f = f1
        
        self.assertEqual(f.lamb, sum_lamb)
        numpy.testing.assert_allclose(f.data, sum_data)
        #print(f.cutoff_time, sum_cutoff)
        self.assertEqual(f.cutoff_time, sum_cutoff)
        self.assertEqual(f.temperature, sum_temp)
        
        #self.assertFalse(f.is_analytical())
        self.assertTrue(f._is_composed)
        self.assertFalse(f._is_empty)        

        # test if inplace addition really happend
        self.assertEqual(fs, f1)
        
        #
        # Loops 
        #
        with energy_units("1/cm"):
            f1 = CorrelationFunction(t, params1)
            f2 = CorrelationFunction(t, params1)
        
        f1_data = f1.data.copy()
        for i in range(5):
            f1 += f1
        
        with energy_units("1/cm"):
            self.assertEqual(f1.lamb, 
              32.0*Manager().convert_energy_2_internal_u(params1["reorg"]))
            
        numpy.testing.assert_allclose(f1.data, 32.0*f1_data)
        self.assertEqual(f1.temperature, params1["T"])  
        
        #self.assertFalse(f1.is_analytical())
        self.assertTrue(f1._is_composed)
        self.assertFalse(f1._is_empty)            
        
        with energy_units("1/cm"):
            f1 = CorrelationFunction(t, params1)
            for i in range(5):
                f1 += f2

            self.assertEqual(f1.lamb, 
                6.0*Manager().convert_energy_2_internal_u(params1["reorg"]))
        numpy.testing.assert_allclose(f1.data, 6.0*f1_data)
        self.assertEqual(f1.temperature, params1["T"])  
        
        #self.assertFalse(f1.is_analytical())
        self.assertTrue(f1._is_composed)
        self.assertFalse(f1._is_empty)  

        #self.assertEqual(f1.params["ftype"],"Value-defined")           
        
        
    def test_reorganization_energy_consistence(self):
        """(CorrelationFunction) Checking that reorganization energy is represented consistently
        
        """
        
        t = TimeAxis(0.0, 2000, 1.0)
        
        params1 = dict(ftype="OverdampedBrownian",
                       reorg = 30.0,
                       cortime = 100.0,
                       T = 300.0)
        params2 = dict(ftype="OverdampedBrownian",
                       reorg = 40.0,
                       cortime = 100.0,
                       T = 300.0)
        params3 = dict(ftype="OverdampedBrownian",
                       reorg = 15.0,
                       cortime = 200.0,
                       T = 300.0)
        params4 = dict(ftype="OverdampedBrownian",
                       reorg = 10.0,
                       cortime = 50.0,
                       T = 300.0)     
        with energy_units("1/cm"):
            f1 = CorrelationFunction(t, params1)
            f2 = CorrelationFunction(t, params2)
            f3 = CorrelationFunction(t, params3)
            f4 = CorrelationFunction(t, params4)     
        
        l1 = f1.measure_reorganization_energy()
        l2 = f1.lamb
        #print(l1, l2, abs(l1-l2)/(l1+l2))
        
        l1 = f3.measure_reorganization_energy()
        l2 = f3.lamb
        #print(l1, l2, abs(l1-l2)/(l1+l2))
        
        self.assertTrue(f1.reorganization_energy_consistent())
        self.assertTrue(f2.reorganization_energy_consistent())
        self.assertTrue(f3.reorganization_energy_consistent())
        self.assertTrue(f4.reorganization_energy_consistent())        
        
        for i in range(5):
            f1 += f1        
        
        self.assertTrue(f1.reorganization_energy_consistent())

        for i in range(5):
            f3 += f2
            
        self.assertTrue(f3.reorganization_energy_consistent())
        
    def test_of_correlation_function_as_Saveable(self):
        """(CorrelationFunction) Testing of saving """
        
        t = TimeAxis(0.0, 1000, 1.0)
        params1 = dict(ftype="OverdampedBrownian",
                       reorg = 30.0,
                       cortime = 100.0,
                       T = 300.0)
        params2 = dict(ftype="OverdampedBrownian",
                       reorg = 40.0,
                       cortime = 100.0,
                       T = 300.0)
        
        with energy_units("1/cm"):
            f1 = CorrelationFunction(t, params1)
            f2 = CorrelationFunction(t, params2)

        import tempfile
        with tempfile.TemporaryFile() as f:
        #with h5py.File("test_file_1",driver="core", 
        #                   backing_store=False) as f:
            
            f1.save(f)#, test=True)
            f.seek(0)
            f1_loaded = CorrelationFunction()
            f1_loaded = f1_loaded.load(f) #, test=True)
            
        with tempfile.TemporaryFile() as f:
        #with h5py.File("test_file_2",driver="core", 
        #                   backing_store=False) as f:
            
            f2.save(f) #, test=True)
            f.seek(0)
            f2_loaded = CorrelationFunction()
            f2_loaded = f2_loaded.load(f) #, test=True)
            
            
        numpy.testing.assert_array_equal(f1.data, f1_loaded.data)
        numpy.testing.assert_array_equal(f1.axis.data, f1_loaded.axis.data)        
 
        numpy.testing.assert_array_equal(f2.data, f2_loaded.data)
        numpy.testing.assert_array_equal(f2.axis.data, f2_loaded.axis.data)        
       
            
            