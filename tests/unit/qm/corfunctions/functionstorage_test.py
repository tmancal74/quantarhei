# -*- coding: utf-8 -*-

import unittest
import numpy.testing as npt

#import quantarhei as qr

from quantarhei import TimeAxis
from quantarhei import FunctionStorage

class TestFunctionStorage(unittest.TestCase):
    """Tests the FunctionStorage functionality
    
    
    """
    
    def test_function_storage_1(self):
        """Testing FunctionStorage test 1
        
        """


        t1 = TimeAxis(0.0, 100, 1.0)
        t3 = TimeAxis(0.0, 110, 1.0)
        
        # FIXME: Check the dimensionality of the storage 
        # - it must be consistent with the self.configuration
        gg = FunctionStorage(7, (t1,t3), show_config=False, config=None)

        def goft_gamma(t):
            return t/100.0
        
        def goft_delta(t):
            return (t/50.0)**2
        
        list_fce = [goft_delta, goft_gamma]
        
        gg.set_goft(1, func=goft_gamma)
        gg.set_goft([4,2], func=goft_delta)
        gg.set_goft(3, func=goft_delta)
        gg.set_goft(12, func=goft_gamma)
        gg.set_goft(0, func=goft_delta)
        
        # FIXME: This does not work
        #gg.set_goft([5,14], func=goft_even)
        
        print("Number of functions:", gg.Nf)
        print("Effective size", gg.effective_size())
        print("Mapping", gg.mapping)
        
        print("Stored functions:")
        print(gg.funcs)
        
        #
        # Creating data
        #
        gg.create_data(reset=dict(t2=5.0))
        
        print("DATA SIZE: ", gg.data_size, "MB")
        
        for kk, fc in enumerate(list_fce):
            print("\nTesting function no:", kk)
            
            print("\n0D, t2")
            t2 = 5.0
            ggt2 = gg[kk,"t2"]
            npt.assert_allclose(ggt2, list_fce[kk](t2))
            print("OK")
            
            print("\n1D, t1")
            ggt1 = gg[kk,"t1"]
            npt.assert_allclose(ggt1, list_fce[kk](t1.data))
            print("OK")
            
            print("\n1D, t3")
            ggt3 = gg[kk, gg.time_mapping["t3"]]
            npt.assert_allclose(ggt3, list_fce[kk](t3.data))
            print("OK")
            
            print("\n1D, t1+t2")
            ggt1t2 = gg[kk,"t1+t2"]
            npt.assert_allclose(ggt1t2, list_fce[kk](t2+t1.data),
                                rtol=1.0e-6)
            print("OK")
            
            print("\n2D, t1+t2+t3")
            ggt1t2t3 = gg[kk,"t1+t2+t3"]
            npt.assert_allclose(ggt1t2t3, 
                                gg[kk,gg.time_mapping["t1+t2+t3"]],
                                rtol=1.0e-6)
            npt.assert_allclose(ggt1t2t3, 
                        list_fce[kk](t2+t1.data[:,None]
                        +t3.data[None,:]).reshape(t1.length*t3.length),
                        rtol=1.0e-6)
            print("OK")


    def test_function_storage_2(self):
        """Testing FunctionStorage test 2
        
        """    
        t1 = TimeAxis(0.0, 100, 1.0)
        t3 = TimeAxis(0.0, 110, 1.0)

        def goft_gamma(t):
            return t/100.0
        
        def goft_delta(t):
            return (t/50.0)**2

        config = {"response_times":{
                    "t1":{"reset":False,"integral":{"s1"},"axis":0},
                    "t2":{"reset":True, "integral":None},
                    "t3":{"reset":False,"integral":{"s3"},"axis":1},
                    "t4":{"reset":False,"integral":None,"axis":2}}}
        
        gg2 = FunctionStorage(7, (t1,t3,t3), show_config=False, config=config)
        
        t2 = 10.0
        
        gg2.set_goft(0, func=goft_gamma)
        gg2.set_goft(1, func=goft_delta)
        gg2.create_data(reset=dict(t2=t2))
        print("DATA SIZE: ", gg2.data_size, "MB")
        list_fce = [goft_gamma, goft_delta]
        
        for kk in range(2):
        
            print("\n0D, t2")
            t2 = 10.0
            ggt2 = gg2[kk,"t2"]
            print("OK")
            
            print("\n1D, t1")
            ggt1 = gg2[kk,"t1"]
            npt.assert_allclose(ggt1, list_fce[kk](t1.data))
            print("OK")
            
            print("\n1D, t3")
            ggt3 = gg2[kk, gg2.time_mapping["t3"]]
            npt.assert_allclose(ggt3, list_fce[kk](t3.data))
            print("OK")
            
            print("\n1D, t1+t2")
            ggt1t2 = gg2[kk,"t1+t2"]
            npt.assert_allclose(ggt1t2, list_fce[kk](t2+t1.data),rtol=1.0e-6)
            print("OK")
            
            print("\n2D, t1+t2+t3")
            ggt1t2t3 = gg2[kk,"t1+t2+t3"]
            npt.assert_allclose(ggt1t2t3, 
                                gg2[kk,gg2.time_mapping["t1+t2+t3"]],
                                rtol=1.0e-5)
            npt.assert_allclose(ggt1t2t3, 
                        list_fce[kk](t2+t1.data[:,None]
                        +t3.data[None,:]).reshape(t1.length*t3.length),
                        rtol=1.0e-6)
            print("OK")
