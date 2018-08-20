# -*- coding: utf-8 -*-

import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.Aggregate class


*******************************************************************************
"""

legacy = False
import tempfile
from quantarhei.core.saveable import Saveable
        

class TSaveable(Saveable):

    def __init__(self):

        self.a = 1.0
        self.text = None
        self._txt = None
        self.b1 = False
        self.b2 = False
        self.b3 = False

        self.dat = None

    def set_a(self, a):
        self.a = a

    def set_str(self, text1, text2):
        self.text = text1
        self._txt = text2

    def set_bool(self, b1, b2, b3):
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3

    def set_data(self, dat):
        self.dat = numpy.array(dat)

    def set_saveable(self, obj):
        self.obj = obj


    def set_liple(self, liple): 
        """Sets a list or tuple
        """
        self.liple = liple
        
    def set_cyclic(self, ts):
        
        self.ts = ts
        self.ts.ts = self


class TestSaveable(unittest.TestCase):
    """Tests for the Saveable class
    
    
    """
    
    def setUp(self):
        
        self.obj1 = TSaveable()
        
        self.obj1.set_a(10.0)
        self.obj1.set_str("Hi, there", "Hello World")
        self.obj1.set_bool(True, False, True)
        
        dat = numpy.zeros((5,5), dtype=numpy.complex128)
        
        dat[3,4] = 1.3455 + 1j*0.3456
        
        self.obj1.set_data(dat)
        
        obj2 = TSaveable()
        obj2.set_str("I am from second tier", "and me, too")
        
        self.obj1.set_saveable(obj2)
        
        self.obj3 = TSaveable()
        self.obj3.set_liple([1, 2, 3])
        
        
        
    def testing_cyclic_reference(self):
        
        a = TSaveable()
        a.a = 1.0
        b = TSaveable()
        b.b = 2.0
        a.set_cyclic(b)
        
        if legacy:
            with h5py.File("test_file_1",driver="core", 
                               backing_store=False) as f:
    
                a.save(f, test=True)
            
                c = TSaveable()
                
                c.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                
                a.save(f)
                f.seek(0)
                
                c = TSaveable()
                c = c.load(f)
        
        self.assertEqual(id(a), id(a.ts.ts))
        self.assertEqual(id(c), id(c.ts.ts))

        lst = [a,b, c]
        
        c.set_liple(lst)
        
        if legacy:
            with h5py.File("test_file_1",driver="core", 
                               backing_store=False) as f:
    
                c.save(f, test=True)
                
                d = TSaveable()
                
                d.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                
                c.save(f)
                f.seek(0)
                
                c = TSaveable()
                d = c.load(f)            

        self.assertEqual(d.liple[0].a, a.a)
        self.assertEqual(d.liple[0].ts.a, b.a)
        self.assertEqual(d.liple[2].a, c.a)
        
         
        
        lst = dict(a=a,b=b, c=c)        
        c.set_liple(lst)
        
        if legacy:
            with h5py.File("test_file_1",driver="core", 
                               backing_store=False) as f:
    
                c.save(f, test=True)
                
                d = TSaveable()
                
                d.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                
                c.save(f)
                f.seek(0)
                
                d = TSaveable()
                d = d.load(f)            

        self.assertEqual(d.liple["a"].a, a.a)
        self.assertEqual(d.liple["a"].ts.a, b.a)
        self.assertEqual(d.liple["c"].a, c.a)
       
        
        
    def test_saving_and_loading_1(self):
        """Testing saving and loading Saveable objects with lists of Saveables
        
        
        """
        
        obj1 = self.obj1
        
        if legacy:
            with h5py.File("test_file_1",driver="core", 
                               backing_store=False) as f:
                
                t1 = TSaveable()
                t2 = TSaveable()
                t3 = TSaveable()
                t3.set_str("ahoj","cau")
                
                t1.set_saveable(t3)
                
                obj1.set_liple([t1,t2])
                obj1.save(f, test=True)
                
                
                obj2 = TSaveable()
                obj2.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                t1 = TSaveable()
                t2 = TSaveable()
                t3 = TSaveable()
                t3.set_str("ahoj","cau")
                
                t1.set_saveable(t3)
                
                obj1.set_liple([t1,t2])
                obj1.save(f)
                
                f.seek(0)
                
                obj2 = TSaveable()
                obj2 = obj2.load(f)
                
            
            
        self.assertEqual(obj1.a,obj2.a)
        self.assertEqual(obj1.text,obj2.text)
        self.assertEqual(obj1._txt,obj2._txt)
        self.assertEqual(obj1.b1,obj2.b1)
        self.assertEqual(obj1.b2,obj2.b2)
        self.assertEqual(obj1.b3,obj2.b3)
        
        numpy.testing.assert_array_equal(obj1.dat,obj2.dat)
        
        self.assertEqual(obj1.obj.text,obj2.obj.text)
        self.assertEqual(obj1.liple[0].obj.text, t3.text)


#        class TClass:
#            
#            def __init__(self):
#                self.x = "Neco"
#                
#                if legacy:     
#                    with h5py.File("test_file_1",driver="core", 
#                                       backing_store=False) as f:
#                        
#                        t1 = TSaveable()
#                        a = TClass()
#            
#                        t1.set_liple([a, "Hello World!"])
#                        t1.save(f, test=True)
#                        
#                        obj2 = TSaveable()
#                        obj2.load(f, test=True)   
#                else:
#                    pass


    def test_saving_and_loading_2(self):
        """Testing saving and loading Saveable objects with dictionaries of Saveables
        
        
        """
        
        obj1 = self.obj1
        
        if legacy:
            with h5py.File("test_file_1",driver="core", 
                               backing_store=False) as f:
                
                t1 = TSaveable()
                t2 = TSaveable()
                t3 = TSaveable()
                t3.set_str("ahoj","cau")
                
                t1.set_saveable(t3)
                
                obj1.set_liple(dict(jedna=t1,dve=t2))
                obj1.save(f, test=True)
                
                
                obj2 = TSaveable()
                obj2.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                t1 = TSaveable()
                t2 = TSaveable()
                t3 = TSaveable()
                t3.set_str("ahoj","cau")
                
                t1.set_saveable(t3)
                
                obj1.set_liple(dict(jedna=t1,dve=t2))
                obj1.save(f)
                
                f.seek(0)
                
                obj2 = TSaveable()
                obj2 = obj2.load(f)
                
            
            
        self.assertEqual(obj1.a,obj2.a)
        self.assertEqual(obj1.text,obj2.text)
        self.assertEqual(obj1._txt,obj2._txt)
        self.assertEqual(obj1.b1,obj2.b1)
        self.assertEqual(obj1.b2,obj2.b2)
        self.assertEqual(obj1.b3,obj2.b3)
        
        numpy.testing.assert_array_equal(obj1.dat,obj2.dat)
        
        self.assertEqual(obj1.obj.text,obj2.obj.text)
        self.assertEqual(obj1.liple["jedna"].obj.text, t3.text)
            

    def test_saving_and_loading_3(self):
        """Testing saving and loading Saveable objects with lists and tuples
        
        
        """ 
           
        obj3 = self.obj3
        
        if legacy:
            with h5py.File("test_file_2",driver="core", 
                               backing_store=False) as f:
                
                obj3.save(f, test=True)
                
                
                obj2 = TSaveable()
                obj2.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                obj3.save(f)
                
                f.seek(0)
                
                obj2 = TSaveable()
                obj2 = obj2.load(f)
            
        liple = obj3.liple
        n = len(liple)
        for i in range(n):
            self.assertEqual(obj3.liple[i],obj2.liple[i])            
            
        # -------------------------------------------------
        obj3.set_liple(["abc","bca", "cab"])
        
        if legacy:
            
            with h5py.File("test_file_2",driver="core", 
                               backing_store=False) as f:
                
                obj3.save(f, test=True)
                
                
                obj2 = TSaveable()
                obj2.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                obj3.save(f)
                
                f.seek(0)
                
                obj2 = TSaveable()
                obj2 = obj2.load(f)

                
        liple = obj3.liple
        n = len(liple)
        for i in range(n):
            self.assertEqual(obj3.liple[i],obj2.liple[i])            

        # -------------------------------------------------
        obj3.set_liple([1,"bca", "cab"])
        
        if legacy:
            with h5py.File("test_file_2",driver="core", 
                               backing_store=False) as f:
                
                obj3.save(f, test=True)
                
                
                obj2 = TSaveable()
                obj2.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                obj3.save(f)
                
                f.seek(0)
                
                obj2 = TSaveable()
                obj2 = obj2.load(f)
            
        liple = obj3.liple
        n = len(liple)
        for i in range(n):
            self.assertEqual(obj3.liple[i],obj2.liple[i])            

        # -------------------------------------------------
        obj3.set_liple([[1,2,3],"bca", "cab"])
        
        if legacy:
            with h5py.File("test_file_2",driver="core", 
                               backing_store=False) as f:
                
                obj3.save(f, test=True)
                
                
                obj2 = TSaveable()
                obj2.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                obj3.save(f)
                
                f.seek(0)
                obj2 = TSaveable()
                obj2 = obj2.load(f)
            
        liple = obj3.liple
        n = len(liple)
        for i in range(n):
            self.assertEqual(obj3.liple[i],obj2.liple[i])
            
            
            
        # -------------------------------------------------
        obj3.set_liple([[1,2,3],("bca",2,("a",234.0)), "cab"])
        
        if legacy:
            with h5py.File("test_file_2",driver="core", 
                               backing_store=False) as f:
                
                obj3.save(f, test=True)
                
                
                obj2 = TSaveable()
                obj2.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                obj3.save(f)
                f.seek(0)
                obj2 = TSaveable()
                obj2 = obj2.load(f)
            
        liple = obj3.liple
        n = len(liple)
        for i in range(n):
            #print(obj3.liple[i],obj2.liple[i])
            self.assertEqual(obj3.liple[i],obj2.liple[i])

        # -------------------------------------------------
        obj3.set_liple(dict(ab="ahoj", cd=2.0))
        
        if legacy:
            with h5py.File("test_file_2",driver="core", 
                               backing_store=False) as f:
                
                obj3.save(f, test=True)
                
                
                obj2 = TSaveable()
                obj2.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                obj3.save(f)
                
                f.seek(0)
                obj2 = TSaveable()
                obj2 = obj2.load(f)
            
        liple = obj3.liple
        for i in liple.keys():
            #print(obj3.liple[i],obj2.liple[i])
            self.assertEqual(obj3.liple[i],obj2.liple[i])



        # -------------------------------------------------
        obj3.set_liple([[1,2,3],("bca",2,("a",234.0)), 
                        dict(ab="ahoj", cd=2.0)])
        
        if legacy:
            with h5py.File("test_file_2",driver="core", 
                               backing_store=False) as f:
                
                obj3.save(f, test=True)
                
                
                obj2 = TSaveable()
                obj2.load(f, test=True)
        else:
            with tempfile.TemporaryFile() as f:
                obj3.save(f)
                
                f.seek(0)
                obj2 = TSaveable()
                obj2 = obj2.load(f)
            
        liple = obj3.liple
        n = len(liple)
        for i in range(n):
            #print(obj3.liple[i],obj2.liple[i])
            self.assertEqual(obj3.liple[i],obj2.liple[i])


        
        #self.assertEqual(1,2)

    def test_saving_and_loading_with_loader(self):
        """Testing loading Saveable objects with load(filename) function
        
        
        """ 
        if legacy:
            from quantarhei.core.saveable import load, read_info
        else:
             from quantarhei.core.parcel import load_parcel, check_parcel
        from quantarhei import Manager
        
        obj3 = self.obj3
        
        if legacy:
            with h5py.File("test_file_2",driver="core", 
                               backing_store=False) as f:
                
                obj3.save(f, test=True)
                
                obj = load(f, test=True)
                
                info = read_info(f)
        else:
            with tempfile.TemporaryFile() as f:
                obj3.save(f)
                
                f.seek(0)
                obj = load_parcel(f)
                
                f.seek(0)
                info = check_parcel(f)
            
        liple = obj3.liple
        n = len(liple)
        for i in range(n):
            self.assertEqual(obj3.liple[i],obj.liple[i])
        
        if legacy:
            self.assertEqual(info["version"],Manager().version)
        else:
            self.assertEqual(info["qrversion"], Manager().version)
            