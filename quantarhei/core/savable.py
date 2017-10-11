# -*- coding: utf-8 -*-
import numpy
import h5py

class Savable:
    """Class implementing object persistence through saving to hdf5 files
    
    
    
    
    """
    
    def save(self, filename):
        """This method should be implemented in classes that inherit from here
        
        """
        raise Exception("'save' method is not implemented")
    
    def save_in(self, loc, name):
        pass
    
    def load(self, filename):
        """This method should be implemented in classes that inherit from here
        
        """        
        raise Exception("'save' method is not implemented")
        
    def load_from(self, loc, name):
        pass
    
    
    #
    # saving methods
    #
    def _save_string_attributes(self, loc, dictionary):
        """
        
        """
        pass
    
    
    def _save_numeric_attributes(self, loc, dictionary):
        """
        
        """        
        pass
    
    def _save_bool_attributes(self, loc, dictionary):
        pass
    
    def _save_array(self, loc, name, data):
        pass
    
    #
    # loading methods
    #   
    def _load_string_attributes(self, loc, attrlist):
        pass
    
    def _load_numeric_attributes(self, loc, attrlist):
        pass
    
    def _load_bool_attributes(self, loc, attrlist):
        pass
    
    def _load_array(self, loc, name):
        pass
    
     
        
    
class TestObject(Savable):
    """
    
    

    
    """
    def __init__(self):
        
        self.name = "Ahoj"
        self.class_name = "TestObject"
        self.N = 23
        
        
    def save(self, filename):
        
        # strings
        strs = dict(name=self.name,
                    class_name=self.class_name) 
        # integers
        ints = dict(N=self.N)
        
        with h5py.File(filename,"w") as file:

            #
            # Saving all types
            #
            self.save_string_attributes(file, strs)
            self.save_numeric_attributes(file, ints)
        
        
        
    def load(self, file):
        
        string_attr_names = ["name", "class_name"]
        string_attrs = self.load_string_attributes(file, string_attr_names)
        self.name = string_attrs["name"]
        self.class_name = string_attrs["class_name"]
        
        numeric_attr_names = ["N"]
        numeric_attrs = self.load_numeric_attributes(file, numeric_attr_names)
        self.N = numeric_attrs["N"]
        
                