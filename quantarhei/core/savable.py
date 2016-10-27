# -*- coding: utf-8 -*-
import numpy

class SavableObject:
    """
    
    
    
    
    """
    
    def save_object(self):
        """This method should be implemented in classes that inherit from here
        
        """
        raise Exception("'save_object' method is not implemented")
    
    
    def load_object(self):
        """This method should be implemented in classes that inherit from here
        
        """        
        raise Exception("'save_object' method is not implemented")
        
    
    def save_string_properties(self, file, properties):
        """Creates a numpy array of strings and saves it
        
        """
        
        for key in properties.keys():
            print(key, properties[key])

    def save_integer_properties(self, file, properties):
        """Creates a numpy array of integers and saves it
        
        """        
        for key in properties.keys():
            print(key, properties[key])          
        
        
class TestObject(SavableObject):
    """
    
    

    
    """
    def __init__(self):
        
        self.name = "Ahoj"
        self.class_name = "TestObject"
        self.N = 23
        
        
    def save_object(self, file):
        
        # strings
        strs = dict(name=self.name,
                    class_name=self.class_name) 
        # integers
        ints = dict(N=self.N)
        
        #
        # Saving all types
        #
        self.save_string_properties(file, strs)
        self.save_integer_properties(file, ints)
        
        
    def load_object(self, file):
        pass
    
        #
        # Load all types
        #
    
        # strings
    
        # integers
    
        # ...
        