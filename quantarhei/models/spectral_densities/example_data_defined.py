# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot as plt

from quantarhei.models.spectdens import DatabaseEntry 
from quantarhei.models.spectdens import DataDefinedEntry

from quantarhei import SpectralDensity
from quantarhei import energy_units

class example_data_defined_array(DataDefinedEntry):
    
    direct_implementation = DatabaseEntry.SPECTRAL_DENSITY
        
    def get_data(self):
        
        data = numpy.array([
                
                  0.0,  0.0,
                  1.0,  0.1,
                  2.0,  0.2,
                  3.0,  0.25,
                  4.0,  0.27,
                  5.0,  0.28,
                  6.0,  0.28,
                  7.0,  0.27,
                  8.0,  0.25,
                  9.0,  0.23,
                  10.0, 0.20
                  
                  ])
        
        return data
    
class example_data_defined_string(DataDefinedEntry):
    
    direct_implementation = DatabaseEntry.SPECTRAL_DENSITY
        
    def get_data_string(self):
        
        return """
                
                  0.0,  0.0
                  1.0,  0.1
                  2.0,  0.2
                  3.0,  0.25
                  4.0,  0.27
                  5.0,  0.28
                  6.0,  0.28
                  7.0,  0.27
                  8.0,  0.25
                  9.0,  0.23
                  10.0, 0.20
                  
                """

class example_data_defined_comment(DataDefinedEntry):
    """
                
                  0.0,  0.0
                  1.0,  0.1
                  2.0,  0.2
                  3.0,  0.25
                  4.0,  0.27
                  5.0,  0.28
                  6.0,  0.28
                  7.0,  0.27
                  8.0,  0.25
                  9.0,  0.23
                  10.0, 0.20
                  
    """
    direct_implementation = DatabaseEntry.SPECTRAL_DENSITY
    
    
    
    
       

if __name__ == "__main__":
    
    ex1 = example_data_defined_array()
    
    plt.plot(ex1._rawdata[:,0], ex1._rawdata[:,1])
    plt.show()    
    
    ex2 = example_data_defined_string()
    print(ex2.get_data_string())
    