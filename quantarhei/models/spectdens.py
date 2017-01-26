# -*- coding: utf-8 -*-

from ..qm.corfunctions import SpectralDensity, CorrelationFunction
from ..core.managers import energy_units

class DatabaseEntry:
    """Database entry in the database of spectral densities and correlation functions
    
    """
    SPECTRAL_DENSITY     = 0
    CORRELATION_FUNCTION = 1
    
    def __init__(self):
        
        self.identificator = None
        self.att_ident = []
        self.direct_implementation = None
            
        
    def get_CorrelationFunction(self, temperature):
        pass
    
    def get_SpectralDensity(self):
        pass
    
    def direct_implementation(self):
        return self.direct_implementation
    
    
class SpectralDensityDatabaseEntry(DatabaseEntry):
    """Database entry for spectral density
    
    This class defines a method to obtain correlation function from spectral
    density. Spectral density calculation has to be implemented.
    
    """

    direct_implementation = DatabaseEntry.SPECTRAL_DENSITY
    
    def get_CorrelationFunction(self, temperature):
        """Returns CorrelationFunction based on SpectralDensity
        
        """
        return self.get_SpectralDensity().get_CorrelationFunction(temperature)
    

class CorrelationFunctionDatabaseEntry(DatabaseEntry):
    """Database entry for correlation function
    
    This class defines a method to obtain spectral density from correlation
    function. Spectral density calculation has to be implemented.
    
    """
    direct_implementation = DatabaseEntry.CORRELATION_FUNCTION
    
    def get_SpectralDensity(self):
        """Returns SpectralDensity on CorrelationFunction based
        
        """
        temperature = self.temperature
        return self.get_CorrelationFunction(temperature).get_SpectralDensity()    


class SpectralDensityDB:
    """Database of spectral densities and correlation functions
    
    This is only a temporary solution which holds all correlation
    functions in memory once it is instantiates. In future, only a list of
    identificators should be held and instantiation of database entries 
    will occure only when required.
    
    """
    
    def __init__(self):
        
        from . import spectral_densities as SD
        import pkgutil
        import importlib 
        import inspect
        
        self.entries = {}
        
        count = 0
        
        # walk trough the spectral_densities directory, importing all modules
        package = SD
        for importer, modname, ispkg in pkgutil.walk_packages(
                    path=package.__path__, prefix=package.__name__+".",
                    onerror=lambda x: None):

            # import a found module
            wen = importlib.import_module(modname)

            # find all present classes
            for name, obj in inspect.getmembers(wen):
                if inspect.isclass(obj):
                    
                    # get the class
                    cls = getattr(wen, name)
                
                    base_classes = cls.__bases__

                    #
                    # The classes below have to be subclassed
                    #
                    if ((DatabaseEntry in base_classes) 
                        or (SpectralDensityDatabaseEntry in base_classes)
                        or (CorrelationFunctionDatabaseEntry in base_classes)):
                        
                        # instantiate the class
                        inst = cls()
                    
                        if isinstance(inst, DatabaseEntry):    
                            # print identificator
                            if inst.identificator is not None:
                                print("Registering: ", inst.identificator)
                                count += 1
                                self.entries[inst.identificator] = inst
                                            
                        else:
                            raise Exception("This cannot be!")

        if count != len(self.entries):
            raise Exception()
            
        
    def get_SpectralDensity(self, axis, ident=None):
        """Returns spectral density based on an identificator

        """
       
        if ident in self.entries.keys():
            inst = self.entries[ident]
        else:
            raise Exception("Identificator: "+ident+" not found")
        
        try:
            sd = inst.get_SpectralDensity(axis)
        except:
            raise Exception("Could not obtain spectral density from "+
                            "definition class for : "+ident)
            
        return sd



class CorrelationFunctionDB(SpectralDensityDB):
    """This class is identical to SectralDensityDB
    
    """
    pass
    
