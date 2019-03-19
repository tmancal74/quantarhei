# -*- coding: utf-8 -*-

class DatabaseEntry:
    """Database entry in the database of spectral densities and correlation functions
    
    """
    SPECTRAL_DENSITY     = 0
    CORRELATION_FUNCTION = 1
    direct_implementation = None    
    identificator = None   
    alt_ident = []
    units = "int"
    
    
    
    def __init__(self):
        pass

    def get_CorrelationFunction(self, temperature):
        pass
    
    def get_SpectralDensity(self):
        pass

    
    
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


class DataDefinedEntry(DatabaseEntry):
    
    direct_implementation = DatabaseEntry.SPECTRAL_DENSITY
    
    def __init__(self):
        
        self._have_data = False
        
        #
        # try to get data defined in get_data method
        #
        data = self.get_data()
        
        if data is not None:
            self._rawdata = self.get_data()
            length = len(self._rawdata)
            self._rawdata.shape = (length//2, 2)
            self._have_data = True
        else:
            #
            # try to get data from a string defined in get_data_string method
            #
            dstr = self.get_data_string()
            if dstr is not None:
                pass
            else:
                # 
                # try to get data from the comment
                #
                pass
            
        if self._have_data:

            if self.direct_implementation == DatabaseEntry.SPECTRAL_DENSITY:
                pass
            elif (self.direct_implementation 
                  == DatabaseEntry.CORRELATION_FUNCTION):
                pass
            else:
                raise Exception("Don't know if this is spectral density or "+
                                "correlation function")
                
        else:
            
            raise Exception("Data not defined in any expected way")
            
          
    def get_data(self):
        """Returns spectral density data
        
        Should return a numpy array with two columns of data to construct the
        spectral density
        
        """
        return None
    
    def get_data_string(self):
        """Returns spectral density data as a string 
        
        Should return a string containing two columns of data to construct the
        spectral density
        
        """
        return None
    
            
class SpectralDensityDB:
    """Database of spectral densities and correlation functions
    
    This is only a temporary solution which holds all correlation
    functions in memory once it is instantiates. In future, only a list of
    identificators should be held and instantiation of database entries 
    will occure only when required.
    
    """
    
    def __init__(self, verbose=False):
        
        from . import spectral_densities as SD
        import pkgutil
        import importlib 
        import inspect
        
        if verbose:
            print("Initializing Spectral Density Database")
            
        self.entries = {}
        self.loaded = {}
        self.entry_count = 0
        
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
                
                if (inspect.isclass(obj) and
                    (not ((name == "DataDefinedEntry")
                    or (name == "DatabaseEntry")))):
                    
                    # get the class
                    cls = getattr(wen, name)
                
                    base_classes = cls.__bases__

                    #
                    # Take only subclasses of the following four
                    #
                    allowed_base_classes = [DatabaseEntry,
                                    SpectralDensityDatabaseEntry,
                                    CorrelationFunctionDatabaseEntry,
                                    DataDefinedEntry]
                    
                    can_instantiate = False
                    for bcl in allowed_base_classes:
                        can_instantiate = (can_instantiate or 
                                           (bcl in base_classes))
                        
                    if can_instantiate:
                        
                        # instantiate the class
                        inst = cls()
                    
                        if isinstance(inst, DatabaseEntry):    
                            
                            if inst.identificator is not None:
                                # print identificator
                                if verbose:
                                    print("... registering: ",
                                          inst.identificator)
                                count += 1
                                
                                #
                                # WE DO NOT SAVE THE OBJECT ON THIS FIRST LOAD
                                #
                                #self.entries[inst.identificator] = inst
                                            
                                self.entries[inst.identificator] = (modname,
                                                                    name)
                                self.loaded[inst.identificator] = False
                                
                        else:
                            raise Exception("This cannot be!")

        if count != len(self.entries):
            raise Exception()
        self.entry_count = count
        
        if verbose:
            print(self.get_status_string())
            print("Database initialized and ready!")


    def get_status_string(self):
        """Returns a status report of the database
        
        """

        # count loaded modules and print a report
        kount = 0
        for ld in self.loaded.keys():
            if self.loaded[ld]:
                kount += 1

        out = "Spectral Density Database Status"
        out += "\n - "+str(self.entry_count)+" spectral densities registered"
        out += "\n - "+str(kount)+" loaded"
            
        return out

        
    def get_SpectralDensity(self, axis, ident=None):
        """Returns spectral density based on an identificator

        """
        import importlib
        
        if ident in self.entries.keys():
            
            # when class is not loaded and instantiated, create its instance
            if not self.loaded[ident]:
                imodname, clname = self.entries[ident]
                imod = importlib.import_module(imodname)
                klas = getattr(imod, clname)
                inst = klas()
                self.entries[ident] = inst
                self.loaded[ident] = True
            # just look it up if it was already loaded
            else:
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
    
