# -*- coding: utf-8 -*-

from .units import conversion_facs_frequency 
from .singleton import Singleton

#from .wrappers import deprecated

import os
import json
import numpy

class Manager(metaclass=Singleton):
    """ Main package Manager  
    
    
    This class handles several important package wide tasks:
    
    1) Usage of units across objects storing data
    2) Basis conversion of all registered objects
    3) Calls to proper optimized implementations of numerically heavy
       sections of the calculations
       
       
    Manager is a singleton class, only one instance exists at all times
    and all managing objects have the instance of the Manager. 
    
    Units Management
    ----------------
    Units management is performed for all classes derived from 
    quantarhei.managers.UnitsManaged class.
    
    
    Basis Conversion Management
    ---------------------------
    Units management is performed for all classes derived from 
    quantarhei.managers.BasisManaged class. 

    Basis management works like this: when an class is defined, and its
    property needs to be basis managed, one should use a predefined type
    `basis_managed_array_property`    
    
    
    
    """
    
    # hard wired unit options
    allowed_utypes = ["energy",
                      "frequency",
                      "dipolemoment",
                      "temperature",
                      "time"]
                      
    units = {"energy"       : ["1/cm","2pi/fs","eV","meV","THz"],
             "frequency"    : ["1/cm","2pi/fs","THz"],
             "dipolemoment" : ["Debye"],
             "temperature"  : ["Kelvin","Celsius","1/cm"],
             "time"         : ["as","fs","ps","ns"]}
             
    units_repre = {"Kelvin":"K",
                   "Celsius":"C",
                   "Debye":"D",
                   "1/cm":"1/cm",
                   "THz":"THz",
                   "2pi/fs":"2pi/fs",
                   "meV":"meV"}
    
    def __init__(self):
        
        self.current_units = {}
        
        # main configuration file
        cfile = "~/.quantarhei/quantarhei.json"        
        
        # test the presence of configuration directory
        conf_path = os.path.dirname(cfile)
        self.conf_path = os.path.expanduser(conf_path)
        self.cfile = os.path.expanduser(cfile)
        
        exists = os.path.exists(self.conf_path)
        isdir  = os.path.isdir(self.conf_path)
        if not exists:
            # create directory
            os.mkdir(self.conf_path)
            
            # write default configuration
            self.main_conf = {"units":"units.json",
                              "implementations":"implementations.json"
                              }
            
            # save it
            with open(self.cfile,'w') as f:
                json.dump(self.main_conf,f)
            
        elif (exists and (not isdir)):
            raise Exception("Cannot create configuration directory.")
            
        else:
            # load the main configuration file
            with open(self.cfile,'r') as f:
                self.main_conf = json.load(f)
                
                
        
            
        #
        #  Setting physical units
        #
            
        # internal units are hardwired
        self.internal_units = {"energy":"2pi/fs", "frequency":"2pi/fs",
                              "dipolemoment":"Debye",
                              "temperature":"Kelvin"}
        
        # current units are read from conf file
        if not exists:
            # set hard wired defaults and save them
            self.current_units = {"energy":"2pi/fs", "frequency":"2pi/fs",
                              "dipolemoment":"Debye",
                              "temperature":"Kelvin"}
                              

            # save them
            self.save_units()
            
        else:
            self.load_units()
                
                
        #
        #  Setting implementations
        #
                
        self.implementation_points = {
        
           "secular-standard-Redfield-rates":"redfield.ssRedfieldRateMatrix"
           
           }
        self.all_implementations = {
            
           "redfield.ssRedfieldRateMatrix": {
               '0':"quantarhei.implementations.python",
               '1':"quantarhei.implementations.cython"
               }
           
           }
                
        self.default_implementations = {
            "redfield.ssRedfieldRateMatrix":'0'
            }

        self.optimal_implementations = {
            "redfield.ssRedfieldRateMatrix":'1'
            }
                
        self.current_implementations = {
            "redfield.ssRedfieldRateMatrix":'0'        
            }
                
        if not exists:
            
            # and save them
            self.save_implementations()  
            
        else:
            self.load_implementations()
            
            
        self.change_implementation_at_runtime = False
        
        self.basis_stack = []
        self.basis_stack.append(0)
        self.basis_transformations = []
        self.basis_transformations.append(1)
        self.basis_registered = {}
        
            
    def save_settings(self):

        # main configuration file
        with open(self.cfile,'w') as f:
            json.dump(self.main_conf,f)
                
        # units setting
        self.save_units()
        # implementations setting
        self.save_implementations()
        
        
        
    def save_implementations(self):
            # set the implementations to standard
        implementations = {
            "imp_points":self.implementation_points,
            "all_available":self.all_implementations,
            "default":self.default_implementations,
            "optimal":self.optimal_implementations,
            "current":self.current_implementations
            }
        imp_file = self.main_conf["implementations"]
        imp_file = os.path.join(self.conf_path,imp_file)
        with open(imp_file,'w') as f:
            json.dump(implementations,f)
            
    def load_implementations(self):
        imp_file = self.main_conf["implementations"]
        imp_file = os.path.join(self.conf_path,imp_file)
        with open(imp_file,'r') as f:
            implementations = json.load(f)
            self.implementation_points = implementations["imp_points"]
            self.all_implementations = implementations["all_available"]
            self.default_implementations = implementations["default"]
            self.optimal_implementations = implementations["optimal"]
            self.current_implementations = implementations["current"]
            
            
        
        
    def save_units(self):
        units_file = self.main_conf["units"]
        units_file = os.path.join(self.conf_path,units_file)
        with open(units_file,'w') as f:
            json.dump(self.current_units,f)       
            
    def load_units(self):
        units_file = self.main_conf["units"]
        units_file = os.path.join(self.conf_path,units_file)
        with open(units_file,'r') as f:
            self.current_units = json.load(f)       
    

        
    def unit_repr(self,utype="energy",mode="current"):
        """Returns a string representing the currently used units
        
        
        """        
    
    
        if utype in self.allowed_utypes:
            if mode == "current":
                return self.units_repre[self.current_units[utype]]
            elif mode == "internal":
                return self.units_repre[self.internal_units[utype]]
            else:
                raise Exception("Unknown representation mode")
            
        else:
            raise Exception("Unknown unit type")
            
    def set_current_units(self,utype,units):
        """Sets current units
        
        
        """
        if utype in self.allowed_utypes:
            if units in self.units[utype]:
                self.current_units[utype] = units
            else:
                raise Exception("Unknown units of %s" % utype)
        else:
            raise Exception("Unknown type of units")
            
    def get_current_units(self,utype):
        """
        
        """
        if utype in self.allowed_utypes:
            return self.current_units[utype]        
        else:
            raise Exception("Unknown type of units")
            
           
#    @deprecated
    def cu_energy(self,val,units="1/cm"):
        """Converst to current energy units
        
        """
        if units in self.units["energy"]:
            x = conversion_facs_frequency[units]
            i_val = x*val
            
            cu = self.current_units["energy"] 
            if cu != "2pi/fs":
                y = conversion_facs_frequency[units] 
                return i_val/y
                
            return i_val
            
#    @deprecated       
    def iu_energy(self,val,units="1/cm"):
        """Converst to internal energy units
        
        """
        if units in self.units["energy"]:
            x = conversion_facs_frequency[units]
            i_val = x*val
            return i_val    
            
    def convert_energy_2_internal_u(self,val):
        """Convert energy from currently used units to internal units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        return val*conversion_facs_frequency[self.current_units["energy"]]
        
            
    def convert_energy_2_current_u(self,val):
        """Convert energy from internal units to currently used units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        return val/conversion_facs_frequency[self.current_units["energy"]] 
        

    def convert_frequency_2_internal_u(self,val):
        """Convert frequency from currently used units to internal units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        return val*conversion_facs_frequency[self.current_units["frequency"]]

        
    def convert_frequency_2_current_u(self,val):   
        """Convert frequency from internal units to currently used units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        return val/conversion_facs_frequency[self.current_units["frequency"]] 
        
        
        
            
    def get_implementation_prefix(self,package="",taskname=""):
        #default_imp_prefix = "quantarhei.implementations.python"
        
        pname = package+"."+taskname
        whichone = self.current_implementations[pname]
        imp_prefix = self.all_implementations[pname][str(whichone)]
        
        return imp_prefix


    def get_implementation_points(self):
        return self.implementation_points
        
        
    def get_all_implementations(self):
        return self.all_implementations
        
    def get_all_implementations_of(self,imp):
        imp_id = self.implementation_points[imp]
        return self.all_implementations[imp_id]
                
    def get_current_implementation(self,imp):
        imp_id = self.implementation_points[imp]
        whichone = self.current_implementations[imp_id]
        return self.all_implementations[imp_id][str(whichone)]
        
    def set_current_implementation(self,imp,choice):
        imp_id = self.implementation_points[imp]
        self.current_implementations[imp_id] = choice
        
    def register_implementation(self,imp_point,prefix,asint=None):
        pass
    
    def commit_implementation(self,imp_point,prefix,asint=None):
        pass
    
    
    def get_current_basis(self):
        """Returns the current basis id
        
        """
        l = len(self.basis_stack)
        return self.basis_stack[l-1]
        
    def set_new_basis(self,SS):
        nb = self.get_current_basis() + 1
        self.basis_stack.append(nb)
        self.basis_transformations.append(SS)
        self.basis_registered[nb] = []
        return nb
        
    def transform_to_current_basis(self,operator):
        
        ob = operator.get_current_basis()
        cb = self.get_current_basis()
                
        if ob != cb:
                            
            SS = numpy.diag(numpy.ones(operator._data.shape[0]))
            # find out if current basis of the object is in the stack (i.e. it 
            # was used sometime in the past)
            if ob in self.basis_stack:
                sl = len(self.basis_stack)
                # scroll back over the bases
                for k in range(1,sl):

                    # take the basis transformation to the earlier used basis
                    ZZ = self.basis_transformations[sl-k]

                    # included into the transformation matrix
                    SS = numpy.dot(ZZ,SS)                
                    # if the basis is found, break away from the loop
                    if self.basis_stack[sl-k-1] == ob:
                        break
            else:
                raise Exception("Basis of the object is not on stack.")
            
            operator.transform(SS)
            operator.set_current_basis(cb)
            self.register_with_basis(cb,operator)
        

    def register_with_basis(self,nb,operator):
        self.basis_registered[nb].append(operator)
        
                

class Managed:
    """Base class for managed objects 
    
    
    
    """
    
    manager = Manager()


    
class UnitsManaged(Managed):
    """Base class for objects with management of units
    
    
    """
    
    def convert_energy_2_internal_u(self,val):
        return self.manager.convert_energy_2_internal_u(val)
        
    def convert_energy_2_current_u(self,val):
        return self.manager.convert_energy_2_current_u(val)
        
    def unit_repr(self,utype="energy"):
        return self.manager.unit_repr(utype)
        
        
class BasisManaged(Managed):
    """Base class for objects with managed basis


    """
    _current_basis = Manager().get_current_basis()
    
    
    def get_current_basis(self):
        return self._current_basis
        
    def set_current_basis(self,bb):
        self._current_basis = bb
            
        
        

class units_context_manager:
    """General context manager to manage physical units of values 
    
    
    """
    
    def __init__(self,utype="energy"):
        self.manager = Manager()
        if utype in self.manager.allowed_utypes:
            self.utype = utype
        else:
            raise Exception("Unknown units type")
    
    def __enter__(self):
        pass
    
    def __exit__(self):
        pass


class energy_units(units_context_manager):
    """Context manager for units of energy
    
    
    """
    
    def __init__(self,units):
        super().__init__(utype="energy")
        
        if units in self.manager.units["energy"]:
            self.units = units
        else:
            raise Exception("Unknown energy units")
            
    def __enter__(self):
        # save current energy units
        self.units_backup = self.manager.get_current_units("energy")
        self.manager.set_current_units(self.utype,self.units)
        
    def __exit__(self,ext_ty,exc_val,tb):
        self.manager.set_current_units("energy",self.units_backup)
        
        
class basis_context_manager:
    """General context manager to manage basis 
    
    
    """    
    def __init__(self):
        self.manager = Manager()

    
    def __enter__(self):
        pass
    
    def __exit__(self,ext_ty,exc_val,tb):
        pass


        
class eigenbasis_of(basis_context_manager):
    """Context manager for basis change
    
    
    """
    
    def __init__(self,operator):
        super().__init__()
        self.op = operator
        
        
    def __enter__(self):

        cb = self.manager.get_current_basis()
        ob = self.op.get_current_basis()
        if cb != ob:
            
            self.manager.transform_to_current_basis(self.op)
            
        
        SS = self.op.diagonalize()
        nb = self.manager.set_new_basis(SS)

        self.manager.register_with_basis(nb,self.op)
        self.op.set_current_basis(nb)
        
        
    def __exit__(self,ext_ty,exc_val,tb):

        
        # This is the basis we are leaving
        bb = self.manager.basis_stack.pop()
        # this is the transformation we got here with
        SS = self.manager.basis_transformations.pop()
        # This is the new basis
        bss = len(self.manager.basis_stack)
       
        nb = self.manager.basis_stack[bss-1]
        
        # inverse of the transformation matrix
        S1 = numpy.linalg.inv(SS)     
        
        # transform all registered objects
        operators = self.manager.basis_registered[bb]
        
        if nb != 0:
            # operators registered with the context above this one
            ops_above = self.manager.basis_registered[nb]

        for op in operators:
            op.transform(S1,inv=SS) #op._data = numpy.dot(SS,numpy.dot(op._data,S1))
            op.set_current_basis(nb)
            
            # operators which appeared in this context and where not
            # register in the one above are now registerd
            if nb != 0:
                if op not in ops_above:
                    self.manager.register_with_basis(nb,op)
            
            
        del self.manager.basis_registered[bb]
        
        
        
        
            
        
        