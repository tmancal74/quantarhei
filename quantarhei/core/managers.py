# -*- coding: utf-8 -*-
import os
import json
import pkg_resources

import numpy

from .units import conversion_facs_frequency
from .units import conversion_facs_energy
from .units import conversion_facs_length

from .singleton import Singleton


class Manager(metaclass=Singleton):
    """ Main package Manager


    This class handles several important package wide tasks:

    1) Usage of units across objects storing data
    2) Basis conversion of all registered objects
    3) Calls to proper optimized implementations of numerically heavy
       sections of the calculations


    Manager is a singleton class, only one instance exists at all times
    and all managing objects have the instance of the Manager.
    
    Properies
    ---------
    
    version : string
        contains the package version number
 
    
    allower_utypes : list
        contains a list of unit types which can be controlled by the Manager
        
    units : dictionary
        dictionary of available units for each units type
        
    units_repre : dictionary
        dictionary of abreviations used to represent various units
        
    units_repre_latex : dictionary
        dictionary of latex prepresentations of available units
    
    """
    


    #version = "0.0.11"
    version = pkg_resources.require("quantarhei")[0].version

    # hard wired unit options
    allowed_utypes = ["energy",
                      "frequency",
                      "dipolemoment",
                      "temperature",
                      "time",
                      "length"]

    units = {"energy"       : ["1/fs", "int", "1/cm", "eV", "meV", "THz",
                               "J", "SI", "nm", "Ha", "a.u."],
             "frequency"    : ["1/fs", "int", "1/cm", "THz", "Hz", "SI",
                               "nm", "Ha", "a.u."],
             "dipolemoment" : ["Debye", "a.u"],
             "temperature"  : ["1/fs", "int", "Kelvin", "Celsius",
                               "1/cm", "eV", "meV", "Thz", "SI"],
             "time"         : ["fs", "int", "as", "ps", "ns", "Ms","ms",
                               "s", "SI"],
             "length"       : ["int", "A", "nm", "Bohr", "a.u.", "m", "SI"]}

    units_repre = {"Kelvin":"K",
                   "Celsius":"C",
                   "Debye":"D",
                   "1/cm":"1/cm",
                   "THz":"THz",
                   "eV":"eV",
                   "1/fs":"1/fs",
                   "int":"2pi/fs",
                   "meV":"meV",
                   "nm":"nm",
                   "Ha":"Ha",
                   "a.u.":"a.u."}
                   
    units_repre_latex = {"Kelvin":"K",
                   "Celsius":"C",
                   "Debye":"D",
                   "1/cm":"cm$^-1$",
                   "THz":"THz",
                   "eV":"eV",
                   "1/fs":"fs$^{-1}$",
                   "meV":"meV",
                   "nm":"nm",
                   "Ha":"Ha",
                   "a.u.":"a.u."}                  

    def __init__(self):

        self.current_units = {}

        # main configuration file
        cfile = "~/.quantarhei/quantarhei.json"

        # test the presence of configuration directory
        conf_path = os.path.dirname(cfile)
        self.conf_path = os.path.expanduser(conf_path)
        self.cfile = os.path.expanduser(cfile)

        exists = os.path.exists(self.conf_path)
        isdir = os.path.isdir(self.conf_path)
        if not exists:
            # create directory
            os.mkdir(self.conf_path)

            # write default configuration
            self.main_conf = {"units":"units.json",
                              "implementations":"implementations.json"
                             }

            # save it
            with open(self.cfile, 'w') as f:
                json.dump(self.main_conf, f)

        elif exists and (not isdir):
            raise Exception("Cannot create configuration directory.")

        else:
            # load the main configuration file
            with open(self.cfile, 'r') as f:
                self.main_conf = json.load(f)




        #
        #  Setting physical units
        #

        # internal units are hardwired
        self.internal_units = {"energy":"1/fs", "frequency":"1/fs",
                               "dipolemoment":"Debye",
                               "temperature":"Kelvin", "length":"A"}

        # current units are read from conf file
        if not exists:
            # set hard wired defaults and save them
            self.current_units = {"energy":"1/fs", "frequency":"1/fs",
                                  "dipolemoment":"Debye",
                                  "temperature":"Kelvin", "length":"A"}


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
        
        #
        #  All available implementations
        #
        self.all_implementations = {
           "redfieldrates.ssRedfieldRateMatrix": {
               '0':"quantarhei.implementations.python",
               '1':"quantarhei.implementations.cython"
               }
           }
            
        self.all_implementations["redfieldtensor.ssRedfieldTensor"] = \
            {'0':"quantarhei.implementations.python",
             '1':"quantarhei.implementations.cython"}
            
                
        self.default_implementations = {
            "redfieldrates.ssRedfieldRateMatrix":'0',
            "redfieldtensor.ssRedfieldRateTensor":'0'
            }

        self.optimal_implementations = {
            "redfieldrates.ssRedfieldRateMatrix":'1'
            }
                
        self.current_implementations = {
            "redfieldrates.ssRedfieldRateMatrix":'0',
            "redfieldtensor.ssRedfieldRateTensor":'0'        
            }
        
        self.parallel_implementations = {}
        self.parallel_implementations["redfieldrates."
                                 +"ssRedfieldRateMatrix"] = \
            {'0':"quantarhei.implementations.python.parallel",
             '1':"quantarhei.implementations.cython.parallel"}
                
        if not exists:
            
            # and save them
            self.save_implementations()  
            
        else:
            self.load_implementations()
            
            
        self.change_implementation_at_runtime = True
        
        self.basis_stack = []
        self.basis_stack.append(0)
        self.basis_transformations = []
        self.basis_transformations.append(1)
        self.basis_registered = {}
        
        self.warn_about_basis_change = False
        self.warn_about_basis_changing_objects = False
        
        self.parallel_conf = None
        
        self.save_dict = {}
        
            
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
            
    def unit_repr_latex(self,utype="energy",mode="current"):
        """Returns a string representing the currently used units
        
        
        """        
    
    
        if utype in self.allowed_utypes:
            if mode == "current":
                return self.units_repre_latex[self.current_units[utype]]
            elif mode == "internal":
                return self.units_repre_latex[self.internal_units[utype]]
            else:
                raise Exception("Unknown representation mode")
            
        else:
            raise Exception("Unknown unit type")            
            
            
            
    def set_current_units(self, utype, units):
        """Sets current units
        
        
        """
        self._saved_units = {}
        self._saved_units[utype] = self.get_current_units(utype)
        
        if utype in self.allowed_utypes:
            if units in self.units[utype]:
                self.current_units[utype] = units
            else:
                raise Exception("Unknown units of %s" % utype)
        else:
            raise Exception("Unknown type of units")
        
    def unset_current_units(self, utype):
        """Restores previously saved units of a given type
        
        """
        try:
            cunits = self._saved_units[utype]
        except KeyError:
            raise Exception("Units to restore not found")
            
        if utype in self.allowed_utypes:
            if cunits in self.units[utype]:
                self.current_units[utype] = cunits
            else:
                raise Exception("Unknown units of %s" % utype)
        else:
            raise Exception("Unknown type of units")
        
        
        
    def get_current_units(self, utype):
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
            x = conversion_facs_energy[units]
            i_val = x*val
            
            cu = self.current_units["energy"] 
            if cu != "2pi/fs":
                y = conversion_facs_energy[units] 
                return i_val/y
                
            return i_val
            
#    @deprecated       
    def iu_energy(self,val,units="1/cm"):
        """Converst to internal energy units
        
        """
        if units in self.units["energy"]:
            x = conversion_facs_energy[units]
            i_val = x*val
            return i_val    
            
            
    def convert_energy_2_internal_u(self,val):
        """Convert energy from currently used units to internal units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        units = self.current_units["energy"]
        cfact = conversion_facs_energy[self.current_units["energy"]]
        
        # special handling for nano meters
        if units == "nm":
            # zero is interpretted as zero energy
            try:
                ret = numpy.zeros(val.shape, dtype=val.dtype)
                ret[val!=0.0] = 1.0/val[val!=0]
                return ret/cfact
            except:            
                return (1.0/val)/cfact
            #if val == 0.0:
            #    return 0.0
            #return (1.0/val)/cfact
        else:
            return val*cfact
        
            
    def convert_energy_2_current_u(self,val):
        """Converts energy from internal units to currently used units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        units = self.current_units["energy"]
        cfact = conversion_facs_energy[units]
        
        # special handling for nanometers
        if units == "nm":
            # zero is interpretted as zero energy
            try:
                ret = numpy.zeros(val.shape, dtype=val.dtype)
                ret[val!=0.0] = 1.0/val[val!=0]
                return ret/cfact
            except:            
                return (1.0/val)/cfact
        else:
            return val/cfact 
        

    def convert_frequency_2_internal_u(self,val):
        """Converts frequency from currently used units to internal units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        return val*conversion_facs_frequency[self.current_units["frequency"]]

        
    def convert_frequency_2_current_u(self,val):   
        """Converts frequency from internal units to currently used units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        return val/conversion_facs_frequency[self.current_units["frequency"]] 
        

    def convert_length_2_internal_u(self,val):
        """Converts length from currently used units to internal units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        return val*conversion_facs_length[self.current_units["length"]]        


    def convert_length_2_current_u(self,val):   
        """Converts frequency from internal units to currently used units
        
        Parameters
        ==========

        val : number, array, list, tuple of numbers
            values to convert            
        
        """
        return val/conversion_facs_length[self.current_units["length"]]   
      
            
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
        
    def set_current_implementation(self, imp, choice):
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
        
    def transform_to_current_basis(self, operator):
        """Transforms an operator to the currently used basis
        
        Parameters
        ----------
        
        operator : operator
            Any basis managed operator
            
            
        """

        if operator.is_basis_protected:
            return
            
        ob = operator.get_current_basis()
        cb = self.get_current_basis()

        if self.warn_about_basis_changing_objects:
            print("Object ", operator.__class__,
                  id(operator), " is changing basis from ", ob, " to: ", cb)
                
        if ob != cb:
                            
            SS = numpy.diag(numpy.ones(operator.dim))
            # find out if current basis of the object is in the stack (i.e. it 
            # was used sometime in the past)
            if ob in self.basis_stack:
                sl = len(self.basis_stack)
                # scroll back over the bases
                for k in range(1,sl):

                    # take the basis transformation to the earlier used basis
                    ZZ = self.basis_transformations[sl-k]

                    # included it into the transformation matrix
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
        
        
    def get_DistributedConfiguration(self):
        """
        
        """
        from .parallel import DistributedConfiguration
        
        if self.parallel_conf is None:
            self.parallel_conf = DistributedConfiguration()
        return self.parallel_conf







"""

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



            
class Managed:
    """Base class for managed objects 
    
    
    
    """
    
    manager = Manager()


    
class UnitsManaged(Managed):
    """Base class for objects with management of units
    
    
    """    
    
    def convert_energy_2_internal_u(self, val):
        return self.manager.convert_energy_2_internal_u(val)
        
    def convert_energy_2_current_u(self, val):
        return self.manager.convert_energy_2_current_u(val)
 
    def convert_length_2_internal_u(self, val):
        return self.manager.convert_length_2_internal_u(val)
        
    def convert_length_2_current_u(self, val):
        return self.manager.convert_length_2_current_u(val)    
       
    def unit_repr(self,utype="energy"):
        return self.manager.unit_repr(utype)
        
    def unit_repr_latex(self,utype="energy"):
        return self.manager.unit_repr_latex(utype)
        
        
class EnergyUnitsManaged(Managed):
    
    utype = "energy"
    units = "2pi/fs"    
    
    def convert_2_internal_u(self,val):
        return self.manager.convert_energy_2_internal_u(val)
        
    def convert_2_current_u(self,val):
        return self.manager.convert_energy_2_current_u(val)

    def unit_repr(self):
        return self.manager.unit_repr("energy")

    def unit_repr_latex(self, utype="energy"):
        return self.manager.unit_repr_latex(utype)
    

class LengthUnitsManaged(Managed):
    """Class providing functions for length units conversion
    
    """
    
    utype = "length"
    units = "A"
    
    def convert_2_internal_u(self, val):
        return self.manager.convert_length_2_internal_u(val)
 
    def convert_2_current_u(self,val):
        return self.manager.convert_length_2_current_u(val)
    
    def unit_repr(self):
        return self.manager.unit_repr(self.utype)

    def unit_repr_latex(self):
        return self.manager.unit_repr_latex(self.utype)
    
        
class BasisManaged(Managed):
    """Base class for objects with managed basis


    """
    _current_basis = Manager().get_current_basis()
    is_basis_protected = False
    
    def get_current_basis(self):
        return self._current_basis
        
    def set_current_basis(self,bb):
        self._current_basis = bb
            
    def protect_basis(self):
        self.is_basis_protected = True
        
    def unprotect_basis(self):
        self.is_basis_protected = False
        
        



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
        
        
class frequency_units(energy_units):
    """Context manager for units of frequency
    
    It behaves exactly the same as ``energy_units`` context manager.
    
    """
    pass
        

class length_units(units_context_manager):
    """Context manager for length units
    
    
    """
    
    def __init__(self, units):
        super().__init__(utype="length")
        
        if units in self.manager.units["length"]:
            self.units = units
        else:
            raise Exception("Unknown length units")
            
    def __enter__(self):
        # save current energy units
        self.units_backup = self.manager.get_current_units("length")
        self.manager.set_current_units(self.utype,self.units)
        
    def __exit__(self,ext_ty,exc_val,tb):
        self.manager.set_current_units("length",self.units_backup)


        
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
    
    def __init__(self, operator):
        super().__init__()
        self.op = operator
        
        
    def __enter__(self):

        if self.manager.warn_about_basis_change:
            print("\nQr >>> Entering basis context manager ...")
            
        cb = self.manager.get_current_basis()
        ob = self.op.get_current_basis()
        
        if cb != ob:
            
            self.manager.transform_to_current_basis(self.op)
            
        
        #SS = self.op.diagonalize()
        SS = self.op.get_diagonalization_matrix()
        nb = self.manager.set_new_basis(SS)

        #self.manager.register_with_basis(nb,self.op)
        #self.op.set_current_basis(nb)

        if self.manager.warn_about_basis_change:
            print("\nQr >>>  ... setting context done")        

    
        
    def __exit__(self,ext_ty,exc_val,tb):
  
        if self.manager.warn_about_basis_change:
            print("\nQr >>> Returning from basis context manager. Cleaning ...")  
            
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
            # the operator might have been set to protected mode
            # inside the context
            if not op.is_basis_protected:
                op.transform(S1,inv=SS) 
            op.set_current_basis(nb)
            
            # operators which appeared in this context and where not
            # register in the one above are now registerd
            if nb != 0:
                if op not in ops_above:
                    self.manager.register_with_basis(nb,op)
            
            
        del self.manager.basis_registered[bb]

        if self.manager.warn_about_basis_change:
            print("\nQr >>> ... cleaning done")        
            

def set_current_units(units=None):
    """Sets units globaly without the need for a context manager
    
    """
    manager = Manager()    
    if units is not None:
        # set units using a supplied dictionary
        for utype in units:
            if utype in manager.allowed_utypes:
                un = units[utype]
                # handle the identity of "frequency" and "energy"
                if utype=="frequency":
                    utype="energy"
                    un = units["frequency"]
                    
                manager.set_current_units(utype,un)
            else:
                raise Exception("Unknown units type %s" % utype)

    else:
        # reset units to the default
        for utype in manager.internal_units:
            if utype in manager.allowed_utypes:
                manager.set_current_units(utype,manager.internal_units[utype])
            else:
                raise Exception("Unknown units type %s" % utype)
        
        
        