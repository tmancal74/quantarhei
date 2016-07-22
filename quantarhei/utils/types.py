# -*- coding: utf-8 -*-

from functools import partial
import numpy
import numbers

def array_property(name,shape=None):
    """ Controls the access to an array property of package classes 
    
    
    """
    storage_name = '_'+name
    
    @property
    def prop(self):
        return getattr(self,storage_name)
    
    @prop.setter
    def prop(self,value):
        try:
            vl = check_numpy_array(value)
            if not (shape == None):
                if not (shape == vl.shape):
                    raise TypeError(
                    '{} must be of shape {}'.format(name,shape))  
            setattr(self,storage_name,vl)
        except:
            raise TypeError(
            '{} must be either a list or numpy.array'.format(name))
        
    return prop
    
    
def basis_managed_array_property(name,dtype,shape=None):

    storage_name = '_'+name

    @property
    def prop(self):

        # check if basis id corresponds to the one known by Manager
        cb = self.manager.get_current_basis()
        # get object's current basis
        ob = self.get_current_basis()

                    
        if cb == ob:
            pass 
        else:
            # change basis
            self.manager.transform_to_current_basis(self)
        
        return getattr(self,storage_name)        
        
    
    @prop.setter
    def prop(self,value):
        
        # check if basis id corresponds to the one known by Manager
        cb = self.manager.get_current_basis()
        # get object's current basis
        ob = self.get_current_basis()
            
        if cb == ob:
            pass
        else:
            # change basis
            self.manager.transform_to_current_basis(self)

        try:
            vl = check_numpy_array(value)
            if not (shape == None):
                if not (shape == vl.shape):
                    raise TypeError(
                    '{} must be of shape {}'.format(name,shape))  
            setattr(self,storage_name,vl)
        except:
            raise TypeError(
            '{} must be either a list or numpy.array'.format(name))
        
    return prop        
        


def check_numpy_array(val):
    """ Checks if argument is a numpy array. 
    
    If the argument is a list, it converts it into numpy array. Otherwise,
    error occurs.
    
    """
    if isinstance(val,numpy.ndarray):
        return val
    elif isinstance(val,list):
        try:
            vl = numpy.array(val)
        except:
            raise TypeError('Numerical array is required')
        return numpy.array(vl)
    else:
        raise TypeError('List or numpy.ndarray required')
    
def typed_property(name,dtype):
    storage_name = '_'+name
            
    @property
    def prop(self):
        return getattr(self,storage_name)
    
    @prop.setter
    def prop(self,value):
        if isinstance(value,dtype):
            setattr(self,storage_name,value)
        else:
            raise TypeError('{} must be of type {}'.format(name,dtype))
            
    return prop
    
    
 
def list_of_typed_property(name,dtype):
    storage_name = '_'+name
    
    @property
    def prop(self):
        return getattr(self,storage_name)
        
    @prop.setter
    def prop(self,value):
        if isinstance(value,list):
            for el in value:
                if not isinstance(el,dtype):
                    raise TypeError('{} must contain \
                    values of type {})'.format(name,dtype),dtype)
            
    return prop
    



   
def alt_type_property(name,dtype):
    storage_name = '_'+name
    
    @property
    def prop(self):
        return getattr(self,storage_name)
    
    @prop.setter
    def prop(self,value):
        if isinstance(dtype,list):
            isset = False
            for tps in dtype:
                if isinstance(value,tps):
                    setattr(self,storage_name,value)
                    isset = True
                    break
            if not isset:
                raise TypeError('{} must be of type {}'.format(name,dtype))
        elif isinstance(value,dtype):
            setattr(self,storage_name,value)
        else:
            raise TypeError('{} must be of type {}'.format(name,dtype))
            
    return prop    
    
    

def derived_type(name,dtype_in):
    if isinstance(dtype_in,list):
        tp = partial(alt_type_property,dtype=dtype_in)
    else:
        tp = partial(typed_property,dtype=dtype_in)
    return tp(name)                    
            

      
Float   = partial(typed_property,dtype=numbers.Real)
Integer = partial(typed_property,dtype=numbers.Integral)
Bool    = partial(typed_property,dtype=bool)

BasisManagedComplex = partial(basis_managed_array_property,
                              dtype=numbers.Complex)
BasisManagedReal = partial(basis_managed_array_property,
                              dtype=numbers.Real)

