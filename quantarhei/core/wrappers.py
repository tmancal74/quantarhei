# -*- coding: utf-8 -*-
from functools import wraps
#import os

from .managers import Manager
   
   
def deprecated(func):
    @wraps(func)
    def wrapper(*arg,**kwargs):
        print("function ",func," is deprecated")
        return func(*arg,**kwargs)
    return wrapper
    
    

def prevent_basis_context(func):
    @wraps(func)
    def wrapper(*arg,**kwargs):
        m = Manager()
        if m._in_eigenbasis_of_context and m._enforce_contexts:
            raise Exception("This function MUST NOT be called"+
                            " from within an 'eigenbasis_of' context.")
        return func(*arg,**kwargs)
    return wrapper        


def enforce_basis_context(func):
    @wraps(func)
    def wrapper(*arg,**kwargs):
        m = Manager()
        if (not m._in_eigenbasis_of_context) and m._enforce_contexts:
            raise Exception("This function MUST be called"+
                            " from within an 'eigenbasis_of' context.")
        return func(*arg,**kwargs)
    return wrapper 


def prevent_energy_units_context(func):
    @wraps(func)
    def wrapper(*arg,**kwargs):
        m = Manager()
        if m._in_energy_units_context and m._enforce_contexts:
            raise Exception("This function MUST NOT be called"+
                            " from within an 'energy_units' context.")
        return func(*arg,**kwargs)
    return wrapper    


def enforce_energy_units_context(func):
    @wraps(func)
    def wrapper(*arg,**kwargs):
        m = Manager()
        if not m._in_energy_units_context and m._enforce_contexts:
            raise Exception("This function MUST be called"+
                            " from within an 'energy_units' context.")
        return func(*arg,**kwargs)
    return wrapper  



