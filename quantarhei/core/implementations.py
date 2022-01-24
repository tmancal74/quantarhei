# -*- coding: utf-8 -*-
from functools import wraps
import os
from importlib import import_module

from .managers import Manager


def implementation(package="",
                   taskname="",
                   at_runtime=False,
                   fallback_local=False,
                   always_local=False):
    """Decorator to select numerical implememtation
    
    
    """
    m = Manager()
    
    def decorate_at_runtime(func):
        """Decoration at run time 
        
        The wrapper decides which function to return at runtime.

        """           
        @wraps(func)
        def wrapper(*arg,**kwargs):
            fc = get_function(func, package, taskname,
                              default_local=fallback_local,
                              always_local=always_local)
            return fc(*arg,**kwargs)
                
        return wrapper
        
        
    def decorate_at_loadtime(func):
        """Decoration at load time 
        
        The wrapper decides which function to return when the Manager module
        is loaded, i.e. at the start of the application.

        """           
 
        fc = get_function(func, package, taskname,
                          default_local=fallback_local,
                          always_local=always_local) 

        @wraps(func)
        def wrapper(*arg,**kwargs):
            return fc(*arg,**kwargs)
                
        return wrapper
        
        
    if (at_runtime and m.change_implementation_at_runtime):
        
        return decorate_at_runtime
    
    else:
        
        return decorate_at_loadtime
    
    
    
#
#  Auxiliary function
#    
    
def load_function(lib,fce):
    """Load the module and get the desired function
    
    """
    try: 
        a =  import_module(lib)
    except:
        print("Cannot load module", lib)
        
    if hasattr(a,fce):
        fc = getattr(a,fce)
    else:
        raise Exception("Cannot reach implementation of %s " % fce)


    return fc    

def get_function(func, package, taskname, default_local, always_local):
    """Decide which function to use
    
    
    """
    if always_local:
        return func
        
    m = Manager()
    # default implementation package
    default_imp_prefix = "quantarhei.implementations.python"
            
    # decide which implementation will be used
    imp_prefix = m.get_implementation_prefix(package=package,
                                             taskname=taskname)
    
    # load the package
    try:
        imp_name = imp_prefix + "." + package
        fc = load_function(imp_name, taskname)
        
    except:
    
        try:
            # fall back on pure Python implementation
            if default_local:
                fc = func
            else:
                imp_name = default_imp_prefix + "." + package
                fc = load_function(imp_name,taskname)

            # FIXME: issue a warning
            print("WARNING: import failed, falling back on pure Python")
        except:
            # do not provide implementation, call the decorated function itself
            # FIXME: issue a warning (this is an unwanted result)
            print("WARNING: calling decorated function itself")
            fc = func
            
    return fc
