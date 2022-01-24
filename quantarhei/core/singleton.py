# -*- coding: utf-8 -*-

class Singleton(type):
    """Base type of singletons, such as the main Quantarhei class Manager
    
    
    Recipe "Creating a singleton in Python" from Stack Overflow.
    
    
    Usage:
    ------
    class MyClass(BaseClass, metaclass=Singleton)
    
    """

    _instances = {}
    
    def __call__(cls,*args,**kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton,cls).__call__(*args,
                                                                **kwargs)
        return cls._instances[cls]
        
        
        
