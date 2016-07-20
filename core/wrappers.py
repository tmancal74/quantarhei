# -*- coding: utf-8 -*-
from functools import wraps
import os
   
   
def deprecated(func):
    @wraps(func)
    def wrapper(*arg,**kwargs):
        print("function ",func," is deprecated")
        return func(*arg,**kwargs)
    return wrapper
    
    
    