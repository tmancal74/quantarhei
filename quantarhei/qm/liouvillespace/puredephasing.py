# -*- coding: utf-8 -*-
"""
    Pure dephasing Tensor
    
    
    Class Details
    -------------
    
"""

#
# We create this intentionally as not basis managed. It must be applied
# only while working in the correct basis!
#
#from ...core.managers import BasisManaged

class PureDephasing: #(BasisManaged):
    
    
    def __init__(self, drates=None, dtype="Lorentzian", cutoff_time=None):
        
        self.data=drates
        self.dtype=dtype
        self.cutoff_time=cutoff_time
        
        