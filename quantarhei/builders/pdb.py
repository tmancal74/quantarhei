# -*- coding: utf-8 -*-

class PDBFile:
    
    def __init__(self, fname=None):
        
        if fname is not None:
            
            #load file
            pass
        
        
    def get_Molecules(self, model=None):
        """Returns all molecules corresponding to a given model
        
        """
        
        molecules = []
        
        if model is None:
            return molecules
            
        