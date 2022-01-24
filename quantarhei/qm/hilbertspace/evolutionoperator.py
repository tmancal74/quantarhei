# -*- coding: utf-8 -*-



class EvolutionOperator:
    
    
    def __init__(self, timeaxis, dim=None, hamiltonian=None, data=None):
        
        self.timeaxis = timeaxis
        self.hamiltonian = hamiltonian
        
        pass
    
    
    def apply(opvec):
        """Apply the evolution operator to an operator or state vector
        
        Returns
        -------
        StateVectorEvolution or OperatorEvolution
        
        """
        pass
    
    def apply_left(opvec):
        """Apply the evolution operator to an operator or vector on the left
        
        Returns
        -------
        StateVectorEvolution or OperatorEvolution
        
        """
        pass
    
    