# -*- coding: utf-8 -*-
"""Test class for SystemBathInteraction


    Provides ready to use classes of the type SystemBathInteraction
    
    
    Class Details
    -------------

"""
from .systembathinteraction import SystemBathInteraction


class TestSystemBathInteraction(SystemBathInteraction):
    """Test class for SystemBathInteraction


    Provides ready to use classes of the type SystemBathInteraction
    
    
    >>> tsbi = TestSystemBathInteraction("trimer-2-lind")
    >>> tsbi.N
    6
    
    >>> tsbi = TestSystemBathInteraction()
    Traceback (most recent call last):
        ...
    Exception: Name of the test must be specified

    """
     
    def __init__(self, name=None):
        
        if name is None:
            raise Exception("Name of the test must be specified")
            
            
        from ...builders.aggregate_test import TestAggregate
        from ...qm.hilbertspace.operators import ProjectionOperator
        
        if name == "dimer-2-env":
            # we get SBI from here
            agg = TestAggregate(name=name)
            agg.build()
    
            super().__init__()
            
            # copy it into the newly created object
            self.__dict__ = agg.sbi.__dict__.copy()

        elif name == "trimer-2-env":
            # we get SBI from here
            agg = TestAggregate(name=name)
            agg.build()
    
            super().__init__()
            
            # copy it into the newly created object
            self.__dict__ = agg.sbi.__dict__.copy()        

        elif name == "dimer-2-lind":
            agg = TestAggregate(name="dimer-2")
            agg.build()
            
            N = agg.get_Hamiltonian().dim
            
            P1 = ProjectionOperator(1, 2, dim=N)
            P2 = ProjectionOperator(2, 1, dim=N)
            
            sys_ops = [P1, P2]
            rates = [1.0/100.0, 1.0/200]
            
            super().__init__(sys_operators=sys_ops, rates=rates)
            
            
        elif name == "trimer-2-lind":
            agg = TestAggregate(name="trimer-2")
            agg.build()
            
            N = agg.get_Hamiltonian().dim
            
            P1 = ProjectionOperator(1, 2, dim=N)
            P2 = ProjectionOperator(2, 1, dim=N)
            P3 = ProjectionOperator(2, 3, dim=N)
            P4 = ProjectionOperator(3, 2, dim=N)
            P5 = ProjectionOperator(1, 3, dim=N)
            P6 = ProjectionOperator(3, 1, dim=N)
            
            sys_ops = [P1, P2, P3, P4, P5, P6]
            rates = [1.0/100.0, 1.0/200, 1.0/100.0, 
                     1.0/300.0, 1.0/200, 1.0/500.0]
            
            super().__init__(sys_operators=sys_ops, rates=rates)
            
        elif name == "dimer-2-lorentz":
            
            N = agg.get_Hamiltonian().dim
            
            P1 = ProjectionOperator(1, 2, dim=N)
            P2 = ProjectionOperator(2, 1, dim=N)
            
            sys_ops = [P1, P2]
            rates = [1.0/100.0, 1.0/200]
            ctimes = []
            
            super().__init__(sys_operators=sys_ops, rates=rates)            
            

        
            
        