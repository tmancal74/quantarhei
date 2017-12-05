# -*- coding: utf-8 -*-
"""

    Quantarhei package 
    ------------------
    
    Author: Tomas Mancal
    Email: tmancal74@gmail.com
    
    
    Aggregate module. Contains definition of the Aggregate class
    
    
    Inheritance
    -----------
    
    The dependency of the classes is the following
    
        AggregateBase
              |
              |
              V
        AggregateSpectroscopy
              |
              |  Adds Liouville pathway generation
              V
        AggregateExcitonAnalysis
              |
              |  Adds analysis of excitons
              V
        Aggregate
        
              
"""
from .aggregate_excitonanalysis import AggregateExcitonAnalysis


class Aggregate(AggregateExcitonAnalysis):
    pass