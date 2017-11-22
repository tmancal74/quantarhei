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
        Aggregate
        
              
"""
from .aggregate_spectroscopy import AggregateSpectroscopy


class Aggregate(AggregateSpectroscopy):
    pass