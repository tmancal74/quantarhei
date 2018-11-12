# -*- coding: utf-8 -*-
"""  
    
    
    This is the class representing tightly organized molecular aggregates such
    as photosynthetic antenna and light-harvesting complexes in
    Quantarhei. Appart from represeting data, this class also provides a
    simplified interface to much of Quantarhei's functionality, such as 
    calculation of spectra and dynamics. In order to make the core more 
    organized, the class `Aggregate` is the tip of series of mutually 
    inheriting classes. They start with AggregateBase, a class which implements
    some of the core functionality and add functionality in classes like 
    `AggregateSpectroscopy`, `AggregateExcitonAnalysis` etc.    

    Inheritance in Aggregate class
    ------------------------------
    
    The dependency of the classes is the following
    
        AggregateBase : 
            basic functionality of the Aggregate
                            
        AggregateSpectroscopy : 
            adds Liouville pathway generation

        AggregateExcitonAnalysis :
            adds analysis of excitons

        AggregatePureDephasing :
            adds calculation of effective pure dephasing rates
            
        Aggregate :
            wraps everything up
      

    Class Details
    -------------
              
"""
from .aggregate_pdeph import AggregatePureDephasing


class Aggregate(AggregatePureDephasing):
    
    """
    This clas wraps up the definition of the Aggregate class. It is the end
    of a long series of mutually inheriting classes starting with
    AggregateBase.
           
    """
    
    pass