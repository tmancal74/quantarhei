# -*- coding: utf-8 -*-

import numpy

from ..utils import Float
from ..utils import Integer

class ValueAxis:
    """
    Linear array of values which are used as variables of the numerical
    functions and parameters dependent matrices
    """
    
    step = Float("step")  
    noPoint = Integer("noPoints")
    startValue = Float("startValue")
    endValue = Float("endValue")
    
    def __init__(self,startValue=0.0,noPoints=1,step=1.0):
        if step > 0:
            self.step = step
        else:
            raise Exception("Parameters step has to be > 0")
            
        self.data = numpy.linspace(startValue,\
                      startValue+(noPoints-1)*step,noPoints)
                      
        self.noPoints = noPoints
        self.startValue = startValue        
        self.endValue = self.data[noPoints-1]
        
        
        
    def locate(self,x):
        """
        Returns the index of the lower neigbor of the value x 
        and the remaining distance.
        
        """
        
        """ nearest smaller neighbor index """
        n0 = numpy.int(numpy.floor((x-self.startValue)/self.step))
        
        """ if n0 is within bounds calculate distance 
            from the lower neighbor """
        if (n0 >= 0) and (n0 < self.noPoints):
            dval = x-self.data[n0]
            return n0, dval
        else:
            raise Exception("Value out of bounds")
        
        
        
    def nearest(self,x):
        """
        Returns the index of the nearest neighbor

        """         
        
        """ nearest smaller neighbor index """
        n0 = numpy.int(numpy.floor((x-self.startValue)/self.step))
        
        
        if (n0 >= 0) and (n0 < self.noPoints):
            
            """ if n0 is with bounds calculate difference
                from the lower neighbor"""
            diff1 = numpy.abs(x-self.data[n0])
            
            """ if the upper neighbor is within bounds calculate difference
                from the upper neighbor """
            if n0+1 < self.noPoints:
                diff2 = numpy.abs(x - self.data[n0+1])
            else:
                diff2 = 5*self.step
            
            """ return the closer neighbor """
            if diff1 < diff2:
                return n0
            else:
                if n0+1 < self.noPoints:
                    return n0 + 1
                else:
                    """ if the upper neigbor is out of bounds
                        return the lower one """
                    return n0   
                    
        else:
            raise Exception("Value out of bounds")     
