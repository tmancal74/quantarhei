# -*- coding: utf-8 -*-

import scipy.interpolate
from .valueaxis import ValueAxis 

class DFunction:
    """
    Discrete function - discrete representation of a function
    
    """

    allowed_interp_types = ("linear","spline")    
    
    def __init__(self,x,y):
        
        if isinstance(x,ValueAxis):
            self.valueAxis = x
        else:
            raise Exception("First argument has to be of a ValueAxis type")

        if isinstance(y,numpy.ndarray):
            
            if len(y.shape) == 1:
                if y.shape[0] != self.valueAxis.noPoints:
                    raise Exception("Wrong number of elements in 1D numpy.ndarray")
                """
                Set values of the function
                """
                self.data = y

            else:
                raise Exception("Second argument has to be one-dimensional numpy.ndarray")

        else:
            raise Exception("Second argument has to be one-dimensional numpy.ndarray")
         
        self.__splines_initialized = False         
         
            
    #def plot(self,how="-k"):
    #    plt.plot(self.valueAxis.values,self.values,how)
        
    
    def at(self,x,approx="linear"):
        
        if not approx in self.allowed_interp_types:
            raise Exception("Unknown interpolation type")
            
        if approx == "linear":
            return self.__get_linear_approx(x)
        elif approx == "spline":
            return self.__get_spline_approx(x)

    """
    Implementations of various interpolation types
    """        
    def __get_linear_approx(self,x):
        n,dval = self.valueAxis.locate(x)
        if n+1 >= self.valueAxis.noPoints:
            val = self.data[n] \
            + dval/self.valueAxis.step*(self.data[n]-self.data[n-1])
        else:
            val = self.data[n] \
            + dval/self.valueAxis.step*(self.data[n+1]-self.data[n])
        
        return val
        
    def __get_spline_approx(self,x):   
        if not self.__splines_initialized:
            self.__set_splines()
        return self.__spline_value(x)
            
    def __set_splines(self):
        
        self.__spline = \
               scipy.interpolate.UnivariateSpline(\
                  self.valueAxis.data,self.data,s=0)
        self.__splines_initialized = True
        print("Calculating splines")
        
    def __spline_value(self,x):
        return self.__spline(x)
