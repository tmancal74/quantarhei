# -*- coding: utf-8 -*-
import numpy


class Disorder:
    

        
    def __init__(self, data=None, distribution="Gaussian", dtype="diagonal",
                 seed=None):
        
        if data is None:
            raise Exception("Data not specified")
        
        self.dtype = dtype
        self.distribution = distribution
        self.seed_pool = []
        self.data = data.copy()
        self.shape = self.data.shape

        if seed is not None:
            numpy.random.seed(seed)        
        
    def disorder_update(self, i_dis, H, ignore_first=False):
        """Adds disorder to an excitonic Hamiltonian 
    
        """
        if (i_dis == 0) and ignore_first: 
            return    

        N = H.data.shape[0] - 1
        
        if self.dtype == "diagonal":
   
            if self.distribution == "Gaussian":
            
                sigma = self.width/numpy.sqrt(2.0*numpy.log(2))

                de = numpy.random.normal(0.0, sigma, N)
    
            else:
                
                raise Exception("Unknown distribution")
    
            # Update the Hamiltonian energies
            for i in range(N):
                H._data[1+i,1+i] = self.data[1+i,1+i] + de[i]

        else:
            
            raise Exception("Unknown disorder type")
        
        
    def set_distribution(self, distribution, params):
        
        if distribution == "Gaussian":
            
            self.width = params["width"]
            
        else:
            
            raise Exception("Unknown distribution")
