# -*- coding: utf-8 -*-

import numpy

from .abscalculator import AbsSpectrumCalculator
from .abs2 import AbsSpectrum
from .. import REAL

class MockAbsSpectrumCalculator(AbsSpectrumCalculator):
        
        
    def bootstrap(self,rwa=0.0, pathways=None, verbose=False, 
                  shape="Gaussian"):
        
        self.shape = shape
        
        self.verbose = verbose
        self.rwa = rwa
        self.pathways = pathways

        atype = self.TimeAxis.atype
        self.TimeAxis.atype = 'complete'
        self.oa1 = self.TimeAxis.get_FrequencyAxis() 
        self.oa1.data += self.rwa
        self.oa1.start += self.rwa
        self.TimeAxis.atype = atype
        
        self.tc = 0

        
    def set_width(self, val):
        self.widthx = val


    def set_deph(self, val):
        self.dephx = val


    def set_pathways(self, pathways):
        self.pathways = pathways
        
        
    def calculate(self, raw=False):
        """Calculate the absorption spectrum for all pathways
        
        """
        
        one = AbsSpectrum()
        one.set_axis(self.oa1)
        one.set_data(numpy.zeros(self.TimeAxis.length,
                                dtype=REAL))
        
        k = 0
        for pwy in self.pathways:
            
            print(k)
            data = self.calculate_pathway(pwy, shape=self.shape, raw=raw)
            one.add_data(data)  
            k += 1
            
        return one


    def calculate_pathway(self, pathway, shape="Gaussian", raw=False):
        """Calculate the shape of a Liouville pathway
        
        """
 
        #noe = 1+pathway.order
        
        pref = pathway.pref
        
        cen1 = pathway.frequency[0]
            
        N1 = self.oa1.length
        
        if pathway.widths[1] < 0.0:
            widthx = self.widthx
        else:
            widthx = pathway.widths[1]
                        
        if pathway.dephs[1] < 0.0:
            dephx = self.dephx
        else:
            dephx = pathway.dephs[1]
                    
        data = numpy.zeros(N1, dtype=REAL)
        
        if shape == "Gaussian":
            o1 = self.oa1.data                    
                    
            data[:] = pref*numpy.exp(-((o1-cen1)/widthx)**2) \
                      /(numpy.sqrt(numpy.pi)*widthx)
            
            if not raw:
                data = o1*data
                    
        
        elif shape == "Lorentzian":
            o1 = self.oa1.data                    
                    
            data[:] = pref*(dephx/numpy.pi)/((o1-cen1)**2 + dephx**2)

            if not raw:
                data = o1*data
            
        else:
            raise Exception("Unknown line shape: "+shape)   
        
        return data

        
            
        
        