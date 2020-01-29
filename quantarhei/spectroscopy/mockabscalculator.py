# -*- coding: utf-8 -*-

import numpy

from .abscalculator import AbsSpectrumCalculator
from . import lineshapes 
from .abs2 import AbsSpectrum
from .. import REAL
from ..spectroscopy.labsetup import LabSetup
from ..builders.aggregates import Aggregate

class MockAbsSpectrumCalculator(AbsSpectrumCalculator):
        
        
    def bootstrap(self,rwa=0.0, pathways=None, lab=None, 
                  shape="Gaussian", verbose=False):
        """Prepare for the spectrum calculation
        
        
        """
        self.shape = shape
        
        self.verbose = verbose
        self.rwa = rwa
        self.pathways = pathways

        # TimeAxis is set in the constructor (__init__ method)
        atype = self.TimeAxis.atype
        self.TimeAxis.atype = 'complete'
        self.oa1 = self.TimeAxis.get_FrequencyAxis() 
        self.oa1.data += self.rwa
        self.oa1.start += self.rwa
        self.TimeAxis.atype = atype
        
        self.tc = 0

        if lab is None:
            lab = LabSetup()
            
        if isinstance(self.system, Aggregate):
            self.system.diagonalize()
            
        if pathways is None:
            # FIXME: this needs to be made in a more transparent manner  
            rho0 = self.system.get_DensityMatrix(condition_type="thermal",
                                                 temperature=0.0)
            ham = self.system.get_Hamiltonian()
            pthways = self.system.liouville_pathways_1(lab=lab, ham=ham,
                                                       etol=1.0e-5,
                                                       verbose=0) 
            self.set_pathways(pthways)
        
    # FIXME: energy/frequency units
    def set_width(self, val):
        """Set the spectral width as FWHM
        
        
        """
        self.widthx = val


    def set_deph(self, val):
        """Set dephasing rate as 1/(dephasing time)
        
        
        """        
        self.dephx = val


    def set_pathways(self, pathways):
        """Set Liouville pathways to be used for spectral calculation
        
        """
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
            
            data = self.calculate_pathway(pwy, shape=self.shape, raw=raw)
            one.add_data(data)  
            k += 1
            
        return one


    def calculate_pathway(self, pathway, shape="Gaussian", raw=False):
        """Calculate the shape of a Liouville pathway
        
        
        Parameters
        ----------
        
        pathway : Liouville pathways object
            A Liouville pathway which spectrum will be calculated
            
        shape : str
            One of the "Gaussian", "Lorentzian" and "Voigt" line shapes
            
        raw : bool
            If set to True, no frequency dependent prefactor will be overlayed
            over the returned spectrum
            
        """ 
        
        pref = pathway.pref
        cen1 = pathway.frequency[0]
        N1 = self.oa1.length
        
        if pathway.widths[1] <= 0.0:
            widthx = self.widthx
        else:
            widthx = pathway.widths[1]
        
        if pathway.dephs[1] < 0.0:
            dephx = self.dephx
        else:
            dephx = pathway.dephs[1]
                    
        data = numpy.zeros(N1, dtype=REAL)
        o1 = self.oa1.data 
        
        if shape == "Gaussian":
                              
            data[:] = pref*lineshapes.gaussian(o1, cen1, widthx)
                      #numpy.sqrt(numpy.log(2.0)/numpy.pi)\
                      #*numpy.exp(-numpy.log(2.0)*((o1-cen1)/widthx)**2) \
                      #/widthx
            
        elif shape == "Lorentzian":                   
                   
            if (dephx is None) or (dephx <= 0.0):
                if pathway.widths[1] < 0.0:
                    widthx = self.widthx
                else:
                    widthx = pathway.widths[1]
                
                # by width we specify FWHM, but gamma is HWHM
                dephx = 2.0*widthx  
                
            data[:] = pref*lineshapes.lorentzian(o1, cen1, dephx)
                      #(dephx/numpy.pi)/((o1-cen1)**2 + dephx**2)
                
        elif shape == "Voigt":
            
            #z = (o1 - cen1 + 1j*dephx)*numpy.sqrt(numpy.log(2.0))/widthx
            
            data[:] = pref*lineshapes.voigt(o1, cen1, widthx, dephx)
                      #numpy.sqrt(numpy.log(2.0))*\
                      #numpy.real(special.wofz(z)) \
                      #/(numpy.sqrt(numpy.pi)*widthx)
            
        else:
            raise Exception("Unknown line shape: "+shape)   

        if not raw:
            data = o1*data
                
        return data

        
            
        
        