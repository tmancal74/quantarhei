# -*- coding: utf-8 -*-

import numpy

from .twodcalculator import TwoDSpectrumCalculator
from .twodcontainer import TwoDSpectrumContainer
from .twod2 import TwoDSpectrum
from ..core.units import convert
from .. import COMPLEX
from .lineshapes import gaussian2D
from .lineshapes import lorentzian2D
from ..core.managers import Manager
from ..core.managers import energy_units

class MockTwoDSpectrumCalculator(TwoDSpectrumCalculator):
    """Calculator of the 2D spectrum from LiouvillePathway objects
    
    
    This class is used to represent LiouvillePatjway objects. Lineshape is
    Gaussian 
    
    """

    def __init__(self, t1axis, t2axis, t3axis,temp=None):
        #t2axis = TimeAxis()
        super().__init__(t1axis, t2axis, t3axis)
        self.widthx = convert(300, "1/cm", "int")
        self.widthy = convert(300, "1/cm", "int")
        self.dephx = convert(300, "1/cm", "int")
        self.dephy = convert(300, "1/cm", "int")   
        self.temp = temp # Temperature

        
    def bootstrap(self,rwa=0.0, pathways=None, verbose=False, 
                  shape="Gaussian", all_positive=False):
        
        self.shape = shape
        self.all_positive = all_positive
        
        self.verbose = verbose
        self.rwa = Manager().convert_energy_2_internal_u(rwa)
        self.pathways = pathways

        with energy_units("int"):
            atype = self.t1axis.atype
            self.t1axis.atype = 'complete'
            self.oa1 = self.t1axis.get_FrequencyAxis() 
            self.oa1.data += self.rwa
            self.oa1.start += self.rwa
            self.t1axis.atype = atype
            
            atype = self.t3axis.atype
            self.t3axis.atype = 'complete'
            self.oa3 = self.t3axis.get_FrequencyAxis() 
            self.oa3.data += self.rwa
            self.oa3.start += self.rwa
            self.t3axis.atype = atype        
        
        self.tc = 0
            

    def set_width(self, val):
        m = Manager()
        self.widthx = m.convert_energy_2_internal_u(val)
        self.widthy = m.convert_energy_2_internal_u(val)
        
    def set_deph(self, val):
        m = Manager()
        self.dephx = m.convert_energy_2_internal_u(val)
        self.dephy = m.convert_energy_2_internal_u(val)


    def set_pathways(self, pathways):
        self.pathways = pathways
        
        
    def calculate_next(self):

        sone = self.calculate_one(self.tc)
        #print(self.tc, sone)
        self.tc += 1
        return sone

    def set_next(self, tc):
        self.tc = tc
        
    def calculate_one(self, tc):
        """Calculate the 2D spectrum for all pathways
        
        """
        
        onetwod = TwoDSpectrum()
        onetwod.set_axis_1(self.oa1)
        onetwod.set_axis_3(self.oa3)
        onetwod.set_resolution("signals")

        k = 0        
        for pwy in self.pathways:
            
            data = self.calculate_pathway(pwy, shape=self.shape)

            if pwy.pathway_type == "R":
                onetwod._add_data(data, dtype="REPH")
            elif pwy.pathway_type == "NR":
                onetwod._add_data(data, dtype="NONR")
            else:
                raise Exception("Unknown pathway type")

            k += 1
        
        if k == 0:
            pwy = None
            data = self.calculate_pathway(pwy, shape=self.shape)
            onetwod._add_data(data, dtype="REPH")
            print("Warning: calculating empty 2D spectrum")

        onetwod.set_t2(self.t2axis.data[tc])    
            
        return onetwod


    def calculate(self):
        """Calculate the 2D spectrum for all pathways
        
        """
        
        onetwod = TwoDSpectrum()
        onetwod.set_axis_1(self.oa1)
        onetwod.set_axis_3(self.oa3)
        onetwod.set_resolution("signals")
        
        k = 0
        for pwy in self.pathways:
            
            data = self.calculate_pathway(pwy, shape=self.shape)
            
            if pwy.pathway_type == "R":
                onetwod._add_data(data, dtype="REPH")
            elif pwy.pathway_type == "NR":
                onetwod._add_data(data, dtype="NONR")
            else:
                raise Exception("Unknown pathway type")
                
            k += 1
        
        if k == 0:
            pwy = None
            data = self.calculate_pathway(pwy, shape=self.shape)
            onetwod._add_data(data, dtype="REPH")

        onetwod.set_t2(0.0)    
            
        return onetwod


    def calculate_all_system(self, sys, H, eUt, lab, show_progress=False):
        """Calculates all 2D spectra for a system and evolution superoperator
        
        """
        tcont = TwoDSpectrumContainer(t2axis=self.t2axis)
        
        kk = 1
        Nk = self.t2axis.length
        for T2 in self.t2axis.data:
            
            if show_progress:
                print(" - calculating", kk, "of", Nk, "at t2 =", T2, "fs")
            
            twod1 = self.calculate_one_system(T2, sys, H, eUt, lab)
        
            tcont.set_spectrum(twod1, tag=T2)
            kk += 1
            
        return tcont


    def calculate_one_system(self, t2, sys, H, eUt, lab, pways=None):
        """Returns 2D spectrum at t2 for a system and evolution superoperator
        
        """
        try:
            Uin = eUt.at(t2)
        except:
            Uin = eUt
    
        # FIXME: this needs to be set differently, and it mu
        # density matrix stored into sys.rho0
        #temp = sys.get_temperature()
        if self.temp is None:
            rho0 = sys.get_DensityMatrix(condition_type="thermal",
                                     temperature=0.0)
        else:
            rho0 = sys.get_DensityMatrix(condition_type="thermal",
                                     temperature=self.temp)
        
        # if the Hamiltonian is larger than eUt, we will calculate ESA
        has_ESA = True
        H1 = sys.get_Hamiltonian()
#        if H1.dim == eUt.dim:
#            has_ESA = False
        
        # get Liouville pathways
        if has_ESA:
            pws = sys.liouville_pathways_3T(ptype=("R1g", "R2g", "R3g",
                                                   "R4g", "R1f*", "R2f*",
                                                   "R1gE", "R2gE"),
                                                   eUt=Uin, ham=H, t2=t2,
                                                   lab=lab)
        else:
            pws = sys.liouville_pathways_3T(ptype=("R1g", "R2g", "R3g",
                                                   "R4g", "R1gE", "R2gE"),
                                                   eUt=Uin, ham=H, t2=t2,
                                                   lab=lab)

        self.set_pathways(pws)
        
        if pways is not None:
            pways[str(t2)] = pws
            
        twod1 = self.calculate_next()
        
        return twod1
        

    def calculate_pathway(self, pathway, shape="Gaussian"):
        """Calculate the shape of a Liouville pathway
        
        """
        
        # we can calculate empty pathway
        if pathway is None:
            N1 = self.oa1.length
            N3 = self.oa3.length            
            reph2D = numpy.zeros((N1, N3), dtype=COMPLEX)
            return reph2D
        
        # FIXME: remove the old version sometime soon
        oldv = False
        
        noe = 1+pathway.order+pathway.relax_order 
        
        cen1 = pathway.frequency[0]
        cen3 = pathway.frequency[noe-2]
        if self.all_positive:
            pref = numpy.abs(pathway.pref)
        else:
            pref = pathway.pref
            
        N1 = self.oa1.length
        N3 = self.oa3.length
        
        if pathway.widths[1] < 0.0:
            widthx = self.widthx
        else:
            widthx = pathway.widths[1]
            
        if pathway.widths[3] < 0.0:
            widthy = self.widthy
        else:
            widthy = pathway.widths[3]
            
        if pathway.dephs[1] < 0.0:
            dephx = self.dephx
        else:
            dephx = pathway.dephs[1]
            
        if pathway.widths[3] < 0.0:
            dephy = self.dephy
        else:
            dephy = pathway.dephs[3]
        
        if (widthx < 1e-5) or (widthy < 1e-5):
            print("Shape, widthx, widthy:",shape, widthx, widthy)
        
        if pathway.pathway_type == "R":

            reph2D = numpy.zeros((N1, N3), dtype=COMPLEX)
            
            if shape == "Gaussian":
                oo3 = self.oa3.data[:]
                
                if oldv:
                    for i1 in range(N1):
                        o1 = -self.oa1.data[i1]                    
                        
                        reph2D[:, i1] = \
                        pref*numpy.exp(-((o1-cen1)/widthx)**2)\
                            *numpy.exp(-((oo3-cen3)/widthy)**2)\
                            /(numpy.pi*widthx*widthy)
            
                else:

                    oo1 = -self.oa1.data[:]

                    reph2D = pref*gaussian2D(oo1, cen1, widthx,
                                                 oo3, cen3, widthy)

                    
            elif shape == "Lorentzian":
                oo3 = self.oa3.data[:]
                
                if oldv:
                    for i1 in range(N1):
                        o1 = -self.oa1.data[i1]                    
                        
                        reph2D[:, i1] = \
                        pref*((dephx/numpy.pi)/((o1-cen1)**2 + dephx**2))\
                            *((dephy/numpy.pi)/((oo3-cen3)**2 + dephy**2))
                    
                else:

                    oo1 = -self.oa1.data[:]
                    
                    reph2D = pref*lorentzian2D(oo1, cen1, dephx,
                                               oo3, cen3, dephy)
                    
            else:
                raise Exception("Unknown line shape: "+shape)   
            
            return reph2D
            
        elif pathway.pathway_type == "NR":
           
            nonr2D = numpy.zeros((N1, N3), dtype=COMPLEX)
            
            if shape == "Gaussian":
                oo3 = self.oa3.data[:]
                
                if oldv:
                    for i1 in range(N1):
                        o1 = self.oa1.data[i1]                    
                    
                        nonr2D[:, i1] = \
                        pref*numpy.exp(-((o1-cen1)/widthx)**2)\
                            *numpy.exp(-((oo3-cen3)/widthy)**2)\
                                /(numpy.pi*widthx*widthy)
                            
                else:

                    oo1 = self.oa1.data[:]                    
                    
                    nonr2D = pref*gaussian2D(oo1, cen1, widthx,
                                               oo3, cen3, widthy)
                        
            elif shape == "Lorentzian":
                oo3 = self.oa3.data[:]
                
                if oldv:
                    for i1 in range(N1):
                        o1 = self.oa1.data[i1]                    
                        
                        nonr2D[:, i1] = \
                        pref*((dephx/numpy.pi)/((o1-cen1)**2 + dephx**2))\
                                *((dephy/numpy.pi)/((oo3-cen3)**2 + dephy**2))
                                
                else:

                    oo1 = self.oa1.data[:]
                    
                    nonr2D = pref*lorentzian2D(oo1, cen1, dephx,
                                               oo3, cen3, dephy)

            else:
                raise Exception("Unknown line shape: "+shape)
            
            return nonr2D
