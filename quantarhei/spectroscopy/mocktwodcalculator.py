# -*- coding: utf-8 -*-

import numpy

from .twodcalculator import TwoDResponseCalculator
from .twodcontainer import TwoDResponseContainer
from .pathwayanalyzer import LiouvillePathwayAnalyzer
from .twod2 import TwoDResponse
from ..core.units import convert
from .. import COMPLEX
from .. import signal_REPH, signal_NONR
from .lineshapes import gaussian2D
from .lineshapes import lorentzian2D
from ..core.managers import Manager
from ..core.managers import energy_units

class MockTwoDResponseCalculator(TwoDResponseCalculator):
    """Calculator of the third order non-linear response 
    
    
    This class is used to represent LiouvillePathway objects. Lineshape is
    Gaussian 
    
    """

    def __init__(self, t1axis, t2axis, t3axis):
        super().__init__(t1axis, t2axis, t3axis)
        self.widthx = convert(300, "1/cm", "int")
        self.widthy = convert(300, "1/cm", "int")
        self.dephx = convert(300, "1/cm", "int")
        self.dephy = convert(300, "1/cm", "int")        

        
    def bootstrap(self,rwa=0.0, pathways=None, verbose=False, 
                  shape="Gaussian"):
        
        self.shape = shape
        
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
        """Calculate next spectrum and increment t2 time
        
        """
        sone = self.calculate_one(self.tc)
        self.tc += 1
        return sone

    def set_next(self, tc):
        """Set next current t2 time
        
        """
        self.tc = tc
        
    def calculate_one(self, tc):
        """Calculate the 2D spectrum for all pathways
        
        """
        
        onetwod = TwoDResponse()
        onetwod.set_axis_1(self.oa1)
        onetwod.set_axis_3(self.oa3)
        onetwod.set_resolution("signals")
        
        # First we fill it with zeros
        pwy = None
        data = self.calculate_pathway(pwy, shape=self.shape)
        onetwod._add_data(data, dtype=signal_REPH)
        onetwod._add_data(data, dtype=signal_NONR)

        if self.pathways is not None:        
            for pwy in self.pathways:
                
                data = self.calculate_pathway(pwy, shape=self.shape)
    
                if pwy.pathway_type == "R":
                    onetwod._add_data(data, dtype=signal_REPH)
                elif pwy.pathway_type == "NR":
                    onetwod._add_data(data, dtype=signal_NONR)
                else:
                    raise Exception("Unknown pathway type")

        onetwod.set_t2(self.t2axis.data[tc])    
            
        return onetwod


    def calculate(self):
        """Calculate the 2D spectrum for all pathways
        
        """
        
        onetwod = TwoDResponse()
        onetwod.set_axis_1(self.oa1)
        onetwod.set_axis_3(self.oa3)
        onetwod.set_resolution("signals")
        
        k = 0
        for pwy in self.pathways:
            
            data = self.calculate_pathway(pwy, shape=self.shape)
            
            if pwy.pathway_type == "R":
                onetwod._add_data(data, dtype=signal_REPH)
            elif pwy.pathway_type == "NR":
                onetwod._add_data(data, dtype=signal_NONR)
            else:
                raise Exception("Unknown pathway type")
                
            k += 1
        
        if k == 0:
            pwy = None
            data = self.calculate_pathway(pwy, shape=self.shape)
            onetwod._add_data(data, dtype=signal_REPH)

        onetwod.set_t2(0.0)    
            
        return onetwod


    def calculate_all_system(self, sys, eUt, lab, 
                             selection=None, show_progress=False, dtol=0.0001):
        """Calculates all 2D spectra for a system and evolution superoperator
        
        """
        self.tc = 0
        
        tcont = TwoDResponseContainer(t2axis=self.t2axis)
        
        kk = 1
        Nk = self.t2axis.length
        for T2 in self.t2axis.data:
            
            if show_progress:
                print(" - calculating", kk, "of", Nk, "at t2 =", T2, "fs")
            
            pways = dict()
            twod1 = self.calculate_one_system(T2, sys, eUt, lab,
                                              selection=selection, pways=pways,
                                              dtol=dtol)
        
            if show_progress:
                try:
                    print("   "+str(len(pways[str(T2)])), "pathways used")
                except:
                    pass
                
            if twod1 is not None:
                tcont.set_spectrum(twod1, tag=T2)
            else:
                print("No pathways found; no spectrum calculated")
                
            kk += 1
            
        return tcont


    def calculate_one_system(self, t2, sys, eUt, lab, 
                             selection=None, pways=None, dtol=1.0e-12):
        """Returns 2D spectrum at t2 for a system and evolution superoperator
        
        """
        try:
            Uin = eUt.at(t2)
        except:
            Uin = eUt
            
        H = eUt.get_Hamiltonian()
    
        # FIXME: this needs to be set differently, and it mu
        rho0 = sys.get_DensityMatrix(condition_type="thermal",
                                     temperature=0.0)
        
        # if the Hamiltonian is larger than eUt, we will calculate ESA
        has_ESA = True
        H1 = sys.get_Hamiltonian()
        if H1.dim == eUt.dim:
            has_ESA = False
        
        #print(has_ESA)
        # get Liouville pathways
        if has_ESA:
            pws = sys.liouville_pathways_3T(ptype=("R1g", "R2g", "R3g",
                                                   "R4g", "R1f*", "R2f*"),
                                                   eUt=Uin, ham=H, t2=t2,
                                                   lab=lab, dtol=dtol)
        else:
            pws = sys.liouville_pathways_3T(ptype=("R1g", "R2g", "R3g",
                                                   "R4g"),
                                                   eUt=Uin, ham=H, t2=t2,
                                                   lab=lab, dtol=dtol)
            
        #print("Number of pathways:", len(pws))
            
        if selection is not None:
            
            anl = LiouvillePathwayAnalyzer()
            anl.pathways = pws
            pws = anl.order_by_amplitude(replace=False)
            
            for rule in selection:
                
                if rule[0] == "omega2":
                    interval = rule[1]
                    anl.pathways = pws
                    pws = anl.select_omega2(interval, replace=False)
                    
                if rule[0] == "order":
                    anl.pathways = pws
                    pws = anl.order_by_amplitude(replace=False)
                    
                if rule[0] == "number":
                    N = rule[1]
                    if len(pws) > N:
                        pws = pws[0:N-1]
                        
                if rule[0] == "percent":
                    mx = pws[0].pref*rule[1]/100.0
                    print(mx, numpy.abs(mx))
                    anl.pathways = pws
                    pws = anl.select_amplitude_GT(numpy.abs(mx))
                    
                if rule[0] == "secular":
                    pass
                
                if rule[0] == "incoherent":
                    pws = []
                    for ii in range(len(anl.pathways)):
                        pw = anl.pathways[ii]
                        states = pw.get_states()
                        pair = states[2]
                        if pair[0] == pair[1]:
                            #print(ii, pair)
                            pws.append(pw)
                        
                    
            
        self.set_pathways(pws)
        
        if pways is not None:
            pways[str(t2)] = pws
            
        (n2, err) = self.t2axis.locate(t2)
        
        twod1 = self.calculate_one(n2)
        
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
        
        #print(shape, widthx, widthy)
        
        prefk = 4.0*numpy.log(2.0)
        
        if pathway.pathway_type == "R":

            reph2D = numpy.zeros((N1, N3), dtype=COMPLEX)
            
            
            if shape == "Gaussian":
                oo3 = self.oa3.data[:]
                
                if oldv:
                    for i1 in range(N1):
                        o1 = -self.oa1.data[i1]                    
                        
                        reph2D[:, i1] = \
                            prefk*pref*numpy.exp(
                                    -prefk*(((o1-cen1)/widthx)**2
                                           +((oo3-cen3)/widthy)**2))\
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
