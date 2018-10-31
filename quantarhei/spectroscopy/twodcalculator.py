# -*- coding: utf-8 -*-

import numpy
import time

from ..utils import derived_type
from ..builders.aggregates import Aggregate
from ..builders.molecules import Molecule
from ..core.time import TimeAxis
from ..core.managers import eigenbasis_of
from ..qm.propagators.poppropagator import PopulationPropagator
from .twod2 import TwoDSpectrum


try:

    import aceto.nr3td as nr3td            
    from aceto.lab_settings import lab_settings
    from aceto.band_system import band_system            
    _have_aceto = True
    
except:
    #
    # FIXME: There should be an optional warning and a fall back onto
    # quantarhei.implementations.aceto module
    #
    #raise Exception("Aceto not available")
    #from ..implementations.aceto import nr3td            
    _have_aceto = False 


class TwoDSpectrumCalculator:
    """Calculator of the 2D spectrum
    
    
    Enables setting up parameters of 2D spectrum calculation for later
    evaluation. The method `calculate` returns TwoDSpectrumContainer
    with a 2D spectrum.
    
    Parameters
    ----------
    
    
    """

    t1axis = derived_type("t1axis",TimeAxis)
    t2axis = derived_type("t2axis",TimeAxis)
    t3axis = derived_type("t3axis",TimeAxis)
    
    system = derived_type("system",[Molecule,Aggregate])
    
    def __init__(self, t1axis, t2axis, t3axis,
                 system=None,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None,
                 effective_hamiltonian=None):
            
            
        self.t1axis = t1axis
        self.t2axis = t2axis
        self.t3axis = t3axis
        
        #FIXME: check the compatibility of the axes 
        
        if system is not None:
            self.system = system
        
        #FIXME: properties to be protected
        self.dynamics = dynamics
        
        # unprotected properties
        self.data = None
        
        self._relaxation_tensor = None
        self._rate_matrix = None
        self._relaxation_hamiltonian = None
        self._has_relaxation_tensor = False
        if relaxation_tensor is not None:
            self._relaxation_tensor = relaxation_tensor
            self._has_relaxation_tensor = True
        if effective_hamiltonian is not None:
            self._relaxation_hamiltonian = effective_hamiltonian
        if rate_matrix is not None:
            self._rate_matrix = rate_matrix
            self._has_rate_matrix = True
            
        #self._have_aceto = False
        
        # after bootstrap information
        self.sys = None
        self.lab = None
        self.t1s = None
        self.t3s = None
        self.rmin = None
        self.rwa = None
        self.oa1 = None
        self.oa3 = None
        self.Uee = None
        self.Uc0 = None
        
        self.tc = 0
        
       
    def _vprint(self, string):
        """Prints a string if the self.verbose attribute is True
        
        """
        if self.verbose:
            print(string)
            
    def bootstrap(self,rwa=0.0, lab=None, verbose=False):
        """Sets up the environment for 2D calculation
        
        """


        self.verbose = verbose
    
    
        if True:
            
            # calculate 2D spectrum using aceto library

            ###############################################################################
            #
            # Create band_system from quantarhei classes
            #
            ###############################################################################
            
            if isinstance(self.system, Aggregate):
            
                pass
            
            else:
                
                raise Exception("Molecule 2D not implememted")
                
            agg = self.system
            
            #
            # hamiltonian and transition dipole moment operators
            #
            H = agg.get_Hamiltonian()
            D = agg.get_TransitionDipoleMoment()
            
            #
            # Construct band_system object
            #
            Nb = 3
            Ns = numpy.zeros(Nb, dtype=numpy.int)
            Ns[0] = 1
            Ns[1] = agg.nmono
            Ns[2] = Ns[1]*(Ns[1]-1)/2
            self.sys = band_system(Nb, Ns)
            
            
            #
            # Set energies
            #
            en = numpy.zeros(self.sys.Ne, dtype=numpy.float64)
            #if True:
            with eigenbasis_of(H):
                for i in range(self.sys.Ne):
                    en[i] = H.data[i,i]
                self.sys.set_energies(en)
            
                #
                # Set transition dipole moments
                #
                dge_wr = D.data[0:Ns[0],Ns[0]:Ns[0]+Ns[1],:]
                def_wr = D.data[Ns[0]:Ns[0]+Ns[1],
                                (Ns[0]+Ns[1]):(Ns[0]+Ns[1]+Ns[2]),:]
            
                dge = numpy.zeros((3,Ns[0],Ns[1]), dtype=numpy.float64)
                deff = numpy.zeros((3,Ns[1],Ns[2]), dtype=numpy.float64)
                
                for i in range(3):
                    dge[i,:,:] = dge_wr[:,:,i]
                    deff[i,:,:] = def_wr[:,:,i]
                self.sys.set_dipoles(0,1,dge)
                self.sys.set_dipoles(1,2,deff)
            
            
            #
            # Relaxation rates
            #
            KK = agg.get_RedfieldRateMatrix()
            
            # relaxation rate in single exciton band
            Kr = KK.data[Ns[0]:Ns[0]+Ns[1],Ns[0]:Ns[0]+Ns[1]] #*10.0
            #print(1.0/Kr)
            
            self.sys.init_dephasing_rates()
            self.sys.set_relaxation_rates(1,Kr)
            
            
            #
            # Lineshape functions
            #
            sbi = agg.get_SystemBathInteraction()
            cfm = sbi.CC
            cfm.create_double_integral()
            
            
            #
            # Transformation matrices
            #
            SS = H.diagonalize()
            SS1 = SS[1:Ns[1]+1,1:Ns[1]+1]
            SS2 = SS[Ns[1]+1:,Ns[1]+1:]
            H.undiagonalize()
            
            self.sys.set_gofts(cfm._gofts)    # line shape functions
            self.sys.set_sitep(cfm.cpointer)  # pointer to sites
            self.sys.set_transcoef(1,SS1)  # matrix of transformation coefficients  
            self.sys.set_transcoef(2,SS2)  # matrix of transformation coefficients  

            #
            # Finding population evolution matrix
            #
            prop = PopulationPropagator(self.t1axis, Kr)
      #      Uee, Uc0 = prop.get_PropagationMatrix(self.t2axis,
      #                                            corrections=True)
            self.Uee, cor = prop.get_PropagationMatrix(self.t2axis,
                                                  corrections=3)

            # FIXME: Order of transfer is set by hand here - needs to be moved
            # to some reasonable place
            
            #Ucor = Uee
            self.Uc0 = cor[0]
            
            #for ko in range(No+1):
            #    print("Subtracting ", ko)
            #    Ucor[:,:,tc] -= cor[ko]

            #
            # define lab settings
            #
            if lab is None:
                self.lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
                X = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
                self.lab.set_laser_polarizations(X,X,X,X)
            else:
                self.lab = lab
            
            #
            # Other parameters
            #
            #dt = self.t1axis.step
            self.rmin = 0.0001
            self.t1s = self.t1axis.data 
            self.t3s = self.t3axis.data
            self.rwa = rwa



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
            
        else:
            
            raise Exception("So far, no 2D outside aceto")
            
        self.tc = 0
            
    def calculate_next(self):

        sone = self.calculate_one(self.tc)
        self.tc += 1
        return sone
    
        
    def calculate_one(self, tc):

        tt2 = self.t2axis.data[tc]        
        Nr1 = self.t1axis.length
        Nr3 = self.t3axis.length        
        #
        # Initialize response storage
        #
        resp_r = numpy.zeros((Nr1, Nr3), 
                             dtype=numpy.complex128, order='F')
        resp_n = numpy.zeros((Nr1, Nr3), 
                             dtype=numpy.complex128, order='F')

        # FIXME: on which axis we should be looking for it2 ??? 
        (it2, err) = self.t1axis.locate(tt2) 
        self._vprint("t2 = "+str(tt2)+"fs (it2 = "+str(it2)+")")
        #tc = it2
    
        #
        # calcute response
        #
        self._vprint("calculating response: ")

        t1 = time.time()
        
        self._vprint(" - ground state bleach")
        # GSB
        nr3td.nr3_r3g(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r) 
        nr3td.nr3_r4g(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
    
        self._vprint(" - stimulated emission")
        # SE
        nr3td.nr3_r1g(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
        nr3td.nr3_r2g(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r)
        
        self._vprint(" - excited state absorption")
        # ESA
        nr3td.nr3_r1fs(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r)
        nr3td.nr3_r2fs(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
        
        # Transfer
        
        Utr = self.Uee[:,:,self.tc]-self.Uc0[:,:,self.tc] #-Uc1[:,:,tc]-Uc2[:,:,tc]
        self.sys.set_population_propagation_matrix(Utr) 
        
        self._vprint(" - stimulated emission with transfer")    
        # SE
        nr3td.nr3_r1g_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
        nr3td.nr3_r2g_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r)
        
#                # This contributes only when No > 0
#                nr3td.nr3_r2g_trN(lab, sys, No, it2, t1s, t3s, rwa, rmin, resp_r)
#                
    
        self._vprint(" - excited state absorption with transfer") 
        # ESA
        nr3td.nr3_r1fs_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r)
        nr3td.nr3_r2fs_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
        
        
        t2 = time.time()
        self._vprint("... calculated in "+str(t2-t1)+" sec")


        #
        # Calculate corresponding 2D spectrum
        #
        
        ftresp = numpy.fft.fft(resp_r,axis=1)
        ftresp = numpy.fft.ifft(ftresp,axis=0)
        reph2D = numpy.fft.fftshift(ftresp)
        
        ftresp = numpy.fft.ifft(resp_n,axis=1)
        ftresp = numpy.fft.ifft(ftresp,axis=0)*ftresp.shape[1]
        nonr2D = numpy.fft.fftshift(ftresp)


        onetwod = TwoDSpectrum()
        onetwod.set_axis_1(self.oa1)
        onetwod.set_axis_3(self.oa3)
        onetwod.set_data(reph2D, dtype="Reph")
        onetwod.set_data(nonr2D, dtype="Nonr")
        
        onetwod.set_t2(self.t2axis.data[tc])
        
        
        return onetwod
                
                
            
    def calculate(self):
        """Returns 2D spectrum
        
        Calculates and returns TwoDSpectrumContainer containing 2D spectrum
        based on the parameters specified in this object.
        
        
        """            
        from .twodcontainer import TwoDSpectrumContainer
                   
        if _have_aceto:

            twods = TwoDSpectrumContainer(self.t2axis)
            
            teetoos = self.t2axis.data
            for tt2 in teetoos:

                onetwod = self.calculate_next()
                twods.set_spectrum(onetwod)   
            
            return twods
        
        else:
            
            # fall back on quantarhei's own implementation
        
            ret = TwoDSpectrumContainer()
            
        
        return ret
    
    

