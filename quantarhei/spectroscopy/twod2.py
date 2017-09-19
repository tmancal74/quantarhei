# -*- coding: utf-8 -*-
"""Two-dimensional Fourier Transform Spectrum


"""
import matplotlib.pyplot as plt
import numpy

from ..core.time import TimeAxis
from ..builders.aggregates import Aggregate
from ..builders.molecules import Molecule
from ..core.managers import eigenbasis_of
from ..qm.propagators.poppropagator import PopulationPropagator 

from ..utils import derived_type

import time

class DFunction2:
    """Descrete two-dimensional function 
    
    """
    def __init__(self, x=None, y=None, z=None):
        pass
    
    def save(self, filename):
        pass
    
    def load(self, filename):
        pass
    
    def plot(self,**kwargs):
        pass    


class TwoDSpectrumBase(DFunction2):
    """Basic container for two-dimensional spectrum
    
    
    """
    
    def __init__(self):
        super().__init__()
        self.data = None
        self.xaxis = None
        self.yaxis = None
            
        self.reph2D = None
        self.nonr2D = None


    def set_axis_1(self, axis):
        self.xaxis = axis
        
    def set_axis_3(self, axis):
        self.yaxis = axis
        
    def set_data(self, data, dtype="Tot"):
        if dtype == "Tot":
            self.data = data
            
        elif dtype == "Reph":
            
            self.reph2D = data
            
        elif dtype == "Nonr":
            
            self.nonr2D = data
            
        else:
            
            raise Exception("Unknow type of data: "+dtype)
        
    def save(self, filename):
        super().save(filename)
    
    def load(self, filename):
        super().load(filename)
    
    
class TwoDSpectrumContainer(TwoDSpectrumBase):
    
    def __init__(self):
        pass
    
    
    def plot(self, axis=None, part="ReTot", vmax=None, cbmax = None):
        
        if part == "ReTot":
            # Real part of the total spectrum
            spect2D = numpy.real(self.reph2D) + numpy.real(self.nonr2D)
         
        else:
            raise Exception("Undefined part of the spectrum: "+part)
         
            
        if axis is not None:    
            w1_min = axis[0]
            w1_max = axis[1]
            w3_min = axis[2]
            w3_max = axis[3]

            (i1_min, dist) = self.xaxis.locate(w1_min)
            (i1_max, dist) = self.xaxis.locate(w1_max)

            (i3_min, dist) = self.yaxis.locate(w3_min)
            (i3_max, dist) = self.yaxis.locate(w3_max)   
            
        else:
            i1_min = 0
            i1_max = self.xaxis.length
            i3_min = 0
            i3_max = self.yaxis.length
            
    
        #
        # Plotting with given units on axes
        #
  
        realout = spect2D[i1_min:i1_max,i3_min:i3_max]
    
        Ncontour = 100
        fig, ax = plt.subplots(1,1)
        cm = ax.contourf(self.xaxis.data[i1_min:i1_max],
                     self.yaxis.data[i3_min:i3_max],
                     realout, Ncontour, vmax=vmax)  
        if cbmax is None:
            fig.colorbar(cm)
        else:
            fig.colorbar(cbmax)
            
        return cm
        
        
    
    def devide_by(self, val):
        """Devides the total spectrum by a value
        
        """
        self.reph2D = self.reph2D/val
        self.nonr2D = self.nonr2D/val
        
    
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
            
        self._have_aceto = False
       
    def _vprint(self, string):
        """Prints a string if the self.verbose attribute is True
        
        """
        if self.verbose:
            print(string)
        
    def calculate(self, rwa=0.0, lab=None, verbose=False):
        """Returns 2D spectrum
        
        Calculates and returns TwoDSpectrumContainer containing 2D spectrum
        based on the parameters specified in this object.
        
        
        """
        self.verbose = verbose
        
        try:
            
            import aceto.nr3td as nr3td 
            from aceto.lab_settings import lab_settings
            from aceto.band_system import band_system            
            self._have_aceto = True
            
        except:
            #
            # FIXME: There should be an optional warning and a fall back onto
            # quantarhei.implementations.aceto module
            #
            raise Exception("Aceto not available")
            
            from ..implementations.aceto import nr3td
            
            self._have_aceto = False
    
    
        if self._have_aceto:
            
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
            sys = band_system(Nb, Ns)
            
            
            #
            # Set energies
            #
            en = numpy.zeros(sys.Ne, dtype=numpy.float64)
            #if True:
            with eigenbasis_of(H):
                for i in range(sys.Ne):
                    en[i] = H.data[i,i]
                sys.set_energies(en)
            
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
                sys.set_dipoles(0,1,dge)
                sys.set_dipoles(1,2,deff)
            
            
            #
            # Relaxation rates
            #
            KK = agg.get_RedfieldRateMatrix()
            
            # relaxation rate in single exciton band
            Kr = KK.data[Ns[0]:Ns[0]+Ns[1],Ns[0]:Ns[0]+Ns[1]] #*10.0
            #print(1.0/Kr)
            
            sys.init_dephasing_rates()
            sys.set_relaxation_rates(1,Kr)
            
            
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
            
            sys.set_gofts(cfm._gofts)    # line shape functions
            sys.set_sitep(cfm.cpointer)  # pointer to sites
            sys.set_transcoef(1,SS1)  # matrix of transformation coefficients  
            sys.set_transcoef(2,SS2)  # matrix of transformation coefficients  

            #
            # Finding population evolution matrix
            #
            prop = PopulationPropagator(self.t1axis, Kr)
      #      Uee, Uc0 = prop.get_PropagationMatrix(self.t2axis,
      #                                            corrections=True)
            Uee, cor = prop.get_PropagationMatrix(self.t2axis,
                                                  corrections=3)

            # FIXME: Order of transfer is set by hand here - needs to be moved
            # to some reasonable place
            No = 0
            
            #Ucor = Uee
            Uc0 = cor[0]
            
            #for ko in range(No+1):
            #    print("Subtracting ", ko)
            #    Ucor[:,:,tc] -= cor[ko]

            #
            # define lab settings
            #
            if lab is None:
                lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
                X = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
                lab.set_laser_polarizations(X,X,X,X)
            
            #
            # Other parameters
            #
            #dt = self.t1axis.step
            t1s = self.t1axis.data 
            t3s = self.t3axis.data 

            Nr1 = self.t1axis.length
            Nr3 = self.t3axis.length

            atype = self.t1axis.atype
            self.t1axis.atype = 'complete'
            oa1 = self.t1axis.get_FrequencyAxis() 
            oa1.data += rwa
            oa1.start += rwa
            self.t1axis.atype = atype
            
            atype = self.t3axis.atype
            self.t3axis.atype = 'complete'
            oa3 = self.t3axis.get_FrequencyAxis() 
            oa3.data += rwa
            oa3.start += rwa
            self.t3axis.atype = atype
            
            tc = 0
            twods = []
            
            teetoos = self.t2axis.data
            for tt2 in teetoos:

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
            
                #
                # calcute response
                #
                self._vprint("calculating response: ")
                rmin = 0.0001
                t1 = time.time()
                
                self._vprint(" - ground state bleach")
                # GSB
                nr3td.nr3_r3g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r) 
                nr3td.nr3_r4g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
            
                self._vprint(" - stimulated emission")
                # SE
                nr3td.nr3_r1g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
                nr3td.nr3_r2g(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
                
                self._vprint(" - excited state absorption")
                # ESA
                nr3td.nr3_r1fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
                nr3td.nr3_r2fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
                
                # Transfer
                    
                #sys.set_population_propagation_matrix(Ucor[:,:,tc]) 
                    
                sys.set_population_propagation_matrix(Uee[:,:,tc]-Uc0[:,:,tc]) #-Uc1[:,:,tc]-Uc2[:,:,tc])
                #sys.set_population_propagation_matrix(Uee[:,:,tc])
                
                self._vprint(" - stimulated emission with transfer")    
                # SE
                nr3td.nr3_r1g_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
                nr3td.nr3_r2g_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
                
#                # This contributes only when No > 0
#                nr3td.nr3_r2g_trN(lab, sys, No, it2, t1s, t3s, rwa, rmin, resp_r)
#                
            
                self._vprint(" - excited state absorption with transfer") 
                # ESA
                nr3td.nr3_r1fs_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_r)
                nr3td.nr3_r2fs_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp_n)
                
                
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
        
        
                onetwod = TwoDSpectrumContainer()
                onetwod.set_axis_1(oa1)
                onetwod.set_axis_3(oa3)
                onetwod.set_data(reph2D, dtype="Reph")
                onetwod.set_data(nonr2D, dtype="Nonr")
                
                twods.append(onetwod)
                tc += 1    
            
            return twods
        
        else:
            
            # fall bakc on quantarhei's own implementation
        
            ret = TwoDSpectrumContainer()
            
        
        return ret
    
    
    
        
        