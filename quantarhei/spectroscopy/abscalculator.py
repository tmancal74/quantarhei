# -*- coding: utf-8 -*-
"""
    Linear absorption spectrum 
    
    Linear absorption spectrum of a molecule or an aggregate of molecules.
    
    
"""
import numpy
import scipy

from ..utils import derived_type
from ..builders import Molecule 
from ..builders import Aggregate
from ..core.time import TimeAxis
from ..core.frequency import FrequencyAxis

from ..core.managers import energy_units
from ..core.managers import EnergyUnitsManaged
from ..core.time import TimeDependent

from .abs2 import AbsSpectrum

class AbsSpectrumCalculator(EnergyUnitsManaged):
    """Linear absorption spectrum 
    
    Parameters
    ----------
    timeaxis : quantarhei.TimeAxis
        TimeAxis on which the calculation will be performed
        
    system : quantarhei.Molecule or quantathei.Aggregate
        System for which the absorption spectrum will be calculated.
        
        
    """

    TimeAxis = derived_type("TimeAxis",TimeAxis)
    system = derived_type("system",[Molecule,Aggregate])
    
    def __init__(self, timeaxis,
                 system=None,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None,
                 effective_hamiltonian=None):
        
        # protected properties
        self.TimeAxis = timeaxis
        self.system = system
        
        #FIXME: properties to be protected
        self.dynamics = dynamics
        
        # unprotected properties
        #self.data = None
        
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
            
        self.rwa = 0.0

     
    def bootstrap(self,rwa=0.0, lab=None):
        """
        
        """
        self.rwa = self.convert_2_internal_u(rwa)
        with energy_units("int"):
            # sets the frequency axis for plottig
            self.frequencyAxis = self.TimeAxis.get_FrequencyAxis()
            self.frequencyAxis.data += self.rwa     
        
        #if isinstance(self.system, Aggregate):
        #    self.system.diagonalize()
                    
        
    def calculate(self, raw=False):
        """ Calculates the absorption spectrum 
        
        
        """
        
        with energy_units("int"):
            if self.system is not None:
                if isinstance(self.system,Molecule):
                    #self._calculate_Molecule(rwa)      
                    spect = self._calculate_monomer(raw=raw)
                elif isinstance(self.system, Aggregate):
                    spect = self._calculate_aggregate( 
                                              relaxation_tensor=
                                              self._relaxation_tensor,
                                              rate_matrix=
                                              self._rate_matrix,
                                              relaxation_hamiltonian=
                                              self._relaxation_hamiltonian,
                                              raw=raw)
            else:
                raise Exception("System to calculate spectrum for not defined")
        
        return spect
    

    def _calculateMolecule(self,rwa):
        
        if self.system._has_system_bath_coupling:
            raise Exception("Not yet implemented")
        else: 
            # calculating stick spectra  

            stick_width = 1.0/0.1
            
            
    def _c2g(self,timeaxis,coft):
        """ Converts correlation function to lineshape function
        
        Explicit numerical double integration of the correlation
        function to form a lineshape function.

        Parameters
        ----------

        timeaxis : cu.oqs.time.TimeAxis
            TimeAxis of the correlation function
            
        coft : complex numpy array
            Values of correlation function given at points specified
            in the TimeAxis object
            
        
        """
        
        ta = timeaxis
        rr = numpy.real(coft)
        ri = numpy.imag(coft)
        sr = scipy.interpolate.UnivariateSpline(ta.data,
                            rr,s=0).antiderivative()(ta.data)
        sr = scipy.interpolate.UnivariateSpline(ta.data,
                            sr,s=0).antiderivative()(ta.data)
        si = scipy.interpolate.UnivariateSpline(ta.data,
                            ri,s=0).antiderivative()(ta.data)
        si = scipy.interpolate.UnivariateSpline(ta.data,
                            si,s=0).antiderivative()(ta.data)
        gt = sr + 1j*si
        return gt
        
    def one_transition_spectrum(self,tr):
        """ Calculates spectrum of one transition
        
        
        """
        ta = tr["ta"] # TimeAxis
        dd = tr["dd"] # transition dipole moment
        om = tr["om"] # frequency - rwa
        gg = tr["gg"] # natural broadening (constant or time dependent)
        
        if self.system._has_system_bath_coupling:
            ct = tr["ct"] # correlation function
        
            # convert correlation function to lineshape function
            gt = self._c2g(ta,ct.data)
            # calculate time dependent response
            at = numpy.exp(-gt -1j*om*ta.data)
        else:
            # calculate time dependent response
            at = numpy.exp(-1j*om*ta.data) 
        
        if len(gg) == 1:
            gam = gg[0]
            rt = numpy.exp(gam*ta.data)
            at *= rt
            #print("Constant: ", rt[20], len(at))
        else:
            rt = numpy.exp((gg)*ta.data)          
            at *= rt
            #print("Time dependent: len = ", rt[20], len(rt))
            
        # Fourier transform the result
        ft = dd*numpy.fft.hfft(at)*ta.step
        ft = numpy.fft.fftshift(ft)
        # invert the order because hfft is a transform with -i
        ft = numpy.flipud(ft)   
        # cut the center of the spectrum
        Nt = ta.length #len(ta.data)        
        return ft[Nt//2:Nt+Nt//2]

        
    def _excitonic_coft(self,SS,AG,n):
        """ Returns energy gap correlation function data of an exciton state 
        
        """
        
        # FIXME: works only for 2 level molecules
        
        c0 = AG.monomers[0].get_egcf((0,1))
        Nt = len(c0)
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        
        ct = numpy.zeros((Nt),dtype=numpy.complex128)
        Na = AG.nmono
        for kk in range(Na):
            
            #nkk = AG.monomers[kk].egcf_mapping[0]
            
            for ll in range(Na):
            
                #nll = AG.monomers[ll].egcf_mapping[0]
                
                ct += ((SS[kk+1,n+1]**2)*(SS[ll+1,n+1]**2)*cfm.get_coft(kk,ll))
                #*AG.egcf_matrix.get_coft(nkk,nll))
            
        return ct


        
    def _calculate_monomer(self, raw=False):
        """ Calculates the absorption spectrum of a monomer 
        
        
        """
        ta = self.TimeAxis
        # transition frequency
        om = self.system.elenergies[1]-self.system.elenergies[0]
        # transition dipole moment
        dm = self.system.dmoments[0,1,:]
        # dipole^2
        dd = numpy.dot(dm,dm)
        # natural life-time from the dipole moment
        gama = [-1.0/self.system.get_electronic_natural_lifetime(1)]
        
        if self.system._has_system_bath_coupling:
            # correlation function
            ct = self.system.get_egcf((0,1))            
            tr = {"ta":ta,"dd":dd,"om":om-self.rwa,"ct":ct,"gg":gama}
        else:
            tr = {"ta":ta,"dd":dd,"om":om-self.rwa,"gg":gama}

        # calculates the one transition of the monomer        
        data = numpy.real(self.one_transition_spectrum(tr))
        

        # we only want to retain the upper half of the spectrum
        Nt = len(self.frequencyAxis.data)//2        
        do = self.frequencyAxis.data[1]-self.frequencyAxis.data[0]
        st = self.frequencyAxis.data[Nt//2]
        # we represent the Frequency axis anew
        axis = FrequencyAxis(st,Nt,do)

        # multiply the spectrum by frequency (compulsory prefactor)
        if not raw:
            data = axis.data*data

        
        spect = AbsSpectrum(axis=axis, data=data)
        
        return spect
        
        
    def _calculate_aggregate(self, relaxation_tensor=None,
                             relaxation_hamiltonian=None, rate_matrix=None,
                             raw=False):
        """ Calculates the absorption spectrum of a molecular aggregate
        
        
        
        """
        ta = self.TimeAxis
        
        # Hamiltonian of the system
        if relaxation_hamiltonian is None:
            HH = self.system.get_Hamiltonian()
        else:
            HH = relaxation_hamiltonian
            

        SS = HH.diagonalize() # transformed into eigenbasis

        
        # Transition dipole moment operator
        DD = self.system.get_TransitionDipoleMoment()
        # transformed into the basis of Hamiltonian eigenstates
        DD.transform(SS)         

        # TimeAxis
        tr = {"ta":ta}
        
        if relaxation_tensor is not None:
            RR = relaxation_tensor
            RR.transform(SS)
            gg = []            
            if isinstance(RR, TimeDependent):
                for ii in range(HH.dim):
                    gg.append(RR.data[:,ii,ii,ii,ii])
            else:
                for ii in range(HH.dim):
                    gg.append([RR.data[ii,ii,ii,ii]])
            tr["gg"] = gg[1]
        elif rate_matrix is not None:
            RR = rate_matrix  # rate matrix is in excitonic basis
            gg = []
            if isinstance(RR, TimeDependent):
                for ii in range(HH.dim):
                    gg.append(RR.data[:,ii,ii])
            else:
                for ii in range(HH.dim):
                    gg.append([RR.data[ii,ii]])
            tr["gg"] = gg[1]
        else:
            tr["gg"] = [0.0]
        
        # get square of transition dipole moment here    #print(H_RC)
        #tr.append(DD.dipole_strength(0,1))
        tr["dd"] = DD.dipole_strength(0,1)
        # first transition energy
        #tr.append(HH.data[1,1]-HH.data[0,0]-rwa)
        tr["om"] = HH.data[1,1]-HH.data[0,0]-self.rwa
        # get a transformed ct here
        ct = self._excitonic_coft(SS,self.system,0)
        #tr.append(ct)
        tr["ct"] = ct
        self.system._has_system_bath_coupling = True
        
        #
        # Calculates spectrum of a single transition
        #
        data = numpy.real(self.one_transition_spectrum(tr))
        
        for ii in range(2,HH.dim):
            if relaxation_tensor is not None:
                tr["gg"] = gg[ii]
            else:
                tr["gg"] = [0.0]
            #tr[1] = DD.dipole_strength(0,ii) # update transition dipole moment
            tr["dd"] = DD.dipole_strength(0,ii)
            #tr[2] = HH.data[ii,ii]-HH.data[0,0]-rwa
            tr["om"] = HH.data[ii,ii]-HH.data[0,0]-self.rwa
            #tr[3] = self._excitonic_coft(SS,self.system,ii-1) # update ct here
            tr["ct"] = self._excitonic_coft(SS,self.system,ii-1)
            
            #
            # Calculates spectrum of a single transition
            #
            data += numpy.real(self.one_transition_spectrum(tr))


        # we only want to retain the upper half of the spectrum
        Nt = len(self.frequencyAxis.data)//2        
        do = self.frequencyAxis.data[1]-self.frequencyAxis.data[0]
        st = self.frequencyAxis.data[Nt//2]
        # we represent the Frequency axis anew
        axis = FrequencyAxis(st,Nt,do)
        
        # multiply the spectrum by frequency (compulsory prefactor)
        if not raw:
            data = axis.data*data
        
        # transform all quantities back
        S1 = numpy.linalg.inv(SS)
        HH.transform(S1)
        DD.transform(S1)
        
        if relaxation_tensor is not None:
            RR.transform(S1)

        spect = AbsSpectrum(axis=axis, data=data)
        
        return spect        

                   
