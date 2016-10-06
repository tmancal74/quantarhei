# -*- coding: utf-8 -*-

import numpy
import scipy

from ..utils import derived_type
from ..builders import Molecule 
from ..builders import Aggregate
from ..core.time import TimeAxis
from ..core.frequency import FrequencyAxis
from ..core.dfunction import DFunction

from ..core.managers import energy_units
from ..core.managers import EnergyUnitsManaged

class AbsSpect(DFunction,EnergyUnitsManaged):
    """Linear absorption spectrum 
    
    Linear absorption spectrum of a molecule or an aggregate of molecules.
    
    Parameters
    ----------
    timeaxis : quantarhei.TimeAxis
        TimeAxis on which the calculation will be performed
        
    system : quantarhei.Molecule or quantathei.Aggregate
        System for which the absorption spectrum will be calculated.
        
        
    """

    TimeAxis = derived_type("TimeAxis",TimeAxis)
    system = derived_type("system",[Molecule,Aggregate])
    
    def __init__(self,timeaxis,system=None,dynamics="secular"):
        # protected properties
        self.TimeAxis = timeaxis
        self.system = system
        
        #FIXME: properties to be protected
        self.dynamics = dynamics
        
        # unprotected properties
        self.data = None
        
    def calculate(self,rwa=0.0):
        """ Calculates the absorption spectrum 
        
        
        """
        #rwa = rwa*cm2int
        
        rwa = self.convert_2_internal_u(rwa)
        
        with energy_units("int"):
            if self.system is not None:
                if isinstance(self.system,Molecule):
                    #self._calculate_Molecule(rwa)      
                    self._calculate_monomer(rwa)
                elif isinstance(self.system,Aggregate):
                    self._calculate_aggregate(rwa)
            else:
                raise Exception("System to calculate spectrum for not defined")
        


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
        gg = tr["gg"] # natural broadening
        if self.system._has_system_bath_coupling:
            ct = tr["ct"] # correlation function
        
            # convert correlation function to lineshape function
            gt = self._c2g(ta,ct.data)
            # calculate time dependent response
            at = numpy.exp(-gt -1j*om*ta.data - gg*ta.data)
        else:
            # calculate time dependent response
            at = numpy.exp(-1j*om*ta.data - gg*ta.data) 
        
        # Fourier transform the result
        ft = dd*numpy.fft.hfft(at)*ta.step
        ft = numpy.fft.fftshift(ft)
        # invert the order because hfft is a transform with -i
        ft = numpy.flipud(ft)   
        # cut the center of the spectrum
        Nt = ta.length #len(ta.data)        
        return ft[Nt//2:Nt+Nt//2]
        
    def _frequency(self,dt):
        """ Calculates the frequency axis corresponding to TimeAxis
        
        
        """
        Nt = self.TimeAxis.length
        return numpy.pi*numpy.fft.fftshift(
              numpy.fft.fftfreq(Nt,d=dt))

        
    def _excitonic_goft(self,SS,AG,n):
        """ Returns energy gap correlation function of an exciton state 
        
        """
        
        # FIXME: works only for 2 level molecules
        
        c0 = AG.monomers[0].get_egcf((0,1))
        Nt = len(c0)
        
        ct = numpy.zeros((Nt),dtype=numpy.complex128)
        Na = AG.nmono
        for kk in range(Na):
            nkk = AG.monomers[kk].egcf_mapping[0]
            for ll in range(Na):
                nll = AG.monomers[ll].egcf_mapping[0]
                # FIXME: This can be obtained directly from AG ???
                ct += ((SS[kk+1,n+1]**2)*(SS[ll+1,n+1]**2)
                *AG.egcf_matrix.get_coft(nkk,nll))
            
        return ct


        
    def _calculate_monomer(self,rwa):
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
        gama = 1.0/self.system.get_electronic_natural_lifetime(1)
        
        if self.system._has_system_bath_coupling:
            # correlation function
            ct = self.system.get_egcf((0,1))            
            tr = {"ta":ta,"dd":dd,"om":om-rwa,"ct":ct,"gg":gama}
        else:
            tr = {"ta":ta,"dd":dd,"om":om-rwa,"gg":gama}

        # calculates the one transition of the monomer        
        self.data = numpy.real(self.one_transition_spectrum(tr))
        
        # sets the frequency axis for plottig
        self.frequencyAxis = self.TimeAxis.get_FrequencyAxis()
        self.frequencyAxis.data += rwa
        # we only want to retain the upper half of the spectrum
        Nt = len(self.frequencyAxis.data)//2        
        do = self.frequencyAxis.data[1]-self.frequencyAxis.data[0]
        st = self.frequencyAxis.data[Nt//2]
        # we represent the Frequency axis anew
        self.axis = FrequencyAxis(st,Nt,do)
        
        #self.frequency = self.frequencyAxis.data #self._frequency(ta.dt) + rwa
        self.frequency = self._frequency(ta.step) + rwa
        
        
    def _calculate_aggregate(self,rwa):
        """ Calculates the absorption spectrum of a molecular aggregate
        
        
        
        """
        ta = self.TimeAxis
        
        # Hamiltonian of the system
        HH = self.system.get_Hamiltonian()
        SS = HH.diagonalize() # transformed into eigenbasis

        
        # Transition dipole moment operator
        DD = self.system.get_TransitionDipoleMoment()
        # transformed into the basis of Hamiltonian eigenstates
        DD.transform(SS) 
        
        tr = {"ta":ta}
        tr["gg"] = 0.0
        # get square of transition dipole moment here
        #tr.append(DD.dipole_strength(0,1))
        tr["dd"] = DD.dipole_strength(0,1)
        # first transition energy
        #tr.append(HH.data[1,1]-HH.data[0,0]-rwa)
        tr["om"] = HH.data[1,1]-HH.data[0,0]-rwa
        # get a transformed ct here
        ct = self._excitonic_goft(SS,self.system,0)
        #tr.append(ct)
        tr["ct"] = ct
        self.system._has_system_bath_coupling = True
        self.data = numpy.real(self.one_transition_spectrum(tr))
        
        for ii in range(2,HH.dim):
            #tr[1] = DD.dipole_strength(0,ii) # update transition dipole moment
            tr["dd"] = DD.dipole_strength(0,ii)
            #tr[2] = HH.data[ii,ii]-HH.data[0,0]-rwa
            tr["om"] = HH.data[ii,ii]-HH.data[0,0]-rwa
            #tr[3] = self._excitonic_goft(SS,self.system,ii-1) # update ct here
            tr["ct"] = self._excitonic_goft(SS,self.system,ii-1)
            self.data += numpy.real(self.one_transition_spectrum(tr))

        # sets the frequency axis for plottig
        self.frequencyAxis = self.TimeAxis.get_FrequencyAxis()
        self.frequencyAxis.data += rwa
        # we only want to retain the upper half of the spectrum
        Nt = len(self.frequencyAxis.data)//2        
        do = self.frequencyAxis.data[1]-self.frequencyAxis.data[0]
        st = self.frequencyAxis.data[Nt//2]
        # we represent the Frequency axis anew
        self.axis = FrequencyAxis(st,Nt,do)
        
        self.frequency = self._frequency(ta.step) + rwa
        
        
    def normalize2(self,norm=1.0):
        mx = numpy.max(self.data)
        self.data = self.data/mx
        
        
    def plot(self,**kwargs):
        """ Plotting absorption spectrum using the DFunction plot method
        
        """
        if "ylabel" not in kwargs:
            ylabel = r'$\alpha(\omega)$ [a.u.]'
            kwargs["ylabel"] = ylabel
            
        super(AbsSpect,self).plot(**kwargs)
        
        
                    
        