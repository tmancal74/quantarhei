# -*- coding: utf-8 -*-

import numpy
import scipy

from ..utils import derived_type
from ..builders import Molecule 
from ..builders import Aggregate
from ..core.time import TimeAxis

class AbsSpect:
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
        if self.system is not None:
            if isinstance(self.system,Molecule):
                #self._calculate_Molecule(rwa)      
                self._calculate_monomer(rwa)
            elif isinstance(self.system,Aggregate):
                self._calculate_aggregate(rwa)
        else:
            raise Exception("System to calculate spectrum for is not defined")
        


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
        sr = scipy.interpolate.UnivariateSpline(ta.time,
                            rr,s=0).antiderivative()(ta.time)
        sr = scipy.interpolate.UnivariateSpline(ta.time,
                            sr,s=0).antiderivative()(ta.time)
        si = scipy.interpolate.UnivariateSpline(ta.time,
                            ri,s=0).antiderivative()(ta.time)
        si = scipy.interpolate.UnivariateSpline(ta.time,
                            si,s=0).antiderivative()(ta.time)
        gt = sr + 1j*si
        return gt
        
    def one_transition_spectrum(self,tr):
        """ Calculates spectrum of one transition
        
        
        """
        ta = tr[0] # TimeAxis
        dd = tr[1] # transition dipole moment
        om = tr[2] # frequency - rwa
        ct = tr[3] # correlation function
        
        # convert correlation function to lineshape function
        gt = self._c2g(ta,ct.data)
        # calculate time dependent response
        at = numpy.exp(-gt -1j*om*ta.time) 
        # Fourier transform the result
        ft = dd*numpy.fft.hfft(at)*ta.dt
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
        Na = AG.nmono #monomers.shape[0]
        for kk in range(Na):
            nkk = AG.monomers[kk].egcf_mapping[0]
            for ll in range(Na):
                nll = AG.monomers[ll].egcf_mapping[0]
                ct += ((SS[kk+1,n+1]**2)*(SS[ll+1,n+1]**2)
                *AG.egcf_matrix.get_coft(nkk,nll))
            
        return ct


        
    def _calculate_monomer(self,rwa):
        """ Calculates the absorption spectrum of a monomer 
        
        
        """
        ta = self.TimeAxis
        # correlation function
        ct = self.system.get_egcf((0,1))
        # transition frequency
        om = self.system.elenergies[1]-self.system.elenergies[0]
        # transition dipole moment
        dm = self.system.dmoments[0,1,:]
        # dipole^2
        dd = numpy.dot(dm,dm)
        
        tr = list()
        tr.append(ta)
        tr.append(dd)
        tr.append(om-rwa)
        tr.append(ct)

        # calculates the one transition of the monomer        
        self.data = numpy.real(self.one_transition_spectrum(tr))
        # sets the frequency axis for plottig
        self.frequency = self._frequency(ta.dt) + rwa
        
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
        
        tr = list()
        tr.append(ta) 
        # get square of transition dipole moment here
        tr.append(DD.dipole_strength(0,1))  
        # first transition energy
        tr.append(HH.data[1,1]-HH.data[0,0]-rwa)
        # get a transformed ct here
        ct = self._excitonic_goft(SS,self.system,0)
        tr.append(ct)
        self.data = numpy.real(self.one_transition_spectrum(tr))
        for ii in range(2,HH.dim):
            tr[1] = DD.dipole_strength(0,ii-1) # update transition dipole moment
            tr[2] = HH.data[ii,ii]-HH.data[0,0]-rwa
            tr[3] = self._excitonic_goft(SS,self.system,ii-1) # update ct here
            self.data += numpy.real(self.one_transition_spectrum(tr))
        
        self.frequency = self._frequency(ta.dt) + rwa
        
        
    def normalize2(self,norm=1.0):
        mx = numpy.max(self.data)
        self.data = self.data/mx
                    
        