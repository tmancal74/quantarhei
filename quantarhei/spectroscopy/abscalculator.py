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
from ..builders import OpenSystem
from ..core.time import TimeAxis
from ..core.frequency import FrequencyAxis
from ..core.managers import energy_units
from ..core.managers import EnergyUnitsManaged
from ..core.time import TimeDependent

from ..core.wrappers import prevent_basis_context

from ..qm.hilbertspace.operators import ReducedDensityMatrix
from .abs2 import AbsSpectrum
from .. import COMPLEX, REAL


class AbsSpectrumCalculator(EnergyUnitsManaged):
    """Linear absorption spectrum 
    
    Parameters
    ----------
    timeaxis : quantarhei.TimeAxis
        TimeAxis on which the calculation will be performed
        
    system : quantarhei.Molecule or quantathei.Aggregate
        System for which the absorption spectrum will be calculated.
        
        
    Examples
    --------

    Calcutor has to be created using a Molecule, Aggregate or OpenSystem

    >>> time = TimeAxis(0.0, 1000, 1.0)
    >>> absc = AbsSpectrumCalculator(time)
    Traceback (most recent call last):
        ...
    TypeError: system must be of type [<class 'quantarhei.builders.molecules.Molecule'>, <class 'quantarhei.builders.aggregates.Aggregate'>, <class 'quantarhei.builders.opensystem.OpenSystem'>]
    
    >>> with energy_units("1/cm"):
    ...     mol = Molecule([0.0, 10000.0])
    >>> absc = AbsSpectrumCalculator(time, mol)
    
    The method bootstrap() must be called before calculate()
    
    >>> abs = absc.calculate() 
    Traceback (most recent call last):
        ...    
    Exception: Calculator must be bootstrapped first: call bootstrap() method of this object.
    
    
    Molecule does not provide RWA information automatically
    
    >>> HH = mol.get_Hamiltonian()
    >>> HH.has_rwa
    False
    
    Setting the RWA frequency can be therefore done explicitely
    
    >>> absc.bootstrap(rwa=1.0)
    >>> print(absc.rwa)
    1.0
    
    When RWA is specified for the Hamiltonian
    
    >>> HH.set_rwa([0, 1])
    
    setting RWA explicitely in bootstrap() is ignored. The value from the 
    Molecule is used.
    
    >>> absc.bootstrap(rwa=1.1)
    >>> print("%5.5f" % absc.rwa)
    1.88365

    RWA can be set by the Molecule
    
    >>> with energy_units("1/cm"):
    ...     mol = Molecule([0.0, 10000.0])
    >>> absc = AbsSpectrumCalculator(time, mol)
    >>> mol.set_electronic_rwa([0, 1])
    >>> absc.bootstrap()
    
    """

    TimeAxis = derived_type("TimeAxis",TimeAxis)
    system = derived_type("system",[Molecule,Aggregate,OpenSystem])
    
    def __init__(self, timeaxis,
                 system=None,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None,
                 effective_hamiltonian=None):
        
        # protected properties
        self.TimeAxis = timeaxis
        self.system = system
        
        if dynamics in ["secular", "non-secular"]:
            self.dynamics = dynamics
        else:
            raise Exception("Unknown type of dynamics: "+str(dynamics))
                    
        # unprotected properties
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
        self.prop_has_rwa = False
        
        self.bootstrapped = False

     
    def bootstrap(self, rwa=0.0, prop=None, lab=None):
        """This function sets some additional information before calculation
        
        
        >>> time = TimeAxis(0.0, 1000, 1.0)
        >>> with energy_units("1/cm"):
        ...     mol = Molecule([0.0, 10000.0])
        >>> absc = AbsSpectrumCalculator(time, mol)  
        
        >>> absc.bootstrap()
        Traceback (most recent call last):
            ...
        Exception: RWA not set by system nor explicitely.
        """
        
        self.prop = prop
        
        HH = self.system.get_Hamiltonian()
        if HH.has_rwa:
            # if HH has RWA, the 'rwa' argument is ignored
            HR = HH.get_RWA_skeleton()
            Ng = HH.rwa_indices[0]
            Ne = HH.rwa_indices[1]
            self.rwa = self.convert_2_internal_u(HR[Ne]-HR[Ng])
            self.prop_has_rwa = True  # Propagator is in RWA
        else:
            if rwa > 0.0:
                self.rwa = self.convert_2_internal_u(rwa)
            else:
                raise Exception("RWA not set by system nor explicitely.")

        #self.rwa = self.convert_2_internal_u(rwa)
            
        with energy_units("int"):
            # sets the frequency axis for plottig
            self.frequencyAxis = self.TimeAxis.get_FrequencyAxis()
            
            # it must be shifted by rwa value
            self.frequencyAxis.data += self.rwa     
            
        self.bootstrapped = True

    @prevent_basis_context                        
    def calculate(self, raw=False, from_dynamics=False, alt=False):
        """ Calculates the absorption spectrum 
        
        
        """
        
        if not self.bootstrapped:
            raise Exception("Calculator must be bootstrapped first: "+
                            "call bootstrap() method of this object.")
        
        with energy_units("int"):
            
            if self.system is not None:
                
                if from_dynamics:
                    
                    # alt = True is for testing only
                    spect = self._calculate_abs_from_dynamics(raw=raw,
                                                              alt=alt)
                        
                elif isinstance(self.system, Aggregate):
                    spect = self._calculate_aggregate( 
                                            relaxation_tensor=
                                            self._relaxation_tensor,
                                            rate_matrix=
                                            self._rate_matrix,
                                            relaxation_hamiltonian=
                                            self._relaxation_hamiltonian,
                                            raw=raw)
                    
                elif isinstance(self.system, (Molecule, OpenSystem)):
                    #self._calculate_Molecule(rwa) 
                    spect = self._calculate_monomer(raw=raw)

            else:
                raise Exception("System to calculate spectrum for not defined")
        
        return spect

        
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
            gt = _c2g(ta,ct.data)
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
        
        c0 = AG.monomers[0].get_transition_environment((0,1)).data #get_egcf((0,1))
        Nt = len(c0)
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        
        ct = numpy.zeros((Nt),dtype=numpy.complex128)
        
        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        for el1 in elst:
            for el2 in elst:
                if cfm.cpointer[el1-1,el2-1] == 0:
                    continue
                coft = cfm.get_coft(el1-1,el2-1) 
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        ct += ((SS[kk,n]**2)*(SS[ll,n]**2)*coft)
            
        return ct

        
    def _calculate_monomer(self, raw=False):
        """ Calculates the absorption spectrum of a monomer 
        
        
        """
        ta = self.TimeAxis

        # loop over first band transitions
        for kk in range(self.system.Nb[1]):

            k_st = self.system.Nb[0] + kk
            # transition frequency (assuming just one ground state)
            om = self.system.elenergies[k_st]-self.system.elenergies[0]
            # transition dipole moment
            dm = self.system.dmoments[0,k_st,:]
            # dipole^2
            dd = numpy.dot(dm,dm)
            # natural life-time from the dipole moment
            gama = [-1.0/self.system.get_electronic_natural_lifetime(k_st)]
            
            if self.system._has_system_bath_coupling:
                # correlation function
                ct = self.system.get_transition_environment((0,k_st)).data #get_egcf((0,1))            
                tr = {"ta":ta,"dd":dd,"om":om-self.rwa,"ct":ct,"gg":gama}
            else:
                tr = {"ta":ta,"dd":dd,"om":om-self.rwa,"gg":gama}

            # calculates the one transition of the monomer        
            data_loc = numpy.real(self.one_transition_spectrum(tr))

            if kk == 0:
                data = data_loc
            else:
                data += data_loc
        

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

    
    def _calculate_abs_from_dynamics(self, raw=False, alt=False):
        """Calculates the absorption spectrum of a molecule from its dynamics
        
        """

        #
        # Frequency axis
        #
        # we only want to retain the upper half of the spectrum
        Nt = len(self.frequencyAxis.data)//2        
        do = self.frequencyAxis.data[1]-self.frequencyAxis.data[0]
        st = self.frequencyAxis.data[Nt//2]
        # we represent the Frequency axis anew
        axis = FrequencyAxis(st,Nt,do)
        
        # the spectrum will be calculate for this system
        system = self.system
        HH = system.get_Hamiltonian()
        if not HH.has_rwa:
            raise Exception("Hamiltonian has to define"+
                            " Rotating Wave Approximation")


        #with energy_units("1/cm"):
        #    print("Skeleton")
        #    print(numpy.diag(HH.get_RWA_skeleton()))
        #    print(HH.get_RWA_data())

        # transition dipole moment
        DD = system.get_TransitionDipoleMoment()
        
        rhoeq = system.get_thermal_ReducedDensityMatrix()

        # the propagator to propagate optical coherences
        prop = self.prop
        # FIXME: time axis must be consistent with other usage of the calculator
        time = prop.TimeAxis
        
        secular = False
                 
        if alt:
            
            at = _spect_from_dyn(time, HH, DD, prop, rhoeq, secular)
            
        else:
            
            at = _spect_from_dyn_single(time, HH, DD, prop, rhoeq, secular)


        #
        # Fourier transform of the time-dependent result
        #
        ft = numpy.fft.hfft(at)*time.step
        ft = numpy.fft.fftshift(ft)
        # invert the order because hfft is a transform with -i
        ft = numpy.flipud(ft)   
        # cut the center of the spectrum
        Nt = time.length #len(ta.data)        
        data = ft[Nt//2:Nt+Nt//2]
       
        #
        # multiply the response by \omega
        #
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
        ct = self._excitonic_coft(SS,self.system,1)
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
            tr["ct"] = self._excitonic_coft(SS,self.system,ii)
            
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


def _spect_from_dyn(time, HH, DD, prop, rhoeq, secular=False):
    """Calculation of the first order signal field.
    

    Parameters
    ----------
    time : TimeAxis
        Time axis on which the spectrum is calculated
    HH : Hamiltonian
        Hamiltonian of the system
    DD : TransitionDipoleMoment
        Transition dipole moment operator
    prop : ReducedDensityMatrixPropagator
        Propagator of the density matrix
    rhoeq : ReducedDensityMatrix
        Equilibrium reduced density matrix 
    secular : bool, optional
        Should secular approximation be used? The default is False.

    Returns
    -------
    at : numpy.array, COMPLEX
        Time dependent signal field.

    """
    # we will loop over transitions
    rhoi = ReducedDensityMatrix(dim=HH.dim)
    #
    # Time dependent data
    #
    at = numpy.zeros(time.length, dtype=COMPLEX)
    # transitions to loop over
    # ig represents all states in the ground state block
    for ig in range(HH.rwa_indices[1]):
        if len(HH.rwa_indices) <= 2:
            Nfin = HH.dim
        else:
            Nfin = HH.rwa_indices[2]
        for ie in range(HH.rwa_indices[1], Nfin):
            rhoi.data[:,:] = 0.0
            rhoi.data[ie,ig] = rhoeq.data[ig,ig]
            dvec_1 = DD.data[ie,ig,:]
            dd = numpy.dot(dvec_1,dvec_1)  
            
            rhot = prop.propagate(rhoi)
            
            if secular:
                at += (dd/3.0)*rhot.data[:,ie,ig]
            else:
                #
                # non-secular loops over all coherences
                #
                for jg in range(HH.rwa_indices[1]):
                    for je in range(HH.rwa_indices[1], Nfin):        
                        dvec_2 = DD.data[je,jg,:]
                        d12 = numpy.dot(dvec_2,dvec_1)
                        at += (d12/3.0)*rhot.data[:,je,jg]  
                        
        return at


def _spect_from_dyn_single(time, HH, DD, prop, rhoeq, secular=False):
    """Calculation of the first order signal field.
    
    Calculation in a single propagation

    Parameters
    ----------
    time : TimeAxis
        Time axis on which the spectrum is calculated
    HH : Hamiltonian
        Hamiltonian of the system
    DD : TransitionDipoleMoment
        Transition dipole moment operator
    prop : ReducedDensityMatrixPropagator
        Propagator of the density matrix
    rhoeq : ReducedDensityMatrix
        Equilibrium reduced density matrix 
    secular : bool, optional
        Should secular approximation be used? The default is False.

    Returns
    -------
    at : numpy.array, COMPLEX
        Time dependent signal field.

    """    
    rhoi = ReducedDensityMatrix(dim=HH.dim)
    #
    # Time dependent data
    #
    at = numpy.zeros(time.length, dtype=COMPLEX)
        
    #deff = numpy.sqrt(numpy.einsum("ijn,ijn->ij", DD.data,DD.data,
    #                               dtype=REAL))
 
    _show = False
    if _show:
        import matplotlib.pyplot as plt

    for kk in range(3):
        deff = DD.data[:,:,kk]
        # excitation by an effective dipole
        rhoi.data = (numpy.dot(deff,rhoeq.data)+numpy.dot(rhoeq.data,deff))/3.0

        rhot = prop.propagate(rhoi)

        
        if _show:
            for ig in range(HH.rwa_indices[1]):
                print("ig=", ig)
                for ll in range(rhot.data.shape[2]):
                    plt.plot(time.data, rhot.data[:,ig,ll])
                plt.title("kk="+str(kk))
                plt.show()


        for ig in range(HH.rwa_indices[1]):
            at += numpy.einsum("j,kj->k",deff[ig,:],rhot.data[:,:,ig])
        

    return at


def _c2g(timeaxis, coft):
    """ Converts correlation function to lineshape function
    
    Explicit numerical double integration of the correlation
    function to form a lineshape function.

    Parameters
    ----------

    timeaxis : TimeAxis
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
