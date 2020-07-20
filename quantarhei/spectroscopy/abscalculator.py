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
from ..core.managers import eigenbasis_of
from ..core.units import convert

from .linear_spectra import LinSpectrum

class LinSpectrumCalculator(EnergyUnitsManaged):
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
        self._gass_lineshape = False
        if relaxation_tensor is not None:
            self._relaxation_tensor = relaxation_tensor
            self._has_relaxation_tensor = True
        if effective_hamiltonian is not None:
            self._relaxation_hamiltonian = effective_hamiltonian
        if rate_matrix is not None:
            self._rate_matrix = rate_matrix
            self._has_rate_matrix = True
            
        self.rwa = 0.0

     
    def bootstrap(self,rwa=0.0, lab=None, HWHH = None, axis=None):
        """
        
        """
        self.rwa = self.convert_2_internal_u(rwa)
        with energy_units("int"):
            # sets the frequency axis for plottig
            self.frequencyAxis = self.TimeAxis.get_FrequencyAxis()
            self.frequencyAxis.data += self.rwa     
        
        if HWHH is not None:  
            self._gass_lineshape = True
            self.HWHH = self.convert_2_internal_u(HWHH)
        
        if axis is not None:
            if axis=="x":
                self.ld_axis = numpy.array([1.0,0.0,0.0],dtype="f8")   
            elif axis=="y":
                self.axis = numpy.array([0.0,1.0,0.0],dtype="f8")
            elif axis=="z":
                self.ld_axis = numpy.array([0.0,0.0,1.0],dtype="f8")
            else:
                self.ld_axis = axis/numpy.linalg.norm(axis)
        else:
            self.ld_axis = numpy.array([0.0,0.0,1.0],dtype="f8")
        
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
    
    def _equilibrium_excit_populations(self, AG, temperature=300,
                                 relaxation_hamiltonian=None):    
        if relaxation_hamiltonian:
            H = relaxation_hamiltonian
        else:
            H = AG.get_Hamiltonian()
        with eigenbasis_of(H):
            rho0 = AG.get_DensityMatrix(condition_type="thermal_excited_state",
                             relaxation_theory_limit="weak_coupling",
                             temperature=temperature,
                             relaxation_hamiltonian=relaxation_hamiltonian)
        return rho0
        
    def one_transition_spectrum_abs(self,tr):
        """ Calculates spectrum of one transition
        
        
        """
        

        ta = tr["ta"] # TimeAxis
        dd = tr["dd"] # transition dipole strength
        om = tr["om"] # frequency - rwa
        gg = tr["gg"] # natural broadening (constant or time dependent)
        
        # CD and fluorescence can be calculated in this step
        # TODO if rotatory strength defined calculate also circular dichroism spectra
        # TOOD calculate fluorescence spectra (for fluorescence there should be a switch because it should be calculated only for the first transition) 
        
        
        if self.system._has_system_bath_coupling:
#            ct = tr["ct"] # correlation function
        
            # convert correlation function to lineshape function
            #gt = self._c2g(ta,ct.data)
            gt = tr["gt"]
            # calculate time dependent response
            at = numpy.exp(-gt -1j*om*ta.data)
        else:
            # calculate time dependent response
            at = numpy.exp(-1j*om*ta.data) 
#        plt.figure()
#        plt.title("Absorption")
#        plt.plot(ta.data,numpy.real(at))
#        plt.plot(ta.data,numpy.imag(at))
    
                    
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
    
    def one_transition_spectrum_ld(self,tr):
        """ Calculates spectrum of one transition
        
        
        """
        

        ta = tr["ta"] # TimeAxis
        ld = tr["ld"] # linear dichroism strength
        om = tr["om"] # frequency - rwa
        gg = tr["gg"] # natural broadening (constant or time dependent)
        
        if self.system._has_system_bath_coupling:
#            ct = tr["ct"] # correlation function
        
            # convert correlation function to lineshape function
            #gt = self._c2g(ta,ct.data)
            gt = tr["gt"]
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
        ft = ld*numpy.fft.hfft(at)*ta.step
        ft = numpy.fft.fftshift(ft)
        # invert the order because hfft is a transform with -i
        ft = numpy.flipud(ft)   
        # cut the center of the spectrum
        Nt = ta.length #len(ta.data)        
        return ft[Nt//2:Nt+Nt//2]
        
    def one_transition_spectrum_gauss(self,tr):
        """ Calculates spectrum of one transition using gaussian broadening
        of the stick spectra. The definition is the same as for  exat tools: 
        https://doi.org/10.1002/jcc.25118
        
        """
        
        
        fa = tr["fa"]     # Frequency axis
        HWHH = tr["HWHH"] # Half width at the half hight (maximum)
        dd = tr["dd"]     # transition dipole strength
        rr = tr["rr"]     # transition dipole strength
        ld = tr["ld"]     # linear dichroism strength
        om = tr["om"]+self.rwa     # frequency
        
        # LineShape = lambda p, x: (x/(p[1]*np.sqrt(2*m.pi))*np.exp(-0.5*((x-p[0])/p[1])**2))
        # broad = broad/np.sqrt(2*np.log(2))
        sigma = HWHH/numpy.sqrt(2*numpy.log(2))
        
        # x = ta.data
        
        data = (fa.data/(sigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-0.5*((fa.data-om)/sigma)**2))
        data_abs = dd*data
        data_CD = rr*data
        data_LD = ld*data
        
        return data_abs,data_CD, data_LD
        
    
    def one_transition_spectrum_fluor(self,tr):
        """ Calculates spectrum of one transition
        
        
        """
        

        ta = tr["ta"] # TimeAxis
        dd = tr["dd"] # transition dipole strength
        om = tr["om"] # frequency - rwa
        gg = tr["gg"] # natural broadening (constant or time dependent)
        
        # CD and fluorescence can be calculated in this step
        # TODO if rotatory strength defined calculate also circular dichroism spectra
        # TOOD calculate fluorescence spectra (for fluorescence there should be a switch because it should be calculated only for the first transition) 
        
        
        if self.system._has_system_bath_coupling:
#            ct = tr["ct"] # correlation function
            re = tr["re"] # reorganisation energy
            
            # convert correlation function to lineshape function
            #gt = self._c2g(ta,ct.data)
            gt = tr["gt"]
            # calculate time dependent response
            at = numpy.exp(-numpy.conjugate(gt) -1j*om*ta.data + 2j*re*ta.data)
        else:
            # calculate time dependent response
            at = numpy.exp(-1j*om*ta.data) 
#        plt.figure()
#        plt.title("Absorption")
#        plt.plot(ta.data,numpy.real(at))
#        plt.plot(ta.data,numpy.imag(at))
    
                    
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
    
    def one_transition_spectrum_cd(self,tr):
        """ Calculates spectrum of one transition
        
        
        """
        

        ta = tr["ta"] # TimeAxis
        rr = tr["rr"] # transition dipole strength
        om = tr["om"] # frequency - rwa
        gg = tr["gg"] # natural broadening (constant or time dependent)
        
        # CD and fluorescence can be calculated in this step
        # TODO if rotatory strength defined calculate also circular dichroism spectra
        # TOOD calculate fluorescence spectra (for fluorescence there should be a switch because it should be calculated only for the first transition) 
        
        
        if self.system._has_system_bath_coupling:
#            ct = tr["ct"] # correlation function
        
            # convert correlation function to lineshape function
            #gt = self._c2g(ta,ct.data)
            gt = tr["gt"]
            # calculate time dependent response
            at = numpy.exp(-gt -1j*om*ta.data)
        else:
            # calculate time dependent response
            at = numpy.exp(-1j*om*ta.data) 
#        plt.figure()
#        plt.title("Absorption")
#        plt.plot(ta.data,numpy.real(at))
#        plt.plot(ta.data,numpy.imag(at))
    
                    
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
        ft = rr*numpy.fft.hfft(at)*ta.step
        ft = numpy.fft.fftshift(ft)
        # invert the order because hfft is a transform with -i
        ft = numpy.flipud(ft)   
        # cut the center of the spectrum
        Nt = ta.length #len(ta.data)        
        return ft[Nt//2:Nt+Nt//2]

        
    def _excitonic_coft_old(self,SS,AG,n):
        """ Returns energy gap correlation function data of an exciton state 
        
        """
        
        # FIXME: works only for 2 level molecules
        
        c0 = AG.monomers[0].get_egcf((0,1))
        Nt = len(c0)
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        
        # get number of monomeric basis states
        Na = 0
        for monomer in AG.monomers:
            Na += monomer.nel-1
        
        ct = numpy.zeros((Nt),dtype=numpy.complex128)
        #Na = AG.nmono
        for kk in range(Na):
            
            #nkk = AG.monomers[kk].egcf_mapping[0]
            
            for ll in range(Na):
            
                #nll = AG.monomers[ll].egcf_mapping[0]
                
                ct += ((SS[kk+1,n+1]**2)*(SS[ll+1,n+1]**2)*cfm.get_coft(kk,ll))
                #*AG.egcf_matrix.get_coft(nkk,nll))
            
        return ct
    
    def _excitonic_coft(self,SS,AG,n):
        """ Returns energy gap correlation function data of an exciton state n
        
        """
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
            
        c0 = AG.monomers[0].get_egcf((0,1))
        Nt = len(c0)
    
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
    
    def _excitonic_coft_all(self,SS,AG):
        """ Returns energy gap correlation function data of an exciton state n
        
        """
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
            
        c0 = AG.monomers[0].get_egcf((0,1))
        Nt = len(c0)
    
        Nst = AG.HamOp.dim
        ct = numpy.zeros((Nst,Nt),dtype=numpy.complex128)

        # electronic states corresponding to single excited states
        import time
        timecount = 0
        elst = numpy.where(AG.which_band == 1)[0]
        start = time.time()
        for el1 in elst:
            for el2 in elst:
                coft = cfm.get_coft(el1-1,el2-1)
                start2 = time.time()
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        ct[:,:] += numpy.dot(
                         numpy.expand_dims((SS[kk,:]**2)*(SS[ll,:]**2),axis=1),
                         numpy.expand_dims(coft,axis=0))
                stop2 = time.time()
                timecount += stop2 - start2
        stop = time.time()
        print(stop-start,stop-start - timecount)
        return ct

    def _excitonic_reorg_energy(self, SS, AG, n):
        """ Returns the reorganisation energy of an exciton state
        """
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        
        rg = 0.0
        
        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        for el1 in elst:
            reorg = cfm.get_reorganization_energy(el1-1,el1-1)
            for kk in AG.vibindices[el1]:
                rg += ((SS[kk,n]**2)*(SS[kk,n]**2)*reorg)
        return rg    
    
    
    def _excitonic_rotatory_strength(self,SS,AG,energy,n):
        # Initialize rotatory strength
        Rot_n = 0
        
        for ii in range(AG.Ntot):
            for jj in range(ii+1,AG.Ntot):
                Rot_n += SS[ii,n]*SS[jj,n]*AG.RR[ii,jj]
#                print(ii,jj,SS[ii,n],SS[jj,n],AG.RR[ii,jj])
#                
#        
#        # electronic states corresponding to single excited states
#        elst = numpy.where(AG.which_band == 1)[0]
#        for el1 in elst:
#            # get monomer number
#            mon_indx = numpy.nonzero(AG.elsigs[el1])[0][0]
#            mon1 = AG.monomers[mon_indx]
#            Ri = mon1.position
#            for el2 in elst:
#                if el2>el1:
#                    mon_indx = numpy.nonzero(AG.elsigs[el2])[0][0]
#                    mon2 = AG.monomers[mon_indx]
#                    Rj = mon2.position
#                    
#                    for ii in AG.vibindices[el1]:
#                        di = AG.vibdipoles[0,ii]
#                        for jj in AG.vibindices[el2]:
#                            dj = AG.vibdipoles[0,jj]
        
        # Scale by excitation energy:
        Rot_n *= energy[n]
#        print(Rot_n,energy,n)
        return Rot_n           

    def _excitonic_rotatory_strength_fullv(self,SS,AG,energy,n):
        # Initialize rotatory strength
        Rot_n = 0
#        Rot_nm = 0
#        
#        DD_vel = self.system.DD[0].copy()
#        for ii in range(1:DD_vel.shape[0]):
#            DD_vel[ii] *= -self.system.HH[ii,ii]
#        
        for ii in range(AG.Ntot):
            for jj in range(AG.Ntot):
                Rot_n += SS[ii,n]*SS[jj,n]*(AG.RRv[ii,jj]+AG.RRm[ii,jj])
#        
#        for ii in range(AG.Ntot):
#            for jj in range(AG.Ntot):
#                Rot_nm += SS[ii,n]*SS[jj,n]*AG.RRm[ii,jj]
        
        return Rot_n#,Rot_nm
        
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
            gt = self._c2g(ta,ct.data)
            tr = {"ta":ta,"dd":dd,"om":om-self.rwa,"ct":ct,"gt":gt,"gg":gama}
        else:
            tr = {"ta":ta,"dd":dd,"om":om-self.rwa,"gg":gama}

        # calculates the one transition of the monomer        
        data = numpy.real(self.one_transition_spectrum_abs(tr))
        

        # we only want to retain the upper half of the spectrum
        Nt = len(self.frequencyAxis.data)//2        
        do = self.frequencyAxis.data[1]-self.frequencyAxis.data[0]
        st = self.frequencyAxis.data[Nt//2]
        # we represent the Frequency axis anew
        axis = FrequencyAxis(st,Nt,do)

        # multiply the spectrum by frequency (compulsory prefactor)
        if not raw:
            data = axis.data*data

        
        spect = LinSpectrum(axis=axis, data=data)
        
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
        energy = numpy.diag(HH.data)
        
        # Transition dipole moment operator
        DD = self.system.get_TransitionDipoleMoment()
        # transformed into the basis of Hamiltonian eigenstates
        DD.transform(SS)         

        # Frequency axis
        # we only want to retain the upper half of the spectrum
        Nt = len(self.frequencyAxis.data)//2        
        do = self.frequencyAxis.data[1]-self.frequencyAxis.data[0]
        st = self.frequencyAxis.data[Nt//2]
        # we represent the Frequency axis anew
        axis = FrequencyAxis(st,Nt,do) 
        
        # TimeAxis
        tr = {"ta":ta, "fa":axis}
        
        # If gaussian groadening is specified calculate also spectra
        if self._gass_lineshape:
            tr["HWHH"] = self.HWHH
        
        
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
                    gg.append(RR.data[:,ii,ii]/2.0)
            else:
                for ii in range(HH.dim):
                    gg.append([RR.data[ii,ii]/2.0])
            tr["gg"] = gg[1]
        else:
            
            tr["gg"] = [0.0]
        
        # get square of transition dipole moment here
        tr["dd"] = DD.dipole_strength(0,1)
        # first transition energy
        tr["om"] = HH.data[1,1]-HH.data[0,0]-self.rwa
        # get a transformed ct here
        ct = self._excitonic_coft(SS,self.system,1)
        tr["ct"] = ct
        # calculate g(t)   
        tr["gt"] = self._c2g(tr["ta"],tr["ct"].data)
        # get reorganization energy
        tr["re"] = self._excitonic_reorg_energy(SS,self.system,1)
        # get rotatory strength
        tr["rr"] = self._excitonic_rotatory_strength_fullv(SS,self.system,energy,1)
        #tr["rr"] = self._excitonic_rotatory_strength(SS,self.system,energy,1)
#        print(1,convert(tr["rr"],"int","1/cm")*numpy.pi*1e-4)
        dip = DD.data[0,1]
        tr["ld"] = (3*(numpy.dot(dip,self.ld_axis))**2 - tr["dd"]) # *3/2
        #tr["ld"] = 3*(DD.  # *3/2
        

        self.system._has_system_bath_coupling = True
        
        temperature = self.system.sbi.get_temperature()
        rho_eq_exct = self._equilibrium_excit_populations(self.system,
                                               temperature=temperature)
        
#        print(tr["ct"])
#        print(max(tr["ct"]))
#        print(self.one_transition_spectrum(tr))
#        print(max(self.one_transition_spectrum(tr)))
        #
        # Calculates spectrum of a single transition
        #
        data = numpy.real(self.one_transition_spectrum_abs(tr))
        data_fl = numpy.real(rho_eq_exct.data[1, 1]*self.one_transition_spectrum_fluor(tr))
        data_cd = numpy.real(self.one_transition_spectrum_cd(tr))
        data_ld = numpy.real(self.one_transition_spectrum_ld(tr))
        if self._gass_lineshape:
            data_gauss,data_cd_gauss,data_ld_gauss = numpy.real(self.one_transition_spectrum_gauss(tr))
#        print("Population mine:",rho_eq_exct.data[1, 1])

        # FOR THE VIBRONIC SYSTEM THE SPECTRA HAVE TO BE SUMED THROUGH THE GROUND STATES (VIBRONIC)
        for jj in range(1,min(1,self.system.Nb[0])): # sum over the ground states
            tr["dd"] = DD.dipole_strength(jj,1)
            # first transition energy
            tr["om"] = HH.data[1,1]-HH.data[jj,jj]-self.rwa
            data_fl += numpy.real(rho_eq_exct.data[1, 1]*self.one_transition_spectrum_fluor(tr))
            
        
        for ii in range(2,HH.dim):
            if relaxation_tensor is not None or rate_matrix is not None:
                tr["gg"] = gg[ii]
            else:
                tr["gg"] = [0.0]
#            print(tr["gg"])
            #tr[1] = DD.dipole_strength(0,ii) # update transition dipole moment
            tr["dd"] = DD.dipole_strength(0,ii)
            #tr[2] = HH.data[ii,ii]-HH.data[0,0]-rwa
            tr["om"] = HH.data[ii,ii]-HH.data[0,0]-self.rwa
            #tr[3] = self._excitonic_coft(SS,self.system,ii-1) # update ct here
            tr["ct"] = self._excitonic_coft(SS,self.system,ii)
            tr["gt"] = self._c2g(tr["ta"],tr["ct"].data)
            tr["re"] = self._excitonic_reorg_energy(SS,self.system,ii)
            tr["rr"] = self._excitonic_rotatory_strength_fullv(SS,self.system,energy,ii)
            dip = DD.data[0,ii]
            tr["ld"] = (3*(numpy.dot(dip,self.ld_axis))**2 - tr["dd"]) # *3/2
            #tr["rr"] = self._excitonic_rotatory_strength(SS,self.system,energy,ii)
#            print(ii,convert(HH.data[ii,ii]-HH.data[0,0]-tr["re"],"int","1/cm"),convert(HH.data[ii,ii]-HH.data[0,0]-2*tr["re"],"int","1/cm"),convert(tr["re"],"int","1/cm"))
#            print(ii,convert(tr["rr"],"int","1/cm")*numpy.pi*1e-4)
            # conversion factor is convert rotatory strength to inverse centimeters and multiply *numpy.pi*1e-4


            
            #
            # Calculates spectrum of a single transition
            #
            data += numpy.real(self.one_transition_spectrum_abs(tr))
            data_fl += numpy.real(rho_eq_exct.data[ii, ii]*self.one_transition_spectrum_fluor(tr))
            data_cd += numpy.real(self.one_transition_spectrum_cd(tr))
            data_ld += numpy.real(self.one_transition_spectrum_ld(tr))
            if self._gass_lineshape:
                data_gauss_tmp,data_cd_gauss_tmp, data_ld_gauss_tmp = numpy.real(self.one_transition_spectrum_gauss(tr))
                data_gauss += data_gauss_tmp
                data_cd_gauss += data_cd_gauss_tmp
                data_ld_gauss += data_ld_gauss_tmp
#            print("Population mine:",rho_eq_exct.data[ii, ii])

            # FOR THE VIBRONIC SYSTEM THE SPECTRA HAVE TO BE SUMED THROUGH THE GROUND STATES (VIBRONIC)
            for jj in range(1,min(ii,self.system.Nb[0])): # sum over the ground states
                tr["dd"] = DD.dipole_strength(jj,ii)
                # first transition energy
                tr["om"] = HH.data[ii,ii]-HH.data[jj,jj]-self.rwa
                data_fl += numpy.real(rho_eq_exct.data[ii, ii]*self.one_transition_spectrum_fluor(tr))

        
        
        # multiply the spectrum by frequency (compulsory prefactor)
        if not raw:
            data = axis.data*data
            data_fl = axis.data*data_fl
            data_cd =  axis.data*data_cd
            data_ld = axis.data*data_ld
            #data_gauss = axis.data*data_gauss
            #data_cd_gauss = axis.data*data_cd_gauss
        
        # transform all quantities back
        S1 = numpy.linalg.inv(SS)
        HH.transform(S1)
        DD.transform(S1)
        
        if relaxation_tensor is not None:
            RR.transform(S1)

        abs_spect = LinSpectrum(axis=axis, data=data)
        fluor_spect = LinSpectrum(axis=axis, data=data_fl)
        CD_spect = LinSpectrum(axis=axis, data=data_cd)
        LD_spect = LinSpectrum(axis=axis, data=data_ld)
        if self._gass_lineshape:
            abs_spect_gauss = LinSpectrum(axis=axis, data=data_gauss)
            CD_spect_gauss = LinSpectrum(axis=axis, data=data_cd_gauss)
            LD_spect_gauss = LinSpectrum(axis=axis, data=data_ld_gauss)
        
            return {"abs": abs_spect, "fluor": fluor_spect, "CD":  CD_spect,
                    "LD": LD_spect, "LD gauss": LD_spect_gauss,
                    "abs gauss": abs_spect_gauss, "CD gauss": CD_spect_gauss} 
        else:
            return {"abs": abs_spect, "fluor": fluor_spect, "CD":  CD_spect,
                    "LD":  LD_spect}    

                   
class AbsSpectrumCalculator(LinSpectrumCalculator):
    def __init__(self, timeaxis,
                 system=None,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None,
                 effective_hamiltonian=None):

        super().__init__(timeaxis,
                 system=system,
                 dynamics=dynamics,
                 relaxation_tensor=relaxation_tensor,
                 rate_matrix=rate_matrix,
                 effective_hamiltonian=effective_hamiltonian)
    
    def calculate(self, raw=False):
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
                                              raw=raw)["abs"]
            else:
                raise Exception("System to calculate spectrum for not defined")
        
        return spect