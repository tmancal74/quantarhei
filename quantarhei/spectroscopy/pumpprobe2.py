# -*- coding: utf-8 -*-

import numpy
import scipy

from .twod2 import TwoDSpectrum
from .twodcontainer import TwoDSpectrumContainer
from .twodcalculator import TwoDResponseCalculator as TwoDSpectrumCalculator
from .mocktwodcalculator import MockTwoDResponseCalculator as MockTwoDSpectrumCalculator
from ..core.dfunction import DFunction
from ..core.managers import Manager
from ..core.managers import energy_units
from quantarhei import convert
from .. import COMPLEX

from ..utils import derived_type
from ..builders.aggregates import Aggregate
from ..builders.molecules import Molecule
from ..core.time import TimeAxis
from ..core.frequency import FrequencyAxis


import matplotlib.pyplot as plt


class PumpProbeSpectrum(DFunction):
    """Class representing pump-probe spectra
    
    
    
    """
    
    #data = None
    
    def __init__(self):
        self.t2 = -1.0
        self.data = None
        
    def set_axis(self, axis):
        self.xaxis = axis
        self.axis = self.xaxis
        
    def set_data(self, data):
        self.data = data

    def get_PumpProbeSpectrum(self):
        """Returns self
        
        This method is here to override the one inherented from TwoDSpectrum
        
        """
        return self


    def set_t2(self, t2):
        """Sets the t2 (waiting time) of the spectrum
        
        
        """
        self.t2 = t2


    def get_t2(self):
        """Returns the t2 (waiting time) of the spectrum
        
        """
        return self.t2
    
    
    def _add_data(self,data):
        if self.data is None:
            self.set_data(data)
        else:
            if self.data.size != len(data):
                raise IOError("Added data length does not match the current" +
                              "one")
            self.data += data
    
    #FIXME: Add function _add_data (if data None = set_data, else add)



class PumpProbeSpectrumContainer(TwoDSpectrumContainer):
    """Container for a set of pump-probe spectra
    
    """
    def __init__(self, t2axis=None):
        
        self.t2axis = t2axis
        self.spectra = {}
        
    def plot(self):
        
        plt.clf()
        spctr = self.get_spectra()
        for sp in spctr:
            plt.plot(sp.xaxis.data,sp.data)
            
    def set_spectrum(self, spec, tag):
        self.spectra[tag] = spec
        
    def amax(self):
        mxs = []
        for s in self.get_spectra():       
            spect = numpy.real(s.data)   
            mx = numpy.amax(spect)
            mxs.append(mx)
        return numpy.amax(numpy.array(mxs))   
     
    def amin(self):
        mxs = []
        for s in self.get_spectra():       
            spect = numpy.real(s.data)   
            mx = numpy.amin(spect)
            mxs.append(mx)
        return numpy.amin(numpy.array(mxs))   


    def plot2D(self, axis = None, units = "nm", zero_centered = True, lines = None):
        
        t2ax = self.t2axis
        freqax = self.spectra[t2ax.data[0]].axis
        # use only positive frequency axis
        min_ind = numpy.min(numpy.where(freqax.data > 0.0))
        with energy_units(units):
            # Prepare array of pump-probe spectra
            ppspec2D  = numpy.zeros((t2ax.length,freqax.length))
            count = 0
            for T2 in t2ax.data:
                ppspec2D[count] += self.spectra[T2].data
                count += 1
        
            # plot pump-probe spectra
            X,Y =  numpy.meshgrid(t2ax.data,freqax.data[min_ind:])
            fig, ax = plt.subplots(figsize=(18,9))
            
            if zero_centered:
                p = ax.contourf(X, Y, ppspec2D.T[min_ind:], 60, cmap=plt.cm.jet, 
                                vmin=-abs(ppspec2D).max(), vmax=abs(ppspec2D).max()) #, vmin=ppspec2D.min(), vmax=ppspec2D.max())  
            else:
                p = ax.contourf(X, Y, ppspec2D.T[min_ind:], 60, cmap=plt.cm.jet, 
                                vmin=ppspec2D.min(), vmax=ppspec2D.max())
            
            if isinstance(lines,list):
                for line_freq in lines:
                    plt.plot(t2ax.data,numpy.ones(t2ax.length)*line_freq,"w",linewidth=2)
            
            if axis is not None:
                ax.axes.set_ylim(axis[1][0],axis[1][1])   
                ax.axes.set_xlim(axis[0][0],axis[0][1])
            else:
                ax.axes.set_ylim(400,750)   
                ax.axes.set_xlim(0,500)
            plt.xlabel("Time [fs]")
            plt.ylabel("Wavelength [" + units + "]")
            
            ax.legend()
            plt.colorbar(p,ax=ax)
            fig.savefig('PP_2D_spectra.png', format='png', dpi=1200)
            
    def plot_slices(self,freqs,expRes=None,units="nm"):
        
        # Initialize frequency cuts of PP spectra
        spectra_freq = numpy.zeros((len(freqs),self.t2axis.length),dtype="f8")
        count = 0
        for T2 in self.t2axis.data:
            sp = self.spectra[T2]
            sp._has_imag = False
            sp._set_splines()
            freq_int = convert(freqs,units,"int")
            spectra_freq[:,count] = self.spectra[T2].at(freq_int)
            count+=1
        
        with energy_units(units):
            fig = plt.figure(figsize=(18,9))
            Nsp = len(freqs)
            plt_num = numpy.arange(Nsp) + 1
            plt_num = plt_num.reshape((Nsp//3,3)).T.reshape(Nsp)
            plt.xlabel("Delay time [ps]")
            plt.ylabel("Intensity [arb. u.]")
            for ii in range(Nsp):
                plt.subplot(Nsp//3,3,plt_num[ii])
                plt.title(str(numpy.round(freqs[ii]))+" "+units)
                plt.plot(self.t2axis.data/1000,spectra_freq[ii])
                if expRes is not None:
                    PP_exp_freq = expRes["frequency"]
                    PP_exp_spec = expRes["intensity"]
                    plt.plot(PP_exp_freq[ii],PP_exp_spec[ii])
                    
                plt.xscale('symlog', linthreshx= (0.1))
            plt.subplots_adjust(hspace=0.6, wspace=0.6)
            plt.show()
            fig.savefig('PP_slices_spectra.png', format='png', dpi=1200)
        return spectra_freq
            
        
    def make_movie(self, filename, axis=None,
                   cmap=None, vmax=None, vmin=None,
                   frate=20, dpi=100, start=None, end=None,
                   show_states=None, progressbar=False):
        
        import matplotlib.pyplot as plt
        import matplotlib.animation as manimation
        
        FFMpegWriter = manimation.writers["ffmpeg"]

        metadata = dict(title="Test Movie", artist='Matplotlib',
                comment='Movie support!')
        writer = FFMpegWriter(fps=frate, metadata=metadata)
        
        fig = plt.figure() 
        
        spctr = self.get_spectra()
        l = len(spctr)
        last_t2 = spctr[l-1].get_t2()
        first_t2 = spctr[0].get_t2()
        
        if vmax is None:
            mx = self.amax()
        else:
            mx = vmax
            
        if vmin is None:
            mn = self.amin()
        else:
            mn = vmin
            
        mxabs = max(numpy.abs(mx), numpy.abs(mn))
        mx = mx+0.05*mxabs
        mn = mn-0.05*mxabs
        
        if start is None:
            start = first_t2
        if end is None:
            end = last_t2
        
                
        with writer.saving(fig, filename, dpi):  
            k = 0
            # Initial call to print 0% progress
            sp2write = self.get_spectra(start=start, end=end)
            l = len(sp2write)
            if progressbar:
                self._printProgressBar(0, l, prefix = 'Progress:',
                                       suffix = 'Complete', length = 50)
            for sp in self.get_spectra(start=start, end=end):
                # FIXME: this does not work as it should yet
                sp.plot(show=False, fig=fig, axis=axis,
                        label="T="+str(sp.get_t2())+"fs") #, vmax=mx, vmin=mn,
                          #)
                writer.grab_frame()
                if progressbar:
                    self._printProgressBar(k + 1, l, prefix = 'Progress:',
                                           suffix = 'Complete', length = 50)
                
                k += 1


class PumpProbeSpectrumCalculator():
    
    t2axis = derived_type("t2axis",TimeAxis)
    t3axis = derived_type("t3axis",TimeAxis)
    
    system = derived_type("system",[Molecule,Aggregate])
    
    def __init__(self, t2axis, t3axis,
                 system=None,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None,
                 effective_hamiltonian=None):
            
            
        self.t2axis = t2axis
        self.t3axis = t3axis
        
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
        self.goft_matrix = None
        
        # after bootstrap information
        self.lab = None
        self.t3s = None
        self.pathways = None
        self.rwa = None
        
        self.verbose = False
        
        self.tc = 0
    
    def bootstrap(self, rwa=0.0, pathways=None, lab=None, verbose=False):
        """Sets up the environment for pump-probe calculation
        
        """


        self.verbose = verbose
            
        if isinstance(self.system, Aggregate):
        
            pass
        
        else:
            
            raise Exception("Molecule pump-probe not implememted")
            
        
        self.verbose = verbose
        self.rwa = Manager().convert_energy_2_internal_u(rwa)
        self.pathways = pathways
        
        with energy_units("int"):
            #atype = self.t3axis.atype
            #self.t3axis.atype = 'complete'
            #self.oa3 = self.t3axis.get_FrequencyAxis() 
            #self.oa3.data += self.rwa
            #self.oa3.start += self.rwa
            #self.t3axis.atype = atype  
            
            # we only want to retain the upper half of the spectrum
            freq = self.t3axis.get_FrequencyAxis()
            freq.data += self.rwa 
            Nt = len(freq.data)//2        
            do = freq.data[1]-freq.data[0]
            st = freq.data[Nt//2]
            # we represent the Frequency axis anew
            self.oa3 = FrequencyAxis(st,Nt,do)      
        
        self.tc = 0
        self.lab = lab
        
    
    def set_pathways(self, pathways):
        self.pathways = pathways
        
    def calculate_all_system(self, sys, eUt, lab, show_progress=False):
        """Calculates all 2D spectra for a system and evolution superoperator
        
        """
        tcont = PumpProbeSpectrumContainer(t2axis=self.t2axis)
        
        kk = 1
        Nk = self.t2axis.length
        
        printProgressBar(0,Nk,prefix="     - Progress:",suffix= "Complete", length = 50)
        for T2 in self.t2axis.data:
            
            if show_progress:
                print(" - calculating", kk, "of", Nk, "at t2 =", T2, "fs")
            
            ppspec1 = self.calculate_one_system(T2, sys, eUt, lab)
            
            tcont.set_spectrum(ppspec1, tag=T2)
            
            printProgressBar(kk,Nk,prefix="     - Progress:",suffix= "Complete", length = 50)
            kk += 1
        
        return tcont
    
    def calculate_one_system(self, t2, sys, eUt, lab, pways=None):
        """Returns pump-probe spectrum at t2 for a system and evolution 
        superoperator
        
        """
        try:
            #print(t2)
            Uin = eUt.at(t2)
            #print(Uin.data.shape)
        except:
            Uin = eUt
            #print("False")
            
        temp = sys.get_temperature()
        # FIXME: Delete zero temperature
        #temp = 0.0
        rho0 = sys.get_DensityMatrix(condition_type="thermal",
                                     temperature=temp)
        
        # if the Hamiltonian is larger than eUt, we will calculate ESA
        has_ESA = True
        H = self.system.get_Hamiltonian()
        
        # get Liouville pathways
        if has_ESA:
            pws = sys.liouville_pathways_3T(ptype=("R1g", "R2g", "R3g",
                                                   "R4g", "R1f*", "R2f*",
                                                   "R1f*E","R1f*E"), 
                                                   #"R1gE", "R2gE"),
                                                   eUt=eUt, ham=H, t2=t2,
                                                   lab=lab)
        else:
            pws = sys.liouville_pathways_3T(ptype=("R1g", "R2g", "R3g",
                                                   "R4g", "R1f*E","R1f*E"), 
                                                   eUt=eUt, ham=H, t2=t2,
                                                   lab=lab)

        self.set_pathways(pws)
        
        SS  = self.system.SS.copy()
        self.goft_matrix = self._SE_excitonic_gofts(SS,self.system, tau = 0.0)
        
        if pways is not None:
            pways[str(t2)] = pws
            
        pprobe1 = self.calculate_next(t2)
        
        return pprobe1
    
    def calculate_next(self,t2):

        sone = self.calculate_one(self.tc,t2)
        #print(self.tc, sone)
        self.tc += 1
        return sone
    
    def calculate_one(self, tc, t2):
        """Calculate the 2D spectrum for all pathways
        
        """
        #import time
        
        onepp = PumpProbeSpectrum()
        onepp.set_axis(self.oa3)

        #start = time.time()
        data = self.calculate_pathways(self.pathways,t2)
        #end = time.time()
        #print("Single spectra calculation:", end - start)
        onepp._add_data(data)

        onepp.set_t2(self.t2axis.data[tc])    
        
        return onepp
    
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

    def _excitonic_coft(self,SS,AG,n):
        """ Returns energy gap correlation function data of an exciton state n
        
        """
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        
        # get number of monomeric basis states
        Na = 0
        for monomer in AG.monomers:
            Na += monomer.nel-1
    
        ct = numpy.zeros((self.t3axis.length),dtype=numpy.complex128)

        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        for el1 in elst:
            for el2 in elst:
                coft = DFunction(cfm.timeAxis,cfm.get_coft(el1,el2))
                ct3 = coft.at(self.t3axis.data)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        ct += ((SS[kk,n]**2)*(SS[ll,n]**2)*ct3)
#        coft = DFunction(cfm.timeAxis,cfm.get_coft(2,2))
#        ct3 = coft.at(self.t3axis.data)
#        print(numpy.isclose(ct,ct3).all())
        return ct
    
    def _SE_excitonic_cofts(self,SS,AG,tau = 0):
        """ Returns energy gap correlation function data of an exciton state n
        
        """
        
        c0 = AG.monomers[0].get_egcf((0,1))
        Nt = len(c0)
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        ctimeAxis = cfm.timeAxis
        
        if self.t3axis.max + tau > ctimeAxis.max:
            raise IOError("Correlation function should be defined on interval"+ 
                          " (0, t2_max+t3_max).")
        
        # get number of monomeric basis states
        Na = 0
        for monomer in AG.monomers:
            Na += monomer.nel-1
        
        ct3 = numpy.zeros((AG.Ntot,AG.Ntot,self.t3axis.length),dtype=numpy.complex128)
        ct3tau = numpy.zeros((AG.Ntot,AG.Ntot,self.t3axis.length),dtype=numpy.complex128)
        
        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        #print(elst)
        for el1 in elst:
            for el2 in elst:
                # get_coft starts from the excited state (ground not included in indexes)
                coft = DFunction(cfm.timeAxis,cfm.get_coft(el1-1,el2-1))
                coft_t3 = coft.at(self.t3axis.data)
                coft_t3tau = coft.at(self.t3axis.data+tau)
                # FIXME: both at time t3 and t3+tau
                
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
                            for m in range(n, AG.Nb[0]+AG.Nb[1]):
                                ct3[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft_t3)
                                ct3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft_t3tau)
        for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
            for m in range(n+1, AG.Nb[0]+AG.Nb[1]):
                ct3[m,n] += ct3[n,m]
                ct3tau[m,n] += ct3tau[n,m]
        
        # electronic states corresponding to double excited states
        elstd = numpy.where(AG.which_band == 2)[0]
        for el1 in elstd:
            for el2 in elstd:
                #print(el1,el2)
                coft = DFunction(cfm.timeAxis,cfm.get_coft(el1-1,el2-1))
                coft_t3 = coft.at(self.t3axis.data)
                coft_t3tau = coft.at(self.t3axis.data+tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                            for m in range(n, numpy.sum(AG.Nb[0:3])):
                                ct3[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft_t3)
                                ct3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft_t3tau)
        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
            for m in range(n+1, numpy.sum(AG.Nb[0:3])):
                ct3[m,n] += ct3[n,m]
                ct3tau[m,n] += ct3tau[n,m]
        
        # Mixed sigle-double excited state correlation functions
        for el1 in elst:
            for el2 in elstd:
                equal = numpy.where(AG.elsigs[el1] == AG.elsigs[el2])[0]
                if equal.size == 1:
                    coft = DFunction(cfm.timeAxis,cfm.get_coft(el1-1,el1-1))
                    coft_t3 = coft.at(self.t3axis.data)
                    coft_t3tau = coft.at(self.t3axis.data+tau)
                    for kk in AG.vibindices[el1]:
                        for ll in AG.vibindices[el2]:
                            for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
                                for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                                    ct3[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft_t3)
                                    ct3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft_t3tau)
        for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
            for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                ct3[m,n] += ct3[n,m]
                ct3tau[m,n] += ct3tau[n,m]
            
        return ct3,ct3tau
    
    
    def _SE_excitonic_cofts_test(self,SS,AG,tau = 0):
        """ Returns energy gap correlation function data of an exciton state n
        
        """
        
        c0 = AG.monomers[0].get_egcf((0,1))
        Nt = len(c0)
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        ctimeAxis = cfm.timeAxis
        
        if self.t3axis.max + tau > ctimeAxis.max:
            raise IOError("Correlation function should be defined on interval"+ 
                          " (0, t2_max+t3_max).")
        
        # get number of monomeric basis states
        Na = 0
        for monomer in AG.monomers:
            Na += monomer.nel-1
        
        ct = numpy.zeros((AG.Ntot,AG.Ntot,cfm.timeAxis.length),dtype=numpy.complex128)
        gt3 = numpy.zeros((AG.Ntot,AG.Ntot,self.t3axis.length),dtype=numpy.complex128)
        gt3tau = numpy.zeros((AG.Ntot,AG.Ntot,self.t3axis.length),dtype=numpy.complex128)
        
        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        #print(elst)
        for el1 in elst:
            for el2 in elst:
                # get_coft starts from the excited state (ground not included in indexes)
                coft = cfm.get_coft(el1-1,el2-1)
                goft = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,coft))
                goft_t3 = goft.at(self.t3axis.data)
                goft_t3tau = goft.at(self.t3axis.data + tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
                            for m in range(n, AG.Nb[0]+AG.Nb[1]):
                                #ct[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft)
                                gt3[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3)
                                gt3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3tau)
        for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
            for m in range(n+1, AG.Nb[0]+AG.Nb[1]):
                #ct[m,n] += ct[n,m]
                gt3[m,n] += gt3[n,m]
                gt3tau[m,n] += gt3tau[n,m]
        
        # electronic states corresponding to double excited states
        elstd = numpy.where(AG.which_band == 2)[0]
        for el1 in elstd:
            for el2 in elstd:
                coft = cfm.get_coft(el1-1,el2-1) # DFunction(cfm.timeAxis,)
                goft = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,coft))
                goft_t3 = goft.at(self.t3axis.data)
                goft_t3tau = goft.at(self.t3axis.data + tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                            for m in range(n, numpy.sum(AG.Nb[0:3])):
                                #ct[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft)
                                gt3[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3)
                                gt3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3tau)
        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
            for m in range(n+1, numpy.sum(AG.Nb[0:3])):
                #ct[m,n] += ct[n,m]
                gt3[m,n] += gt3[n,m]
                gt3tau[m,n] += gt3tau[n,m]
        
        # Mixed sigle-double excited state correlation functions
        for el1 in elst:
            for el2 in elstd:
                equal = numpy.where(AG.elsigs[el1] == AG.elsigs[el2])[0]
                if equal.size == 1:
                    coft = cfm.get_coft(el1-1,el1-1) # DFunction(cfm.timeAxis,)
                    goft = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,coft))
                    goft_t3 = goft.at(self.t3axis.data)
                    goft_t3tau = goft.at(self.t3axis.data + tau)
                    for kk in AG.vibindices[el1]:
                        for ll in AG.vibindices[el2]:
                            for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
                                for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                                    #ct[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft)
                                    gt3[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3)
                                    gt3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3tau)
        for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
            for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                #ct[m,n] += ct[n,m]
                gt3[m,n] += gt3[n,m]
                gt3tau[m,n] += gt3tau[n,m]
                
        if 0:        
            #check if gts correctly defined
            for ii in range(AG.Nb[0],AG.Ntot):
                for jj in range(AG.Nb[0],AG.Ntot):
                    goft = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,ct[ii,jj]))
                    if not numpy.isclose(gt3[ii,jj],goft.at(self.t3axis.data)).all():
                        print(ii,jj,False)
                        print(gt3[ii,jj] - goft.at(self.t3axis.data))
            
        return gt3,gt3tau
    
    def _SE_excitonic_gofts(self,SS,AG,tau = 0):
        """ Returns energy gap correlation function data of an exciton state n
        
        """

        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        ctimeAxis = cfm.timeAxis
        
        if self.t3axis.max + tau > ctimeAxis.max:
            raise IOError("Correlation function should be defined on interval"+ 
                          " (0, t2_max+t3_max).")
        
        # get number of monomeric basis states
        Na = 0
        for monomer in AG.monomers:
            Na += monomer.nel-1
        
        gt3tau = numpy.zeros((AG.Ntot,AG.Ntot,self.t3axis.length),dtype=numpy.complex128)
        
        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        #print(elst)
        for el1 in elst:
            for el2 in elst:
                # get_coft starts from the excited state (ground not included in indexes)
                coft = cfm.get_coft(el1-1,el2-1)
                goft = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,coft))
                goft_t3tau = goft.at(self.t3axis.data + tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
                            for m in range(n, AG.Nb[0]+AG.Nb[1]):
                                gt3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3tau)
        for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
            for m in range(n+1, AG.Nb[0]+AG.Nb[1]):
                gt3tau[m,n] += gt3tau[n,m]
        
        # electronic states corresponding to double excited states
        elstd = numpy.where(AG.which_band == 2)[0]
        for el1 in elstd:
            for el2 in elstd:
                coft = cfm.get_coft(el1-1,el2-1) # DFunction(cfm.timeAxis,)
                goft = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,coft))
                goft_t3tau = goft.at(self.t3axis.data + tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                            for m in range(n, numpy.sum(AG.Nb[0:3])):
                                gt3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3tau)
        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
            for m in range(n+1, numpy.sum(AG.Nb[0:3])):
                gt3tau[m,n] += gt3tau[n,m]
        
        # Mixed sigle-double excited state correlation functions
        for el1 in elst:
            for el2 in elstd:
                equal = numpy.where(AG.elsigs[el1] == AG.elsigs[el2])[0]
                if equal.size == 1:
                    coft = cfm.get_coft(el1-1,el1-1) # DFunction(cfm.timeAxis,)
                    goft = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,coft))
                    goft_t3tau = goft.at(self.t3axis.data + tau)
                    for kk in AG.vibindices[el1]:
                        for ll in AG.vibindices[el2]:
                            for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
                                for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                                    gt3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3tau)
        for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
            for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                gt3tau[m,n] += gt3tau[n,m]

            
        return gt3tau
    
    def calculate_pathways(self, pathways, tau):
        """Calculate the shape of a Liouville pathway
        
        """
        import time
        
        # we can calculate empty pathway
        if pathways is None:
            N3 = self.oa3.length            
            ppspec = numpy.zeros(N3, dtype=COMPLEX)
            return ppspec
            
        N3 = self.oa3.length
        
        SS  = self.system.SS.copy()
        
        # precalculate single excited state correlation functions
#        start = time.time()
#        ct3,ct3tau = self._SE_excitonic_cofts(SS,self.system, tau = tau)
#        gt3s,gt3tau = self._SE_excitonic_cofts_test(SS,self.system,tau = tau)
        if self.goft_matrix is not None:
            gt3s = self.goft_matrix
        else:
            gt3s = self._SE_excitonic_gofts(SS,self.system, tau = 0.0)
        gt3tau = self._SE_excitonic_gofts(SS,self.system, tau = tau)
#        end = time.time()
#        print("Calculation of coft:", end - start)
        
#        start = time.time()
#        # convert correlation functions to lineshape functions
#        ngs = self.system.Nb[0]
#        nes = self.system.Nb[1]
#        nfs = self.system.Nb[2]
#        gt3s = numpy.zeros((ngs+nes+nfs,ngs+nes+nfs,self.t3axis.length),dtype=numpy.complex128)
#        for ii in range(ngs, ngs + nes + nfs): #self.system.Ntot): # self.system.Nb[0]+self.system.Nb[1]):
#            for jj in range(ii, ngs + nes + nfs):# self.system.Ntot): # self.system.Nb[0]+self.system.Nb[1]):
#                gt3s[ii,jj,:] = self._c2g(self.t3axis,ct3[ii,jj])
#                if ii!=jj:
#                    gt3s[jj,ii,:] = gt3s[ii,jj,:]
#
#        end = time.time()
#        print("Conversion coft 2 goft:", end - start)
        
#        start = time.time()
        ppspec = numpy.zeros(self.t3axis.length,dtype=numpy.complex128)
        for pwy in self.pathways:
            _is_relax = False
            _is_jump = True
            if len(pwy.relaxations) == 1:
                _is_relax = True
                if pwy.relaxations[0][0] == pwy.relaxations[0][1]:
                    _is_jump = False
            
            om = pwy.frequency[-2] - self.rwa
            pref = pwy.pref
            omtau = pwy.frequency[1]
            
            if pwy.pathway_name not in ["R1f*", "R2f*"]:
#                pass
                if _is_relax and _is_jump: # Pathways with relaxation
                    n = pwy.states[-2][0]
                    #print(om,n,pwy.states[2],pwy.relaxations[0],pwy.pathway_name)
                    ppspec += pref*numpy.exp(-gt3s[n,n] - 1j*om*self.t3axis.data)
                
                else:   # Pathways without relaxation                    
                    if pwy.pathway_name in ["R1g","R2g"]: # Stimulated emission
                        state = pwy.states[1]
                        ft = - 1j*om*self.t3axis.data - 1j*omtau*tau
#                        coft = (ct3tau[state[0],state[0]] 
#                                 - numpy.conj(ct3tau[state[0],state[1]])
#                                 + numpy.conj(ct3[state[1],state[0]]))
#                        goft = self._c2g(self.t3axis,coft)
                        gt1 = gt3tau[state[0],state[0]] - numpy.conj(gt3tau[state[0],state[1]])
                        gt2 = numpy.conj(gt3tau[state[1],state[1],0]) - gt3tau[state[0],state[1],0]
                        gt3 = numpy.conj(gt3s[state[1],state[0]])
                        ft -= gt1 + gt2 + gt3
                        
#                        gt1 = self._c2g(self.t3axis,ct3tau[state[0],state[0]] 
#                                    - numpy.conj(ct3tau[state[0],state[1]]))
#                        gt2 = self._c2g(self.t3axis,numpy.conj(ct3tau[state[1],state[1]])
#                                    - ct3tau[state[0],state[1]])
#                        gt3 = numpy.conj(gt3s[state[1],state[0]])
##                        print(numpy.isclose(goft,gt1+gt3).all())
#                        ft -= gt1 + gt2[0] + gt3
##                        
##                        ft -= goft + gt2[0]
                        ppspec += pref*numpy.exp(ft)
                        
                    elif pwy.pathway_name in ["R3g","R4g"]: # Ground state bleach 
                        ft = - 1j*om*self.t3axis.data - 1j*omtau*tau
                        state = pwy.states[-2]
                        n = max(state)
                        ft -= gt3s[n,n]
                        ppspec += pref*numpy.exp(ft)
            else:
                if _is_relax and _is_jump: # Pathways with relaxation
                    # Ecxited state absorption
                    n = pwy.states[-2][0]
                    m = pwy.states[-2][1]

                    ft = - 1j*om*self.t3axis.data
                    ft -= gt3s[n,n] + numpy.conj(gt3s[m,m])
                    ft += gt3s[n,m] + numpy.conj(gt3s[m,n])
                    ppspec += pref*numpy.exp(ft)
                    
                else:   # Pathways without relaxation
                    # Ecxited state absorption
                    pass
                    Fl = pwy.states[-2][0]
                    Ek = pwy.states[-2][1]
                    Ej = pwy.states[1][0]
                    
                    ft = - 1j*om*self.t3axis.data - 1j*omtau*tau
                    
                    if 0:   # Safer (smaller numerical errors), but slower
                        gt1 = self._c2g(self.t3axis,
                            ct3[Fl,Fl] - ct3[Fl,Ek] + ct3[Ej,Ek] - ct3[Ej,Fl])
                    else:   # Faster, but might cause numerical errors.
                        gt1 = gt3s[Fl,Fl] - gt3s[Fl,Ek] + gt3s[Ej,Ek] - gt3s[Ej,Fl]
                        
                    gt2 = (gt3tau[Ej,Ej,0] - numpy.conj(gt3tau[Ej,Ek,0])
                                    - gt3tau[Fl,Ej,0] + numpy.conj(gt3tau[Fl,Ek,0]))
                                    
                    gt3 = (numpy.conj(gt3tau[Ek,Ek]) - gt3tau[Ek,Ej]
                                    + gt3tau[Fl,Ej] - numpy.conj(gt3tau[Fl,Ek]))
                    ft -= gt1 + gt2 + gt3
                    
#                    gt2 = self._c2g(self.t3axis,
#                                    ct3tau[Ej,Ej] - numpy.conj(ct3tau[Ej,Ek])
#                                    - ct3tau[Fl,Ej] + numpy.conj(ct3tau[Fl,Ek]))
#                                    
#                    gt3 = self._c2g(self.t3axis,
#                                    numpy.conj(ct3tau[Ek,Ek]) - ct3tau[Ek,Ej]
#                                    + ct3tau[Fl,Ej] - numpy.conj(ct3tau[Fl,Ek]))
#
#
#                    ft -= gt1 + gt2[0] + gt3
##                    ft2 = -gt3s[Fl,Fl,:] - numpy.conj(gt3s[Ek,Ek,:]) - 1j*om*self.t3axis.data - 1j*omtau*tau
##                    print(tau,convert(pwy.frequency[-2],"int","nm"),convert(pwy.frequency[-2],"int","1/cm"),numpy.isclose(ft,ft2).all())
                    ppspec += pref*numpy.exp(ft)

#        end = time.time()
#        print("Calculation of the pathways:", end - start)

        ppspec = -ppspec
        
        
        # Fourier transform the result
        ft = numpy.fft.hfft(ppspec)*self.t3axis.step
        ft = numpy.fft.fftshift(ft)
        # invert the order because hfft is a transform with -i
        ft = numpy.flipud(ft)   
        # cut the center of the spectrum
        Nt = self.t3axis.length #len(ta.data)        
        
        data = numpy.real(ft[Nt//2:Nt+Nt//2])

        data = self.oa3.data*data
        
        return data

def calculate_from_2D(twod):
    """Calculates pump-probe spectrum from 2D spectrum
    
    Calculates pump-probe spectrum from 2D spectrum
    using projection theorem
    
    Parameters
    ----------
    
    twod: TwoDSpectrum
        2D spectrum from which pump-probe will be calculated
        
    """
    
    pp = PumpProbeSpectrum()
    
    # waiting time
    t2 = twod.get_t2()
    pp.set_t2(t2)
    
    # time step for integration
    xaxis = twod.xaxis
    dw = xaxis.step
    
    # real part of the total 2D spectrum
    tddata = numpy.real(twod.d__data)
    
    # integration over omega_1 axis
    ppdata = -numpy.sum(tddata,axis=1)*dw
    
    # setting pump-probe data
    pp.set_data(ppdata)
    pp.set_axis(twod.yaxis)
    
    return pp
    


class MockPumpProbeSpectrumCalculator(MockTwoDSpectrumCalculator):
    """Effective line shape pump-probe spectrum calculator
    
    
    """
    
    
    def calculate_all_system(self, sys, H, eUt, lab):
        """Calculates all spectra corresponding to a specified t2axis
        
        """
        
        temporary_fix = True
        
        if temporary_fix:
            # calculation via 2D spectrum
            tcont = super().calculate_all_system(sys, H, eUt, lab)
            return tcont.get_PumpProbeSpectrumContainer()
      
    
    def calculate_one_system(self, t2, sys, H, eUt, lab):
        """Calculates one spectru corresponding to a specified t2 time
        
        """
        
        temporary_fix = True
        
        if temporary_fix:
            # calculation via 2D spectrum
            tcont = super().calculate_one_system(t2, sys, H, eUt, lab)
            return tcont.get_PumpProbeSpectrum()


def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'): # █ = U-219 
    """
    Call n a loop to create terminal progress bar

    Parameters
    --------  
    iteration: int (Required)
        Current iteration 
    total: int (Required)
        Total interactions
    prefix: string (optional)
        Prefix string
    suffix: string (optional)
        Suffix string
    decimal: int (optional)
        positive number of decimals in percent complete
    length: int (optional)
        Character length of bar
    fill: str (optional)
        Fill bar character
    """
    
    percent = ("{:0." + str(decimals) + "f}").format(100* (iteration/float(total)))
    filledlength = int(length * iteration // total)
    bar = fill * filledlength + "-" * (length - filledlength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # print New line after completition
    if iteration == total:
        print()