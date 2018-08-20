# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 14:26:16 2018

@author: Johan
"""

# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    fluor module
    
    This module contains classes to support calculation of linear fluorescence
    spectra.

"""
#import h5py
import numpy
import scipy
import matplotlib.pyplot as plt

#from scipy.optimize import minimize, leastsq, curve_fit

from ..utils import derived_type
from ..builders import Molecule 
from ..builders import Aggregate
from ..core.time import TimeAxis
from ..core.frequency import FrequencyAxis
from ..core.dfunction import DFunction

from ..core.managers import energy_units
from ..core.managers import EnergyUnitsManaged
from ..core.managers import eigenbasis_of
from ..core.time import TimeDependent
from ..core.units import cm2int

from ..core.saveable import Saveable

class FluorSpectrumBase(DFunction, EnergyUnitsManaged):
    """Provides basic container for fluorescence spectrum
    
    """
    
    def __init__(self, axis=None, data=None):
        super().__init__()
        self.axis = axis
        self.data = data
        
    def set_axis(self, axis):
        """Sets axis atribute
        
        Parameters
        ----------
        
        axis : FrequencyAxis object
            Frequency axis object. This object has managed energy units
            
        """
        self.axis = axis
        
    def set_data(self, data):
        """Sets data atribute
        
        Parameters
        ----------
        
        data : array like object (numpy array)
            Sets the data of the fluorescence spectrum
            
        """
        self.data = data
        
    def set_by_interpolation(self, x, y, xaxis="frequency"):
        
        from scipy import interpolate
        
        if xaxis == "frequency":
            
            om = self.convert_2_internal_u(x)
            
        elif xaxis == "wavelength":
            # convert to internal (nano meters) units of wavelength
            
            
            # convert to energy (internal units)
            # to cm
            om = 1.0e-7*x
            # to 1/cm
            om = 1.0/om
            # to 1/fs
            om = om*cm2int
          
        if om[1] > om[2]:
            # reverse order
            om = numpy.flip(om,0)
            y = numpy.flip(y,0)
            
        # equidistant points on the x-axis
        omin = numpy.amin(om)
        omax = numpy.amax(om)
        length = om.shape[0]
        step = (omax-omin)/length
        
        # new frequency axis
        waxis = FrequencyAxis(omin, length, step)
        
        # spline interpolation 
        tck = interpolate.splrep(om, y, s=0)
        ynew = interpolate.splev(waxis.data, tck, der=0)
        
        # setting the axis and data
        self.axis = waxis
        self.data = ynew
        
    
    def clear_data(self):
        """Sets spectrum data to zero
        
        """
        shp = self.data.shape
        self.data = numpy.zeros(shp, dtype=numpy.float64)

    def normalize2(self,norm=1.0):
        """Normalizes spectrum to a given value
        
        """
        mx = numpy.max(self.data)
        self.data = norm*self.data/mx

    def normalize(self):
        """Normalization to one
        
        """
        self.normalize2(norm=1.0)
        
    def subtract(self, val):
        """Subtracts a value from the spectrum to shift its base line
        
        """
        self.data -= val
        

    def add_to_data(self, spect):
        """Performs addition on the data.
        
        Expects a compatible object holding fluorescence spectrum
        and adds its data to the present fluorescence spectrum.
        
        Parameters
        ----------
        
        spect : spectrum containing object
            This object should have a compatible axis and some data
        
        """

        
        if self.axis is None:
            self.axis = spect.axis.copy()
            
        if not numpy.allclose(spect.axis.data, self.axis.data):
            numpy.savetxt("spect_data_wrong.dat", spect.axis.data)
            numpy.savetxt("self_data_wrong.dat", self.axis.data)
            raise Exception("Incompatible axis")
            
        if self.data is None:
            self.data = numpy.zeros(len(spect.data),
                                    dtype=spect.axis.data.dtype)
        
        self.data += spect.data
        
        
    def load_data(self, filename, ext=None, replace=False):
        """Load the spectrum from a file
        
        Uses the load method of the DFunction class to load the fluorescence
        spectrum from a file. It sets the axis type to 'frequency', otherwise
        no changes to the inherited method are applied.
        
        Parameters
        ----------
        
        """
        super().load_data(filename, ext=ext, axis='frequency', replace=replace)

    #save method is inherited from DFunction 
    
        
        
    def plot(self, **kwargs):
        """ Plotting fluorescence spectrum using the DFunction plot method
        
        """
        if "ylabel" not in kwargs:
            ylabel = r'$\alpha(\omega)$ [a.u.]'
            kwargs["ylabel"] = ylabel
            
        fig = super().plot(**kwargs)
        if fig is not None:
            return fig


        
    def gaussian_fit(self, N=1, guess=None, plot=False, Nsvf=251):
        from scipy.signal import savgol_filter
        from scipy.interpolate import UnivariateSpline
        """Performs a Gaussian fit of the spectrum based on an initial guess
        
        
        Parameters
        ----------
        
        Nsvf : int
            Length of the Savitzky-Golay filter window (odd integer)
            
            
        """
        x = self.axis.data
        y = self.data
        
        if guess is None:
            
            raise Exception("Guess is required at this time")
            # FIXME: create a reasonable guess
            guess = [1.0, 11000.0, 300.0, 0.2,
                     11800, 400, 0.2, 12500, 300]
            
            #
            # Find local maxima and guess their parameters
            #

            # Fit with a given number of Gaussian functions
            
            if not self._splines_initialized:
                self._set_splines()
            
            # get first derivative and smooth it
            der = self._spline_r.derivative()
            y1 = der(x)
            y1sm = savgol_filter(y1,Nsvf,polyorder=3)
        
            # get second derivative and smooth it
            y1sm_spl_der = UnivariateSpline(x,y1sm,s=0).derivative()(x)
            y2sm = savgol_filter(y1sm_spl_der,Nsvf,polyorder=3)
        
            # find positions of optima by looking for zeros of y1sm
        
        
            # isolate maxima by looking at the value of y2sm
        

            #plt.plot(x, der(x))
            #plt.plot(x, y1sm)
            plt.plot(x, y2sm)
            plt.show()
        
        
        
        def funcf(x, *p):
            return _n_gaussians(x, N, *p)
        
        # minimize, leastsq,
        from scipy.optimize import curve_fit            
        popt, pcov = curve_fit(funcf, x, y, p0=guess)
        
        if plot:
        
            plt.plot(x,y)
            plt.plot(x,_n_gaussians(x, N, *popt))
            for i in range(N):
                a = popt[3*i]
                print(i, a)
                b = popt[3*i+1]
                c = popt[3*i+2]
                y = _gaussian(x, a, b, c)
                plt.plot(x, y,'-r')
            plt.show()
        
        # FIXME: Create a readable report
        
        return popt, pcov
        
#    def convert_to_energy(self, eaxis, units):
#        """
#        
#        """
#        
#        if units == "nm":
#            x = self.axis.data
#            y = self.data
#            
#            # to cm
#            x = 1.0e-7*x
#            # to 1/cm
#            x = 1.0/x
#            # to rad/fs
#            x = x*cm2int
#            
#            xn = numpy.zeros(x.shape, dtype=x.dtype)
#            yn = numpy.zeros(y.shape, dtype=y.dtype) 
#            
#            for i in range(len(x)):
#                xn[i] = x[len(x)-i-1]
#                yn[i] = y[len(x)-i-1]
#                
#            # spline it
#            
#            # evaluate at points if eaxis
#
            
def _gaussian(x, height, center, fwhm, offset=0.0):
    """Gaussian function with a possible offset
    
    
    Parameters
    ----------
    
    x : float array
        values to calculate Gaussian function at
        
    height : float
        height of the Gaussian at maximum
        
    center : float
        position of maximum
        
    fwhm : float
        full width at half maximum of the Gaussian function
        
    offset : float
        the value at infinity; effectively an offset on the y-axis
        
    
    """
    
    return height*numpy.exp(-(((x - center)**2)*4.0*numpy.log(2.0))/
                            (fwhm**2)) + offset   


def _n_gaussians(x, N, *params):
    """Sum of N Gaussian functions plus an offset from zero

    Parameters
    ----------
    
    x : float
        values to calculate Gaussians function at        

    N : int
        number of Gaussians
        
    params : floats
        3*N + 1 parameters corresponding to height, center, fwhm  for each 
        Gaussian and one value of offset
        
    """
    n = len(params)
    k = n//3
    
    if (k*3 == n) and (k == N):
        
        res = 0.0
        pp = numpy.zeros(3)
        for i in range(k):
            pp[0:3] = params[3*i:3*i+3]
            #pp[3] = 0.0
            arg = tuple(pp)
            res += _gaussian(x, *arg)
        res += params[n-1] # last parameter is an offset
        return res
            
    else:
        raise Exception("Inconsistend number of parameters")        


class FluorSpectrum(FluorSpectrumBase):
    """Class representing fluorescence spectrum
    
    
    """
    pass
    
#    def _non_data_attributes(self):
#        
#        return {"object_tag":str(self),
#                "axis_real":[self.axis.start, 
#                             self.axis.step,
#                             self.axis.time_start],
#                "axis_int":[self.axis.length],
#                "axis_str":[self.axis.atype]}
#    
#    def save(self, filename):
#        """Save the whole FluorSpectrum object
#        
#        Attributes saved:
#            axis
#            data
#            
#        """
#        # FIXME: convert into hdf5
#        with energy_units("int"):
#            att = self._non_data_attributes()
#            numpy.savez(filename,
#                    axis_real=att["axis_real"],
#                    axis_int=att["axis_int"],
#                    axis_str=att["axis_str"],
#                    axis_data=self.axis.data,
#                    data=self.data)    
#    
#    def load(self, filename):
#        """Loads the FluorSpectrum object from file 
#        
#        """
#        # FIXME: convert into hdf5        
#        infile = numpy.load(filename)
#        
#        start = infile['axis_real'][0]
#        step = infile['axis_real'][1]
#        time_start = infile['axis_real'][2]
#        atype = infile['axis_str'][0]
#        length = infile['axis_int'][0]
#        
#        with energy_units('int'):
#            faxis = FrequencyAxis(start, length, step, atype=atype,
#                              time_start=time_start)
#        
#        data = infile['data']
#        
#        self.set_axis(faxis)
#        self.set_data(data)

  
    
    
class FluorSpectrumContainer(Saveable):
    
    def __init__(self, axis=None):

        self.axis = axis
        self.count = 0
        self.spectra = {}
        
    def set_axis(self, axis):
        self.axis = axis
        
    def set_spectrum(self, spect, tag=None):
        """Stores fluorescence spectrum 
        
        Checks compatibility of its frequency axis
        
        """
        frq = spect.axis
        
        if self.axis is None:
            self.axis = frq
            
        if self.axis.is_equal_to(frq):
            if tag is None:
                tag1 = str(self.count)
            else:
                tag1 = str(tag)
            self.spectra[tag1] = spect
            self.count += 1
        else:
            raise Exception("Incompatible time axis (equal axis required)")
            
            
    def get_spectrum(self, tag):
        """Returns spectrum corresponing to time t2
        
        Checks if the time t2 is present in the t2axis
        
        """  
        if not isinstance(tag, str):
            tag = str(tag)
            
        if tag in self.spectra.keys():
            return self.spectra[tag]     
        else:
            raise Exception("Unknown spectrum")

        
    def get_spectra(self):
        """Returns a list or tuple of the calculated spectra
        
        """
        
        ven = [value for (key, value) in sorted(self.spectra.items())]
        return ven        
    
#    def save(self, filename):
#        """Save the whole set of spectra into an hdf5 file
#
#            
#        """
#        
#        with h5py.File(filename,"w") as f:
#            sps = f.create_group("spectra")
#            sps.attrs.create("count",self.count)
#            dt = h5py.special_dtype(vlen=str)
#            with energy_units("int"):
#                for ks in self.spectra.keys():
#                    # spectrum object to be saved
#                    sp = self.spectra[ks]
#                    # directory name
#                    dtname = "spectrum_"+str(ks)
#                    # creating directory to save this object
#                    spsub = sps.create_group(dtname)
#                    spsub.attrs.create("tag",ks)
#                    # fluorescence data to be saved
#                    data2save = spsub.create_dataset("data",data=sp.data)
#                    # axis data get some more attributes
#                    atr = sp._non_data_attributes()
#                    data2save.attrs.create("axis_real",atr["axis_real"])
#                    data2save.attrs.create("axis_int",atr["axis_int"])
#                    data2save.attrs.create("axis_str",atr["axis_str"],dtype=dt)
#    
#    
#    
#    def load(self, filename):
#        """Loads the whole set of spectra from an hdf5 file
#        
#        
#        """
#        with h5py.File(filename,"r") as f:
#            sps = f["spectra"]
#            expected_count = sps.attrs["count"]
#            with energy_units("int"):
#                for grp in sps.keys():
#                    spectrum = sps[grp]
#                    tag = spectrum.attrs["tag"]
#                    data = numpy.array(spectrum["data"])
#                    a_real = spectrum["data"].attrs["axis_real"]
#                    a_int = spectrum["data"].attrs["axis_int"]
#                    a_str = spectrum["data"].attrs["axis_str"]
#    
#                    start = a_real[0]
#                    step = a_real[1]
#                    time_start = a_real[2]
#                    atype = a_str[0]
#                    length = a_int[0]
#                    
#                    with energy_units('int'):
#                        faxis = FrequencyAxis(start, length, step, atype=atype,
#                                          time_start=time_start)
#                    aaa = FluorSpectrum(axis=faxis, data=data)
#                    if self.count == 0:
#                        self.set_axis(aaa.axis)
#                    self.set_spectrum(aaa,tag=tag)
#                
#            if self.count != expected_count:
#                raise Exception("Incorrect number of spectra read")
    
    
class FluorSpectrumCalculator(EnergyUnitsManaged):
    """Linear fluorescence spectrum 
    
    Linear fluorescence spectrum of a molecule or an aggregate of molecules.
    
    Parameters
    ----------
    timeaxis : quantarhei.TimeAxis
        TimeAxis on which the calculation will be performed
        
    system : quantarhei.Molecule or quantathei.Aggregate
        System for which the fluorescence spectrum will be calculated.
        
        
    """

    TimeAxis = derived_type("TimeAxis",TimeAxis)
    system = derived_type("system",[Molecule,Aggregate])
    
    def __init__(self, timeaxis,
                 system=None,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None,
                 effective_hamiltonian=None,
                 temperature=300):
        
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
        self.temperature = temperature
        self.rwa = 0.0
        
    def bootstrap(self,rwa=0.0):
        """
        
        """
        self.rwa = self.convert_2_internal_u(rwa)
        with energy_units("int"):
            # sets the frequency axis for plottig
            self.frequencyAxis = self.TimeAxis.get_FrequencyAxis()
            self.frequencyAxis.data += self.rwa     
        
        
        
    def calculate(self):
        """ Calculates the fluorescence spectrum 
        
        
        """
        
        with energy_units("int"):
            if self.system is not None:
                if isinstance(self.system,Molecule):
                    #self._calculate_Molecule(rwa)      
                    spect = self._calculate_monomer()
                elif isinstance(self.system,Aggregate):
                    spect = self._calculate_aggregate( 
                                              relaxation_tensor=
                                              self._relaxation_tensor,
                                              rate_matrix=
                                              self._rate_matrix,
                                              relaxation_hamiltonian=
                                              self._relaxation_hamiltonian)
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
            re = tr["re"] # reorganisation energy

            # convert correlation function to lineshape function
            gt = self._c2g(ta,ct.data)
            # calculate time dependent response
            at = numpy.exp(-numpy.conjugate(gt) -1j*om*ta.data + 2j*re*ta.data)
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

    def _equilibrium_populations(self, AG, temperature=300,
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
    
    def _excitonic_reorg_energy(self, SS, AG, n):
        """ Returns the reorganisation energy of an exciton state
        """
    
        c0 = AG.monomers[0].get_egcf((0,1))
        Nt = len(c0)
        
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        
        reorg = numpy.zeros((Nt),dtype=numpy.complex128)
        Na = AG.nmono
        for kk in range(Na):
            reorg += ((SS[kk+1,n+1]**2)*(SS[kk+1,n+1]**2)*cfm.get_reorganization_energy(kk,kk))
        
        return reorg    
        
    def _calculate_monomer(self):
        """ Calculates the fluorescence spectrum of a monomer 
        
        
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
            sbi = self.system.get_SystemBathInteraction(ta)
            re = sbi.CC.get_reorganization_energy(0,0)
            tr = {"ta":ta,"dd":dd,"om":om-self.rwa,"ct":ct,"gg":gama, "re":re}
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
        
        spect = FluorSpectrum(axis=axis, data=data)
        
        return spect
        
        
    def _calculate_aggregate(self, relaxation_tensor=None,
                             relaxation_hamiltonian=None, rate_matrix=None):
        """ Calculates the fluorescence spectrum of a molecular aggregate
        
        
        
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
        
        re = self._excitonic_reorg_energy(SS,self.system,0)
        tr["re"] = re
        rho_eq = self._equilibrium_populations(self.system,
                                               temperature=self.temperature)
        #
        # Calculates spectrum of a single transition
        #
        data = rho_eq.data[1, 1]*numpy.real(self.one_transition_spectrum(tr))
        
        for ii in range(2,HH.dim):
            if relaxation_tensor is not None:
                tr["gg"] = gg[ii]
            elif rate_matrix is not None:
                tr["gg"] = gg[ii]
            else:
                tr["gg"] = [0.0]
            #tr[1] = DD.dipole_strength(0,ii) # update transition dipole moment
            tr["dd"] = DD.dipole_strength(0,ii)
            #tr[2] = HH.data[ii,ii]-HH.data[0,0]-rwa
            tr["om"] = HH.data[ii,ii]-HH.data[0,0]-self.rwa
            #tr[3] = self._excitonic_coft(SS,self.system,ii-1) # update ct here
            tr["ct"] = self._excitonic_coft(SS,self.system,ii-1)
            tr["re"] = self._excitonic_reorg_energy(SS,self.system,ii-1)            

            
            #
            # Calculates spectrum of a single transition
            #
            data += rho_eq.data[ii, ii]*numpy.real(self.one_transition_spectrum(tr))


        # we only want to retain the upper half of the spectrum
        Nt = len(self.frequencyAxis.data)//2        
        do = self.frequencyAxis.data[1]-self.frequencyAxis.data[0]
        st = self.frequencyAxis.data[Nt//2]
        # we represent the Frequency axis anew
        axis = FrequencyAxis(st,Nt,do)
        
        
        # transform all quantities back
        S1 = numpy.linalg.inv(SS)
        HH.transform(S1)
        DD.transform(S1)
        
        if relaxation_tensor is not None:
            RR.transform(S1)

        spect = FluorSpectrum(axis=axis, data=data)
        
        return spect        

                   
