# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    abs module
    
    This module contains classes to support calculation of linear absorption
    spectra.

"""
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
from ..core.time import TimeDependent
from ..core.units import cm2int

class AbsSpectrumBase(DFunction, EnergyUnitsManaged):
    """Provides basic container for absorption spectrum
    
    
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
            Sets the data of the absorption spectrum
            
        """
        self.data = data
    
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

    def add_to_data(self, spect):
        """Performs addition on the data.
        
        Expects a compatible object holding absorption spectrum
        and adds its data to the present absorption spectrum.
        
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
        
        
    def load(self, filename, ext=None, replace=False):
        """Load the spectrum from a file
        
        Uses the load method of the DFunction class to load the absorption
        spectrum from a file. It sets the axis type to 'frequency', otherwise
        no changes to the inherited method are applied.
        
        Parameters
        ----------
        
        """
        super().load(filename, ext=ext, axis='frequency', replace=replace)

    #save method is inherited from DFunction    
        
    def plot(self,**kwargs):
        """ Plotting absorption spectrum using the DFunction plot method
        
        """
        if "ylabel" not in kwargs:
            ylabel = r'$\alpha(\omega)$ [a.u.]'
            kwargs["ylabel"] = ylabel
            
        super().plot(**kwargs)


        
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
        
    #
    # Temporary methods for compatibility and testing
    #        
    def get_AbsSpectContainer(self, timeaxis):
        
        asc = AbsSpectContainer(timeaxis)
        asc.set_data(self.data)
        asc.axis = self.axis
        
    def convert_to_energy(self, eaxis, units):
        """
        
        """
        
        if units == "nm":
            x = self.axis.data
            y = self.data
            
            # to cm
            x = 1.0e-7*x
            # to 1/cm
            x = 1.0/x
            # to rad/fs
            x = x*cm2int
            
            xn = numpy.zeros(x.shape, dtype=x.dtype)
            yn = numpy.zeros(y.shape, dtype=y.dtype) 
            
            for i in range(len(x)):
                xn[i] = x[len(x)-i-1]
                yn[i] = y[len(x)-i-1]
                
            # spline it
            
            # evaluate at points if eaxis

            
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
         
    
            
     
class AbsSpectrumDifference(EnergyUnitsManaged):
    """Difference between two absorption spectra and its minimuzation
    
    Describes and handles difference between two absorption spectra. Provides
    minimization facility provided that the second spectrum is specified
    by a function.
    
    Parameters
    ----------
    
    target : AbsSpectrum
        Absorption spectrum to be compared to
        
    optfce : function or AbsSpectrum
        If a function is specified it is expected that one will try to optimize
        its parameters to fit the target spectrum. If an absorption spectrum 
        is submitted, one can calculate the difference between the target and
        submitted spectra. Attempted to call `minimize` method will fail.
    
    Methods
    -------

    difference

    minimize
     
        
    """
    def __init__(self, target=None, optfce=None, bounds=(0.0,1.0), tol=1.0e-6 ):
        if target is None:
            raise Exception("Target spectrum has to be specified")
        if optfce is None:
            raise Exception("Second spectrum or a function to be optimized"
                            + " must be specified")
            
        self.target = target
        if callable(optfce):
            self.optfce = optfce
            self.secabs = None
            self._can_minimize = True
        else:
            self.optfce = None
            self.secabs = optfce
            self._can_minimize = False
            
        self.bounds = []
        self.bounds.append(self.convert_2_internal_u(bounds[0]))
        self.bounds.append(self.convert_2_internal_u(bounds[1]))

        with energy_units("int"):
            nl, val = target.axis.locate(self.bounds[0])
            nu, val = target.axis.locate(self.bounds[1])
        
        self.nl = nl
        self.nu = nu
        self.x = target.axis._data[nl:nu]
        self.tol = tol
        self.opt_result = None
        
        self.difftype = "square"
        
    def difference(self, par=None):
        """Calculates difference between spectra
        
        Calculates difference between the target spectrum and the spectrum
        calculated from submitted parameters
        
        Parameters
        ----------
        
        par : list or array (optional)
            parameters of the function 
        
        """
        target = self.target.data[self.nl:self.nu]
        if self._can_minimize:
            if par is None:
                raise Exception("Function parameters must be specified "+
                                "to calculate difference")
            secabs = self.optfce(par)
            sdat = numpy.zeros(len(self.x), dtype=numpy.float64)
            #
            # get values at the correct points by fitting (using the fact that
            # absorption spectrum is a DFunction )
            #
            i = 0

            for xi in self.x:
                sdat[i] = secabs.at(xi)
                i += 1
            secabs = sdat    
        else:
            secabs = self.secabs #.data[self.nl:self.nu]
            sdat = numpy.zeros(len(self.x), dtype=numpy.float64)
            #
            # get values at the correct points by fitting (using the fact that
            # absorption spectrum is a DFunction )
            #
            i = 0
            for xi in self.x:
                sdat[i] = secabs.at(xi)
                i += 1
            secabs = sdat    

            
        #
        # FIXME: Add some more types of differences
        #
            
        if self.difftype == "square":
            diff = 1000.0*numpy.sum(numpy.abs((target-secabs)**2/
                                   (self.x[len(self.x)-1]-self.x[0])))
        elif self.difftype == "measure":
            diff = 0.0
        else:
            raise Exception("Unknown differene type")
            
        print("DIFF: ", diff)
        
        return diff
        
    def minimize(self, init_params, method):
        """Minimizes the submitted function and returns optimal parameters
        
        """
        if self._can_minimize:
            from scipy.optimize import minimize
            self.opt_result = minimize(self.difference, init_params,
                                       method=method, tol=self.tol,
                                       options=dict(disp=True))
            return self.opt_result.x
        else:
            raise Exception("Cannot perform minimization, "+
                            "no function suplied")
        
        
class AbsSpectContainer(DFunction, EnergyUnitsManaged):
    """This class contains a single absorption spectrum
    
    Contains absorption spectrum and enables some manipulations on it.
    
    """
    TimeAxis = derived_type("TimeAxis",TimeAxis)
    
    def __init__(self, timeaxis):
        self.TimeAxis = timeaxis
        self.data = None
        self.axis = None
             
    def set_data(self, data):
        self.data = data

    def clear_data(self):
        shp = self.data.shape
        self.data = numpy.zeros(shp, dtype=numpy.float64)

    def normalize2(self,norm=1.0):
        mx = numpy.max(self.data)
        self.data = self.data/mx

    def add_to_data(self, spect):

        if not numpy.allclose(spect.TimeAxis.data, self.TimeAxis.data):
            raise Exception("Incompatible TimeAxis")
        
        if self.data is None:
            self._copy_internals(spect)
        else:
            self.data += spect.data

            
    def _copy_internals(self, spect):
        self.data = numpy.zeros(self.TimeAxis.length, dtype=numpy.float64)
        self.data = spect.data.copy()
        self.axis = spect.axis

    def _frequency(self,dt):
        """ Calculates the frequency axis corresponding to TimeAxis
        
        
        """
        Nt = self.TimeAxis.length
        return numpy.pi*numpy.fft.fftshift(
              numpy.fft.fftfreq(Nt,d=dt))       
                
        
    def plot(self,**kwargs):
        """ Plotting absorption spectrum using the DFunction plot method
        
        """
        if "ylabel" not in kwargs:
            ylabel = r'$\alpha(\omega)$ [a.u.]'
            kwargs["ylabel"] = ylabel
            
        super(AbsSpectContainer,self).plot(**kwargs)
        

    #
    # Temporary methods for compatibility and testing
    #
    def get_AbsSpectrumBase(self):
        return AbsSpectrumBase(axis=self.axis, data=self.data)
        
        
class AbsSpect(AbsSpectContainer):
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
                    self._calculate_aggregate(rwa, 
                                              relaxation_tensor=
                                              self._relaxation_tensor,
                                              rate_matrix=
                                              self._rate_matrix,
                                              relaxation_hamiltonian=
                                              self._relaxation_hamiltonian)
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
            rt = numpy.exp(gg*ta.data)          
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
        gama = [-1.0/self.system.get_electronic_natural_lifetime(1)]
        
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
        
        
    def _calculate_aggregate(self, rwa, relaxation_tensor=None,
                             relaxation_hamiltonian=None, rate_matrix=None):
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
        tr["om"] = HH.data[1,1]-HH.data[0,0]-rwa
        # get a transformed ct here
        ct = self._excitonic_coft(SS,self.system,0)
        #tr.append(ct)
        tr["ct"] = ct
        self.system._has_system_bath_coupling = True
        
        #
        # Calculates spectrum of a single transition
        #
        self.data = numpy.real(self.one_transition_spectrum(tr))
        
        for ii in range(2,HH.dim):
            if relaxation_tensor is not None:
                tr["gg"] = gg[ii]
            else:
                tr["gg"] = [0.0]
            #tr[1] = DD.dipole_strength(0,ii) # update transition dipole moment
            tr["dd"] = DD.dipole_strength(0,ii)
            #tr[2] = HH.data[ii,ii]-HH.data[0,0]-rwa
            tr["om"] = HH.data[ii,ii]-HH.data[0,0]-rwa
            #tr[3] = self._excitonic_coft(SS,self.system,ii-1) # update ct here
            tr["ct"] = self._excitonic_coft(SS,self.system,ii-1)
            
            #
            # Calculates spectrum of a single transition
            #
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
        
        # transform all quantities back
        S1 = numpy.linalg.inv(SS)
        HH.transform(S1)
        DD.transform(S1)
        
        if relaxation_tensor is not None:
            RR.transform(S1)
        

                    
        