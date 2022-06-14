# -*- coding: utf-8 -*-
"""Class representing evolution superoperator

    This class represents evolution superoperator at discrete times. 
    
    
    
    Examples
    --------
    
    Evolution superoperator can be
    created from a known Hamiltonian and relaxation tensor
    
    >>> import quantarhei as qr
    >>> # We use the predefined Aggregate - dimer of two-level systems
    >>> # with an environment
    >>> tagg = qr.TestAggregate(name="dimer-2-env")
    >>> tagg.set_coupling_by_dipole_dipole()
    >>> tagg.build()
    >>> # we get the existing predefined time axis for the calculations
    >>> time = tagg.get_SystemBathInteraction().TimeAxis
    >>> print(time.start, time.length, time.step)
    0 1000 1.0
    
    We create relaxation tensor and obtain renormalized Hamiltonian from 
    standard Redfield theory (`stR`)
    
    >>> RR, HH = tagg.get_RelaxationTensor(time, relaxation_theory="stR")
    
    We initialize evolution superoperator
    
    >>> eSO = EvolutionSuperOperator(time, ham=HH, relt=RR)
    >>> eSO.calculate()
    
    We create initial condition. Dimension of the problem is 3 (ground state
    and two excited states of the dimer)
    
    >>> dim = HH.dim
    >>> print(dim)
    3
    
    >>> rho = qr.ReducedDensityMatrix(dim=dim)
    >>> rho.data[2,2] = 1.0
    
    Now we calculate density matrix evolution using a propagator

    >>> prop = tagg.get_ReducedDensityMatrixPropagator(time, 
    ...                                                relaxation_theory="stR")
    >>> rhot = prop.propagate(rho)
    
    Comparing the results at one point can be done like this:
    
    >>> rhot_eSO = eSO.apply(100.0, rho)
    >>> dat_100_eSO = rhot_eSO.data
    >>> dat_100     = rhot.data[100,:,:]
    >>> numpy.allclose(dat_100_eSO, dat_100)
    True
    
    To see the full aggreement, we may plot all the times calculated with
    both methods. We plot the dynamics calculated using the evolution
    superoperator only every 20 fs using asterisks "*"
        
    .. plot::
        
        import quantarhei as qr
        import matplotlib.pyplot as plt
        import numpy
        
        # We use the predefined Aggregate - dimer of two-level systems
        # with an environment
        tagg = qr.TestAggregate(name="dimer-2-env")
        tagg.set_coupling_by_dipole_dipole()
        tagg.build()
        
        # we get the existing predefined time axis for the calculations
        time = tagg.get_SystemBathInteraction().TimeAxis
        
        RR, HH = tagg.get_RelaxationTensor(time, relaxation_theory="stR")
        
        eSO = qr.qm.EvolutionSuperOperator(time, ham=HH, relt=RR)
        eSO.calculate()
        
        dim = HH.dim
        
        rho = qr.ReducedDensityMatrix(dim=dim)
        rho.data[2,2] = 1.0
        
        prop = tagg.get_ReducedDensityMatrixPropagator(time,
                                                       relaxation_theory="stR")
        rhot = prop.propagate(rho)
        
        rhot.plot(coherences=False, show=False)
        
        time_skip = numpy.arange(0.0, 1000, 20.0)
        rhot_skip = numpy.zeros((50,3,3))
        k = 0
        for tm in time_skip:
            rhot_skip[k,:,:] = numpy.real(eSO.apply(tm, rho).data)
            k += 1
        
        plt.plot(time_skip,numpy.real(rhot_skip[:,1,1]),"*g")
        plt.plot(time_skip,numpy.real(rhot_skip[:,2,2]),"*g")
        
        plt.show()
        
        
    We may want to calculate the evolotion superoperator with a step larger
    than the one which we specified for evaluation of bath correlation
    functions (In this example we use predefined  SystemBathInteraction object
    which holds this information). Our time axis is too dense for our needs.
    We specify a less dense one
    
    >>> time2 = qr.TimeAxis(0.0, 1000, 50.0)
    
    This one has the step of 50 fs. We define an evolution superoperator
    with this time:
        
    >>> eSO2 = qr.qm.EvolutionSuperOperator(time2, HH, RR)
    
    Now, to obtain the same results as before, we need to set a time step
    of the propagation to 1 fs as before. This is done by setting a "dense"
    time step with is N times shorter than the one specified in the time 
    axis. In our case N = 50
    
    >>> eSO2.set_dense_dt(50)
    >>> eSO2.calculate()
    >>> rhot_eSO = eSO.apply(100.0, rho)
    >>> dat_100_eSO = rhot_eSO.data
    >>> numpy.allclose(dat_100_eSO, dat_100)
    True
    
    We can calculate a similar picture as before, but now with an evolution
    superoperator calculated only every 50 fs.
    
    .. plot::
        
        import quantarhei as qr
        import matplotlib.pyplot as plt
        import numpy
        
        # We use the predefined Aggregate - dimer of two-level systems
        # with an environment
        tagg = qr.TestAggregate(name="dimer-2-env")
        tagg.set_coupling_by_dipole_dipole()
        tagg.build()
        
        # we get the existing predefined time axis for the calculations
        time = tagg.get_SystemBathInteraction().TimeAxis
        time2 = qr.TimeAxis(0.0, int(time.length/50), 50*time.step)
        
        RR, HH = tagg.get_RelaxationTensor(time, relaxation_theory="stR")
        
        eSO = qr.qm.EvolutionSuperOperator(time2, ham=HH, relt=RR)
        eSO.set_dense_dt(50)
        eSO.calculate()
        
        dim = HH.dim
        
        rho = qr.ReducedDensityMatrix(dim=dim)
        rho.data[2,2] = 1.0
        
        prop = tagg.get_ReducedDensityMatrixPropagator(time,
                                                       relaxation_theory="stR")
        rhot = prop.propagate(rho)
        
        rhot.plot(coherences=False, show=False)
        
        
        rhot_skip = numpy.zeros((time2.length,3,3))
        k = 0
        for tm in time2.data:
            rhot_skip[k,:,:] = numpy.real(eSO.apply(tm, rho).data)
            k += 1
        
        plt.plot(time2.data,numpy.real(rhot_skip[:,1,1]),"*g")
        plt.plot(time2.data,numpy.real(rhot_skip[:,2,2]),"*g")
        
        plt.show()
                                                                                                                  

    Class Details
    -------------

"""

# standard library imports
import time
import numbers

# dependencies imports
import numpy

# quantarhei imports
from ..propagators.rdmpropagator import ReducedDensityMatrixPropagator
from ..propagators.dmevolution import ReducedDensityMatrixEvolution
from ..hilbertspace.operators import ReducedDensityMatrix
from ...core.time import TimeAxis
from ...core.saveable import Saveable

from .superoperator import SuperOperator
from ...core.time import TimeDependent
from ... import COMPLEX
import matplotlib.pyplot as plt

import quantarhei as qr

#from ...utils.types import BasisManagedComplexArray

class EvolutionSuperOperator(SuperOperator, TimeDependent, Saveable):
    """Class representing evolution superoperator
    
    
    
    Parameters
    ----------
    
    time : TimeAxis
        TimeAxis obejct specifying the time points of the superoperator
    
    ham: Hamiltonian
        Hamiltonian of the system
        
    relt: relaxation tensor
        Relaxation tensor of the system
        
    mode: str 
        Mode of information storage. It can be "all" which means that all
        points of the evolution are stored, or it can be "jit" = just in (one)
        time. In the "jit" mode, only the "present" time of the evolution
        operator is stored.
    
    """
    
    def __init__(self, time=None, ham=None, relt=None, pdeph=None, mode="all"):
        super().__init__()
        
        self.time = time
        self.ham = ham
        self.relt = relt
        self.mode = mode
        self.pdeph = pdeph
        
        try:
            self.dim = ham.dim
        except:
            self.dim = 1
        
        self.dense_time = None
        self.set_dense_dt(1)
            
        if (self.time is not None) and (self.mode == "all"):
            
            self.data = numpy.zeros((self.time.length, self.dim, self.dim,
                                     self.dim, self.dim),
                                     dtype=qr.COMPLEX)
            #
            # zero time value (unity superoperator)
            #
            dim = self.dim
            for i in range(dim):
                for j in range(dim):
                    self.data[0,i,j,i,j] = 1.0
                    
        else:
            
            self.data = numpy.zeros((self.dim, self.dim, self.dim, self.dim),
                                dtype=qr.COMPLEX) 
            #
            # zero time value (unity superoperator)
            #
            dim = self.dim
            for i in range(dim):
                for j in range(dim):
                    self.data[i,j,i,j] = 1.0
                
        self.now = 0
        

    def get_Hamiltonian(self):
        """Returns the Hamiltonian associated with thise evolution
        
        """
        return self.ham

    def set_dense_dt(self, Nt):
        """Set a denser time axis for calculations between two points of the superoperator
        
        Parameters
        ----------
        
        Nt : int
            Number of steps between two points of the superoperator to be used
            for numerical propagation
            
        """
        self.dense_time = TimeAxis(0.0, Nt+1, self.time.step/Nt)
        
        
    def update_dense_time(self, i):
        """Update the start time of the dense_time
        
        
        """
        start = self.time.data[i]
        
        self.dense_time.start = start
        self.dense_time = TimeAxis(self.dense_time.start, 
                                   self.dense_time.length,
                                   self.dense_time.step) 
        
        
    def set_PureDephasing(self, pdeph):
        """Sets the PureDephasing object for the dynamic calculation
        
        """
        self.pdeph = pdeph        


    def has_PureDephasing(self):
        """Return True if the EvolutionSuperOperator has pure dephasing
        
        """
        if self.pdeph is None:
            return False
        else:
            return True


    def _initialize_data(self, save=False):
        """Initializes EvolutionSuperOperator data
        
        """
        
        #
        # Create new data
        #      

        if (self.mode == "all") or save:          

            # if we are supposed to save all time steps
            Nt = self.time.length                    
            self.data = numpy.zeros((Nt, self.dim, self.dim,
                                         self.dim, self.dim),
                                         dtype=qr.COMPLEX)
            #
            # zero time value (unity superoperator)
            #
            dim = self.dim
            for i in range(dim):
                for j in range(dim):
                    self.data[0,i,j,i,j] = 1.0
                
        elif self.mode == "jit":
            
            # if we need to keep only the last state
            if self.dim != self.data.shape[0]:
                self.data = numpy.zeros((self.dim, self.dim,
                                         self.dim, self.dim),
                                         dtype=qr.COMPLEX)

            #
            # zero time value (unity superoperator)
            #
            dim = self.dim
            for i in range(dim):
                for j in range(dim):
                    self.data[i,j,i,j] = 1.0
        
        

    def calculate(self, show_progress=False):
        """Calculates the data of the evolution superoperator
        
        
        Parameters
        ----------
        
        show_progress : bool
            When set True, reports on its progress and elapsed time
        
        
        """
        if self.mode != "all":
            raise Exception("This method (calculate()) can be used only"+
                            " with mode='all'")
        Nt = self.time.length
        
        self._initialize_data()
            
        if show_progress:
            print("Calculating evolution superoperator ")
        self._init_progress()
        
        #
        # Let us propagate from t0 to t0+dt*(self.dense_time.length-1)
        #
        
        if (self.pdeph is not None) and (self.pdeph.dtype == "Gaussian"):
            
            #
            # We calculate every interval completely
            #
            for ti in range(1, Nt):
                
                t0 = self.time.data[ti-1] # initial time of the propagation

                Ut1 = self._elemental_step_TimeDependent(t0)
                
                self.data[ti,:,:,:,:] = \
                    numpy.tensordot(Ut1, self.data[ti-1,:,:,:,:])
                    
                if show_progress:
                    print("propagation: ", ti, "of", Nt)            
                
        else:
            
            #
            # We calculate the first step of the first interval
            #
            t0 = 0.0

            self.data[1,:,:,:,:] = self._one_step_with_dense_TimeIndep(t0,
                                                    self.dense_time.length,
                                                    self.dense_time.step,
                                                    Nt,
                                                    show_progress) 
            
            
            #
            # repeat propagation over the longer interval
            #
            self._calculate_remainig_using_first_interval(Nt,
                                                          show_progress)

        if show_progress:
            print("...done")
            
            
    def _elemental_step_TimeIndep(self, t0, dens_dt, Nt, show_progress=False):
        """Single elemental step of propagation with the dense time step
        
        """
        dim = self.ham.dim        
        one_step_time = TimeAxis(t0, 2, self.dense_time.step)
        prop = ReducedDensityMatrixPropagator(one_step_time, self.ham, 
                                              RTensor=self.relt, 
                                              PDeph=self.pdeph)
        rhonm0 = ReducedDensityMatrix(dim=dim)
        Ut1 = numpy.zeros((dim, dim, dim, dim), dtype=COMPLEX)
        for n in range(dim):
            if show_progress:
                self._progress(Nt, dim, 0, n, 0)
            for m in range(dim):
                rhonm0.data[n,m] = 1.0
                rhot = prop.propagate(rhonm0)
                Ut1[:,:,n,m] = rhot.data[1,:,:]
                rhonm0.data[n,m] = 0.0
                
        return Ut1


    def _elemental_step_TimeDependent(self, t0):
        """Single step of propagation with the dense time step 
        
        assuming time dependent relaxation tensor
        
        """
        
        dim = self.dim
        one_step_time = TimeAxis(t0, self.dense_time.length,
                                 self.dense_time.step)
        
        prop = ReducedDensityMatrixPropagator(one_step_time, self.ham,
                                              RTensor=self.relt, 
                                              PDeph=self.pdeph)
        rhonm0 = ReducedDensityMatrix(dim=dim)
        Ut1 = numpy.zeros((dim, dim, dim, dim), dtype=COMPLEX)
        for n in range(dim):
            for m in range(dim):
                rhonm0.data[n,m] = 1.0
                rhot = prop.propagate(rhonm0)
                Ut1[:,:,n,m] = rhot.data[one_step_time.length-1,:,:]
                rhonm0.data[n,m] = 0.0
                #self.now += 1
        return Ut1


    def _one_step_with_dense_TimeIndep(self, t0, Ndense, dens_dt, Nt,
                                     show_progress=False):
        """One step of propagation over the standard time step
        
        This step is componsed of Ndense time steps
        
        """
        Ut1 = self._elemental_step_TimeIndep(t0, dens_dt, Nt,
                                           show_progress=show_progress)
        #
        # propagation to the end of the first interval
        #
        Udt = numpy.zeros(Ut1.shape, dtype=COMPLEX)
        Udt[:,:,:,:] = Ut1[:,:,:,:]
        for ti in range(2, self.dense_time.length):
            Udt = numpy.tensordot(Ut1, Udt)
        return Udt

        
    def _calculate_remainig_using_first_interval(self, Nt,
                                                 show_progress=False):
        """Calculate the rest of the superoperator with known first interval
        
        
        """
        Udt = self.data[1,:,:,:,:]
        
        for ti in range(2, Nt):
            if show_progress:
                print("Self propagation: ", ti, "of", Nt)            
            self.data[ti,:,:,:,:] = \
                numpy.tensordot(Udt, self.data[ti-1,:,:,:,:])        


    def calculate_next(self, save=False):
        """Calculates one point of data of the superopetor
        
        """
        
        if self.mode != "jit":
            raise Exception("This method (calculate_next()) can be used only"+
                            " with mode='jit'")
            
        Nt = self.time.length

        if (self.pdeph is not None) and (self.pdeph.dtype == "Gaussian"):


            if self.now == 0:

                #
                # Create new data
                #    
                self._initialize_data(save=save)
                            
            #
            # We calculate every interval completely
            #
            ti = self.now + 1

            t0 = self.time.data[ti-1] # initial time of the propagation
            
            Ut1 = self._elemental_step_TimeDependent(t0)
            
            if save:
                self.data[ti, :,:,:,:] = \
                numpy.tensordot(Ut1, self.data[ti-1,:,:,:,:])
            else:
                self.data[:,:,:,:] = \
                    numpy.tensordot(Ut1, self.data[:,:,:,:])
              
            self.now += 1
                  
        else:
            
            #
            # Propagation in the first interval
            #
            if self.now == 0:
                #
                # Create new data
                #  
                self._initialize_data(save=save)

                t0 = 0.0

                self.Udt = self._one_step_with_dense_TimeIndep(t0,
                                                    self.dense_time.length,
                                                    self.dense_time.step, Nt,
                                                    show_progress=False) 
                
                if save:
                    self.data[1,:,:,:,:] = self.Udt[:,:,:,:]
                else:
                    self.data[:,:,:,:] = self.Udt[:,:,:,:]
                
                self.now += 1
            
            #
            # Propagations of the later intervals
            #
            else:
                
                ti = self.now + 1
                
                if save:
                    self.data[ti, :,:,:,:] = \
                        numpy.tensordot(self.Udt, self.data[ti-1,:,:,:,:])
                else:
                    self.data[:,:,:,:] = \
                        numpy.tensordot(self.Udt, self.data[:,:,:,:])
                
                self.now += 1


    def at(self, time=None):
        """Retruns evolution superoperator tensor at a given time
        
        
        Parameters
        ----------
        
        time : float, None
            Time (in fs) at which the tensor should be returned. If time is
            None, the whole data object is returned
            
            
        """

        if time is not None:
            ti, dt = self.time.locate(time)

            return SuperOperator(data=self.data[ti, :, :, :, :])
        else:
            return SuperOperator(data=self.data)

          
    def apply(self, time, target, copy=True):
        """Applies the evolution superoperator at a given time
        
        
        Parameters
        ----------
        
        time : float, array (list, tupple) of floats or TimeAxis
            Time(s) at which the evolution superoperator should be applied
            
        target : DensityMatrix, ReducedDensityMatrix
            Operator which the evolution superoperator should be applied
            
        copy : bool
            If True, the target object is copied and new value is assigned
            to its `data` attribute. If False, we assign the new values to
            the target object itself and no copying occurs.


        """

        if isinstance(time, numbers.Real):
            
            #
            # Apply at a single point in time and return ReducedDensityMatrix
            #
            ti, dt = self.time.locate(time)
            if copy:
                import copy
                oper_ven = copy.copy(target)
                oper_ven.data = numpy.tensordot(self.data[ti, :, :, :, :],
                                                target.data)
                return oper_ven
            else:
                target.data = numpy.tensordot(self.data[ti, :, :, :, :],
                                              target.data)
                return target
            
        else:
            # Probably evaluation at more than one time
            # here, `copy` parameter is irrelevant
            
            if isinstance(time, str) or (id(time) == id(self.time)):
                
                #
                # Either we say time="all" or time= exactly the time axis
                # of the evolution superoperator
                #
                # we apply the tensor at all its points
                #
                
                if isinstance(time, str):
                    if time != "all":
                        raise Exception("When argument time is a string, "+
                                        "it must be equal to 'all'")

                rhot = ReducedDensityMatrixEvolution(timeaxis=self.time,
                                                     rhoi=target)
                k_i = 0
                for tt in time.data:
                    rhot.data[k_i,:,:] = \
                    numpy.tensordot(self.data[k_i,:,:,:,:],
                                    target.data)
                    k_i += 1
                
                return rhot

            elif isinstance(time, (list, numpy.array, tuple, TimeAxis)):
                
                #
                # we apply at points specified by TimeAxis
                #
                if isinstance(time, TimeAxis):
                    ntime = time
                else:
                    length = len(time)
                    dt = time[1]-time[0]
                    t0 = time[0]
                    ntime = TimeAxis(t0, length, dt)
                
                rhot = ReducedDensityMatrixEvolution(timeaxis=ntime,
                                                     rhoi=target)
                
                k_i = 0
                for tt in ntime.data:
                    Ut = self.at(tt)
                    rhot.data[k_i,:,:] = numpy.tensordot(Ut.data, target.data)
                    k_i += 1
                    
                return rhot
            
            else:
                raise Exception("Invalid argument: time")


    def plot_element(self, elem, part="REAL", show=True):
        """Plots a selected element of the evolution superoperator
        
        
        Parameters
        ----------
        
        elem : tuple
            A tuple of indices determing the element of the superoperator
            
        """
        
        shape = self.data.shape
        tl = self.time.length
        if (len(elem) == len(shape)-1) and (len(elem) == 4):

            if tl == shape[0]:
                if part == "REAL":
                    dat1 = numpy.real(self.data[:,elem[0],elem[1],elem[2],elem[3]])
                    dat2 = None
                elif part == "IMAG":
                    dat1 = numpy.imag(self.data[:,elem[0],elem[1],elem[2],elem[3]])
                    dat2 = None
                elif part == "BOTH":
                    dat1 = numpy.real(self.data[:,elem[0],elem[1],elem[2],elem[3]])
                    dat2 = numpy.imag(self.data[:,elem[0],elem[1],elem[2],elem[3]])
                else:
                    raise Exception("Unknown data part: "+part)
                    
                plt.plot(self.time.data, dat1)
                if dat2 is not None:
                    plt.plot(self.time.data, dat2)
                    
                if show:
                    plt.show()
        
        else:
            print("Nothing to plot")
            
          
    def get_element_fft(self, elem, window=None):
        """Returns a DFunction with the FFT of the element evolution
        
        """

        if window is None:
            winfce = qr.DFunction(self.time, 
                               numpy.ones(self.time.length, dtype=qr.REAL))
        else:
            winfce = window
            
        dat = self.data[:,elem[0],elem[1],elem[2],elem[3]]
        
        fdat = numpy.fft.ifft(dat*winfce.data)
        fdat = numpy.fft.fftshift(fdat)
        
        time = self.time
        freq = time.get_FrequencyAxis()
        
        ffce = qr.DFunction(freq, fdat)
        
        return ffce
    

    # FIXME: In principle, we can define Greens function and return it
    def get_fft(self, window=None, subtract_last=True):
        """Returns Fourier transform of the whole evolution superoperator
        
        
        Parameters
        ----------
        
        window: DFunction or numpy array
            Windowing function by which the data are multiplied
            
        subtract_last: bool
            If True, the value at the last available time is subtracted 
            from all times
        
        """

        if window is None:
            winfce = qr.DFunction(self.time, 
                               numpy.ones(self.time.length, dtype=qr.REAL))
        else:
            winfce = window
            
        dat = self.data
        if subtract_last:
            dat = dat - dat[-1,:,:,:,:]
        
        # FIXME: Implement windowing
        fdat = numpy.fft.ifft(dat, axis=0) #*winfce.data)
        fdat = numpy.fft.fftshift(fdat, axes=0)
        
        time = self.time
        freq = time.get_FrequencyAxis()
        
        return fdat, freq        
        

    #
    # Calculation `progressbar`
    #

    def _estimate_remaining_loops(self, Nt, dim, ti, n, m):
        return (dim-n)*dim 
        
    
    def _init_progress(self):
        self.ccount = 0
        self.strtime = time.time()
        self.oldtime = self.strtime
        self.remlast = 0
        
        
    def _progress(self, Nt, dim, ti, n, m):
        
        curtime = time.time()
        #dt = curtime - self.oldtime
        Dt = curtime - self.strtime
        self.oldtime = curtime
        self.ccount += 1
         
        dt_past = Dt/(self.ccount*dim)
        
        remloops = self._estimate_remaining_loops(Nt, dim, ti, n, m)
        
        remtime = int(remloops*(dt_past))
        
        txt1 = "Propagation cycle "+str(ti)+" of "+str(Nt)+ \
               " : ("+str(n)+" of "+str(dim)+")"
        if True:
        #if remtime <= self.remlast:
            txt2 = " : remaining time "+str(remtime)+" sec"
        #else:
        #    txt2 = " : remaining time ... calculating"
        tlen1 = len(txt1)
        tlen2 = len(txt2)
        
        rem = 65 - tlen1 - tlen2
        txt = "\r"+txt1+(" "*rem)+txt2
        print(txt, end="\r")
        
        self.remlast = remtime


    def __str__(self):
        out  = "\nquantarhei.EvolutionSuperOperator object"
        out += "\n========================================"
#        out += "\nunits of energy %s" % self.unit_repr()
#        out += "\nRotating Wave Approximation (RWA) enabled : "\
#            +str(self.has_rwa)
#        if self.has_rwa:
#            out += "\nNumber of blocks : "+str(self.Nblocks)
#            out += "\nBlock average energies:"
#            for k in range(self.Nblocks):
#                out += "\n "+str(k)+" : "\
#                +str(self.rwa_energies[self.rwa_indices[k]])
        #out += "\ndata = \n"
        #out += str(self.data)
        return out                        