# -*- coding: utf-8 -*-
"""
*******************************************************************************

    REDUCED DENSITY MATRIX PROPAGATOR

*******************************************************************************

Propagator of the Reduced density matrix dynamics

Created on Mon Apr 13 18:27:07 2015

author: Tomas Mancal, Charles University

"""

import numpy
import numpy.linalg

from ..hilbertspace.hamiltonian import Hamiltonian
from ..hilbertspace.operators import Operator
from ...core.time import TimeAxis
from ...core.time import TimeDependent
from ...core.saveable import Saveable
from ..liouvillespace.redfieldtensor import RelaxationTensor
from ..hilbertspace.operators import ReducedDensityMatrix, DensityMatrix
from .dmevolution import ReducedDensityMatrixEvolution
from ...core.matrixdata import MatrixData
from ...core.managers import Manager

import quantarhei as qr

import matplotlib.pyplot as plt

    
_show_debug = False
def debug(msg):
    
    if _show_debug:
        print(msg)

class ReducedDensityMatrixPropagator(MatrixData, Saveable): 
    """
    
    Reduced Density Matrix Propagator calculates the evolution of the
    reduced density matrix based on the Hamiltonian and optionally
    a relaxation tensor. Relaxation tensor may be constant or time
    dependent. 
    
    """
    
    def __init__(self, timeaxis=None, Ham=None, RTensor=None, Iterm=None,
                 Efield=None, Trdip=None, PDeph=None):
        """
        
        Creates a Reduced Density Matrix propagator which can propagate
        a given initial density matrix with the given Hamiltonian and 
        relaxation tensor.Time axis of the propagation is specied by 
        the second argument.
        
        RWA : tuple
            A tuple of block starting indices. It is assumed that if RWA is
            to be used, that there is at least a ground state and a single
            block of excited state starting from the Nth energy level. In this
            case the tuple reads (0, N). The tuple always starts with 0. In
            case there are two excited blocks, it reads (0, N1, N2) where
            N2 > N1 > 0.
        
        
        Examples
        --------
        
        The constructor accepts only numpy.ndarray object, so the
        following code will fail, because it submits normal Python arrays.         
        
        >>> pr = ReducedDensityMatrixPropagator([[0.0, 0.0],
        ...                                     [0.0, 1.0]],[0,1,2,3,4,5])
        Traceback (most recent call last):
            ...
        Exception: TimeAxis expected here.
        
        The correct way to construct the propagator is the following:

        >>> h = numpy.array([[0.0, 0.0],[0.0,1.0]])
        >>> HH = Hamiltonian(data=h)
        >>> times = TimeAxis(0,1000,1)
        >>> pr = ReducedDensityMatrixPropagator(times,HH)
        
        """
        self.has_Trdip = False
        self.has_Efield = False
        self.has_PDeph = False
        self.has_RTensor = False
        self.has_Iterm = False
        self.has_RWA = False
        self.has_EField = False
        
        if not ((timeaxis is None) and (Ham is None)):
            
            #
            # Hamiltonian and TimeAxis are required
            #
            if isinstance(timeaxis, TimeAxis):
                self.TimeAxis = timeaxis
            else:
                raise Exception("TimeAxis expected here.")

            if isinstance(Ham, Hamiltonian):
                self.Hamiltonian = Ham
            elif Ham is None:
                raise Exception("Hamiltonian is required.")
            else:
                raise Exception("Hamiltonian represented by a wrong type.")
            
            #
            # RelaxationTensor is not requited
            #
            if isinstance(RTensor, RelaxationTensor):
                self.RelaxationTensor = RTensor
                self.has_RTensor = True
                self.has_relaxation = True
            elif RTensor is None:
                self.has_RTensor = False
                self.has_relaxation = False
            else:
                raise Exception("RelaxationTensor or None expected here.")
            
            #
            # TransitionDipoleMoment is not required
            #
            if Trdip is not None:            
                if isinstance(Trdip, Operator):
                    self.Trdip = Trdip
                    self.has_Trdip = True
                else:
                    raise Exception("Operator or None expected here.")

            #
            # Driving field is not required
            #    
            if Efield is not None:
                if isinstance(Efield, numpy.ndarray):
                    self.Efield = Efield
                    self.has_Efield = True
                    self.has_EField = False
                else: 
                    self.EField = Efield
                    self.has_EField = True
                    self.has_Efield = False                    
            
            #
            # Pure dephasing also counts as relaxation; not required
            #
            if PDeph is not None:
                self.PDeph = PDeph
                self.has_PDeph = True
                self.has_relaxation = True
            
            #
            # Initial term is not required; has to come with relaxation
            #
            if Iterm is not None:
                if self.has_relaxation:
                    self.Iterm = Iterm
                    self.has_Iterm = True
                else:
                    raise Exception("RelaxationTensor has to be set first.")
            
            
            self.Odt = self.TimeAxis.data[1]-self.TimeAxis.data[0]
            self.dt = self.Odt
            self.Nref = 1
            
            self.Nt = self.TimeAxis.data.shape[0]
            
            N = self.Hamiltonian.data.shape[0]
            self.N = N
            self.data = numpy.zeros((self.Nt,N,N),dtype=numpy.complex64)
        
            self.verbose = Manager().log_conf.verbose
            
        else:
            
            raise Exception("TimeAxis and Hamiltonian are required to"+
                            " initialize the propagator.")

        
    def setDtRefinement(self, Nref):
        """
        The TimeAxis object specifies at what times the propagation
        should be stored. We can tell the propagator to use finer
        time step for the calculation by setting the refinement. The
        refinement is an integer by which the TimeAxis time step should
        be devided to get the finer time step. In the code below, we
        have dt = 10 in the TimeAxis, but we want to calculate with
        dt = 1
        
        >>> HH = Hamiltonian(data=numpy.array([[0.0, 0.0],[0.0,1.0]]))
        >>> times = TimeAxis(0,1000,10.0)
        >>> pr = ReducedDensityMatrixPropagator(times, HH)
        >>> pr.setDtRefinement(10)
        
        """
        self.Nref = Nref
        self.dt = self.Odt/self.Nref
        
        
    def propagate(self, rhoi, method="short-exp", mdata=None, Nref=1):
        """
        
        >>> T0   = 0
        >>> Tmax = 100
        >>> dt   = 1
        >>> Nref = 30


        >>> initial_dm = [[1.0, 0.0, 0.0],
        ...               [0.0, 0.0, 0.0],
        ...               [0.0, 0.0, 0.0]]

        >>> Ham_matrix = [[0.0, 0.1, 0.0],
        ...                [0.1, 0.0, 0.1],
        ...                [0.0, 0.1, 0.1]]        
        
        >>> HH = Hamiltonian(data=Ham_matrix)
        >>> times = TimeAxis(T0,Tmax,dt)
        
        >>> pr = ReducedDensityMatrixPropagator(times, HH)
        >>> pr.setDtRefinement(Nref)

        >>> rhoi = ReducedDensityMatrix(data=numpy.array(initial_dm))  
        >>> rhot=pr.propagate(rhoi)
        
        Refinement of the time-step can also be set when calling 
        the propagate() method.
        
        >>> pr = ReducedDensityMatrixPropagator(times, HH) 
        >>> rhot=pr.propagate(rhoi, Nref=Nref) 
        
        First argument must be ReducedDensityMatrix
        
        >>> rhot = pr.propagate(initial_dm)
        Traceback (most recent call last):
            ...
        Exception: First argument has be of the ReducedDensityMatrix type
        """
        
        if Nref > 1:
            self.setDtRefinement(Nref)
        
        #
        # Testing if the object submitted is density matrix
        #
        if not (isinstance(rhoi, ReducedDensityMatrix) 
             or isinstance(rhoi, DensityMatrix)):
            raise Exception("First argument has be of"+
            " the ReducedDensityMatrix type")
              
        #######################################################################
        #
        #    PROPAGATIONS WITH RELAXATION AND/OR DEPHASING
        #
        #
        #######################################################################
        if self.has_relaxation:

            ###################################################################
            #
            # Time-dependent relaxation tensor
            #
            ###################################################################
            if isinstance(self.RelaxationTensor, TimeDependent):

                ###############################################################
                #
                # Propagation with external field
                #
                ###############################################################
                if (self.has_Efield and self.has_Trdip):
                    
                    if method == "short-exp": 
                        return \
                        self.__propagate_short_exp_with_TD_relaxation_field(\
                        rhoi,L=4)
                    elif method == "short-exp-2":
                        return \
                        self.__propagate_short_exp_with_TD_relaxation_field(\
                        rhoi,L=2)
                    elif method == "short-exp-4":
                        return \
                        self.__propagate_short_exp_with_TD_relaxation_field(\
                        rhoi,L=4)
                    elif method == "short-exp-6":
                        return \
                        self.__propagate_short_exp_with_TD_relaxation_field(\
                        rhoi,L=6)            
                    else:
                        raise Exception("Unknown propagation method: "+method)

                elif (self.has_EField and self.has_Trdip):
                    
                    if method == "short-exp": 
                        return \
                        self.__propagate_short_exp_with_TD_relaxation_EField(\
                        rhoi,L=4)
                    elif method == "short-exp-2":
                        return \
                        self.__propagate_short_exp_with_TD_relaxation_EField(\
                        rhoi,L=2)
                    elif method == "short-exp-4":
                        return \
                        self.__propagate_short_exp_with_TD_relaxation_EField(\
                        rhoi,L=4)
                    elif method == "short-exp-6":
                        return \
                        self.__propagate_short_exp_with_TD_relaxation_EField(\
                        rhoi,L=6)            
                    else:
                        raise Exception("Unknown propagation method: "+method)

                    
                ###############################################################
                #
                # Progation without external field
                #
                ###############################################################                    
                else:

                    if method == "short-exp":
                        return self.__propagate_short_exp_with_TD_relaxation(\
                        rhoi,L=4)
                    elif method == "short-exp-2":
                        return self.__propagate_short_exp_with_TD_relaxation(\
                        rhoi,L=2)
                    elif method == "short-exp-4":
                        return self.__propagate_short_exp_with_TD_relaxation(\
                        rhoi,L=4)
                    elif method == "short-exp-6":
                        return self.__propagate_short_exp_with_TD_relaxation(\
                        rhoi,L=6)            
                    else:
                        raise Exception("Unknown propagation method: "+method)

            ###################################################################
            #
            # Constant relaxation tensor
            #
            ###################################################################
            else: 

                ###############################################################
                #
                # Propagation with external field
                #
                ###############################################################
                if (self.has_Efield and self.has_Trdip):
                    
                    if method == "short-exp":
                        return \
                        self.__propagate_short_exp_with_relaxation_field(
                        rhoi,L=4)
                    elif method == "short-exp-2":
                        return \
                        self.__propagate_short_exp_with_relaxation_field(
                        rhoi,L=2)
                    elif method == "short-exp-4":
                        return \
                        self.__propagate_short_exp_with_relaxation_field(
                        rhoi,L=4)
                    elif method == "short-exp-6":
                        return \
                        self.__propagate_short_exp_with_relaxation_field(
                        rhoi,L=6)            
                    else:
                        raise Exception("Unknown propagation method: "+method)                
                        
                elif (self.has_EField and self.has_Trdip):
                    
                    if method == "short-exp":
                        return \
                        self.__propagate_short_exp_with_relaxation_EField(
                        rhoi,L=4)
                    elif method == "short-exp-2":
                        return \
                        self.__propagate_short_exp_with_relaxation_EField(
                        rhoi,L=2)
                    elif method == "short-exp-4":
                        return \
                        self.__propagate_short_exp_with_relaxation_EField(
                        rhoi,L=4)
                    elif method == "short-exp-6":
                        return \
                        self.__propagate_short_exp_with_relaxation_EField(
                        rhoi,L=6)            
                    else:
                        raise Exception("Unknown propagation method: "+method)
                        
                ###############################################################
                #
                # Progation without external field
                #
                ###############################################################                    
                else:

                    #
                    #  Section used by Time-independent Redfield and similar
                    #

                    if method == "short-exp":
                        return self.__propagate_short_exp_with_relaxation(
                        rhoi,L=4)
                    elif method == "short-exp-2":
                        return self.__propagate_short_exp_with_relaxation(
                        rhoi,L=2)
                    elif method == "short-exp-4":
                        return self.__propagate_short_exp_with_relaxation(
                        rhoi,L=4)
                    elif method == "short-exp-6":
                        return self.__propagate_short_exp_with_relaxation(
                        rhoi,L=6)            

                    else:
                        raise Exception("Unknown propagation method: "+method)   
            
        #######################################################################
        #
        #    PROPAGATIONS WITHOUT RELAXATION
        #
        #
        #######################################################################
        else:

            if (self.has_Efield and self.has_Trdip):   

                if method == "short-exp":
                    return self.__propagate_short_exp_efield(rhoi,L=4)
                elif method == "short-exp-2":
                    return self.__propagate_short_exp_efield(rhoi,L=2)
                elif method == "short-exp-4":
                    return self.__propagate_short_exp_efield(rhoi,L=4)
                elif method == "short-exp-6":
                    return self.__propagate_short_exp_efield(rhoi,L=6)            
    
                else:
                    raise Exception("Unknown propagation method: "+method)                

            elif (self.has_EField and self.has_Trdip):   

                if method == "short-exp":
                    return self.__propagate_short_exp_EField(rhoi,L=4)
                elif method == "short-exp-2":
                    return self.__propagate_short_exp_EField(rhoi,L=2)
                elif method == "short-exp-4":
                    return self.__propagate_short_exp_EField(rhoi,L=4)
                elif method == "short-exp-6":
                    return self.__propagate_short_exp_EField(rhoi,L=6)            
    
                else:
                    raise Exception("Unknown propagation method: "+method)
            else:
                 
                if method == "short-exp":
                    return self.__propagate_short_exp(rhoi,L=4)
                elif method == "short-exp-2":
                    return self.__propagate_short_exp(rhoi,L=2)
                elif method == "short-exp-4":
                    return self.__propagate_short_exp(rhoi,L=4)
                elif method == "short-exp-6":
                    return self.__propagate_short_exp(rhoi,L=6)            
    
                else:
                    raise Exception("Unknown propagation method: "+method)
        
 

    def __propagate_short_exp(self, rhoi, L=4):
        """Short expansion of an exponention to integrate equations of motion
        
        
        Propagation with Hamiltonian only
        
        
        """
        
        debug("(1)")

        (pr, rho1, rho2) = self._INIT_EXP(rhoi)
        
        HH = self._INIT_RWA()
        
        indx = 1
        for ii in self.TimeAxis.data[1:self.Nt]:
            
            for jj in range(0,self.Nref):
                
                for ll in range(1,L+1):
                   
                    rho1 = -_COM(HH, ll, self.dt, rho1)
                             
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pr.data[indx,:,:] = rho2                        
            indx += 1                       
            
        self._CLOSE_RWA(pr)
            
        return pr


    def _INIT_EXP(self, rhoi):
        """

        Parameters
        ----------
        rhoi : TYPE
            DESCRIPTION.

        Returns
        -------
        pr : TYPE
            DESCRIPTION.
        rho1 : TYPE
            DESCRIPTION.
        rho2 : TYPE
            DESCRIPTION.

        """
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        rho1 = rhoi.data
        rho2 = rhoi.data   
        
        return (pr, rho1, rho2)


    def _INIT_RWA(self):
        """Create Hamiltonian matrix respecting RWA settings

        Parameters
        ----------
        Ham : TYPE
            DESCRIPTION.

        Returns
        -------
        HH : TYPE
            DESCRIPTION.

        """
        if self.Hamiltonian.has_rwa:
            HH = self.Hamiltonian.get_RWA_data()
        else:
            HH = self.Hamiltonian.data            
        return HH


    def _CLOSE_RWA(self, pr):
        """Closes the RWA treatment

        Parameters
        ----------
        Ham : TYPE
            DESCRIPTION.
        pr : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True  


    def __propagate_short_exp_efield(self, rhoi,L=4): 

        raise Exception("NOT IMPLEMENTED")
        

    def __propagate_short_exp_EField(self, rhoi,L=4):
        
        raise Exception("NOT IMPLEMENTED")
                
 
    def __propagate_short_exp_with_relaxation(self, rhoi, L=4):
        """Integration by short exponentional expansion
        
        Integration by expanding exponential to Lth order. Time independent
        relaxation tensor
              
        """
        
        if self.RelaxationTensor.as_operators:
            return self.__propagate_short_exp_with_rel_operators(rhoi, L=L)

        
        debug("(2)")
        
        (pr, rho1, rho2) = self._INIT_EXP(rhoi)
        
        HH = self._INIT_RWA()
            
        RR = self.RelaxationTensor.data

        #
        #  Propagation with aditional pure dephasing
        #
        if self.has_PDeph:
            
            self._BOOT_DEPH()
            
            IR = 0.0
            indx = 1
            for ii in range(1, self.Nt): 

                # time at the beginning of the step
                tNt = self.TimeAxis.data[indx-1]  

                # if self.has_Iterm:
                #     IR = self.RelaxationTensor.Iterm[indx,:,:] 
                IR = self._GET_IR(indx)
                
                for jj in range(0, self.Nref):
                    
                    tt = tNt + jj*self.dt  # time right now 
 
                    for ll in range(1, L+1):

                        rhoY = -_COM(HH, ll, self.dt, rho1)
                        _TTI(rhoY, RR, IR, ll, self.dt, rho1)
                        
                        rho1 = rhoY                        
                        rho2 = rho2 + rho1
                        
                    rho2 = self._APPLY_DEPH(tt, rho2)
                        
                    rho1 = rho2    
                    
                pr.data[indx,:,:] = rho2 
                indx += 1   
        #
        #  Standard propagation with Hamiltonian and relaxation
        #        
        else:
            
            #IR = 0.0
            indx = 1
            for ii in range(1, self.Nt): 
                
                # if self.has_Iterm:
                #     IR = self.RelaxationTensor.Iterm[indx,:,:]
                IR = self._GET_IR(indx)

                for jj in range(0, self.Nref):
                    
                    for ll in range(1, L+1):
                        
                        rhoY = -_COM(HH, ll, self.dt, rho1)
                        _TTI(rhoY, RR, IR, ll, self.dt, rho1)
                        
                        rho1 = rhoY                                 
                        rho2 = rho2 + rho1
                    rho1 = rho2    
                    
                pr.data[indx,:,:] = rho2 
                indx += 1   

        self._CLOSE_RWA(pr)    
        
        return pr  


    def _GET_IR(self, indx):
        """Returns in value of the initial term at a give time

        Parameters
        ----------
        indx : int
            Index representing the time.

        Returns
        -------
        IR : numpy.array, 2 indices
            The value of the initial term for a given time.

        """
        if self.has_Iterm:
            IR = self.RelaxationTensor.Iterm[indx,:,:] 
        else:
            IR = 0.0
        return IR


    def _BOOT_DEPH(self):
        """
        

        Returns
        -------
        None.

        """
        if self.PDeph.dtype == "Lorentzian":
            self.expo = numpy.exp(-self.PDeph.data*self.dt)
            self.t0 = 0.0
        elif self.PDeph.dtype == "Gaussian":
            self.expo = numpy.exp(-self.PDeph.data*(self.dt**2)/2.0)
            self.t0 = self.PDeph.data*self.dt


    def _APPLY_DEPH(self, tt, rho2):
        """
        

        Parameters
        ----------
        tt : TYPE
            DESCRIPTION.
        rho2 : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return rho2*self.expo*numpy.exp(-self.t0*tt)
      
                
    def __propagate_short_exp_with_rel_operators(self, rhoi, L=4):
        """Integration by short exponentional expansion
        
        Integration by expanding exponential to Lth order. 
              
            
        """
        
        debug("(3)")

        (pr, rho1, rho2) = self._INIT_EXP(rhoi)
        # pr = ReducedDensityMatrixEvolution(self.TimeAxis, rhoi)
        
        # rho1 = rhoi.data
        # rho2 = rhoi.data
        
        #
        # RWA is applied here
        #
        # if self.Hamiltonian.has_rwa:
        #     HH = self.Hamiltonian.get_RWA_data() #data  - self.HOmega
        # else:
        #     HH = self.Hamiltonian.data
        HH = self._INIT_RWA()
        
        qr.log_detail("PROPAGATION (short exponential with "+
                     "relaxation in operator form): order ", L, 
                     verbose=self.verbose)
        qr.log_detail("Using complex numpy implementation")
        
        Km = self.RelaxationTensor.Km # real
        Lm = self.RelaxationTensor.Lm # complex
        Ld = self.RelaxationTensor.Ld # complex - get by transposition
        Kd = numpy.zeros(Km.shape, dtype=numpy.float64)
        Nm = Km.shape[0]
        for m in range(Nm):
            Kd[m, :, :] = numpy.transpose(Km[m, :, :])
            
        indx = 1

        levs = [qr.LOG_QUICK] #, 8]
        verb = qr.loglevels2bool(levs)

        # after each step we apply pure dephasing (if present)
        if self.has_PDeph:
            
            # if self.PDeph.dtype == "Lorentzian":
            #     expo = numpy.exp(-self.PDeph.data*self.dt)
            #     t0 = 0.0
            # elif self.PDeph.dtype == "Gaussian":
            #     expo = numpy.exp(-self.PDeph.data*(self.dt**2)/2.0)
            #     t0 = self.PDeph.data*self.dt
            self._BOOT_DEPH()
            
            # loop over time
            for ii in range(1, self.Nt):
                qr.printlog(" time step ", ii, "of", self.Nt, 
                            verbose=verb[0], loglevel=levs[0])
                
                # time at the beginning of the step
                tNt = self.TimeAxis.data[indx-1]  
                #print("tNt = ", tNt)
                
                # steps in between saving the results
                for jj in range(0, self.Nref):
                    
                    tt = tNt + jj*self.dt  # time right now 
                    
                    # L interations to get short exponential expansion
                    for ll in range(1, L+1):
                        
                        rhoY = -_COM(HH, ll, self.dt, rho1)
                        # rhoY =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) 
                        #                         - numpy.dot(rho1,HH))
                        
                        _OTI(rhoY, Km, Kd, Lm, Ld, ll, self.dt, rho1)
                        # for mm in range(Nm):
                            
                        #    rhoY += (self.dt/ll)*(
                        #     numpy.dot(Km[mm,:,:],numpy.dot(rho1, Ld[mm,:,:]))
                        #    +numpy.dot(Lm[mm,:,:],numpy.dot(rho1, Kd[mm,:,:]))
                        #    -numpy.dot(numpy.dot(Kd[mm,:,:],Lm[mm,:,:]), rho1)
                        #    -numpy.dot(rho1, numpy.dot(Ld[mm,:,:],Km[mm,:,:]))
                        #    )
                                 
                        rho1 = rhoY #+ rhoX
                        
                        rho2 = rho2 + rho1
                       
                    # pure dephasing is added here                        
                    #rho2 = rho2*expo*numpy.exp(-t0*tt)
                    rho2 = self._APPLY_DEPH(tt, rho2)
                        
                    rho1 = rho2    
                
                pr.data[indx,:,:] = rho2 
                indx += 1
            
        # no extra dephasing
        else:
            
             # loop over time
            for ii in range(1, self.Nt):
                qr.printlog(" time step ", ii, "of", self.Nt, 
                            verbose=verb[0], loglevel=levs[0])
                
                # steps in between saving the results
                for jj in range(0, self.Nref):
                    
                    # L interations to get short exponential expansion
                    for ll in range(1, L+1):
                        
                        rhoY = -_COM(HH, ll, self.dt, rho1)
                        # rhoY =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) 
                        #                          - numpy.dot(rho1,HH))
                        
                        _OTI(rhoY, Km, Kd, Lm, Ld, ll, self.dt, rho1)
                        # for mm in range(Nm):
                            
                        #    rhoY += (self.dt/ll)*(
                        #     numpy.dot(Km[mm,:,:],numpy.dot(rho1, Ld[mm,:,:]))
                        #    +numpy.dot(Lm[mm,:,:],numpy.dot(rho1, Kd[mm,:,:]))
                        #    -numpy.dot(numpy.dot(Kd[mm,:,:],Lm[mm,:,:]), rho1)
                        #    -numpy.dot(rho1, numpy.dot(Ld[mm,:,:],Km[mm,:,:]))
                        #    )
                                 
                        rho1 = rhoY #+ rhoX
                        
                        rho2 = rho2 + rho1
                    
                    rho1 = rho2    
                
                pr.data[indx,:,:] = rho2 
                indx += 1
           
             
        qr.log_detail("...DONE")

        #if self.Hamiltonian.has_rwa:
        #    pr.is_in_rwa = True
        self._CLOSE_RWA(pr)
            
        return pr



        

    def __propagate_short_exp_with_TD_relaxation(self,rhoi,L=4):
        """
              Short exp integration
              
              This is the propagator used with the Time-dependent relaxation
              tensors
              
              
        """
    

        if self.RelaxationTensor.as_operators:
            return self.__propagate_short_exp_with_TDrel_operators(rhoi,L=L)
                                                                       
        debug("(4)")
        
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        
        rho1 = rhoi.data
        rho2 = rhoi.data

        #
        # RWA is applied here
        #
        if self.Hamiltonian.has_rwa:
            HH = self.Hamiltonian.get_RWA_data()
        else:
            HH = self.Hamiltonian.data
       
        #
        # set cut-off index by the tensor cut-off time
        #
        if self.RelaxationTensor._has_cutoff_time:
            cutoff_indx = \
            self.TimeAxis.nearest(self.RelaxationTensor.cutoff_time)
        else:
            sbi = self.RelaxationTensor.SystemBathInteraction
            cutoff_indx = sbi.TimeAxis.length
            
        indx = 1
        indxR = 1
        
        # check if the relaxation tensor requires inhomogeneous term
        
        try:
            self.has_Iterm = self.RelaxationTensor.has_Iterm
        except:
            print(self, "This tensor does not support has_Iterm")
            self.has_Iterm = False
            
        if self.has_Iterm:
            self.RelaxationTensor.initial_term(rhoi)
            
            
        #
        # Propagate with time-dependent relaxation tensor and inhomogeneous
        # term if there is one.
        #
        sysstep = self.RelaxationTensor.SystemBathInteraction.TimeAxis.step
        Nref_max = round(self.TimeAxis.step/sysstep)
        Nref_req = self.Nref
        
        if Nref_max % Nref_req == 0:
            
            stride = Nref_max//Nref_req
            
        else:
            time = self.TimeAxis
            print("Relaxation tensor has a timestep  :", sysstep, "fs")
            print("Propagation timestep is:", time.step, "fs")
            print("ERROR:", Nref_req,"steps of", sysstep,"fs do not fit"+
                  " neatly into a step of", time.step,"fs")
            raise Exception("Incompatible number of refinement steps")
            
        IR = 0.0 
        dt = sysstep*stride
        for ii in self.TimeAxis.data[1:self.Nt]:
            
            for jj in range(self.Nref):
                
                
                RR = self.RelaxationTensor.data[indxR,:,:,:,:]
                if self.has_Iterm:
                    IR = self.RelaxationTensor.Iterm[indxR,:,:]                           
                
                for ll in range(1,L+1):
                    
                    rhoY = -_COM(HH, ll, dt, rho1)
                    _TTI(rhoY, RR, IR, ll, dt, rho1)
                    # rho1 =  (dt/ll)*(numpy.tensordot(RR,rho1) \
                    #    - 1j*(numpy.dot(HH,rho1)- numpy.dot(rho1,HH)) + IR)
                    rho1 = rhoY      
                    rho2 = rho2 + rho1
                rho1 = rho2    

                if indxR < cutoff_indx - 1:                      
                    indxR += stride
                else:
                    indxR = cutoff_indx

                
            pr.data[indx,:,:] = rho2
            
            #
            # We respect the tensor cut-off
            #
            indx += 1
  
                
        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True
            
        return pr     


    def __propagate_short_exp_with_TDrel_operators(self, rhoi, L=4):
        """
            Short exp integration with time-dependent relaxation tensor
              
            The tensor is in operator form.
            
            
        """
        import time
        debug("(5)")
        pr = ReducedDensityMatrixEvolution(self.TimeAxis, rhoi)
        
        rho1 = rhoi.data
        rho2 = rhoi.data
        
        #HH = self.Hamiltonian.data  
        #
        # RWA is applied here
        #
        # FIXME: RWA has to be applied to ralaxation tensor, too!!!!
        #
        if self.Hamiltonian.has_rwa:
            HH = self.Hamiltonian.get_RWA_data() #data  - self.HOmega
        else:
            HH = self.Hamiltonian.data
            
        if self.RelaxationTensor._has_cutoff_time:
            cutoff_indx = \
            self.TimeAxis.nearest(self.RelaxationTensor.cutoff_time)
        else:
            cutoff_indx = self.TimeAxis.length

        Km = self.RelaxationTensor.Km
        Kd = numpy.zeros(Km.shape, dtype=numpy.float64)
        Nm = Km.shape[0]
        for m in range(Nm):
            Kd[m, :, :] = numpy.transpose(Km[m, :, :])
                        
        indx = 1
        indxR = 1
        print("TEST")
        t1 = time.time()
        for ii in range(1, self.Nt): 

            Lm = self.RelaxationTensor.Lm[indxR,:,:,:]
            Ld = self.RelaxationTensor.Ld[indxR,:,:,:]
       
            for jj in range(0, self.Nref):
                
                for ll in range(1, L+1):
                    
                    rhoY =  - _COM(HH, ll, self.dt,rho1) 
                    
                    #(1j*self.dt/ll)*(numpy.dot(HH,rho1) 
                    #                         - numpy.dot(rho1,HH))
                    
                    _OTI(rhoY, Km, Kd, Lm, Ld, ll, self.dt, rho1)
                    
                    # for mm in range(Nm):
                        
                    #     rhoY += (self.dt/ll)*(
                    #     numpy.dot(Km[mm,:,:],numpy.dot(rho1, Ld[mm,:,:]))
                    #     +numpy.dot(Lm[mm,:,:],numpy.dot(rho1, Kd[mm,:,:]))
                    #     -numpy.dot(numpy.dot(Kd[mm,:,:],Lm[mm,:,:]), rho1)
                    #     -numpy.dot(rho1, numpy.dot(Ld[mm,:,:],Km[mm,:,:]))
                    #     )
                             
                    rho1 = rhoY # + rhoX
                    
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pr.data[indx,:,:] = rho2 
            indx += 1             
            if indxR < cutoff_indx-1:                      
                indxR += 1             
        t2 = time.time()
        print("Duration:", t2-t1, "s")
        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True

            
        return pr
        
        
    def __propagate_short_exp_with_TD_relaxation_field(self,rhoi,L=4):
        """Short exp integration of the density matrix with external driving
        
        
        
        """
        if self.RelaxationTensor.as_operators:
            return self.__propagate_short_exp_TDrelax_field_oper(rhoi, L=L)

            
        debug("(6)")    
        
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        
        rho1 = rhoi.data
        rho2 = rhoi.data

        #
        # RWA is applied here
        #
        if self.Hamiltonian.has_rwa:
            HH = self.Hamiltonian.get_RWA_data()
        else:
            HH = self.Hamiltonian.data

        #
        # We do not have an information on polarization - we take X as default
        #
        MU = self.Trdip.data[:,:,0]
        
        EE = self.Efield
        
        #
        # set cut-off index by the tensor cut-off time
        #
        if self.RelaxationTensor._has_cutoff_time:
            cutoff_indx = \
            self.TimeAxis.nearest(self.RelaxationTensor.cutoff_time)
        else:
            sbi = self.RelaxationTensor.SystemBathInteraction
            cutoff_indx = sbi.TimeAxis.length
            
        indx = 1
        indxR = 1
        
        # check if the relaxation tensor requires inhomogeneous term
        
        try:
            self.has_Iterm = self.RelaxationTensor.has_Iterm
        except:
            print(self, "This tensor does not support has_Iterm")
            self.has_Iterm = False
            
        if self.has_Iterm:
            self.RelaxationTensor.initial_term(rhoi)
            
            
        #
        # Propagate with time-dependent relaxation tensor and inhomogeneous
        # term if there is one.
        #
        sysstep = self.RelaxationTensor.SystemBathInteraction.TimeAxis.step
        Nref_max = round(self.TimeAxis.step/sysstep)
        Nref_req = self.Nref
        
        if Nref_max % Nref_req == 0:
            
            stride = Nref_max//Nref_req
            
        else:
            time = self.TimeAxis
            print("Relaxation tensor has a timestep  :", sysstep, "fs")
            print("Propagation timestep is:", time.step, "fs")
            print("ERROR:", Nref_req,"steps of", sysstep,"fs do not fit"+
                  " neatly into a step of", time.step,"fs")
            raise Exception("Incompatible number of refinement steps")
            
        IR = 0.0 
        dt = sysstep*stride
        for ii in self.TimeAxis.data[1:self.Nt]:
            
            for jj in range(self.Nref):
                
                
                RR = self.RelaxationTensor.data[indxR,:,:,:,:]
                if self.has_Iterm:
                    IR = self.RelaxationTensor.Iterm[indxR,:,:]                           
                MuE = MU*EE[indx]
                
                for ll in range(1,L+1):
                    
                    rhoY = -_COM(HH, ll, dt, rho1)
                    _TTI(rhoY, RR, IR, ll, dt, rho1)
                    rhoY += _COM(MuE, ll, dt, rho1)
                    # rho1 = (dt/ll)*(numpy.tensordot(RR,rho1) 
                    #      - 1j*(numpy.dot(HH,rho1)-numpy.dot(rho1,HH)) + IR
                    #      + 1j*(numpy.dot(MU,rho1)-numpy.dot(rho1,MU))*EE[indx])
                    rho1 = rhoY        
                    rho2 = rho2 + rho1
                    
                rho1 = rho2    

                if indxR < cutoff_indx - 1:                      
                    indxR += stride
                else:
                    indxR = cutoff_indx

                
            pr.data[indx,:,:] = rho2
            
            #
            # We respect the tensor cut-off
            #
            indx += 1
  
                
        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True
            
        return pr     



    def __propagate_short_exp_TDrelax_field_oper(self, rhoi, L=4):
        """

        """
        debug("(7)")    



    def __propagate_short_exp_with_TD_relaxation_EField(self,rhoi,L=4):
        """Short exp integration of the density matrix with external driving
        
        
        
        """
        if self.RelaxationTensor.as_operators:
            return self.__propagate_short_exp_TDrelax_EField_oper(rhoi, L=L)

            
        debug("(8)")    
        
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        
        rho1 = rhoi.data
        rho2 = rhoi.data

        #
        # RWA is applied here
        #
        if self.Hamiltonian.has_rwa:
            HH = self.Hamiltonian.get_RWA_data()
        else:
            HH = self.Hamiltonian.data

        #
        # We do not have an information on polarization - we take X as default
        #
        MU = self.Trdip.data[:,:,0]
        
        EE = self.EField.field
        
        #
        # set cut-off index by the tensor cut-off time
        #
        if self.RelaxationTensor._has_cutoff_time:
            cutoff_indx = \
            self.TimeAxis.nearest(self.RelaxationTensor.cutoff_time)
        else:
            sbi = self.RelaxationTensor.SystemBathInteraction
            cutoff_indx = sbi.TimeAxis.length
            
        indx = 1
        indxR = 1
        
        # check if the relaxation tensor requires inhomogeneous term
        
        try:
            self.has_Iterm = self.RelaxationTensor.has_Iterm
        except:
            print(self, "This tensor does not support has_Iterm")
            self.has_Iterm = False
            
        if self.has_Iterm:
            self.RelaxationTensor.initial_term(rhoi)
            
            
        #
        # Propagate with time-dependent relaxation tensor and inhomogeneous
        # term if there is one.
        #
        sysstep = self.RelaxationTensor.SystemBathInteraction.TimeAxis.step
        Nref_max = round(self.TimeAxis.step/sysstep)
        Nref_req = self.Nref
        
        if Nref_max % Nref_req == 0:
            
            stride = Nref_max//Nref_req
            
        else:
            time = self.TimeAxis
            print("Relaxation tensor has a timestep  :", sysstep, "fs")
            print("Propagation timestep is:", time.step, "fs")
            print("ERROR:", Nref_req,"steps of", sysstep,"fs do not fit"+
                  " neatly into a step of", time.step,"fs")
            raise Exception("Incompatible number of refinement steps")
            
        IR = 0.0 
        dt = sysstep*stride
        for ii in self.TimeAxis.data[1:self.Nt]:
            
            for jj in range(self.Nref):
                
                
                RR = self.RelaxationTensor.data[indxR,:,:,:,:]
                if self.has_Iterm:
                    IR = self.RelaxationTensor.Iterm[indxR,:,:]                           
                MuE = MU*EE[indx]
                
                for ll in range(1,L+1):
                    
                    rhoY = -_COM(HH, ll, dt, rho1)
                    _TTI(rhoY, RR, IR, ll, dt, rho1)
                    rhoY += _COM(MuE, ll, dt, rho1) 
                    # rho1 = (dt/ll)*(numpy.tensordot(RR,rho1) 
                    #      - 1j*(numpy.dot(HH,rho1)-numpy.dot(rho1,HH)) + IR
                    #      + 1j*(numpy.dot(MU,rho1)-numpy.dot(rho1,MU))*EE[indx])
                    rho1 = rhoY       
                    rho2 = rho2 + rho1
                rho1 = rho2    

                if indxR < cutoff_indx - 1:                      
                    indxR += stride
                else:
                    indxR = cutoff_indx

                
            pr.data[indx,:,:] = rho2
            
            #
            # We respect the tensor cut-off
            #
            indx += 1
  
                
        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True
            
        return pr     



    def __propagate_short_exp_TDrelax_EField_oper(self, rhoi, L=4):
        """

        """
        debug("(9)")    




    def __propagate_short_exp_with_relaxation_field(self,rhoi,L=4):
        """Short exponential integration with relaxation and array field
        
        Propagates the system defined by Hamiltonian and the relaxation
        tensor under the influence of an electric field defined as 
        a real-valued array. The relaxation tensor must be independent
        of time.
        
        It is implicitly assumed that the field is X-polarized.
        
        Rotating wabe approximation is not considered.
        
        """
        if self.RelaxationTensor.as_operators:
            return \
            self.__propagate_short_exp_with_relaxation_field_oper(rhoi, L=L)
        
        debug("(10)")
        
        # we forbid the refinement of the time step
        if self.Nref != 1:
            raise Exception("Cannot propagation with refined time-step.")
            
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        
        rho1 = rhoi.data
        rho2 = rhoi.data
        
        HH = self.Hamiltonian.data        
        RR = self.RelaxationTensor.data   
        
        #
        # We do not have an information on polarization - we take X as default
        #
        MU = self.Trdip.data[:,:,0]
        
        EE = self.Efield    
        
        IR = 0.0
        indx = 1
        for tt in self.TimeAxis.data[1:self.Nt]:
            
            #for jj in range(0,self.Nref):
            if self.has_Iterm:
                IR = self.RelaxationTensor.Iterm[indx,:,:]
            MuE = MU*EE[indx]
            
            for ll in range(1,L+1):
                
                rhoY =  -_COM(HH, ll, self.dt, rho1)
                _TTI(rhoY,RR, IR, ll, self.dt, rho1)
                rhoY += _COM(MuE, ll, self.dt, rho1)
                       # - (1j*self.dt/ll)*(numpy.dot(HH,rho1) \
                       #   - numpy.dot(rho1,HH) ) \
                       # + (self.dt/ll)*numpy.tensordot(RR,rho1) \
                       #  + (1j*self.dt/ll)*( numpy.dot(MU,rho1) \
                       #   - numpy.dot(rho1,MU) )*EE[indx]                             
                rho1 = rhoY        
                rho2 = rho2 + rho1
            rho1 = rho2   
            
            # the for jj loop would end here
                
            pr.data[indx,:,:] = rho2                        
            indx += 1             
            
        return pr


    def __propagate_short_exp_with_relaxation_field_oper(self, rhoi, L=4):
        
        debug("(12)")


    def __propagate_short_exp_with_relaxation_EField(self,rhoi,L=4):
        """Short exponential integration with relaxation and EField like object

        Propagates the system defined by Hamiltonian and the relaxation
        tensor under the influence of an electric field defined by
        the LabField or EField object. The relaxation tensor must be 
        independent of time.
        
        The EField carries the information on the field polarization, and 
        this information is used. 
        
        Rotating wabe approximation is considered if defined by the Hamiltonian
        object.

        """
        if self.RelaxationTensor.as_operators:
            return \
            self.__propagate_short_exp_with_relaxation_EField_oper(rhoi, L=L)
        
        debug("(11)") 
        # we forbid the refinement of the time step
        if self.Nref != 1:
            raise Exception("Cannot propagation with refined time-step.")
            
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
            
        rho1 = rhoi.data
        rho2 = rhoi.data
 
        #
        # RWA is applied here
        #
        if self.Hamiltonian.has_rwa:
            
            HH = self.Hamiltonian.get_RWA_data() 
            
            # get the field corresponding to RWA
            om = self.Hamiltonian.rwa_energies[self.Hamiltonian.rwa_indices[1]]
            
            try:
                Nfields = len(self.EField)
            except:
                Nfields = 1
            
            
            if Nfields > 1:
                
                # rotating wave frequency is set to all of the fields globally
                self.EField[0].set_rwa(om)

            else:
                self.EField.set_rwa(om)
            
            
            # the two complex components of the field
            if Nfields > 1:
                Epls = []
                Emin = []
                for kk in range(Nfields):
                    Epls.append(self.EField[kk].field_p)
                    Emin.append(self.EField[kk].field_m)                
            else:
                Epls = self.EField.field_p
                Emin = self.EField.field_m
 
            if Nfields > 1:
                self.EField[0].restore_rwa()
            else:            
                self.EField.restore_rwa()
            
            # upper and lower triagle
            N = self.Hamiltonian.dim
            Mu = numpy.zeros((N,N), dtype=qr.REAL)
            for ii in range(N):
                for jj in range(ii+1,N):
                    Mu[ii,jj] = 1.0
            Ml = numpy.transpose(Mu)
                   
            
        else:
            
            HH = self.Hamiltonian.data 
            EField = self.EField.field
        
        RR = self.RelaxationTensor.data
        
        if Nfields > 1:

            MU = []
            for kk in range(Nfields):
                pol = self.EField[kk].pol

                MU.append(numpy.dot(self.Trdip.data[:,:,:], pol))
        else:
            pol = self.EField.pol
            MU = numpy.dot(self.Trdip.data[:,:,:], pol)
        
        #
        # Propagation
        #
        IR = 0.0
        indx = 1
        for ii in self.TimeAxis.data[1:self.Nt]:

            if self.has_Iterm:
                IR = self.RelaxationTensor.Iterm[indx,:,:] 
                
            if self.Hamiltonian.has_rwa:
                if Nfields > 1:
                    MuE = []
                    for kk in range(Nfields):
                        MuE.append(MU[kk]*(Ml*Epls[kk][indx]+Mu*Emin[kk][indx]))
                else:
                    MuE = [MU*(Ml*Epls[indx]+Mu*Emin[indx])]
            else:
                if Nfields > 1:
                    MuE = []
                    for kk in range(Nfields):
                        MuE.append(MU[kk]*EField[kk][indx])
                else:
                    MuE = [MU*EField[indx]]
            
            #for jj in range(0,self.Nref):
            for ll in range(1,L+1):
                
                rhoY = -_COM(HH, ll, self.dt, rho1)
                _TTI(rhoY, RR, IR, ll, self.dt, rho1)
                
                for jj in range(Nfields):
                    rhoY += _COM(MuE[jj], ll, self.dt, rho1)
                             
                rho1 = rhoY
                
                rho2 = rho2 + rho1
            rho1 = rho2  
                
            # for jj loop would end here
                
            pr.data[indx,:,:] = rho2                        
            indx += 1
             
        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True   
            
        return pr
    

    def __propagate_short_exp_with_relaxation_EField_oper(self, rhoi, L=4):
        
        debug("(13")
    
    
def _OTI(rhoY, Km, Kd, Lm, Ld, ll, dt, rho1):
    """Operator form of the Time Indepedent relaxation tensor

    Parameters
    ----------
    rhoY : TYPE
        DESCRIPTION.
    Km : TYPE
        DESCRIPTION.
    Kd : TYPE
        DESCRIPTION.
    Lm : TYPE
        DESCRIPTION.
    Ld : TYPE
        DESCRIPTION.
    ll : TYPE
        DESCRIPTION.
    dt : TYPE
        DESCRIPTION.
    rho1 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    
    """
    
    Nm = Km.shape[0]
    for mm in range(Nm):
        
       rhoY += (dt/ll)*(
        numpy.dot(Km[mm,:,:],numpy.dot(rho1, Ld[mm,:,:]))
       +numpy.dot(Lm[mm,:,:],numpy.dot(rho1, Kd[mm,:,:]))
       -numpy.dot(numpy.dot(Kd[mm,:,:],Lm[mm,:,:]), rho1)
       -numpy.dot(rho1, numpy.dot(Ld[mm,:,:],Km[mm,:,:]))
       )

            
def _COM(HH, ll, dt, rho1):
    """Commutator with the a given operator (e.g. Hamiltonian)
    
    Parameters
    ----------
    rhoY : TYPE
        DESCRIPTION.
    HH : TYPE
        DESCRIPTION.
    ll : TYPE
        DESCRIPTION.
    dt : TYPE
        DESCRIPTION.
    rho1 : TYPE
        DESCRIPTION.

    Returns
    -------
    ret : numpy.array, 2 indices
        The resulting commutaror
    
    
    """
    ret = (1j*dt/ll)*(numpy.dot(HH,rho1) - numpy.dot(rho1,HH))
    return ret


def _TTI(rhoY, RR, IR, ll, dt, rho1):
    """ Tensor form of the Time Indendent relaxation tensor

    Parameters
    ----------
    rhoY : numpy.array, 2 indices
        Input and output density matrix 
    RR : numpy.array, 4 indices
        Relaxation tensor
    ll : int
        Index counting the steps in exponential expansion
    dt : float
        Time step
    rho1 : numpy.array, 2 indices
        State before the time step.

    Returns
    -------
    None.

    """
    rhoY += (dt/ll)*(numpy.tensordot(RR,rho1) + IR)

