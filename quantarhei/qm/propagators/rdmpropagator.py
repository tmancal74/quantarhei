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

#import scipy.integrate
import numpy.linalg

import matplotlib.pyplot as plt

#import cu.oqs.cython.propagators as prop

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

    

class ReducedDensityMatrixPropagator(MatrixData, Saveable): 
    """
    
    Reduced Density Matrix Propagator calculates the evolution of the
    reduced density matrix based on the Hamiltonian and optionally
    a relaxation tensor. Relaxation tensor may be constant or time
    dependent. 
    
    """
    
    def __init__(self, timeaxis=None, Ham=None, RTensor=None,
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
        
        >>> pr = RDMPropagator([[0.0, 0.0],[0.0, 1.0]],[0,1,2,3,4,5])
        Traceback (most recent call last):
            ...
        Exception
        
        The correct way to construct the propagator is the following:

        >>> h = numpy.array([[0.0, 0.0],[0.0,1.0]])
        >>> HH = Hamiltonian(data=h)
        >>> times = TimeAxis(0,1000,1)
        >>> pr = RDMPropagator(HH,times)
        
        """
        self.has_Trdip = False
        self.has_Efield = False
        self.has_PDeph = False
        self.has_RTensor = False
        self.has_RWA = False
        self.has_EField = False
        
        if not ((timeaxis is None) and (Ham is None)):
            
            #
            # Hamiltonian and TimeAxis are required
            #
            if isinstance(Ham,Hamiltonian):
                self.Hamiltonian = Ham
            else:
                raise Exception
            
            if isinstance(timeaxis,TimeAxis):
                self.TimeAxis = timeaxis
            else:
                raise Exception
            
            #
            # RelaxationTensor is not requited
            #
            if isinstance(RTensor,RelaxationTensor):
                self.RelaxationTensor = RTensor
                self.has_RTensor = True
                self.has_relaxation = True
            elif RTensor is None:
                self.has_RTensor = False
                self.has_relaxation = False
            else:
                raise Exception
    
            if Trdip is not None:            
                if isinstance(Trdip,Operator):
                    self.Trdip = Trdip
                    self.has_Trdip = True
                else:
                    raise Exception
    
            if Efield is not None:
                if isinstance(Efield,numpy.ndarray):
                    self.Efield = Efield
                    self.has_Efield = True
                    self.has_EField = False
                else: 
                    self.EField = Efield
                    self.has_EField = True
                    self.has_Efield = False                    


            
            
            #
            # RWA
            #
            #ham = self.Hamiltonian
            #if ham.has_rwa:
                
                #self.RWA = self.Hamiltonian.rwa_indices
                #self.RWU = numpy.zeros(self.RWA.shape, dtype=self.RWA.dtype)
                
                #HH = ham.data
                #shape = HH.shape
                #HOmega = numpy.zeros(shape, dtype=qr.REAL)
                #for ii in range(shape[0]):
                #    HOmega[ii,ii] = ham.rwa_energies[ii]
                                    
                #self.HOmega = self.ham.get_RWA_skeleton()
                
                #print(self.RWA)
                #print(self.RWU)
            
            #
            # Pure dephasing also counts as relaxation
            #
            if PDeph is not None:
                self.PDeph = PDeph
                self.has_PDeph = True
                self.has_relaxation = True
            
            
            self.Odt = self.TimeAxis.data[1]-self.TimeAxis.data[0]
            self.dt = self.Odt
            self.Nref = 1
            
            self.Nt = self.TimeAxis.data.shape[0]
            
            N = self.Hamiltonian.data.shape[0]
            self.N = N
            self.data = numpy.zeros((self.Nt,N,N),dtype=numpy.complex64)
            self.propagation_name = ""
        
            self.verbose = Manager().log_conf.verbose

        
    def setDtRefinement(self, Nref):
        """
        The TimeAxis object specifies at what times the propagation
        should be stored. We can tell the propagator to use finer
        time step for the calculation by setting the refinement. The
        refinement is an integer by which the TimeAxis time step should
        be devided to get the finer time step. In the code below, we
        have dt = 10 in the TimeAxis, but we want to calculate with
        dt = 1
        
        >>> HH = numpy.array([[0.0, 0.0],[0.0,1.0]])
        >>> times = numpy.linspace(0,1000,10)
        >>> pr = RDMPropagator(HH,times)
        >>> pr.setDtRefinement(10)
        
        """
        self.Nref = Nref
        self.dt = self.Odt/self.Nref
        
        
    def propagate(self, rhoi, method="short-exp", mdata=None, name=""):
        """
        
        >>> T0   = 0
        >>> Tmax = 100
        >>> dt   = 1
        >>> Nref = 30


        >>> initial_dm = [[1.0, 0.0, 0.0],
        ...               [0.0, 0.0, 0.0],
        ...               [0.0, 0.0, 0.0]]

        >>> Hamiltonian = [[0.0, 0.1, 0.0],
        ...                [0.1, 0.0, 0.1],
        ...                [0.0, 0.1, 0.1]]        
        
        >>> HH = numpy.array(Hamiltonian)
        >>> times = numpy.linspace(T0,Tmax,(Tmax-T0)/dt+1)
        
        >>> pr = RDMPropagator(HH,times)
        >>> pr.setDtRefinement(Nref)

        >>> rhoi = numpy.array(initial_dm)  

        >>> pr.propagate(rhoi,method="primitive")
        
        
        """
        
        #
        # FIXME: Remove this
        #
        self.propagation_name = name
        
        #
        # Testing if the object submitted is density matrix
        #
        if not (isinstance(rhoi, ReducedDensityMatrix) 
             or isinstance(rhoi, DensityMatrix)):
            raise Exception("First argument has be of"+
            "the ReducedDensityMatrix type")
              
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

                    #
                    # FIXME: These methods are untested
                    #
                    elif method == "primitive":
                        return self.__propagate_primitive_with_relaxation(rhoi)
                    elif method == "Runge-Kutta":
                        return self.__propagate_Runge_Kutta(rhoi)
                    elif method == "diagonalization":
                        return self.__propagate_diagonalization(rhoi)

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

                raise Exception("NOT IMPLEMENTED")

            else:
                 
                if method == "short-exp":
                    return self.__propagate_short_exp(rhoi,L=4)
                elif method == "short-exp-2":
                    return self.__propagate_short_exp(rhoi,L=2)
                elif method == "short-exp-4":
                    return self.__propagate_short_exp(rhoi,L=4)
                elif method == "short-exp-6":
                    return self.__propagate_short_exp(rhoi,L=6)            
    
                #
                # FIXME: These methods are not tested
                #
                elif method == "primitive":
                    return self.__propagate_primitive(rhoi)
                elif method == "Runge-Kutta":
                    return self.__propagate_Runge_Kutta(rhoi)
                elif method == "diagonalization":
                    return self.__propagate_diagonalization(rhoi)
    
                else:
                    raise Exception("Unknown propagation method: "+method)
        
            
        
    def __propagate_primitive(self, rhoi):
        """Primitive integration of equantion of motion
        
        This is of no other than perhaps pedagogical value
        
        """
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)

        rhoPrim = rhoi.data
        HH = self.Hamiltonian.data  
        
        
        indx = 0
        for ii in self.TimeAxis.time: 

            for jj in range(0,self.Nref):   
                
                drho  = -1j*(  numpy.dot(HH,rhoPrim) \
                             - numpy.dot(rhoPrim,HH) )
                             
                rhoPrim = rhoPrim + drho*self.dt

            pr.data[indx,:,:] = rhoPrim            
            
            indx += 1            
            
        return pr

        
    def __propagate_primitive_with_relaxation(self, rhoi):
        """Primitive integration of equantion of motion with relaxation
        
        This is of no other than perhaps pedagogical value
        
        """
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
       
        rhoPrim = rhoi.data
        HH = self.Hamiltonian.data        
        RR = self.RelaxationTensor.data        
        
        indx = 0
        for ii in self.TimeAxis.data: 

            for jj in range(0,self.Nref):   
                
                drho  = -1j*(  numpy.dot(HH,rhoPrim) \
                             - numpy.dot(rhoPrim,HH) ) \
                        + numpy.tensordot(RR,rhoPrim)
                             
                rhoPrim = rhoPrim + drho*self.dt

            pr.data[indx,:,:] = rhoPrim            
            
            indx += 1            
            
        return pr
        
        
    def __propagate_Runge_Kutta(self, rhoi):
        """Runge-Kutta integration of equation of motion
              
        NOT IMPLEMENTED
        
        """
        indx = 0
        for ii in self.timeaxis: 

            self.rho[:,:,indx] = rhoi            
            
            indx += 1  
 

    def __propagate_short_exp(self, rhoi, L=4):
        """Short expansion of an exponention to integrate equations of motion
        
        
        """
        
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        rho1 = rhoi.data
        rho2 = rhoi.data
        
        #HH = self.Hamiltonian.data
        #
        # RWA is applied here
        #
        if self.Hamiltonian.has_rwa:
            HH = self.Hamiltonian.get_RWA_data() #data  - self.HOmega
        else:
            HH = self.Hamiltonian.data
        
        indx = 1
        for ii in self.TimeAxis.data[1:self.Nt]:
            
            for jj in range(0,self.Nref):
                
                for ll in range(1,L+1):
                    
                    rho1 = -1j*(self.dt/ll)*(numpy.dot(HH,rho1) \
                             - numpy.dot(rho1,HH) )
                             
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pr.data[indx,:,:] = rho2                        
            indx += 1                       
            
        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True
            
        return pr
 
 
    def __propagate_short_exp_with_relaxation(self, rhoi, L=4):
        """Integration by short exponentional expansion
        
        Integration by expanding exponential to Lth order
              
              
        """
        
        try:
            if self.RelaxationTensor.as_operators:
                return self.__propagate_short_exp_with_rel_operators(rhoi, L=L)
        except:
            raise Exception("Operator propagation failed")
        
        
        pr = ReducedDensityMatrixEvolution(self.TimeAxis, rhoi,
                                           name=self.propagation_name)
        
        rho1 = rhoi.data
        rho2 = rhoi.data
        
        #HH = self.Hamiltonian.data  
        #
        # RWA is applied here
        #
        if self.Hamiltonian.has_rwa:
            HH = self.Hamiltonian.get_RWA_data() #data  - self.HOmega
        else:
            HH = self.Hamiltonian.data
            
        RR = self.RelaxationTensor.data

        if self.has_PDeph:
            
            if self.PDeph.dtype == "Lorentzian":
                expo = numpy.exp(-self.PDeph.data*self.dt)
                t0 = 0.0
            elif self.PDeph.dtype == "Gaussian":
                expo = numpy.exp(-self.PDeph.data*(self.dt**2)/2.0)
                t0 = self.PDeph.data*self.dt

            
            indx = 1
            for ii in range(1, self.Nt): 

                # time at the beginning of the step
                tNt = self.TimeAxis.data[indx-1]  
                
                for jj in range(0, self.Nref):
                    
                    tt = tNt + jj*self.dt  # time right now 
 
                    for ll in range(1, L+1):
                        
                        rho1 =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) 
                                                 - numpy.dot(rho1,HH)) \
                               + (self.dt/ll)*numpy.tensordot(RR,rho1)
                                 
                        rho2 = rho2 + rho1
                        
                    # pure dephasing is added here                        
                    rho2 = rho2*expo*numpy.exp(-t0*tt)
                        
                    rho1 = rho2    
                    
                pr.data[indx,:,:] = rho2 
                indx += 1   
                
        else:
            
            indx = 1
            for ii in range(1, self.Nt): 
                
                for jj in range(0, self.Nref):
                    
                    for ll in range(1, L+1):
                        
                        rho1 =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) 
                                                 - numpy.dot(rho1,HH)) \
                               + (self.dt/ll)*numpy.tensordot(RR,rho1)
                                 
                                 
                        rho2 = rho2 + rho1
                    rho1 = rho2    
                    
                pr.data[indx,:,:] = rho2 
                indx += 1   
           

        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True
            
        return pr  

        
    def __propagate_short_exp_with_rel_operators(self, rhoi, L=4):
        """Integration by short exponentional expansion
        
        Integration by expanding exponential to Lth order. 
              
            
        """
        mana = Manager()
        save_pytorch = None
        
        legacy = mana.gen_conf.legacy_relaxation
        
        if mana.num_conf.gpu_acceleration:
            save_pytorch = mana.num_conf.enable_pytorch
            mana.num_conf.enable_pytorch = True
            
        if mana.num_conf.enable_pytorch and (not legacy):
            ret =  self._propagate_SExp_RTOp_ReSymK_Re_pytorch(rhoi,
                                        self.Hamiltonian,
                                        self.RelaxationTensor,
                                        self.dt, 
                                        use_gpu=mana.num_conf.gpu_acceleration,
                                        L=L)
            
            if save_pytorch is not None:
                mana.num_conf.enable_pytorch = save_pytorch
                
            return ret
        
        elif not legacy:
            return self._propagate_SExp_RTOp_ReSymK_Re_numpy(rhoi,
                                                 self.Hamiltonian,
                                                 self.RelaxationTensor,
                                                 self.dt, L=L)
        
        #
        # legacy version
        #

        pr = ReducedDensityMatrixEvolution(self.TimeAxis, rhoi,
                                           name=self.propagation_name)
        
        rho1 = rhoi.data
        rho2 = rhoi.data
        
        #
        # RWA is applied here
        #
        if self.Hamiltonian.has_rwa:
            HH = self.Hamiltonian.get_RWA_data() #data  - self.HOmega
        else:
            HH = self.Hamiltonian.data
        
        qr.log_detail("PROPAGATION (short exponential with "+
                     "relaxation in operator form): order ", L, 
                     verbose=self.verbose)
        qr.log_detail("Using complex numpy implementation")
        
        try:
            Km = self.RelaxationTensor.Km # real
            Lm = self.RelaxationTensor.Lm # complex
            Ld = self.RelaxationTensor.Ld # complex - get by transposition
            Kd = numpy.zeros(Km.shape, dtype=numpy.float64)
            Nm = Km.shape[0]
            for m in range(Nm):
                Kd[m, :, :] = numpy.transpose(Km[m, :, :])
        except:
            raise Exception("Tensor is not in operator form")
            
        indx = 1

        levs = [qr.LOG_QUICK] #, 8]
        verb = qr.loglevels2bool(levs)

        # after each step we apply pure dephasing (if present)
        if self.has_PDeph:
            
            if self.PDeph.dtype == "Lorentzian":
                expo = numpy.exp(-self.PDeph.data*self.dt)
                t0 = 0.0
            elif self.PDeph.dtype == "Gaussian":
                expo = numpy.exp(-self.PDeph.data*(self.dt**2)/2.0)
                t0 = self.PDeph.data*self.dt
            
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
                        
                        rhoY =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) 
                                                 - numpy.dot(rho1,HH))
                        
                        #rhoX = numpy.zeros(rho1.shape, dtype=numpy.complex128)
                        for mm in range(Nm):
                            
                           rhoY += (self.dt/ll)*(
                            numpy.dot(Km[mm,:,:],numpy.dot(rho1, Ld[mm,:,:]))
                           +numpy.dot(Lm[mm,:,:],numpy.dot(rho1, Kd[mm,:,:]))
                           -numpy.dot(numpy.dot(Kd[mm,:,:],Lm[mm,:,:]), rho1)
                           -numpy.dot(rho1, numpy.dot(Ld[mm,:,:],Km[mm,:,:]))
                           )
                                 
                        rho1 = rhoY #+ rhoX
                        
                        rho2 = rho2 + rho1
                       
                    # pure dephasing is added here                        
                    rho2 = rho2*expo*numpy.exp(-t0*tt)
                        
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
                        
                        rhoY =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) 
                                                 - numpy.dot(rho1,HH))
                        
                        #rhoX = numpy.zeros(rho1.shape, dtype=numpy.complex128)
                        for mm in range(Nm):
                            
                           rhoY += (self.dt/ll)*(
                            numpy.dot(Km[mm,:,:],numpy.dot(rho1, Ld[mm,:,:]))
                           +numpy.dot(Lm[mm,:,:],numpy.dot(rho1, Kd[mm,:,:]))
                           -numpy.dot(numpy.dot(Kd[mm,:,:],Lm[mm,:,:]), rho1)
                           -numpy.dot(rho1, numpy.dot(Ld[mm,:,:],Km[mm,:,:]))
                           )
                                 
                        rho1 = rhoY #+ rhoX
                        
                        rho2 = rho2 + rho1
                    
                    rho1 = rho2    
                
                pr.data[indx,:,:] = rho2 
                indx += 1
           
             
        qr.log_detail("...DONE")

        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True
            
        return pr


    def _propagate_SExp_RTOp_ReSymK_Re_numpy(self, rhoi, Ham, RT, dt, L=4):
        """Integration by short exponentional expansion
        
        Integration by expanding exponential (_SExp_) to Lth order. 
        This is a numpy (_numpy) implementation with real (_Re_) matrices
        for  a system part of the system-bath interaction operator  ``K``
        in a form of real symmetric operator (ReSymK). The relaxation tensor
        is assumed in form of a set of operators (_RTOp_)
              
            
        """

        Nref = self.Nref
        Nt = self.Nt
        verbose = self.verbose
        timea = self.TimeAxis
        prop_name = self.propagation_name
        
        # no self beyond this point
        
        qr.log_detail("PROPAGATION (short exponential with "+
                    "relaxation in operator form): order ", L, 
                    verbose=verbose)
        qr.log_detail("Using real valued numpy implementation")
        
        pr = ReducedDensityMatrixEvolution(timea, rhoi,
                                           name=prop_name)
        
        rho1_r = numpy.real(rhoi.data)
        rho2_r = numpy.real(rhoi.data)
        rho1_i = numpy.imag(rhoi.data)
        rho2_i = numpy.imag(rhoi.data)
         
        HH = Ham.data
                
        try:
            Km = RT.Km #self.RelaxationTensor.Km # real
            Lm_r = numpy.real(RT.Lm) #self.RelaxationTensor.Lm) # complex
            Lm_i = numpy.imag(RT.Lm) #self.RelaxationTensor.Lm)
            Nm = Km.shape[0]
        except:
            raise Exception("Tensor is not in operator form")
            
        indx = 1
        
        # verbosity inside loops
        levs = [qr.LOG_QUICK] 
        verb = qr.loglevels2bool(levs, verbose=self.verbose)

        # after each step we apply pure dephasing (if present)
        if self.has_PDeph:
        
            # loop over time
            for ii in range(1, Nt):
                qr.printlog("time step ", ii, "of", Nt, 
                            verbose=verb[0], loglevel=levs[0], end="\r")
                
                # steps in between saving the results
                for jj in range(Nref):
                    
                    # L interations to get short exponential expansion
                    for ll in range(1, L+1):
    
                        A = numpy.dot(HH,rho1_i)
                        B = numpy.dot(HH,rho1_r)
                        rhoY_r =  (dt/ll)*(A + numpy.transpose(A))
                        rhoY_i = -(dt/ll)*(B - numpy.transpose(B)) 
                        
                        for mm in range(Nm):
                        
                            a = numpy.dot(Lm_r[mm,:,:], rho1_r)
                            A = a - numpy.transpose(a)
                            b = numpy.dot(Lm_i[mm,:,:], rho1_i)
                            B = b - numpy.transpose(b)
                            c = numpy.dot(Lm_r[mm,:,:], rho1_i)
                            C = -(c + numpy.transpose(c))
                            d = numpy.dot(Lm_i[mm,:,:], rho1_r)
                            D = d + numpy.transpose(d)
                            
                            E = B - A
                            F = C - D
                            
                            A = numpy.dot(Km[mm,:,:], E)
                            B = numpy.dot(Km[mm,:,:], F)
                            rhoY_r += (dt/ll)*(A + numpy.transpose(A))
                            rhoY_i += (dt/ll)*(B - numpy.transpose(B))
                            
                        rho1_r = rhoY_r 
                        rho1_i = rhoY_i
                        
                        rho2_r +=  rho1_r
                        rho2_i +=  rho1_i
                        
                        rho2_r = rho2_r*numpy.exp(-self.PDeph.data*dt)
                        rho2_i = rho2_i*numpy.exp(-self.PDeph.data*dt)
                        
                    rho1_r = rho2_r
                    rho1_i = rho2_i
                    
                pr.data[indx,:,:] = rho2_r + 1j*rho2_i 
                indx += 1             

        # propagatiomn with no extra dephasing
        else:
            
            # loop over time
            for ii in range(1, Nt):
                qr.printlog("time step ", ii, "of", Nt, 
                            verbose=verb[0], loglevel=levs[0], end="\r")
                
                # steps in between saving the results
                for jj in range(Nref):
                    
                    # L interations to get short exponential expansion
                    for ll in range(1, L+1):
    
                        A = numpy.dot(HH,rho1_i)
                        B = numpy.dot(HH,rho1_r)
                        rhoY_r =  (dt/ll)*(A + numpy.transpose(A))
                        rhoY_i = -(dt/ll)*(B - numpy.transpose(B)) 
                        
                        for mm in range(Nm):
                        
                            a = numpy.dot(Lm_r[mm,:,:], rho1_r)
                            A = a - numpy.transpose(a)
                            b = numpy.dot(Lm_i[mm,:,:], rho1_i)
                            B = b - numpy.transpose(b)
                            c = numpy.dot(Lm_r[mm,:,:], rho1_i)
                            C = -(c + numpy.transpose(c))
                            d = numpy.dot(Lm_i[mm,:,:], rho1_r)
                            D = d + numpy.transpose(d)
                            
                            E = B - A
                            F = C - D
                            
                            A = numpy.dot(Km[mm,:,:], E)
                            B = numpy.dot(Km[mm,:,:], F)
                            rhoY_r += (dt/ll)*(A + numpy.transpose(A))
                            rhoY_i += (dt/ll)*(B - numpy.transpose(B))
                            
                        rho1_r = rhoY_r 
                        rho1_i = rhoY_i
                        
                        rho2_r +=  rho1_r
                        rho2_i +=  rho1_i
                        
                    rho1_r = rho2_r
                    rho1_i = rho2_i
                    
                pr.data[indx,:,:] = rho2_r + 1j*rho2_i 
                indx += 1             

        
        qr.log_detail()
        qr.log_detail("...DONE")

        return pr

    def _propagate_SExp_RTOp_ReSymK_Re_pytorch(self, rhoi, Ham, RT, dt,
                                               use_gpu=False, L=4):
        """Integration by short exponentional expansion
        
        Integration by expanding exponential (_SExp_) to Lth order. 
        This is a PyTorch (_pytorch) implementation with real (_Re_) matrices
        for  a system part of the system-bath interaction operator  ``K``
        in a form of real symmetric operator (ReSymK). The relaxation tensor
        is assumed in form of a set of operators (_RTOp_)
              
            
        """

        Nref = self.Nref
        Nt = self.Nt
        verbose = self.verbose
        timea = self.TimeAxis
        prop_name = self.propagation_name
        
        try: 
            import torch
        except:
            raise Exception("PyTorch not installed")
        
        # no self beyond this point
        
        qr.log_detail("PROPAGATION (short exponential with "+
                    "relaxation in operator form): order ", L, 
                    verbose=verbose)
        qr.log_detail("Using pytorch implementation")
        qr.log_detail("Using GPU: ", use_gpu & torch.cuda.is_available())
        
        pr = ReducedDensityMatrixEvolution(timea, rhoi,
                                           name=prop_name)
        
        rho1_r = torch.from_numpy(numpy.real(rhoi.data))
        rho2_r = torch.from_numpy(numpy.real(rhoi.data))
        rho1_i = torch.from_numpy(numpy.imag(rhoi.data))
        rho2_i = torch.from_numpy(numpy.imag(rhoi.data))
         
        HH = torch.from_numpy(Ham.data)
                
        try:
            Km = torch.from_numpy(RT.Km) #self.RelaxationTensor.Km # real
            Lm_r = torch.from_numpy(numpy.real(RT.Lm)) #self.RelaxationTensor.Lm) # complex
            Lm_i = torch.from_numpy(numpy.imag(RT.Lm)) #self.RelaxationTensor.Lm)
            Nm = RT.Km.shape[0]
        except:
            raise Exception("Tensor is not in operator form")
            
        if use_gpu & torch.cuda.is_available():
            rho1_r = rho1_r.cuda()
            rho2_r = rho1_r
            rho1_i = rho1_i.cuda()
            rho2_i = rho1_i
            HH = HH.cuda()
            Km = Km.cuda()
            Lm_r = Lm_r.cuda()
            Lm_i = Lm_i.cuda()
 
        indx = 1
        
        # verbosity inside loops
        levs = [qr.LOG_QUICK] 
        verb = qr.loglevels2bool(levs)
        
        # loop over time
        for ii in range(1, Nt):
            qr.printlog(" time step ", ii, "of", Nt, 
                        verbose=verb[0], loglevel=levs[0])
            
            # steps in between saving the results
            for jj in range(Nref):
                
                # L interations to get short exponential expansion
                for ll in range(1, L+1):

                    A = torch.matmul(HH,rho1_i)
                    B = torch.matmul(HH,rho1_r)
                    rhoY_r = torch.mul(A + torch.transpose(A, 0, 1), dt/ll)
                    rhoY_i = torch.mul(B - torch.transpose(B, 0, 1), -dt/ll) 

                    for mm in range(Nm):
                    
                        a = torch.matmul(Lm_r[mm,:,:], rho1_r)
                        A = a - torch.transpose(a, 0, 1)
                        b = torch.matmul(Lm_i[mm,:,:], rho1_i)
                        B = b - torch.transpose(b, 0, 1)
                        c = torch.matmul(Lm_r[mm,:,:], rho1_i)
                        C = -(c + torch.transpose(c, 0, 1))
                        d = torch.matmul(Lm_i[mm,:,:], rho1_r)
                        D = d + torch.transpose(d, 0, 1)
                        
                        E = B - A
                        F = C - D
                        
                        A = torch.matmul(Km[mm,:,:], E)
                        B = torch.matmul(Km[mm,:,:], F)
                        rhoY_r += torch.mul(A + torch.transpose(A, 0, 1),dt/ll)
                        rhoY_i += torch.mul(B - torch.transpose(B, 0, 1),dt/ll)
 
                    rho1_r = rhoY_r 
                    rho1_i = rhoY_i
                    
                    rho2_r += rho1_r
                    rho2_i += rho1_i
                    
                rho1_r = rho2_r
                rho1_i = rho2_i
            
            if use_gpu & torch.cuda.is_available():
                rho2_sr = rho2_r.cpu()
                rho2_si = rho2_i.cpu()
            else:
                rho2_sr = rho2_r
                rho2_si = rho2_i                
    
            pr.data[indx,:,:] = rho2_sr.numpy() + 1j*rho2_si.numpy() 
            indx += 1             
         
        qr.log_detail("...DONE")
        return pr


    def __propagate_short_exp_with_TD_relaxation(self,rhoi,L=4):
        """
              Short exp integration
        """

        try:
            if self.RelaxationTensor.as_operators:
                return self.__propagate_short_exp_with_TDrel_operators(rhoi, L=L)
        except:
            raise Exception("Operator propagation failed")
        
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        
        rho1 = rhoi.data
        rho2 = rhoi.data

        
        #HH = self.Hamiltonian.data  
        #
        # RWA is applied here
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
            
        indx = 1
        indxR = 1
        for ii in self.TimeAxis.data[1:self.Nt]:

            RR = self.RelaxationTensor.data[indxR,:,:]        
            
            for jj in range(0,self.Nref):
                for ll in range(1,L+1):
                    
                    rho1 =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) \
                             - numpy.dot(rho1,HH) ) \
                           + (self.dt/ll)*numpy.tensordot(RR,rho1)
                             
                             
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pr.data[indx,:,:] = rho2
            
            indx += 1
            if indxR < cutoff_indx-1:                      
                indxR += 1             
            
        return pr     


    def __propagate_short_exp_with_TDrel_operators(self, rhoi, L=4):
        """
              Short exp integration
        """

        pr = ReducedDensityMatrixEvolution(self.TimeAxis, rhoi,
                                           name=self.propagation_name)
        
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

        try:
            Km = self.RelaxationTensor.Km
            Kd = numpy.zeros(Km.shape, dtype=numpy.float64)
            Nm = Km.shape[0]
            for m in range(Nm):
                Kd[m, :, :] = numpy.transpose(Km[m, :, :])
        except:
            raise Exception("Tensor is not in operator form")
                        
        indx = 1
        indxR = 1
        for ii in range(1, self.Nt): 

            Lm = self.RelaxationTensor.Lm[indxR,:,:,:]
            Ld = self.RelaxationTensor.Ld[indxR,:,:,:]
       
            for jj in range(0, self.Nref):
                
                for ll in range(1, L+1):
                    
                    rhoY =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) 
                                             - numpy.dot(rho1,HH))
                    
                    rhoX = numpy.zeros(rho1.shape, dtype=numpy.complex128)
                    for mm in range(Nm):
                        
                       rhoX += (self.dt/ll)*(
                        numpy.dot(Km[mm,:,:],numpy.dot(rho1, Ld[mm,:,:]))
                       +numpy.dot(Lm[mm,:,:],numpy.dot(rho1, Kd[mm,:,:]))
                       -numpy.dot(numpy.dot(Kd[mm,:,:],Lm[mm,:,:]), rho1)
                       -numpy.dot(rho1, numpy.dot(Ld[mm,:,:],Km[mm,:,:]))
                       )
                             
                    rho1 = rhoY + rhoX
                    
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pr.data[indx,:,:] = rho2 
            indx += 1             
            if indxR < cutoff_indx-1:                      
                indxR += 1             
            
        return pr
        
        
    def __propagate_diagonalization(self,rhoi):
        pass
        
        
        
    def __propagate_short_exp_with_TD_relaxation_field(self,rhoi,L=4):
        """Short exp integration of the density matrix with external driving
        
        
        
        """
        try:
            if self.RelaxationTensor.as_operators:
                return self.__propagate_short_exp_with_TD_relaxation_field_operators(rhoi, L=L)
        except:
            raise Exception("Operator propagation failed")
            
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        
        rho1 = rhoi.data
        rho2 = rhoi.data
        
        #HH = self.Hamiltonian.data    
        MU = self.Trdip.data 
        #
        # RWA is applied here
        #
        # FIXME: RWA has to be applied to ralaxation tensor, too!!!!
        #
        if self.Hamiltonian.has_rwa:
            HH = self.Hamiltonian.get_RWA_data() #data  - self.HOmega
        else:
            HH = self.Hamiltonian.data
        
        indx = 1
        for ii in self.TimeAxis.time[1:self.Nt]:

            RR = self.RelaxationTensor.data[indx,:,:]        
            EE = self.Efield[indx]
            
            for jj in range(0,self.Nref):
                for ll in range(1,L+1):
                    
                    rho1 =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) \
                             - numpy.dot(rho1,HH) ) \
                           + (self.dt/ll)*numpy.tensordot(RR,rho1) \
                            + (1j*self.dt/ll)*( numpy.dot(MU,rho1) \
                             - numpy.dot(rho1,MU) )*EE
                             
                             
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pr.data[indx,:,:] = rho2                        
            indx += 1             
            
        return pr         


    def __propagate_short_exp_with_TD_relaxation_field_operators(self, rhoi, L=4):
        """

        """
        pass        


    def __propagate_short_exp_with_relaxation_field(self,rhoi,L=4):
        """
              Short exp integration
        """

            
        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        
        rho1 = rhoi.data
        rho2 = rhoi.data
        
        HH = self.Hamiltonian.data        
        RR = self.RelaxationTensor.data        
        MU = self.Trdip.data
        
        indx = 1
        for ii in self.TimeAxis.data[1:self.Nt]:
            
            EE = self.Efield[indx]
            
            for jj in range(0,self.Nref):
                for ll in range(1,L+1):
                    
                    rho1 =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) \
                             - numpy.dot(rho1,HH) ) \
                           + (self.dt/ll)*numpy.tensordot(RR,rho1) \
                            + (1j*self.dt/ll)*( numpy.dot(MU,rho1) \
                             - numpy.dot(rho1,MU) )*EE                             
                             
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pr.data[indx,:,:] = rho2                        
            indx += 1             
            
        return pr


    def __propagate_short_exp_with_relaxation_EField(self,rhoi,L=4):
        """
              Short exp integration
        """

        pr = ReducedDensityMatrixEvolution(self.TimeAxis,rhoi)
        
        rho1 = rhoi.data
        rho2 = rhoi.data
 
        #
        # RWA is applied here
        #
        if self.Hamiltonian.has_rwa:
            
            HH = self.Hamiltonian.get_RWA_data() #data  - self.HOmega
            
            # get the field corresponding to RWA
            #print(self.Hamiltonian.rwa_energies)
            #print(self.Hamiltonian.rwa_indices)
            om = self.Hamiltonian.rwa_energies[self.Hamiltonian.rwa_indices[1]]
            self.EField.subtract_omega(om)
            
            # the two complex components of the field
            Epls = self.EField.field(sign=1)
            Emin = self.EField.field(sign=-1)
            self.EField.restore_omega()
            
            # upper and lower triagle
            N = self.Hamiltonian.dim
            Mu = numpy.zeros((N,N), dtype=qr.REAL)
            for ii in range(N):
                for jj in range(ii+1,N):
                    Mu[ii,jj] = 1.0
            Ml = numpy.transpose(Mu)
            
        else:
            
            HH = self.Hamiltonian.data 
            EField = self.EField.field()
        
        RR = self.RelaxationTensor.data
        
        pol = self.EField.pol
        MU = numpy.dot(self.Trdip.data[:,:,:], pol)
        
        #
        # Propagation
        #
        
        indx = 1
        for ii in self.TimeAxis.data[1:self.Nt]:
            
            if self.Hamiltonian.has_rwa:
                MuE = MU*(Mu*Epls[indx]+Ml*Emin[indx])
            else:
                MuE = MU*EField[indx]
            
            for jj in range(0,self.Nref):
                for ll in range(1,L+1):
                    
                    rho1 =  - (1j*self.dt/ll)*(numpy.dot(HH,rho1) \
                             - numpy.dot(rho1,HH) ) \
                           + (self.dt/ll)*numpy.tensordot(RR,rho1) \
                            + (1j*self.dt/ll)*( numpy.dot(MuE,rho1) \
                             - numpy.dot(rho1,MuE) )                             
                             
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pr.data[indx,:,:] = rho2                        
            indx += 1
             
        if self.Hamiltonian.has_rwa:
            pr.is_in_rwa = True   
            
        return pr

