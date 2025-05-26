# -*- coding: utf-8 -*-

import numpy
import time
import os

from ..utils import derived_type
from ..builders.aggregates import Aggregate
from ..builders.molecules import Molecule
from ..builders.opensystem import OpenSystem
from ..core.time import TimeAxis
from ..core.managers import eigenbasis_of
from ..qm.propagators.poppropagator import PopulationPropagator
from .twod2 import TwoDResponse
from .. import signal_REPH, signal_NONR

from ..spectroscopy.responses import NonLinearResponse

from .. import REAL

# deprecated class
from ..spectroscopy.responses import LiouvillePathway

import quantarhei as qr

#try:
#
#    import aceto.nr3td as nr3td            
#    from aceto.lab_settings import lab_settings
#    from aceto.band_system import band_system            
#    _have_aceto = True
#    
#except:
    #
    #
    #raise Exception("Aceto not available")
from ..implementations.aceto import nr3td  
from ..implementations.aceto.band_system import band_system 
from ..implementations.aceto.lab_settings import lab_settings         
_have_aceto = False # we assume we have Python implementation


class TwoDResponseCalculator:
    """Calculator of the 2D spectrum
    
    
    Enables setting up parameters of 2D spectrum calculation for later
    evaluation. The method `calculate` returns TwoDSpectrumContainer
    with a 2D spectrum.
    
    Parameters
    ----------
    
    
    """

    t1axis = derived_type("t1axis",TimeAxis)
    t2axis = derived_type("t2axis",TimeAxis)
    t3axis = derived_type("t3axis",TimeAxis)
    
    system = derived_type("system",[Molecule,Aggregate,OpenSystem])
    
    _has_responses = False
    _has_system = False
    
    def __init__(self, t1axis, t2axis, t3axis, system=None, responses=None,
                 dynamics="secular", relaxation_tensor=None, rate_matrix=None,
                 effective_hamiltonian=None):
            
            
        self.t1axis = t1axis
        self.t2axis = t2axis
        self.t3axis = t3axis
        
        #FIXME: check the compatibility of the axes 
        
        if system is not None:
            self.system = system
            self._has_system = True
        else:
            self._has_system = False
            
        if responses is not None:
            self.resp_fcions = responses
            self._has_responses = True
        else:
            self._has_responses = False
        
        #FIXME: properties to be protected
        self.dynamics = dynamics
        
        # unprotected properties
        self.data = None

        self.responses = []
        
        self._relaxation_tensor = None
        self._rate_matrix = None
        self._relaxation_hamiltonian = None
        self._has_relaxation_tensor = False
        self._has_rate_matrix = False
        if relaxation_tensor is not None:
            self._relaxation_tensor = relaxation_tensor
            self._has_relaxation_tensor = True
        if effective_hamiltonian is not None:
            self._relaxation_hamiltonian = effective_hamiltonian
        if rate_matrix is not None:
            self._rate_matrix = rate_matrix
            self._has_rate_matrix = True

        
        #
        # after bootstrap information
        #
        self.sys = None
        self.lab = None
        self.t1s = None
        self.t3s = None
        self.rmin = None
        self.rwa = None
        self.oa1 = None
        self.oa3 = None
        self.Uee = None
        self.Uc0 = None
        
        self.tc = 0
        
       
    def _vprint(self, *args, **kwargs):
        """Prints a string if the self.verbose attribute is True
        
        """
        if self.verbose:
            print(*args, **kwargs)

            
    def bootstrap(self, rwa=0.0, pad=0, lab=None, verbose=False, 
                  write_resp=False, keep_resp=False):
        """Sets up the environment for 2D calculation
        write_resp takes a string, creates a directory with the name of
        the string and saves the respoonses and time axis as a npz file
        
        keep_resp saves the responses as a list of dictionaries. The 
        list goes through the time points in t2.
        
        """

        self.verbose = verbose
        self.pad = pad
        self.write_resp = write_resp
        self.keep_resp = keep_resp
        self.rwa = qr.Manager().convert_energy_2_internal_u(rwa)

        with qr.energy_units("int"):

            if self.write_resp:
                try:
                    os.mkdir(write_resp)
                except OSError:
                    print ("Creation of the directory failed, "+
                           "it either already exists "+
                           "or you didn't give a string")
    
              
            # calculate 2D spectrum using aceto library
              
            ###############################################################
            #
            # Create band_system from quantarhei classes
            #
            ###############################################################
            
            if self._has_system:
            
                if isinstance(self.system, (Aggregate, OpenSystem)):
                
                    pass
                
                else:
                    
                    raise Exception("Molecule 2D not implememted")
                    
                sys = self.system
                sys.diagonalize()
                
                #
                # hamiltonian and transition dipole moment operators
                #
                H = sys.get_Hamiltonian()
                D = sys.get_TransitionDipoleMoment()
                
                #
                # Construct band_system object
                #
                Nb = 3
                Ns = numpy.zeros(Nb, dtype=numpy.int32)
                Ns[0] = sys.Nb[0] #1
                Ns[1] = sys.Nb[1] #agg.nmono
                Ns[2] = sys.Nb[2] #Ns[1]*(Ns[1]-1)/2
                
                _use_aceto = False

                
                #
                # Relaxation rates
                #
                if not self._has_rate_matrix:
                    #try:
                    if False:
                        KK = sys.get_RedfieldRateMatrix()
                    #except:
                    else:

                        # FIXME: This is a quick fix to make a zero rate matrix
                        class hlp:
                            def __init__(self, N):
                                self.data =  numpy.zeros((N,N), dtype=REAL)

                        KK = hlp(Ns[1])
                else:
                    KK = self._rate_matrix
                
                # relaxation rate in single exciton band
                Kr = KK.data[Ns[0]:Ns[0]+Ns[1],Ns[0]:Ns[0]+Ns[1]] #*10.0
                #print(1.0/KK.data)
                
                # FIXME: we need also 2 exciton rates
                #                
                
                #
                # Lineshape functions
                #
                sbi = sys.get_SystemBathInteraction()
                cfm = sbi.CC
                cfm.create_double_integral()
                             

#
#  This section will also be removed - It goes to the new Response class
#        

                #
                # Finding population evolution matrix
                #
                prop = PopulationPropagator(self.t2axis, Kr)
                #Uee, Uc0 = prop.get_PropagationMatrix(self.t2axis,
                #                                     corrections=True)
                self.Uee, cor = prop.get_PropagationMatrix(self.t2axis,
                                                      corrections=3)
                  
                # FIXME: Order of transfer is set by hand here 
                # - needs to be moved to some reasonable place
                
                #Ucor = Uee
                self.Uc0 = cor[0]
                
###############################################################################
             
            #
            # bootstrap responses
            #
            if self._has_responses:
                
                for rsp in self.resp_fcions:
                    rsp.set_rwa(self.rwa)
                    
            elif self._has_system:
                
                # FIXME: create responses from system
                
                pass
                # FIXME: set _has_responses to True after they are calculated
                    
             
            #
            # define lab settings
            #
            if lab is None:
                self.lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
                X = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
                self.lab.set_laser_polarizations(X,X,X,X)
            else:
                self.lab = lab
            
            #
            # Other parameters
            #
            #dt = self.t1axis.step
            self.rmin = 0.0001
            self.t1s = self.t1axis.data 
            self.t3s = self.t3axis.data
             
            atype = self.t1axis.atype
            self.t1axis.atype = 'complete'
            self.oa1 = self.t1axis.get_FrequencyAxis() 
            self.oa1.data += self.rwa
            self.oa1.start += self.rwa
            #print(self.oa1.start, self.oa1.data[0])
            self.t1axis.atype = atype
              
            atype = self.t3axis.atype
            self.t3axis.atype = 'complete'
            self.oa3 = self.t3axis.get_FrequencyAxis() 
            self.oa3.data += self.rwa
            self.oa3.start += self.rwa
            #print(self.oa3.start, self.oa3.data[0])
            self.t3axis.atype = atype
            
            
            self.tc = 0

        
    def reset_t2_time(self):
        """Resets the population time of the calculations
        
        """
        self.tc = 0
        

    def calculate_next(self):
        """ Calculate next population time of a 2D spectrum
        
        """
        sone = self.calculate_one(self.tc)
        self.tc += 1
        return sone
    
        
    def calculate_one(self, tc):
        """ Calculate one population time
        
        
        """
        
        try:
            tt2 = self.t2axis.data[tc]  
        except:
            print("Time axis error:\n"+
                  "  perhaps tc =", tc, " (representing t2 population time) is outside range?")
            print("You can reset automatic calculation along the population time axes by calling:")
            print("> twodcalc.reset_t2_time() ")
            print("where 'twodcalc' is the TwoDResponseCalculator object.")
        
        
        Nt2 = self.t2axis.length
        Nr1 = self.t1axis.length
        Nr3 = self.t3axis.length   
        
        # FIXME: on which axis we should be looking for it2 ??? 
        (it2, err) = self.t2axis.locate(tt2) 
        self._vprint("t2 = "+str(tt2)+"fs (it2 = "+str(it2)
                     +" of "+str(Nt2)+")", end="\r")      
        
        #
        # Initialize response storage
        #
        #if _have_aceto and self._has_system:
        #    order = 'F'
        #else:
        order = 'C'
            
        ntype = numpy.complex128
        
        
        # FIXME:  Fix the axis of time
        # the order of axis is wrong. 2D code works only if it is Nr3, Nr1
        
        resp_r = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_n = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Rgsb = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Ngsb = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)

        resp_Rse = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Nse = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Resa = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Nesa = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Rsewt = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Nsewt = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Resawt = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Nesawt = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        
    
        if False: #_have_aceto and self._has_system:
            #
            # calcute response
            #
            self._vprint("calculating response: ")
    
            t1 = time.time()
    
            self._vprint(" - ground state bleach")
            # GSB
            nr3td.nr3_r3g(self.lab, self.sys, it2, self.t1s, self.t3s, 
                          self.rwa, self.rmin, resp_Rgsb)
            nr3td.nr3_r4g(self.lab, self.sys, it2, self.t1s, self.t3s, 
                          self.rwa, self.rmin, resp_Ngsb)
    
            self._vprint(" - stimulated emission")
            # SE
            nr3td.nr3_r2g(self.lab, self.sys, it2, self.t1s, self.t3s, 
                          self.rwa, self.rmin, resp_Rse)
            nr3td.nr3_r1g(self.lab, self.sys, it2, self.t1s, self.t3s, 
                          self.rwa, self.rmin, resp_Nse)
    
            self._vprint(" - excited state absorption")
            # ESA
            nr3td.nr3_r1fs(self.lab, self.sys, it2, self.t1s, self.t3s, 
                           self.rwa, self.rmin, resp_Resa)
            nr3td.nr3_r2fs(self.lab, self.sys, it2, self.t1s, self.t3s, 
                           self.rwa, self.rmin, resp_Nesa)
            
            #
            # Transfer
            #
            Utr = self.Uee[:,:,self.tc] - self.Uc0[:,:,self.tc] 
                                      # -Uc1[:,:,tc]-Uc2[:,:,tc] etc.
            self.sys.set_population_propagation_matrix(Utr)
            
            self._vprint(" - stimulated emission with transfer")
            # SE
            nr3td.nr3_r2g_trans(self.lab, self.sys, it2, self.t1s, self.t3s,
                                self.rwa, self.rmin, resp_Rsewt)
            nr3td.nr3_r1g_trans(self.lab, self.sys, it2, self.t1s, self.t3s,
                                self.rwa, self.rmin, resp_Nsewt)
    
            ## This contributes only when No > 0
            #nr3td.nr3_r2g_trN(lab, sys, No, it2, t1s, t3s, rwa, rmin, resp_r)
            #
    
            self._vprint(" - excited state absorption with transfer")
            # ESA
            nr3td.nr3_r1fs_trans(self.lab, self.sys, it2, self.t1s, self.t3s,
                                 self.rwa, self.rmin, resp_Resawt)
            nr3td.nr3_r2fs_trans(self.lab, self.sys, it2, self.t1s, self.t3s,
                                 self.rwa, self.rmin, resp_Nesawt)
    
            t2 = time.time()
            self._vprint("... calculated in "+str(t2-t1)+" sec")

            self._has_responses = False  # responses are ignored if we use aceto
            


        if self._has_system and not self._has_responses:
            #
            # Calculating all responses from the system
            #
            self.resp_fcions = []

            # basic pathways
            Nr1g = qr.NonLinearResponse(self.lab, self.system, "R1g",
                                        self.t1axis, self.t2axis, self.t3axis)
            Nr2g = qr.NonLinearResponse(self.lab, self.system, "R2g",
                                        self.t1axis, self.t2axis, self.t3axis)
            Nr3g = qr.NonLinearResponse(self.lab, self.system, "R3g",
                                        self.t1axis, self.t2axis, self.t3axis)
            Nr4g = qr.NonLinearResponse(self.lab, self.system, "R4g",
                                        self.t1axis, self.t2axis, self.t3axis)

            self.resp_fcions.append(Nr1g)
            self.resp_fcions.append(Nr2g)
            self.resp_fcions.append(Nr3g)
            self.resp_fcions.append(Nr4g)

            if self.system.mult > 1:
                # ESA (if mult > 1)
                Nr1f = qr.NonLinearResponse(self.lab, self.system, "R1f",
                                            self.t1axis, self.t2axis, self.t3axis)
                Nr2f = qr.NonLinearResponse(self.lab, self.system, "R2f",
                                            self.t1axis, self.t2axis, self.t3axis)

                self.resp_fcions.append(Nr1f)
                self.resp_fcions.append(Nr2f)            


            # relaxation (if relax neq 0)
            Nr1g_scM0g = qr.NonLinearResponse(self.lab, self.system, 
                                              "R1g_scM0g", self.t1axis, self.t2axis, self.t3axis)
            Nr2g_scM0g = qr.NonLinearResponse(self.lab, self.system, 
                                              "R2g_scM0g", self.t1axis, self.t2axis, self.t3axis)
            
            KK = Nr1g_scM0g.KK
            if numpy.all(numpy.isclose(KK, 0.0, atol=1e-9)):

                pass  # we avoid calculating relaxation if the matrix is zero

            else:

                #print("Including relaxation")
                self.resp_fcions.append(Nr1g_scM0g)
                self.resp_fcions.append(Nr2g_scM0g)            


            if self.system.mult > 1:
                Nr1f_scM0g = qr.NonLinearResponse(self.lab, self.system, 
                                                "R1f_scM0g", self.t1axis, self.t2axis, self.t3axis)
                Nr2f_scM0g = qr.NonLinearResponse(self.lab, self.system, 
                                                "R2f_scM0g", self.t1axis, self.t2axis, self.t3axis)
                Nr1f_scM0e = qr.NonLinearResponse(self.lab, self.system, 
                                                "R1f_scM0e", self.t1axis, self.t2axis, self.t3axis)
                Nr2f_scM0e = qr.NonLinearResponse(self.lab, self.system, 
                                                "R2f_scM0e", self.t1axis, self.t2axis, self.t3axis)

                self.resp_fcions.append(Nr1f_scM0g)
                self.resp_fcions.append(Nr2f_scM0g)
                self.resp_fcions.append(Nr1f_scM0e)
                self.resp_fcions.append(Nr2f_scM0e)


            self._has_responses = True
                


        if self._has_responses:
            #
            # Calculation from predefined non-linear responses
            #
                
            for resp in self.resp_fcions:
                
                if isinstance(resp, NonLinearResponse):
                    
                    if resp.rtype == "R":
                        
                        resp_Rgsb += resp.calculate_matrix(tt2)
    
                    elif resp.rtype == "NR":
                        
                        resp_Ngsb += resp.calculate_matrix(tt2)
                        
                    else:
                        
                        raise Exception("Unknown response type")
                
                
                elif isinstance(resp, LiouvillePathway):
                
                    if resp.rtype == "R":
                        
                        resp_Rgsb += resp.calculate_matrix(self.lab, None, tt2, 
                                                           self.t1s, self.t3s, 
                                                           self.rwa)
    
                    elif resp.rtype == "NR":
                        
                        resp_Ngsb += resp.calculate_matrix(self.lab, None, tt2, 
                                                           self.t1s, self.t3s, 
                                                           self.rwa)
                    else:
                        
                        raise Exception("Unknown response type")
                    
            
        else:
            
            raise Exception("Calculation method not implemented")
            

        # only for Aceto we need the sum 
        #
        # FIXME: discontinue Aceto and remove the sum (and the code above)
        #
        resp_r = resp_Rgsb #+ resp_Rse + resp_Resa + resp_Rsewt + resp_Resawt
        resp_n = resp_Ngsb #+ resp_Nse + resp_Nesa + resp_Nsewt + resp_Nesawt



        #
        # Calculate corresponding 2D spectrum
        #
        onetwod = TwoDResponse()

        # pad is set to 0 by default. If changed in the bootstrap,
        # responses are padded with 0s and the time axis is lengthened
        t13Pad = TimeAxis(self.t1axis.start, self.t1axis.length + self.pad,
                          self.t1axis.step)
        if self.pad > 0:
            self._vprint('padding by - ' + str(self.pad))

            t13Pad.atype = 'complete'
            t13PadFreq = t13Pad.get_FrequencyAxis()
            t13PadFreq.data += self.rwa
            t13PadFreq.start += self.rwa

            onetwod.set_axis_1(t13PadFreq)
            onetwod.set_axis_3(t13PadFreq)

            # Sloping the end of the data down to 0 so there isn't a hard
            # cutoff at the end of the data
            from scipy.signal import windows as sig
            window = 20
            tuc = sig.tukey(window * 2, 1, sym = False)
            for k in range(len(resp_r)):
                resp_r[len(resp_r)-window:,k] *= tuc[window:]
                resp_r[k,len(resp_r)-window:] *= tuc[window:]
                resp_n[len(resp_n)-window:,k] *= tuc[window:]
                resp_n[k,len(resp_n)-window:] *= tuc[window:]

            resp_r = numpy.hstack((resp_r, 
                                   numpy.zeros((resp_r.shape[0], self.pad))))
            resp_r = numpy.vstack((resp_r, 
                                   numpy.zeros((self.pad, resp_r.shape[1]))))
            resp_n = numpy.hstack((resp_n, 
                                   numpy.zeros((resp_n.shape[0], self.pad))))
            resp_n = numpy.vstack((resp_n, 
                                   numpy.zeros((self.pad, resp_n.shape[1]))))

        else:
            onetwod.set_axis_1(self.oa1)
            onetwod.set_axis_3(self.oa3)

        if self.keep_resp:
            resp = {
                'time': self.t1axis.data, 'time_pad': t13Pad.data,
                'rTot': resp_r, 'nTot': resp_n,
                'rGSB': resp_Rgsb, 'nGSB': resp_Ngsb,
                'rSE': resp_Rse, 'nSE': resp_Nse,
                'rESA': resp_Resa, 'nESA': resp_Nesa,
                'rSEWT': resp_Rsewt, 'nSEWT': resp_Nsewt,
                'rESAWT': resp_Resawt, 'nESAWT': resp_Nesawt
                }
            self.responses.append(resp)

        if self.write_resp:
            numpy.savez('./'+self.write_resp+'/respT'+str(int(tt2))+'.npz',
                time = self.t1axis.data, time_pad=t13Pad.data,
                rTot=resp_r, nTot=resp_n,
                rGSB=resp_Rgsb, nGSB=resp_Ngsb,
                rSE=resp_Rse, nSE=resp_Nse,
                rESA=resp_Resa, nESA=resp_Nesa,
                rSEWT=resp_Rsewt, nSEWT=resp_Nsewt,
                rESAWT=resp_Resawt, nESAWT=resp_Nesawt)

        # FIXME: This only applies when 
        resp_r[:,0] = resp_r[:,0]*0.5
        resp_n[:,0] = resp_n[:,0]*0.5
        resp_r[0,:] = resp_r[0,:]*0.5
        resp_n[0,:] = resp_n[0,:]*0.5       
        
        ftresp = numpy.fft.fft(resp_r,axis=1)   # \omega_1
        ftresp = numpy.fft.ifft(ftresp,axis=0)  # \omega_3        
        reph2D = numpy.fft.fftshift(ftresp)
        
        ftresp = numpy.fft.ifft(resp_n,axis=1)*ftresp.shape[1]  # \omega_1
        ftresp = numpy.fft.ifft(ftresp,axis=0)                  # \omega_3
        nonr2D = numpy.fft.fftshift(ftresp)


        onetwod.set_resolution("signals")
        onetwod._add_data(reph2D, dtype=signal_REPH)
        onetwod._add_data(nonr2D, dtype=signal_NONR)

        onetwod.set_t2(self.t2axis.data[tc])

        return onetwod            
                
            
    def calculate(self):
        """Returns 2D spectrum
        
        Calculates and returns TwoDSpectrumContainer containing 2D spectrum
        based on the parameters specified in this object.
        
        
        """            
        # FIXME: we will later use only one branch below
        from .twodcontainer import TwoDResponseContainer
 
                  
        # if _have_aceto and self._has_system:

        #     twods = TwoDResponseContainer(self.t2axis)
            
        #     teetoos = self.t2axis.data
            
        #     kk = 0
        #     Nk = teetoos.shape[0]
            
        #     for tt2 in teetoos:

        #         self._vprint_r(" Calculating t2 =", tt2, "fs (",kk,"of",Nk,")")
        #         onetwod = self.calculate_next()
        #         twods.set_spectrum(onetwod)   
        #         kk += 1
            
        #     return twods
        
        if self._has_responses or self._has_system:
            
            # calculate user defined responses

            twods = TwoDResponseContainer(self.t2axis)
            
            teetoos = self.t2axis.data
 
            kk = 0
            Nk = teetoos.shape[0]
            for tt2 in teetoos:

                onetwod = self.calculate_next()
                twods.set_spectrum(onetwod) 
                kk += 1
            
            return twods
        
        else: 
            
            raise Exception("2D calculation in this mode not implemented.")


    def reset_evaluation_functions(self, fcions):
        """Resets the evaluation functions used by the reseponse functions
        
        
        
        
        """
        
        if self._has_responses:
            
            for rsp, fce in zip(self.resp_fcions, fcions):
                rsp.set_evaluation_function(fce)
                print(rsp)
            
        else:
            raise Exception("Calculatore has no responses defined.")
            