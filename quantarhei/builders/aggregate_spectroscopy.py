# -*- coding: utf-8 -*-
"""
    Class comprising the aggregate methods for support of spectroscopic simulations



    Class Details
    -------------

"""

import numpy

from .aggregate_base import AggregateBase
from ..spectroscopy import diagramatics as diag
from ..core.managers import eigenbasis_of

import quantarhei as qr

class AggregateSpectroscopy(AggregateBase):
    """ Class comprising the aggregate methods for support of spectroscopic simulations


    """
    
    
    ########################################################################
    #
    #   SPECTROSCOPY
    #
    ########################################################################
                       
    def liouville_pathways_3(self, ptype="R3g", dtol=0.01, ptol=1.0e-3, 
                             lab=None, verbose=0):
        """ Generator of Liouville pathways """

        ham = self.get_Hamiltonian()
        return self.liouville_pathways_3T(ptype, dtol=dtol, ptol=ptol, lab=lab,
                                eUt2=qr.qm.SOpUnity(dim=ham.dim),
                                verbose=verbose)
        #
        # Rest is ignored for now (may be valuabel in the future)
        #
        
        pop_tol = ptol
        dip_tol = numpy.sqrt(self.D2_max)*dtol
        
        # Check if the ptype is a tuple
        if not isinstance(ptype, (tuple,list)):
            ptype_tuple = (ptype,)
        else:
            ptype_tuple = ptype
        lst = []
         
        for ptp in ptype_tuple:
        
            if ptp == "R3g":
            
                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                
                k = 0
                l = 0
                for i1g in ngs:
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:

                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break
                            
                            for i3g in ngs:
                                
                                if self.D2[i3g,i2e] < dip_tol:
                                    break
                                
                                for i4e in nes:
                            
                                    if ((self.D2[i4e,i1g] < dip_tol)
                                    and (self.D2[i3g,i4e] < dip_tol)) :
                                        break
                                   
                                    l += 1

                                    #      Diagram R3g
                                    #
                                    #                                     
                                    #      |g_i3> <g_i3|
                                    # <----|-----------|
                                    #      |e_i4> <g_i3|
                                    # ---->|-----------|
                                    #      |g_i1> <g_i3|
                                    #      |-----------|---->
                                    #      |g_i1> <e_i2|
                                    #      |-----------|<----
                                    #      |g_i1> <g_i1|
                                    
                                    try:

                                        lp = \
                                        diag.liouville_pathway("R", i1g,
                                               aggregate=self,
                                               order=3,pname=ptp)         
                                        # |g_i1> <g_i1|
                                        lp.add_transition((i2e,i1g),-1)
                                        # |g_i1> <e_i2|
                                        lp.add_transition((i3g,i2e),-1)
                                        # |g_i1> <g_i3|
                                        lp.add_transition((i4e,i1g),+1)
                                        # |e_i5> <g_i3|
                                        lp.add_transition((i3g,i4e),+1)
                                        # |g_i3> <g_i3|

                                    except:
                                        
                                        
                                        break
                                    
                                    lp.build()
                                    lst.append(lp)
                                    k += 1
            
            if ptp == "R2g":

                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                
                k = 0
                l = 0
                for i1g in ngs:
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break
                            
                            for i3e in nes:
                                
                                if self.D2[i3e,i1g] < dip_tol:
                                    break
                                
                                for i4g in ngs:

                                    if ((self.D2[i4g,i2e] < dip_tol)
                                     or (self.D2[i4g,i3e] < dip_tol)):
                                        break
                                    
                                    l += 1

                                    #      Diagram R2g
                                    #
                                    #                                     
                                    #      |g_i4> <g_i4|
                                    # <----|-----------|
                                    #      |e_i3> <g_i4|
                                    #      |-----------|---->
                                    #      |e_i3> <e_i2|
                                    # ---->|-----------|
                                    #      |g_i1> <e_i2|
                                    #      |-----------|<----
                                    #      |g_i1> <g_i1|
                                    
                                    try:
                                        lp = \
                                        diag.liouville_pathway("R", i1g,
                                                           aggregate=self,
                                                           order=3,pname=ptp,
                                                           popt_band=1)
                                        #      |g_i1> <g_i1|
                                        lp.add_transition((i2e,i1g),-1)
                                        #      |g_i1> <e_i2|
                                        lp.add_transition((i3e,i1g),+1)
                                        #      |e_i3> <e_i2|
                                        lp.add_transition((i4g,i2e),-1)
                                        #      |e_i3> <g_i4|
                                        lp.add_transition((i4g,i3e),+1)
                                        #      |g_i4> <g_i4|

                                    except:
                                        
                                        break
                                    
                                    lp.build()
                                    lst.append(lp)
                                    k += 1
                
             
            if ptp == "R1g":
                
                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                
                #nrg = len(ngs)
                #nre = len(nes) 
                
                #print("Ground state : ", nrg)
                #print("Excited state: ", nre)
                #print("R1g: ",nrg*nre*nre*nrg)
                
                k = 0
                l = 0
                for i1g in ngs:
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break

                            for i3e in nes:

                                if self.D2[i3e,i1g] < dip_tol:
                                    break

                                for i4g in ngs:

                                    if ((self.D2[i4g,i3e] < dip_tol)
                                     or (self.D2[i4g,i2e] < dip_tol)):
                                        break

                                    l += 1

                                    #      Diagram R1g
                                    #
                                    #                                     
                                    #      |g_i4> <g_i4|
                                    # <----|-----------|
                                    #      |e_i2> <g_i4|
                                    #      |-----------|---->
                                    #      |e_i2> <e_i3|
                                    #      |-----------|<----
                                    #      |e_i2> <g_i1|
                                    # ---->|-----------|
                                    #      |g_i1> <g_i1|
                                
                                    try:
                                        lp = \
                                        diag.liouville_pathway("NR",i1g,
                                                           aggregate=self,
                                                           order=3,pname=ptp,
                                                           popt_band=1)
                                        #      |g_i1> <g_i1|                                                           
                                        lp.add_transition((i2e,i1g),+1)
                                        #      |e_i2> <g_i1|        
                                        lp.add_transition((i3e,i1g),-1)
                                        #      |e_i2> <e_i3|
                                        lp.add_transition((i4g,i3e),-1)
                                        #      |e_i2> <g_i4|
                                        lp.add_transition((i4g,i2e),+1)
                                        #      |g_i4> <g_i4|

                                    except:
                                        
                                        break
                                    
                                    lp.build()
                                    lst.append(lp)
                                    k += 1
            
            if ptp == "R4g":
                
                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                
                #nrg = len(ngs)
                #nre = len(nes) 
                
                #print("Ground state : ", nrg)
                #print("Excited state: ", nre)
                #print("R4g: ",nrg*nre*nrg*nrg*nre)
                
                k = 0
                l = 0
                for i1g in ngs:
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break
                            
                            for i3g in ngs:

                                if self.D2[i3g,i2e] < dip_tol:
                                    break
                                
                                for i4e in nes:

                                    if ((self.D2[i4e,i3g] < dip_tol)
                                     or (self.D2[i1g,i4e] < dip_tol)):
                                        break
                                    
                                    l += 1
                                    

                                    #      Diagram R4g
                                    #
                                    #                                     
                                    #      |g_i1> <g_i1|
                                    # <----|-----------|
                                    #      |e_i4> <g_i1|
                                    # ---->|-----------|
                                    #      |g_i3> <g_i1|
                                    # <----|-----------|
                                    #      |e_i2> <g_i1|
                                    # ---->|-----------|
                                    #      |g_i1> <g_i1|

                                    try:
                                        lp = \
                                        diag.liouville_pathway("NR",i1g,
                                                           aggregate=self,
                                                           order=3,pname=ptp)
                                        #      |g_i1> <g_i1|                                                           
                                        lp.add_transition((i2e,i1g),+1)
                                        #      |e_i2> <g_i1|
                                        lp.add_transition((i3g,i2e),+1)
                                        #      |g_i3> <g_i1|
                                        lp.add_transition((i4e,i3g),+1)
                                        #      |e_i4> <g_i1|
                                        lp.add_transition((i1g,i4e),+1)
                                        #      |g_i1> <g_i1|

                                    except:
                                        
                                        break
                                    
                                    lp.build()
                                    lst.append(lp)
                                    k += 1
            
            if ptp == "R1f*":
                
                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                try:
                    nfs = self.get_excitonic_band(band=2)
                except:
                    break
                
#                print(ngs)
#                print(nes)
#                print(nfs)
#                for a in nes:
#                    for b in nfs:
#                        print(a,b," : ",self.D2[a,b],self.D2[b,a])

                k = 0
                l = 0
                for i1g in ngs:
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break
                            
                            for i3e in nes:
                                
                                if self.D2[i3e,i1g] < dip_tol:
                                    break
                                
                                for i4f in nfs:

                                    if ((self.D2[i4f,i3e] < dip_tol)
                                     or (self.D2[i2e,i4f] < dip_tol)):
                                        #print("Breaking")
                                        #print(self.D2[i4f,i3e],self.D2[i2e,i4f])
                                        break
                                    
                                    l += 1
                                    

                                    #      Diagram R4g
                                    #
                                    #                                     
                                    #      |e_i2> <e_i2|
                                    # <----|-----------|
                                    #      |f_i4> <e_i2|
                                    # ---->|-----------|
                                    #      |e_i3> <e_i2|
                                    # ---->|-----------|
                                    #      |g_i1> <e_i2|
                                    #      |-----------|<----
                                    #      |g_i1> <g_i1|

                                    try:

                                        lp = \
                                        diag.liouville_pathway("R",i1g,
                                                           aggregate=self,
                                                           order=3,pname=ptp,
                                                           popt_band=1)
                                        #      |g_i1> <g_i1|                                                           
                                        lp.add_transition((i2e,i1g),-1)
                                        #      |g_i1> <e_i2|
                                        lp.add_transition((i3e,i1g),+1)
                                        #      |e_i3> <e_i2|
                                        lp.add_transition((i4f,i3e),+1)
                                        #      |f_i4> <e_i2|
                                        lp.add_transition((i2e,i4f),+1)
                                        #      |e_i2> <e_i2|

                                    except:
                                        
                                        break
                                    
                                    lp.build()
                                    lst.append(lp)
                                    k += 1
            
            if ptp == "R2f*":
                
                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                
                try:
                    nfs = self.get_excitonic_band(band=2)
                except:
                    break
                

                k = 0
                l = 0
                for i1g in ngs:
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break                                                        
                            
                            for i3e in nes:
                                
                                if self.D2[i3e,i1g] < dip_tol:
                                    break

                                for i4f in nfs:

                                    if ((self.D2[i4f,i2e] < dip_tol)
                                     or (self.D2[i3e,i4f] < dip_tol)):
                                        break
                                    
                                    l += 1
                                    

                                    #      Diagram R4g
                                    #
                                    #                                     
                                    #      |e_i3> <e_i3|
                                    # <----|-----------|
                                    #      |f_i4> <e_i3|
                                    # ---->|-----------|
                                    #      |e_i2> <e_i3|
                                    #      |-----------|<----
                                    #      |e_i2> <g_i1|
                                    # ---->|-----------|
                                    #      |g_i1> <g_i1|

                                    try:

                                        lp = \
                                        diag.liouville_pathway("NR",i1g,
                                                           aggregate=self,
                                                           order=3,pname=ptp,
                                                           popt_band=1)
                                        #      |g_i1> <g_i1|                                                           
                                        lp.add_transition((i2e,i1g),+1)
                                        #      |e_i2> <g_i1|
                                        lp.add_transition((i3e,i1g),-1)
                                        #      |e_i2> <e_i3|
                                        lp.add_transition((i4f,i2e),+1)
                                        #      |f_i4> <e_i3|
                                        lp.add_transition((i3e,i4f),+1)
                                        #      |e_i3> <e_i3|

                                    except:
                                        
                                        break
                                    
                                    lp.build()
                                    lst.append(lp)
                                    k += 1
                    

            if ptp == "R2g->3g":

                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                
                k = 0
                l = 0
                for i1g in ngs:
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break
                            
                            for i3e in nes:
                                
                                if self.D2[i3e,i1g] < dip_tol:
                                    break
                                
                                # relaxation 
                                for i4g in ngs:
                                    for i5g in ngs:
                                
                                        for i6e in nes:

                                            if ((self.D2[i6e,i4g] < dip_tol)
                                            or (self.D2[i5g,i6e] < dip_tol)):
                                                break
                                    
                                            l += 1

                                    #      Diagram R2g_ETICS
                                    #      (Compensates R3g)
                                    #
                                    #                                     
                                    #      |g_i5> <g_i5|
                                    # <----|-----------|
                                    #      |e_i6> <g_i5|
                                    # ---->|-----------|
                                    #      |g_i4> <g_i5|
                                    #      |***********|
                                    #      |e_i3> <e_i2|
                                    # ---->|-----------|
                                    #      |g_i1> <e_i2|
                                    #      |-----------|<----
                                    #      |g_i1> <g_i1|
                                    
                                            if True:
                                            #try:
                                                lp = \
                                                diag.liouville_pathway("R_E",
                                                           i1g,
                                                           aggregate=self,
                                                           order=3,
                                                           relax_order=1,
                                                           pname=ptp)
                                                #      |g_i1> <g_i1|
                                                lp.add_transition((i2e,i1g),-1)
                                                #      |g_i1> <e_i2|
                                                lp.add_transition((i3e,i1g),+1)
                                                #      |e_i3> <e_i2|
                                                lp.add_transfer((i4g,i5g),
                                                                  (i3e,i2e))
                                                #      |g_i4> <g_i5|
                                                lp.add_transition((i6e,i4g),+1)
                                                #      |e_i6> <g_i5|
                                                lp.add_transition((i5g,i6e),+1)
                                                #      |g_i5> <g_i5|

                                            #except:
                                        
                                            #    break
                                    
                                            lp.build()
                                            lst.append(lp)
                                            k += 1

            if ptp == "R1g->4g":

                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                
                k = 0
                l = 0
                for i1g in ngs:
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break
                            
                            for i3e in nes:
                                
                                if self.D2[i3e,i1g] < dip_tol:
                                    break
                                
                                # relaxation 
                                for i4g in ngs:
                                    for i5g in ngs:
                                
                                        for i6e in nes:

                                            if ((self.D2[i6e,i4g] < dip_tol)
                                            or (self.D2[i5g,i6e] < dip_tol)):
                                                break
                                    
                                            l += 1

                                    #      Diagram R2g_ETICS
                                    #      (Compensates R3g)
                                    #
                                    #                                     
                                    #      |g_i5> <g_i5|
                                    # <----|-----------|
                                    #      |e_i6> <g_i5|
                                    # ---->|-----------|
                                    #      |g_i4> <g_i5|
                                    #      |***********|
                                    #      |e_i2> <e_i3|
                                    #      |-----------|<----
                                    #      |e_i2> <g_i1|
                                    # ---->|-----------|
                                    #      |g_i1> <g_i1|
                                    
                                            #if True:
                                            try:
                                                lp = \
                                                diag.liouville_pathway("NR_E",
                                                           i1g,
                                                           aggregate=self,
                                                           order=3,
                                                           relax_order=1,
                                                           pname=ptp)
                                                #      |g_i1> <g_i1|
                                                lp.add_transition((i2e,i1g),+1)
                                                #      |e_i2> <g_i1|
                                                lp.add_transition((i3e,i1g),-1)
                                                #      |e_i2> <e_i3|
                                                lp.add_transfer((i4g,i5g),
                                                                  (i2e,i3e))
                                                #      |g_i4> <g_i5|
                                                lp.add_transition((i6e,i4g),+1)
                                                #      |e_i6> <g_i5|
                                                lp.add_transition((i5g,i6e),+1)
                                                #      |g_i5> <g_i5|

                                            except:
                                        
                                                break
                                    
                                            lp.build()
                                            lst.append(lp)
                                            k += 1

         
        if lab is not None:
            for l in lst:
                l.orientational_averaging(lab)
         
        return lst     
                

    
    def liouville_pathways_3T(self, ptype="R3g", eUt=None, ham=None, t2=0.0,
                              dtol=1.0e-12, ptol=1.0e-3, etol=1.0e-6,
                              verbose=0, lab=None):
        """ Generator of Liouville pathways with energy transfer
        
        
        
        
        Parameters
        ----------
        
        ptype : tuple, list, str
            List of strings or a string representing one or more
            Liouville pathway types that are to be calculated
            
        eUt : EvolutionSuperOperator
            Evolution superoperator representing the energy 
            transfer in the system 
            
        t2 : float
            Waiting time at which the spectrum is calculated
            
        dtol : float
            Minimum acceptable strength of the transition from ground
            to excited state, relative to the maximum dipole strength 
            available in the system
            
        ptol : float
            Minimum acceptable population of the ground state (e.g. states
            not thermally populated are excluded)

        lab : LaboratorySetup
            Object representing laboratory setup - number of pulses, 
            polarization etc.
            
        Returns
        -------
        
        lst : list
            List of LiouvillePathway objects
            
            
        """
        if self._diagonalized:
            if verbose > 0:
                print("Diagonalizing aggregate")
            self.diagonalize()
            if verbose > 0:
                print("..done")
        
        pop_tol = ptol
        dip_tol = numpy.sqrt(self.D2_max)*dtol
        evf_tol = etol
        
        # Check if the ptype is a tuple
        if not isinstance(ptype, (tuple,list)):
            ptype_tuple = (ptype,)
        else:
            ptype_tuple = ptype
        lst = []
        
        if verbose > 0:
            print("Pathways", ptype_tuple)
         
        #
        # data of the evolution superoperator in eigenstate basis
        #
        
        try:
            # either the eUt is a complete evolution superoperator
            eUt2 = eUt.at(t2)
            eUt2_dat = numpy.zeros(eUt2.data.shape, dtype=eUt2.data.dtype)
            HH = eUt.get_Hamiltonian()
            with eigenbasis_of(HH):
                eUt2_dat[:,:,:,:] = eUt2.data
        except:
            # or it is only a super operator at a given time t2
            # in this case 'ham' must be specified
            eUt2 = eUt
            eUt2_dat = numpy.zeros(eUt2.data.shape, dtype=eUt2.data.dtype)
            with eigenbasis_of(ham):
                eUt2_dat[:,:,:,:] = eUt2.data
    
        
        
        
        for ptp in ptype_tuple:
        
            if ptp == "R1g":
                
                generate_R1g(self, lst, eUt2_dat,
                             pop_tol, dip_tol, evf_tol, verbose)

            elif ptp == "R2g":
                
                generate_R2g(self, lst, eUt2_dat,
                             pop_tol, dip_tol, evf_tol, verbose)

            elif ptp == "R3g":
            
                generate_R3g(self, lst, eUt2_dat, pop_tol, dip_tol, verbose)
                
            elif ptp == "R4g":
                
                generate_R4g(self, lst, eUt2_dat, pop_tol, dip_tol, verbose)
                
            
            elif ptp == "R1f*":
                
                generate_R1f(self, lst, eUt2_dat, 
                             pop_tol, dip_tol, evf_tol, verbose)
                
            
            elif ptp == "R2f*":
                
                generate_R2f(self, lst, eUt2_dat, 
                             pop_tol, dip_tol, evf_tol, verbose)
                                   
            else:
                
                raise Exception("Unknown pythway type: "+str(ptp))
                
        
        if lab is not None:
            for l in lst:
                l.orientational_averaging(lab)
         
        return lst     


    def liouville_pathways_1(self, eUt=None, ham=None, dtol=0.01, ptol=1.0e-3,
                             etol=1.0e-6, verbose=0, lab=None):
        """ Generator of the first order Liouville pathways 
        
        
        Generator of the pathways for an absorption spectrum
        calculation.
        
        
        
        Parameters
        ----------
        
            
        eUt : EvolutionSuperOperator
            Evolution superoperator representing the evolution of optical
            coherence in the system 
            
            
        dtol : float
            Minimum acceptable strength of the transition from ground
            to excited state, relative to the maximum dipole strength 
            available in the system
            
        ptol : float
            Minimum acceptable population of the ground state (e.g. states
            not thermally populated are excluded)

        lab : LaboratorySetup
            Object representing laboratory setup - number of pulses, 
            polarization etc.
            
        Returns
        -------
        
        lst : list
            List of LiouvillePathway objects
            
            
        """
        if self._diagonalized:
            if verbose > 0:
                print("Diagonalizing aggregate")
            self.diagonalize()
            if verbose > 0:
                print("..done")
        
        pop_tol = ptol
        dip_tol = numpy.sqrt(self.D2_max)*dtol
        evf_tol = etol
                        
        if eUt is None:
            
            # secular absorption spectrum calculation
            eUt2_dat = None
            sec = True
        
        else:
            
            raise Exception("Not implemented yet")

        lst = []        

        if sec:
            generate_1orderP_sec(self, lst,
                                 pop_tol, dip_tol, verbose)
        else:
            raise Exception("Not implemented yet")                                
        
        if lab is not None:
            for l in lst:
                l.orientational_averaging(lab)
         
        return lst   





        
def generate_R1g(self, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose=0):        
        
    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)
    
    ver = verbose              
    
    if verbose > 0:
        print("Liouville pathway R1g")
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)
        print("Evolution amplitude:  ", evf_tol)
    
    k = 0
    l = 0
    for i1g in ngs:

        if verbose > 0: 
            print("Ground state: ", i1g, "of", len(ngs))
        
        # Only thermally allowed starting states are considered
        if self.rho0[i1g,i1g] > pop_tol:
    
            for i2e in nes:

                if verbose > 1: 
                    print("Excited state: ", i2e, "of", len(nes))
                
                if self.D2[i2e, i1g] > dip_tol:

                    for i3e in nes:
                    
                        if self.D2[i3e, i1g] > dip_tol:
                                        
                            for i2d in nes:
                                for i3d in nes:

                                    evf = eUt2[i2d, i3d, i2e, i3e]
                                    if abs(evf) > evf_tol:

                                        for i4g in ngs:

                                            if ((self.D2[i4g,i3d] > dip_tol)
                                            and (self.D2[i4g,i2d] > dip_tol)):
                                   

                                                l += 1
                                                
                                                lp = _generate_R1g(self, i1g,
                                                                   i2e, i3e,
                                                                   i2d, i3d,
                                                                   i4g, evf,
                                                                   verbose=ver)
                                                
                                                lp.build()
                                                lst.append(lp)
                                                k += 1


def _generate_R1g(self, i1g, i2e, i3e, i2d, i3d, i4g, evf, verbose=0):

    #      Diagram R1g
    #
    #                                     
    #      |g_i4> <g_i4|
    # <----|-----------|
    #      |d_i2> <g_i4|
    #      |-----------|---->
    #      |d_i2> <d_i3|
    #      |***********|
    #      |e_i2> <e_i3|
    #      |-----------|<----
    #      |e_i2> <g_i1|
    # ---->|-----------|
    #      |g_i1> <g_i1|
        
    try:
        if verbose > 5:
            print(" * Generating R1g", i1g, i2e, i3e)
        lp = \
        diag.liouville_pathway("NR",
                   i1g,
                   aggregate=self,
                   order=3, pname="R1g",
                   popt_band=1,
                   relax_order=1)

        # first transition lineshape
        width1 = \
        self.get_transition_width((i2e, i1g))
        deph1 = \
        self.get_transition_dephasing((i2e, 
                   i1g))
        # third transition lineshape
        width3 = \
        self.get_transition_width((i2d, i4g))
        deph3 = \
        self.get_transition_dephasing((i2d,
                   i4g))


        #      |g_i1> <g_i1|                                                           
        lp.add_transition((i2e,i1g),+1,
              interval=1, 
              width=width1, 
              deph=deph1)
        #      |e_i2> <g_i1|        
        lp.add_transition((i3e,i1g),-1)
        #      |e_i2> <e_i3|
        lp.add_transfer(((i2d, i3d)),
             (i2e, i3e))
        lp.set_evolution_factor(evf)
        #      |d_i2> <d_i3|                                                                                
        lp.add_transition((i4g,i3d),-1)
        #      |d_i2> <g_i4|
        lp.add_transition((i4g,i2d),+1,
              interval=3, 
              width=width3, 
              deph=deph3)
        #      |g_i4> <g_i4|

    except:

        raise Exception("Pathway generation failed")
        
    return lp



def generate_R2g(self, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose=0):
    
    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)
    
    if verbose > 0:
        print("Liouville pathway R2g")
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)
        print("Evolution amplitude:  ", evf_tol)
    
    k = 0
    l = 0
    for i1g in ngs:

        if verbose > 0: 
            print("Ground state: ", i1g, "of", len(ngs))

        # Only thermally allowed starting states are considered
        if self.rho0[i1g,i1g] > pop_tol:
    
            for i2e in nes:
                
                if verbose > 1: 
                    print("Excited state: ", i2e, "of", len(nes))
               
                if self.D2[i2e,i1g] > dip_tol:
                
                    for i3e in nes:
                    
                        if self.D2[i3e,i1g] > dip_tol:

                            for i3d in nes:
                                for i2d in nes:
                            
                                    evf = eUt2[i3d, i2d, i3e, i2e]
                                    if abs(evf) > evf_tol:

                    
                                        for i4g in ngs:

                                            if ((self.D2[i4g,i2e] > dip_tol)
                                            and (self.D2[i4g,i3e] > dip_tol)):
                                    
                                
                                                l += 1

                                #      Diagram R2g
                                #
                                #                                     
                                #      |g_i4> <g_i4|
                                # <----|-----------|
                                #      |d_i3> <g_i4|  
                                #      |-----------|---->
                                #      |d_i3> <d_i2|
                                #      |***********|
                                #      |e_i3> <e_i2|
                                # ---->|-----------|
                                #      |g_i1> <e_i2|
                                #      |-----------|<----
                                #      |g_i1> <g_i1|
                                
                                                try:
                                                    if verbose > 5:
                                                        print(" * Generating R2g", i1g, i2e, i3e)
                                                    
                                                    lp = \
                                                    diag.liouville_pathway("R", 
                                                                           i1g,
                                                                aggregate=self,
                                                                order=3,pname="R2g",
                                                                popt_band=1,
                                                                relax_order=1)
                
                                                    # first transition lineshape
                                                    width1 = \
                                            self.get_transition_width((i2e, i1g))
                                                    deph1 = \
                                            self.get_transition_dephasing((i2e, 
                                                                           i1g))
                                                    # third transition lineshape
                                                    width3 = \
                                            self.get_transition_width((i3d, i4g))
                                                    deph3 = \
                                            self.get_transition_dephasing((i3d,
                                                                           i4g))
                                                    
                                        
                                                    #      |g_i1> <g_i1|
                                                    lp.add_transition((i2e,i1g),-1,
                                                                      interval=1, 
                                                                      width=width1, 
                                                                      deph=deph1)
                                                    #      |g_i1> <e_i2|
                                                    lp.add_transition((i3e,i1g),+1)
                                                    #      |e_i3> <e_i2|
                                                    lp.add_transfer(((i3d, i2d)),
                                                                     (i3e, i2e))
                                                    lp.set_evolution_factor(evf)
                                                    lp.add_transition((i4g,i2d),-1)
                                                    #      |e_i3> <g_i4|
                                                    lp.add_transition((i4g,i3d),+1,
                                                                      interval=3, 
                                                                      width=width3, 
                                                                      deph=deph3)
                                                    #      |g_i4> <g_i4|
                
                                                except:
                                                    
                                                    raise Exception()
                                                    break
                                                
                                                lp.build()
                                                lst.append(lp)
                                                k += 1
                    


def generate_R3g(self, lst, eUt2, pop_tol, dip_tol, verbose=0):

    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)
    
    if verbose > 0:
        print("Liouville pathway R3g")
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)
        #print("Evolution amplitude:  ", evf_tol)
        
    
    k = 0
    l = 0
    for i1g in ngs:
        
        if verbose > 0: 
            print("Ground state: ", i1g, "of", len(ngs))
            
        # Only thermally allowed starting states are considered
        if self.rho0[i1g,i1g] > pop_tol:

            for i2e in nes:
                
                if verbose > 1: 
                    print("Excited state: ", i2e, "of", len(nes))
                
                if self.D2[i2e,i1g] > dip_tol:

                    for i3g in ngs:
                    
                        if self.D2[i3g,i2e] > dip_tol:
                            
                            evf = eUt2[i1g, i3g, i1g, i3g]
                    
                            for i4e in nes:
                
                                if ((self.D2[i4e,i1g] > dip_tol)
                                and (self.D2[i3g,i4e] > dip_tol)) :

                       
                                    l += 1

                        #      Diagram R3g
                        #
                        #                                     
                        #      |g_i3> <g_i3|
                        # <----|-----------|
                        #      |e_i4> <g_i3|
                        # ---->|-----------|
                        #      |g_i1> <g_i3|
                        #      |-----------|---->
                        #      |g_i1> <e_i2|
                        #      |-----------|<----
                        #      |g_i1> <g_i1|
                        
                                    try:
                                        if verbose > 5:
                                            print(" * Generating R3g", i1g, i2e)
            
                                        lp = \
                                        diag.liouville_pathway("R", i1g,
                                               aggregate=self,
                                               order=3, pname="R3g")     
            
                                        # first transition lineshape
                                        width1 = \
                                        self.get_transition_width((i2e, i1g))
                                        deph1 = \
                                        self.get_transition_dephasing((i2e, 
                                                                       i1g))
                                        # third transition lineshape
                                        width3 = \
                                        self.get_transition_width((i4e,i3g))
                                        deph3 = \
                                        self.get_transition_dephasing((i4e,
                                                                       i3g))
                                        
                                        # |g_i1> <g_i1|
                                        lp.add_transition((i2e,i1g),-1,
                                                                  interval=1, 
                                                                  width=width1, 
                                                                  deph=deph1)
                                        # |g_i1> <e_i2|
                                        lp.add_transition((i3g,i2e),-1)
                                        # |g_i1> <g_i3|
                                        lp.add_transition((i4e,i1g),+1)
                                        # |e_i5> <g_i3|
                                        lp.add_transition((i3g,i4e),+1,
                                                                  interval=3, 
                                                                  width=width3, 
                                                                  deph=deph3)
                                        # |g_i3> <g_i3|
            
                                        lp.set_evolution_factor(evf)
            
                                    except:

                                        raise Exception("Generation of pathway failed")
                                    
                                    lp.build()
                                    lst.append(lp)
                                    k += 1


def generate_R4g(self, lst, eUt2, pop_tol, dip_tol, verbose=0):
    
    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)
                    
    if verbose > 0:
        print("Liouville pathway R4g")
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)
    
    k = 0
    l = 0
    for i1g in ngs:

        if verbose > 0: 
            print("Ground state: ", i1g, "of", len(ngs))
        
        # Only thermally allowed starting states are considered
        if self.rho0[i1g,i1g] > pop_tol:
    
            for i2e in nes:
                
                if verbose > 1: 
                    print("Excited state: ", i2e, "of", len(nes))
                    
                    #if i2e == 4:
                    #    print("Changing verbosity to 10")
                    #    verbose = 10
                        

                #print(self.D2[i2e,i1g], dip_tol, self.D2[i2e,i1g] > dip_tol)
                if self.D2[i2e,i1g] > dip_tol:

                    for i3g in ngs:

                        if self.D2[i3g,i2e] > dip_tol:
                    
                            evf = eUt2[i1g, i3g, i1g, i3g]

                            for i4e in nes:

                                if ((self.D2[i4e,i3g] > dip_tol)
                                and (self.D2[i1g,i4e] > dip_tol)):
                        
                                    l += 1
                        

                        #      Diagram R4g
                        #
                        #                                     
                        #      |g_i1> <g_i1|
                        # <----|-----------|
                        #      |e_i4> <g_i1|
                        # ---->|-----------|
                        #      |g_i3> <g_i1|
                        # <----|-----------|
                        #      |e_i2> <g_i1|
                        # ---->|-----------|
                        #      |g_i1> <g_i1|

                                    try:
                                        if verbose > 5:
                                            print(" * Generating R4g", i1g, i2e)
                                        lp = \
                                        diag.liouville_pathway("NR",i1g,
                                                           aggregate=self,
                                                           order=3,pname="R4g")
                                                                                                          
                                        # first transition lineshape
                                        width1 = \
                                        self.get_transition_width((i2e, i1g))
                                        deph1 = \
                                        self.get_transition_dephasing((i2e, 
                                                                       i1g))
                                        # third transition lineshape
                                        width3 = \
                                        self.get_transition_width((i4e,i1g))
                                        deph3 = \
                                        self.get_transition_dephasing((i4e,
                                                                       i1g))
                                        
                                        
                                        #      |g_i1> <g_i1|                                                           
                                        lp.add_transition((i2e,i1g),+1,
                                                                  interval=1, 
                                                                  width=width1, 
                                                                  deph=deph1)
                                        #      |e_i2> <g_i1|
                                        lp.add_transition((i3g,i2e),+1)
                                        #      |g_i3> <g_i1|
                                        lp.add_transition((i4e,i3g),+1)
                                        #      |e_i4> <g_i1|
                                        lp.add_transition((i1g,i4e),+1,
                                                                  interval=3, 
                                                                  width=width3, 
                                                                  deph=deph3)
                                        #      |g_i1> <g_i1|
                                        
                                        lp.set_evolution_factor(evf)
            
                                    except:
                                        
                                        break
                                    
                                    lp.build()
                                    lst.append(lp)
                                    k += 1
                #if verbose == 10:
                #    print("////")
                #    qr.stop()                    
                    

def generate_R1f(self, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose=0):
    
    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)
    try:
        nfs = self.get_excitonic_band(band=2)
    except:
        raise Exception("Excited states not available for R1f* pathway"+
                        " generation")
    
    if verbose > 0:
        print("Liouville pathway R1f*")    
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)
        print("Evolution amplitude:  ", evf_tol)
    
    k = 0
    l = 0
    for i1g in ngs:
    
        if verbose > 0: 
            print("Ground state: ", i1g, "of", len(ngs))
        
        # Only thermally allowed starting states are considered
        if self.rho0[i1g,i1g] > pop_tol:
    
            for i2e in nes:
                
                if self.D2[i2e,i1g] > dip_tol:
                
                    for i3e in nes:
                    
                        if self.D2[i3e,i1g] > dip_tol:
                             
                            if verbose > 2:
                                print("Excited state: ", i2e, i3e, "of",
                                      nes[len(nes)-1])
                            for i3d in nes:
                                for i2d in nes:
                            
                                    evf = eUt2[i3d, i2d, i3e, i2e]
                                    if abs(evf) > evf_tol:
                    
                                        for i4f in nfs:
    
                                            if ((self.D2[i4f,i3d] > dip_tol)
                                            and (self.D2[i2d,i4f] > dip_tol)):
                        
                                                l += 1           
    
                        #      Diagram R1f*
                        #
                        #                                     
                        #      |d_i2> <d_i2|
                        # <----|-----------|
                        #      |f_i4> <d_i2|
                        # ---->|-----------|
                        #      |d_i3> <d_i2|
                        #      |***********|
                        #      |e_i3> <e_i2|
                        # ---->|-----------|
                        #      |g_i1> <e_i2|
                        #      |-----------|<----
                        #      |g_i1> <g_i1|
    
                                                try:
                                                    if verbose > 5:
                                                        print(" * Generating R1f*", i1g, i2e)
                                                    lp = \
                                                    diag.liouville_pathway("R",i1g,
                                                               aggregate=self,
                                                               order=3,pname="R1f*",
                                                               popt_band=1,
                                                               relax_order=1)
                                                    
                                                    # first transition lineshape
                                                    width1 = \
                                            self.get_transition_width((i2e, i1g))
                                                    deph1 = \
                                            self.get_transition_dephasing((i2e, 
                                                                           i1g))
                                                    # third transition lineshape
                                                    width3 = \
                                            self.get_transition_width((i4f,i2d))
                                                    deph3 = \
                                            self.get_transition_dephasing((i4f,
                                                                           i2d))
                                                    
                                                    #      |g_i1> <g_i1|                                                           
                                                    lp.add_transition((i2e,i1g),-1,
                                                                      interval=1, 
                                                                      width=width1, 
                                                                      deph=deph1)
                                                    #      |g_i1> <e_i2|
                                                    lp.add_transition((i3e,i1g),+1)
                                                    #      |e_i3> <e_i2|
                                                    lp.add_transfer(((i3d, i2d)),
                                                                     (i3e, i2e))
                                                    lp.set_evolution_factor(evf)
                                                    #      |d_i3> <d_i2|                                        
                                                    lp.add_transition((i4f,i3d),+1)
                                                    #      |f_i4> <d_i2|
                                                    lp.add_transition((i2d,i4f),+1,
                                                                      interval=3, 
                                                                      width=width3, 
                                                                      deph=deph3)
                                                    #      |d_i2> <d_i2|
                    
                                                except:
                                                    
                                                    raise Exception("Construction"+
                                                    "relaxation pathway failed")
                                        
                                                lp.build()
                                                lst.append(lp)
                                                k += 1

def generate_R2f(self, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose=0):
    
    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)
    
    try:
        nfs = self.get_excitonic_band(band=2)
    except:
        raise Exception("Excited states not available for R2f* pathway"+
                        " generation")
    
    if verbose > 0:
        print("Liouville pathway R2f*")
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)
        print("Evolution amplitude:  ", evf_tol)

    k = 0
    l = 0
    for i1g in ngs:

        if verbose > 0: 
            print("Ground state: ", i1g, "of", len(ngs))
        
        # Only thermally allowed starting states are considered
        if self.rho0[i1g,i1g] > pop_tol:
    
            for i2e in nes:
                
                if self.D2[i2e,i1g] > dip_tol:
                
                    for i3e in nes:
                    
                        if self.D2[i3e,i1g] > dip_tol:
                    
                            if verbose > 2:
                                print("Excited state: ", i2e, i3e, "of",
                                      nes[len(nes)-1])
                                
                            for i2d in nes:
                                for i3d in nes:
                            
                                    evf = eUt2[i2d, i3d, i2e, i3e]
                                    if abs(evf) > evf_tol:

                                        for i4f in nfs:

                                            if ((self.D2[i4f,i2d] > dip_tol)
                                            and (self.D2[i3d,i4f] > dip_tol)):
                                
                                                l += 1
                                

                                #      Diagram R2f*
                                #
                                #                                     
                                #      |d_i3> <d_i3|
                                # <----|-----------|
                                #      |f_i4> <d_i3|
                                # ---->|-----------|
                                #      |d_i2> <d_i3|
                                #      |***********|
                                #      |e_i2> <e_i3|
                                #      |-----------|<----
                                #      |e_i2> <g_i1|
                                # ---->|-----------|
                                #      |g_i1> <g_i1|

                                                try:
                                                    if verbose > 5:
                                                        print(" * Generating R1f*", i1g, i2e)
                
                                                    lp = \
                                                    diag.liouville_pathway("NR",
                                                                           i1g,
                                                                aggregate=self,
                                                                order=3,pname="R2f*",
                                                                popt_band=1,
                                                                relax_order=1)
                                                    
                                                    # first transition lineshape
                                                    width1 = \
                                            self.get_transition_width((i2e, i1g))
                                                    deph1 = \
                                            self.get_transition_dephasing((i2e, 
                                                                           i1g))
                                                    # third transition lineshape
                                                    width3 = \
                                            self.get_transition_width((i4f,i3d))
                                                    deph3 = \
                                            self.get_transition_dephasing((i4f,
                                                                           i3d))
                                                    
                                                    #      |g_i1> <g_i1|                                                           
                                                    lp.add_transition((i2e,i1g),+1,
                                                                      interval=1, 
                                                                      width=width1, 
                                                                      deph=deph1)
                                                    #      |e_i2> <g_i1|
                                                    lp.add_transition((i3e,i1g),-1)                                        
                                                    #      |e_i2> <e_i3|
                                                    lp.add_transfer(((i2d, i3d)),
                                                                     (i2e, i3e))
                                                    lp.set_evolution_factor(evf)
                                                    #      |d_i2> <d_i3|                                        
                                                    lp.add_transition((i4f,i2d),+1)
                                                    #      |f_i4> <d_i3|
                                                    lp.add_transition((i3d,i4f),+1,
                                                                      interval=3, 
                                                                      width=width3, 
                                                                      deph=deph3)
                                                    #      |d_i3> <d_i3|
                
                                                except:
                                                    
                                                    break
                                                
                                                lp.build()
                                                lst.append(lp)
                                                k += 1




def generate_1orderP_sec(self, lst, pop_tol, dip_tol, verbose=0):
    
    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)
    
    if verbose > 0:
        print("Liouville pathway of first order")
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)

    k = 0
    l = 0
    for i1g in ngs:

        if verbose > 0: 
            print("Ground state: ", i1g, "of", len(ngs))
        
        # Only thermally allowed starting states are considered
        if self.rho0[i1g,i1g] > pop_tol:
    
            for i2e in nes:
                
                if self.D2[i2e,i1g] > dip_tol:
                
                                
                    l += 1
                                

                    #      Diagram P1
                    #
                    #                                     
                    #      |g_i1> <g_i1|
                    # <----|-----------|
                    #      |e_i2> <g_i1|
                    # ---->|-----------|
                    #      |g_i1> <g_i1|

                    try:
                        if verbose > 5:
                            print(" * Generating P1", i1g, i2e)

                        lp = \
                        diag.liouville_pathway("NR",
                                               i1g,
                                    aggregate=self,
                                    order=1,pname="P1",
                                    popt_band=1,
                                    relax_order=1)
                        
                        # first transition lineshape
                        width1 = \
                            self.get_transition_width((i2e, i1g))
                        deph1 = \
                            self.get_transition_dephasing((i2e, 
                                               i1g))
                        
                        #      |g_i1> <g_i1|                                                           
                        lp.add_transition((i2e,i1g),+1,
                                          interval=1, 
                                          width=width1, 
                                          deph=deph1)
                        #      |e_i2> <g_i1|
                        lp.add_transition((i1g,i2e),+1,
                                          interval=1, 
                                          width=width1, 
                                          deph=deph1)
                        #      |g_i1> <g_i1|

                    except:
                        
                        break
                    
                    lp.build()
                    lst.append(lp)
                    k += 1
