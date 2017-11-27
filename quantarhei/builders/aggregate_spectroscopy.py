# -*- coding: utf-8 -*-

import numpy

from .aggregate_base import AggregateBase
from ..spectroscopy import diagramatics as diag

class AggregateSpectroscopy(AggregateBase):


    ########################################################################
    #
    #   SPECTROSCOPY
    #
    ########################################################################
                       
    def liouville_pathways_3(self, ptype="R3g", dtol=-1.0, ptol=1.0e-3, lab=None):
        """ Generator of Liouville pathways """
        
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
        
    
#    def ETICS_evolution_operator(self,popt,K,proj):
#        """ETICS evolution operator
#        
#        Parameters
#        ----------
#        
#        popt : float
#            Population time
#            
#        K : float
#            Relaxation rate
#            
#        proj : array
#            An array of electronic states representing a projector
#            on the state which
#            gets deexcited in the process of ETICS
#            
#        """
#        
#        # all this uses site basis        
#            
#        
#        # this is the evolution superoperator at time popt
#        U_ETICS = numpy.zeros((self.Nb[0],self.Nb[0],
#                               self.Nb[1],self.Nb[1]),dtype=numpy.complex128)
#        U_EX    = numpy.zeros((self.Nb[1],self.Nb[1]),dtype=numpy.complex128)
#                               
#        
#        # FC factors and diagonalization 
#        Ceg = numpy.zeros((self.Nb[1],self.Nb[0]))
#        Cge = numpy.zeros((self.Nb[0],self.Nb[1]))
#        
#        # Cge
#        if len(proj)>1:
#            raise Exception("Not implemented yet")            
#        a = proj[0]
#        
#        esg = self.get_ElectronicState(self.elsigs[0])
#        ag = 0
#        for vg in esg.vsignatures():
#            vs_g = VibronicState(esg,vg)
#            ae = 0
#            for ie in range(1,self.number_of_electronic_states_in_band(1)+1):
#                ese = self.get_ElectronicState(self.elsigs[ie])
#                for ve in ese.vsignatures():
#                    vs_e = VibronicState(ese,ve)
#                    if True:
#                    #if ie == a:
#                        Cge[ag,ae] = self.fc_factor(vs_g,vs_e)
#                    else:
#                        Cge[ag,ae] = 0.0
#                    ae += 1
#            ag += 1
#        
#        
#        # we use the fact that we are diagonalized
#        S1 = self.S1
#        
#        Cge = numpy.dot(Cge,S1[self.Nb[0]:(self.Nb[0]+self.Nb[1]),
#                               self.Nb[0]:(self.Nb[0]+self.Nb[1])])
#                               
#        # Ceg
#        Ceg = Cge.T 
#        
#        # frequencies between states in the ground and excited state 
#        omegas_g = numpy.zeros((self.Nb[0],self.Nb[0]),dtype=numpy.float64)
#        omegas_e = numpy.zeros((self.Nb[1],self.Nb[1]),dtype=numpy.float64)
#        for a in range(self.Nb[0]):
#            for b in range(self.Nb[0]):
#                omegas_g[a,b] = self.HH[a,a]-self.HH[b,b]
#                
#        for a in range(self.Nb[0],self.Nb[0]+self.Nb[1]):
#            for b in range(self.Nb[0],self.Nb[0]+self.Nb[1]):
#                omegas_e[a-self.Nb[0],b-self.Nb[0]] = self.HH[a,a]-self.HH[b,b]
#                        
#        
#        # some precalculated values
#        eom_me = numpy.exp(-1j*popt*omegas_e)
#        eom_mg = numpy.exp(-1j*popt*omegas_g)
#        eom_pg = numpy.exp( 1j*popt*omegas_g)
#
#
#
#        KK = numpy.zeros((self.Nb[0],self.Nb[0],
#                               self.Nb[1],self.Nb[1]),dtype=numpy.float64)
#
#        for ag in range(self.Nb[0]):
#            for bg in range(self.Nb[0]):
#                for ae in range(self.Nb[1]):
#                    for be in range(self.Nb[1]):
#                        KK[ag,bg,ae,be] = K*Cge[ag,ae]*Ceg[be,bg]
#                        
#        eK = numpy.exp(-0.5*K*popt)   
#
#
#        for ag in range(self.Nb[0]):
#            for bg in range(self.Nb[0]):
#                for ae in range(self.Nb[1]):
#                    for be in range(self.Nb[1]):
#                        
#                        U_ETICS[ag,bg,ae,be] = KK[ag,bg,ae,be] \
#                *(((KK[ag,bg,ae,be]+1j*(omegas_g[ag,bg]-omegas_e[ae,be]))/
#                (KK[ag,bg,ae,be]**2 + (omegas_g[ag,bg]-omegas_e[ae,be])**2))
#                *eom_mg[ag,bg]) \
#                *(1.0-eom_pg[ag,bg]*eom_me[ae,be]*(eK**2))
#                
#                        if (abs(omegas_g[ag,bg] - omegas_e[ae,be])>1.0*cm2int):
#                            U_ETICS[ag,bg,ae,be] = 0.0
#                        #print(">",ag,bg,ae,be,"->",U_ETICS[ag,bg,ae,be])
#
#        for ae in range(self.Nb[1]):
#            for be in range(self.Nb[1]):
#                if ae != be:
#                    U_EX[ae,be] = numpy.exp(-0.5*K*popt)*eom_me[ae,be]
#        for ae in range(self.Nb[1]):
#            U_EX[ae,ae] = numpy.exp(-K*popt)                         
#        
#        # we return a result in exciton basis
#        return U_ETICS, U_EX
        
        

    
    def liouville_pathways_3T(self, ptype="R3g", eUt2=None,
                              dtol=-1.0, ptol=1.0e-3, etol=1.0e-6,
                              verbose=True, lab=None):
        """ Generator of Liouville pathways with energy transfer
        
        
        
        
        Parameters
        ----------
        
        ptype : tuple, list, str
            List of strings or a string representing one or more
            Liouville pathway types that are to be calculated
            
        eUt2 : EvolutionSuperOperator
            Evolution superoperator at time t2 representing the energy 
            transfer in the system 
            
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
        
        pop_tol = ptol
        dip_tol = numpy.sqrt(self.D2_max)*dtol
        evf_tol = etol
        
        # Check if the ptype is a tuple
        if not isinstance(ptype, (tuple,list)):
            ptype_tuple = (ptype,)
        else:
            ptype_tuple = ptype
        lst = []
        
        if verbose:
            print("Pathways", ptype_tuple)
         
        for ptp in ptype_tuple:
        
            if ptp == "R3g":
            
                if verbose:
                    print("Liouville pathway R3g")
                    
                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                
                k = 0
                l = 0
                for i1g in ngs:
                    
                    if verbose: 
                        print("Ground state: ", i1g, "of", len(ngs))
                        
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
                
                if verbose:
                    print("Liouville pathway R2g")
                
                k = 0
                l = 0
                for i1g in ngs:

                    if verbose: 
                        print("Ground state: ", i1g, "of", len(ngs))

                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            #print('\rFirst excitation: %i ' % i2e, end = '\r')
                            if self.D2[i2e,i1g] < dip_tol:
                                break
                            
                            for i3e in nes:
                                
                                if self.D2[i3e,i1g] < dip_tol:
                                    break
                                
                                for i3d in nes:
                                    for i2d in nes:
                                        
                                        evf = eUt2.data[i3d, i2d, i3e, i2e]
                                        if abs(evf) < evf_tol:
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
                                                lp = \
                                                diag.liouville_pathway("R", 
                                                                       i1g,
                                                            aggregate=self,
                                                            order=3,pname=ptp,
                                                            popt_band=1,
                                                            relax_order=1)
                                                
                                                #      |g_i1> <g_i1|
                                                lp.add_transition((i2e,i1g),-1)
                                                #      |g_i1> <e_i2|
                                                lp.add_transition((i3e,i1g),+1)
                                                #      |e_i3> <e_i2|
                                                lp.add_transfer(((i3d, i2d)),
                                                                 (i3e, i2e))
                                                lp.set_evolution_factor(evf)
                                                lp.add_transition((i4g,i2d),-1)
                                                #      |e_i3> <g_i4|
                                                lp.add_transition((i4g,i3d),+1)
                                                #      |g_i4> <g_i4|
        
                                            except:
                                                
                                                raise Exception()
                                                break
                                            
                                            lp.build()
                                            lst.append(lp)
                                            k += 1
                
             
            if ptp == "R1g":
                
                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                                
                if verbose:
                    print("Liouville pathway R1g")
                
                k = 0
                l = 0
                for i1g in ngs:

                    if verbose: 
                        print("Ground state: ", i1g, "of", len(ngs))
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break

                            for i3e in nes:

                                if self.D2[i3e,i1g] < dip_tol:
                                    break

                                for i2d in nes:
                                    for i3d in nes:
                                        
                                        evf = eUt2.data[i2d, i3d, i2e, i3e]
                                        if abs(evf) < evf_tol:
                                            break

                                        for i4g in ngs:

                                            if ((self.D2[i4g,i3d] < dip_tol)
                                             or (self.D2[i4g,i2d] < dip_tol)):
                                                break
        
                                            l += 1
        
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
                                                lp = \
                                                diag.liouville_pathway("NR",
                                                                       i1g,
                                                            aggregate=self,
                                                            order=3,pname=ptp,
                                                            popt_band=1,
                                                            relax_order=1)
                                                #      |g_i1> <g_i1|                                                           
                                                lp.add_transition((i2e,i1g),+1)
                                                #      |e_i2> <g_i1|        
                                                lp.add_transition((i3e,i1g),-1)
                                                #      |e_i2> <e_i3|
                                                lp.add_transfer(((i3d, i2d)),
                                                                 (i2e, i3e))
                                                lp.set_evolution_factor(evf)
                                                #      |d_i2> <d_i3|                                                                                
                                                lp.add_transition((i4g,i3d),-1)
                                                #      |d_i2> <g_i4|
                                                lp.add_transition((i4g,i2d),+1)
                                                #      |g_i4> <g_i4|
        
                                            except:
                                                
                                                break
                                            
                                            lp.build()
                                            lst.append(lp)
                                            k += 1
            
            if ptp == "R4g":
                
                ngs = self.get_electronic_groundstate()
                nes = self.get_excitonic_band(band=1)
                                
                if verbose:
                    print("Liouville pathway R1g")
                
                k = 0
                l = 0
                for i1g in ngs:

                    if verbose: 
                        print("Ground state: ", i1g, "of", len(ngs))
                    
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
                
                if verbose:
                    print("Liouville pathway R1f*")
                
                
                k = 0
                l = 0
                for i1g in ngs:

                    if verbose: 
                        print("Ground state: ", i1g, "of", len(ngs))
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break
                            
                            for i3e in nes:
                                
                                if self.D2[i3e,i1g] < dip_tol:
                                    break
                                
                                if verbose:
                                    print("Excited state: ", i2e, i3e, "of",
                                          nes[len(nes)-1])
                                for i3d in nes:
                                    for i2d in nes:
                                        
                                        evf = eUt2.data[i3d, i2d, i3e, i2e]
                                        if abs(evf) < evf_tol:
                                            break
                                
                                        for i4f in nfs:

                                            if ((self.D2[i4f,i3d] < dip_tol)
                                             or (self.D2[i2d,i4f] < dip_tol)):
                                                break
                                    
                                            l += 1
                                    

                                    #      Diagram R4g
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

                                                lp = \
                                                diag.liouville_pathway("R",i1g,
                                                           aggregate=self,
                                                           order=3,pname=ptp,
                                                           popt_band=1,
                                                           relax_order=1)
                                                #      |g_i1> <g_i1|                                                           
                                                lp.add_transition((i2e,i1g),-1)
                                                #      |g_i1> <e_i2|
                                                lp.add_transition((i3e,i1g),+1)
                                                #      |e_i3> <e_i2|
                                                lp.add_transfer(((i3d, i2d)),
                                                                 (i3e, i2e))
                                                lp.set_evolution_factor(evf)
                                                #      |d_i3> <d_i2|                                        
                                                lp.add_transition((i4f,i3d),+1)
                                                #      |f_i4> <d_i2|
                                                lp.add_transition((i2d,i4f),+1)
                                                #      |d_i2> <d_i2|

                                            except:
                                                raise Exception("Construction"+
                                                "relaxation pathway failed")
                                                #break
                                    
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
                
                if verbose:
                    print("Liouville pathway R2f*")

                k = 0
                l = 0
                for i1g in ngs:

                    if verbose: 
                        print("Ground state: ", i1g, "of", len(ngs))
                    
                    # Only thermally allowed starting states are considered
                    if self.rho0[i1g,i1g] > pop_tol:
                
                        for i2e in nes:
                            
                            if self.D2[i2e,i1g] < dip_tol:
                                break                                                        
                            
                            for i3e in nes:
                                
                                if self.D2[i3e,i1g] < dip_tol:
                                    break
                                
                                if verbose:
                                    print("Excited state: ", i2e, i3e, "of",
                                          nes[len(nes)-1])
                                for i2d in nes:
                                    for i3d in nes:
                                        
                                        evf = eUt2.data[i2d, i3d, i2e, i3e]
                                        if abs(evf) < evf_tol:
                                            break

                                        for i4f in nfs:
        
                                            if ((self.D2[i4f,i2d] < dip_tol)
                                             or (self.D2[i3d,i4f] < dip_tol)):
                                                break
                                            
                                            l += 1
                                            
        
                                            #      Diagram R4g
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
        
                                                lp = \
                                                diag.liouville_pathway("NR",
                                                                       i1g,
                                                            aggregate=self,
                                                            order=3,pname=ptp,
                                                            popt_band=1,
                                                            relax_order=1)
                                                #      |g_i1> <g_i1|                                                           
                                                lp.add_transition((i2e,i1g),+1)
                                                #      |e_i2> <g_i1|
                                                lp.add_transition((i3e,i1g),-1)                                        
                                                #      |e_i2> <e_i3|
                                                lp.add_transfer(((i2d, i3d)),
                                                                 (i2e, i3e))
                                                lp.set_evolution_factor(evf)
                                                #      |d_i2> <d_i3|                                        
                                                lp.add_transition((i4f,i2d),+1)
                                                #      |f_i4> <d_i3|
                                                lp.add_transition((i3d,i4f),+1)
                                                #      |d_i3> <d_i3|
        
                                            except:
                                                
                                                break
                                            
                                            lp.build()
                                            lst.append(lp)
                                            k += 1
                    
        
        if lab is not None:
            for l in lst:
                l.orientational_averaging(lab)
         
        return lst     

        
        
        
