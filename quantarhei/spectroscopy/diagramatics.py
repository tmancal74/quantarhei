# -*- coding: utf-8 -*-
"""

    Liouville pathways and their analysis



"""


import numpy

from ..utils.types import Integer
from ..core.units import cm2int
from ..core.managers import UnitsManaged

class liouville_pathway(UnitsManaged):

    order = Integer("order")
    nint = Integer("nint")
    
    def __init__(self, ptype, sinit, aggregate=False, 
                 order=3, pname="",relax_order=0, popt_band=0):
        """Liouville pathway through a molecular aggregate

    
        Liouville pathway is represented as a sequence of transitions
        on the right and left hand sides of a double sided Feynman diagram
        

        
        Parameters
        ----------
        
        ptype : str {"R", "NR", "DC"}
            type of the pathway 
            
        sinit
            starting state (all pathways start from diagonal density matrix)        
        
        aggregate
            specifies the aggregate object in which the Liouville pathway is taken
            
        order
            order of the pathway. Default value is 3. This is also currently
            the only value allowed.
        
        """
        
        if not aggregate:
            raise Exception("aggregate has to be specified")
            
        # initial state of the pathway
        self.sinit = numpy.zeros(2,dtype=numpy.int16)
        self.sinit[0] = sinit
        self.sinit[1] = sinit
        
        # order of the pathway
        self.order = order 
        
        # order of the pathways in terms of relaxation events
        self.relax_order = relax_order
        
        # list events associated with relaxations and light interations
        self.event = [None]*(1+order+relax_order)
        
        # type of the pathway (rephasing, non-rephasing, double-coherence)
        self.pathway_type = ptype
        
        # pathway name (one of the standard names)
        self.pathway_name = pname

        # aggregate for which the pathway is made
        self.aggregate = aggregate 

        
        # current state of the pathway (during building)
        self.current = numpy.zeros(2,dtype=numpy.int16)
        self.current[0] = sinit
        self.current[0] = sinit

        # index counting the number of interactions with light
        self.nint = 0
        
        # index counting relaxations 
        self.nrel = 0
        
        # index of events in the diagram (interactions and relaxations)
        self.ne = 0

        # light indiced transitions
        self.transitions = numpy.zeros((order+1,2), dtype=numpy.int)
        
        # relaxation induced transitions
        self.relaxations = [None]*relax_order
        
        # sides from which the transitions occurred
        self.sides = numpy.zeros(order+1,dtype=numpy.int16)
        
        self.states = numpy.zeros((1+order+relax_order,2),dtype=numpy.int)
        
        # transition dipole moments associated with the transition
        self.dmoments = numpy.zeros((order+1,3))
        
        # FIXME: we probably do not need the energy property at all
        # energy of the transition
        self.energy = numpy.zeros(order+1)
        
        # frequency of the transition
        self.frequency = numpy.zeros(1+order+relax_order)
        
        # orientational prefactor (negative means that it is not initialized)
        self.pref = -1.0
        
        # band through which the pathway travels at population time
        self.popt_band = popt_band
        
        # quantity used to evaluate orientational averaging
        self.F4n = numpy.zeros(3)
        
        # was the pathway already built?
        self.built = False
        
    def __str__(self):
        """String representation of the Liouville pathway 
        
        
        
        """
        
        k = 3
        # output string - starts empty
        out = ""
        # number of events (light interactions + relaxations)
        # equals to the number of horizotal lines in the diagram
        noe = 1+self.order+self.relax_order
        
        ii = 0 # index of light interactions
        rr = 0 # index of relaxations
        # iterate over horizoltal lines
        for ee in range(noe):
            
            # what is the event type
            if self.event[ee] == "I":
                
                # side on which the interaction occurs
                sd = self.sides[ii]
                # first we start with both sides from self.sinit
                if ii == 0:
                    lx = self.sinit[0]
                    rx = self.sinit[0]
                    out = ("    |%i    %i|  \n") % (lx, rx)
                # now it depends if the interaction is from left or right
                if sd == 1: # interaction from left
                    if ii != 3: # if this is not the last interaction print
                                # also the frequency
                        outr =  ("    |%i    %i|      %r\n" % 
                                                     (self.transitions[ii,0],
                                        rx,numpy.round(self.frequency[ee])))
                     
                    else: # last interaction
                        outr =  ("    |%i    %i|  \n" %
                                                     (self.transitions[ii,0],
                                                     rx)) 
                    # output the arrow
                    outr += ("--->|------|  \n")
                    # and an empty space
                    outr += ("    |      |  \n") 
                    # all this is appended at the end
                    out = outr+out
                    lx = self.transitions[ii,0]
                else:  # interaction from the right
                    if ii != 3: # if this is not the last interaction
                                # print also the frequency
                        outl =  ("    |%i    %i|      %r\n" % 
                                              (lx, self.transitions[ii,0],
                                            numpy.round(self.frequency[ee])))
                    # actually, iteraction from the right as last does not
                    # occur by convention
                                            
                    # output the arrow
                    outl += ("    |------|<---  \n") % self.frequency[ee]
                    # and an empty space
                    outl += ("    |      |  \n")
                    # append everything at the end
                    out = outl+out
                    rx = self.transitions[ii,0]
                    
                ii += 1
                
            elif self.event[ee] == "R":
                lf = self.relaxations[rr][0][0]
                rf = self.relaxations[rr][0][1]
                ene = self.convert_energy_2_current_u(self.frequency[ee])
                outR  = "    |%i    %i|      %r\n" % (lf,rf,
                                            numpy.round(ene))
                outR += "   >|******|< \n"
                outR += "    |      | \n"
                out = outR+out
                
                lx = lf
                rx = rf
                rr += 1               
                
            else:
                raise Exception("Unknown event type")
                    
            k -= 1


        outd = ("\n\nLiouville Pathway %s (type = %s) \n" %
                 (self.pathway_name, self.pathway_type))
        outd += ("Orientational prefactor: %r \n\n" % self.pref)
        
        out = outd+out
        
        return out 
    
    def add_transition(self, transition, side):
        """ Adds a transition to the Liouville pathway. 

        Parameters
        ----------

        transition
            Tuple such as (8,2) specifying the states in the Hamiltonian
            between which the transition proceeds. The order of the initial
            and final states is from right to left. So in this case the 
            starting state is 2 and the final state is 8.
            
        side
            Takes values of +1 and -1 one denoting the side of the diagram
            where the interaction occurs. +1 stands for the left,
            -1 for the right. The sign of the pathway is given by the number
            of interaction on the right.
        
        """

        # final state of the transition
        nf = transition[0]
        # initial state of the transition
        ni = transition[1]
        
        # which side is interacted on? Left or right?
        sd = (abs(side)-side)//2
        text = ["right","left"]        
        
        # check if the transition start is consistent with current state
        # in the diagram
        if (self.current[sd] != ni):
            raise Exception(("Transition on the {} hand"+
            "side of the diagram has to start from state {}").format(text[sd],
                                                        self.current[sd]))
            
        # save the transition associated with this interaction
        self.transitions[self.nint,:] = transition
        # save the side on which the transition occurs
        self.sides[self.nint] = side 
        # save the current 
        self.current[sd] = nf
        
        self.states[self.ne,0] = self.current[0]
        self.states[self.ne,1] = self.current[1]
        
        """
          Some values can be stored locally - they are also
          in the aggregate object
        """
        #self.energy[self.nint] = \
        #       (self.aggregate.HH[nf,nf]/cm2int 
        #       -self.aggregate.HH[ni,ni]/cm2int)
        self.energy[self.nint] = \
               (self.aggregate.HH[nf,nf] 
               -self.aggregate.HH[ni,ni])
               
        self.dmoments[self.nint,:] = \
               self.aggregate.DD[nf,ni,:]
        
        # check if the number of interactions is not larger than the order
        # of the pathway
        if self.nint < self.order:
            nl = self.current[0]
            np = self.current[1]
            el = self.aggregate.HH[nl,nl]
            ep = self.aggregate.HH[np,np]
            
            self.frequency[self.ne] = el - ep
        elif self.nint > self.order:
            etext = ("Number of interactions larger than the order"
            +" of the pathway.")
            raise Exception(etext)

        self.nint += 1
        
        self.event[self.ne] = "I"
        self.ne += 1
        
    def add_transfer(self, fin, sta):
        """Adds an energy transfer compoment to the pathway
        
        
        
        """
        # final state of the transition
        nfl = fin[0]
        nfr = fin[1] 
        # initial state of the transition
        nil = sta[0]
        nir = sta[1]
                
        # check if the transition start is consistent with current state
        # in the diagram
        if ((self.current[0] != nil) or (self.current[1] != nir)):
            raise Exception("Relaxation does not start" 
            + "from the current state")
        
        
        if self.nrel < self.relax_order:            
            self.relaxations[self.nrel] = (fin,sta)
        else:
            raise Exception("Expected number of relaxation events exceeded")
             
        self.current[0] = nfl
        self.current[1] = nfr
        self.states[self.ne,0] = self.current[0]
        self.states[self.ne,1] = self.current[1]
        el = self.aggregate.HH[nfl,nfl]
        ep = self.aggregate.HH[nfr,nfr]
        self.frequency[self.ne] = el - ep       

        self.nrel += 1
        
        self.event[self.ne] = "R"
        self.ne += 1
        
        
        
    def build(self):
        """Building the Liouville pathway internals
         
        Here we calculate all important values which depend on all 
        transitions in the pathway such as:
        
        1) F4n ... vector for spatial averaging
        
        """
        
        d = self.dmoments
        
        self.F4n[0] = numpy.dot(d[3,:],d[2,:])*numpy.dot(d[1,:],d[0,:])
        self.F4n[1] = numpy.dot(d[3,:],d[1,:])*numpy.dot(d[2,:],d[0,:])
        self.F4n[2] = numpy.dot(d[3,:],d[0,:])*numpy.dot(d[2,:],d[1,:]) 
        
        self.sign = numpy.prod(self.sides)
        
        self.built = True
        
        
    def get_current_state(self):
        """
        
        """
        return self.current
        
    def get_transition(self,n):
        """ Returns info on the transition occuring on the n-th interaction
        
        """
        return self.side[n], self.transitions[n]
    
    def orientational_averaging(self,lab):
        """ Orientational averaging
        
        Orientational averaging and temperature dependence of the
        initial density matrix 
        
        """
        # weight in the initial state of the first transition
        n0 = self.transitions[0,1]
        self.pref = self.sign*(numpy.dot(lab.F4eM4,self.F4n)
        *numpy.real(self.aggregate.rho0[n0,n0])) 
        
    
    def get_transition_energy(self,n):
        """ Transition energy of a given interaction with light
        
        """
        return self.energy[n]
        
    def get_interval_frequency(self,n):
        """Returns frequency corresponding to a given response interval
        
        """
        return self.frequency[n]
    
    def get_dmoment(self, n):
        """ Transition dipole moment of a given interaction with light
        
        """
        return self.dmoments[n]
        
    def get_prefactor(self):
        """ Dipole and temperature dependent prefactor 
        
        """
        return self.pref
    
    

class LiouvillePathwayAnalyzer(UnitsManaged):
    """Class providing methods for Liouville pathway analysis
    
    
    Parameters
    ----------
    
    pathways : list
        List of pathways to analyze
    
    """

    def __init__(self, pathways=None):
        self.pathways = pathways


    def set_pathways(self, pathways):
        """Sets pathways to be analyzed
        
        """
        self.pathways = pathways

    
    def max_pref(self, pathways):
        """Return the maximum of pathway prefactors
        
        
        Parameters
        ----------
        
        pathways : list
            List of Liouville pathways
            
            
        Returns
        -------
        
        pmax : float
            Maximum prefactor of the pathways
            
        rec : int
            position of the pathway with the maximum prefactor
        
        """
        
        if pathways is None:
            pthways = self.pathways
        else:
            pthways = pathways
        
        pmax = 0.0
        k = 0
        rec = -1
        for pway in pthways:
            if pway.pref > pmax:
                rec = k
                pmax = pway.pref
            k += 1
            
        return (pmax, rec)


    def select_pref_GT(self, val, pathways=None, replace=True,
                       verbatime=False):
        """Select all pathways with prefactors greater than a value
        
        """
        
        if pathways is None:
            pthways = self.pathways
        else:
            pthways = pathways
            
        selected = []
        for pway in pthways:
            if pway.pref > val:
                selected.append(pway)
        
        if verbatime:
            print("Selected", len(selected), "pathways")
            
        # if the pathways were not specified from argument then return selected
        # to self
        if (pathways is None) and replace:
            self.pathways = selected
            
        return selected

    
    def select_frequency_window(self, window, pathways=None, replace=True, 
                                verbatime=False):
        """Selects pathways with omega_1 and omega_3 in a certain range
        
        """
        if pathways is None:
            pthways = self.pathways
        else:
            pthways = pathways

        om1_low = self.convert_energy_2_internal_u(window[0])
        om1_upp = self.convert_energy_2_internal_u(window[1])
        om3_low = self.convert_energy_2_internal_u(window[2])
        om3_upp = self.convert_energy_2_internal_u(window[3])
        
        
        selected = []
        
        for pway in pthways:
            ne = len(pway.frequency)
            #om1 = numpy.abs(pway.frequency[0])
            om1 = numpy.abs(pway.get_interval_frequency(0))
            #om3 = pway.frequency[ne-2]
            om3 = numpy.abs(pway.get_interval_frequency(ne-2))
            
            if (((om1 >= om1_low) and (om1 <= om1_upp)) and 
                ((om3 >= om3_low) and (om3 <= om3_upp))):
                selected.append(pway) 
                
        if verbatime:
            print("Selected", len(selected), "pathways")

        if (pathways is None) and replace:
            self.pathways = selected
            
        return selected


    def select_omega2(self, interval, pathways=None, replace=True,
                      verbatime=False):
        """Selects pathways with omega_2 in a certain interval
        
        """    
        if pathways is None:
            pthways = self.pathways
        else:
            pthways = pathways

        om2_low = self.convert_energy_2_internal_u(interval[0])
        om2_upp = self.convert_energy_2_internal_u(interval[1])  
        
        selected = []
        
        for pway in pthways:
            ne = len(pway.frequency)
            #om2 = pway.frequency[ne-3]
            om2 = pway.get_interval_frequency(ne-3)
            
            if (om2 >= om2_low) and (om2 <= om2_upp):
                selected.append(pway)
                
        if verbatime:
            print("Selected", len(selected), "pathways")

        if (pathways is None) and replace:
            self.pathways = selected
        
        return selected
    
    
    def order_by_pref(self, pthways):
        """Orders the list of pathways by pathway prefactors
        
        """
        lst = sorted(pthways, key=lambda pway: pway.pref, reverse=True)
        return lst
    
        