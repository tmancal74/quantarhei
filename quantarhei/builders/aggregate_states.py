# -*- coding: utf-8 -*-
import numpy

from ..utils.types import Float
from ..utils.types import Integer


class aggregate_state():
    """ State of an aggregate
    
    
    
    """
    
    energy = Float('energy')
    el_energy = Float('el_energy')
    vib_energy = Float('vib_energy')
    
    def __init__(self, aggregate, state_tuple):
        
        desc = state_tuple[0]
        
        self.el_descriptor = desc[0]
        self.vib_descriptor = desc[1:]
        self.vibrations = state_tuple[1]
        
        # This is temporary - numbers have to be computed
        self.energy = 0.0
        self.el_energy = 0.0
        self.vib_energy = 0.0
        
        self.aggregate = aggregate
        
    def aslist(self):
        return [self.el_descriptor,self.vib_descriptor]


class electronic_state:
    """ Represents electronic state of an aggregate 
   
    Parameters
    ----------
    
    aggregate :
        Parent aggregate of the electronic state
        
    elsignature : tuple of int
        Electronic signature of the state
        
    index : int
        Index of the state in the aggregate 
        
    Properties
    ----------
    
    nmono : int
        Number of monomers in the aggregate
        
    elsignature : tuple of int
        Signature of the electronic state
        
    vibmodes :
        Vibrational (sub)modes of the state
        
    vsiglength : int
        Length of the vibrational signatures (= total number of submodes)
    
    """
    
    nmono = Integer('nmono')
    
    def __init__(self, aggregate, elsignature, index=0):
        
        Nsig = len(elsignature)
        nmono = len(aggregate.monomers)
        if Nsig == nmono:
            self.nmono = nmono
        else:
            raise Exception("Incompatible length of elsignature")
        
        self.elsignature = elsignature
        self.aggregate = aggregate
        
        self.band = 0
        for k in self.elsignature:
            self.band += k 
        
        elst = elsignature
        vb_ls = [] #tuple()
        n = 0
        for mn in aggregate.monomers:
            for a in range(mn.nmod):
                vb_ls.append(mn.get_Mode(a).get_SubMode(elst[n]))
            n += 1
            
        self.vibmodes = vb_ls
        self.vsiglength = len(vb_ls)
        
        self.index = index
        
        #
        # Implementation of vibrational substate generation approximations
        #
        self.vibgen_approximation = None
        
        
    def energy(self, vsig=""):
        """ Returns energy of the state (electronic + vibrational)
        
        """
        en = 0.0
        k = 0
        if not vsig=="":
            if not (len(vsig) == self.vsiglength):
                raise Exception()  
            k = 0
            for nn in self.vibmodes:
                en += vsig[k]*nn.omega
                k += 1
            k = 0
        for nn in self.elsignature:
            en += self.aggregate.monomers[k].elenergies[nn]
            k += 1
            
        return en
    
        
    def _spa_ndindex(self, tup):
        """Generates vibrational signatures in SPA 
        
        """
        
        ret = []
        N = len(tup)
        zer = numpy.zeros(N,dtype=numpy.int)
        ret.append(tuple(zer))
        # single particles
        for i in range(N):
            one = zer.copy()
            one[i] = 1
            ret.append(tuple(one))
        return tuple(ret)
        
        
    def _sppma_ndindex(self, tup):
        """Generates vibrational signatures in SPA 
        
        """
        pass

        
    def vsignatures(self, approx=None):
        """ Generator of the vibrational signatures
        
        Parameters
        ----------
        
        approx : str
            Type of approximation in generation of vibrational
            signatures. Values are 
            
            None  : no approximation, all states generated
            'SPA' : single particle approximation. Generates maximum 1 phonon
                    particle per state over all modes
            'SPPMA' : single particle per mode approximation
            
        """
        vibmax = []
        for sm in self.vibmodes:
            vibmax.append(sm.nmax)
            
        if approx is None:    
            return numpy.ndindex(tuple(vibmax))
        elif approx == 'SPA':
            return self._spa_ndindex(tuple(vibmax))
        elif approx == 'SPPMA':
            return self._sppma_ndindex(tuple(vibmax))
        else:
            raise Exception("Unknown vibrational "+
                            "state generation approximation")
            
    
    def number_of_states(self):
        """Number of vibrational sub-states in the electronic state
        
        
        """
        n = 1
        for a in self.vibmodes:
            n *= a.nmax
        return n
        
    def __str__(self):
        out  = "\nquantarhei.electronic_state object"
        out += "\n=================================="
        out += "\nenergy = %f" % self.energy()
        out += "\nnumber of substates = %i " % self.number_of_states()
        out += ("\nmember of an aggregate of %i molecules" 
                % self.aggregate.nmono)
        return out
    
        
class vibronic_state:
    """Represents a vibronic state of an aggregate
    
    Vibronic state of an aggregate is composed of an electronic state
    and its vibronic signature.
    
    """
    
    def __init__(self, elstate, vsig):
        self.elstate = elstate
        self.vsig = vsig
        
    def energy(self):
        """Returns the energy of the state
        
        """
        return self.elstate.energy(self.vsig)
        
        
    def signature(self):
        """Returns a signature of the state
        
        """
        return (self.elstate.elsignature, self.vsig)
    
    
    def __str__(self):
        out  = "\nquantarhei.vibronic_state object"
        out += "\n=================================="
        out += "\nenergy = %f" % self.energy()
        out += ("\nmember of an aggregate of %i molecules" 
                % self.elstate.aggregate.nmono)
        return out        