# -*- coding: utf-8 -*-
import numpy

from ..utils.types import Float
from ..utils.types import Integer

from ..core.managers import UnitsManaged, energy_units


class ElectronicState(UnitsManaged):
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
    
    def __init__(self, aggregate, elsignature, index=None):
        
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
        
        # if index is specified, use it, otherwise try to search for it
        if index is not None:
            self.index = index
        else:
            self.index = self.aggregate.elsigs.index(self.elsignature) 
        
        #
        # Implementation of vibrational substate generation approximations
        #
        #
        self.vibgen_approximation = None
        
        
    def get_signature(self):
        """Returns electronic signature of this state
        
        """
        return self.elsignature


    def energy(self, vsig=None):
        """ Returns energy of the state (electronic + vibrational)
        
        """
        en = 0.0

        if vsig is not None:
            if not (len(vsig) == self.vsiglength):
                raise Exception()  
            k = 0
            for nn in self.vibmodes:
                en += vsig[k]*self.convert_energy_2_current_u(nn.omega)
                k += 1

        k = 0
        for nn in self.elsignature:
            en += \
            self.convert_energy_2_current_u(
                    self.aggregate.monomers[k].elenergies[nn])
            k += 1
            
        return en
    
        
    def vibenergy(self, vsig=None):
        """ Returns vibrational energy of a state
        
        """
        
        en = 0.0
        
        if vsig is not None:
            
            if not (len(vsig) == self.vsiglength):
                raise Exception()  
           
            k = 0
            for nn in self.vibmodes:
                en += vsig[k]*self.convert_energy_2_current_u(nn.omega)
                k += 1
            
        return en
        

    def _spa_ndindex(self, tup, ecut=None):
        """Generates vibrational signatures in SPA 
        
        """
        if ecut is not None:
            vec = self.convert_energy_2_internal_u(ecut)
            
        ret = []
        N = len(tup)
        zer = numpy.zeros(N,dtype=numpy.int)
        ret.append(tuple(zer))
        # single particles
        
        if ecut is not None:
            with energy_units("int"):    
                for i in range(N):
                    one = zer.copy()
                    one[i] = 1
                    sig = tuple(one)   
                    en = self.vibenergy(vsig=sig) 
                    if en <= vec:
                        ret.append(sig)

        else:
            for i in range(N):
                one = zer.copy()
                one[i] = 1
                sig = tuple(one)
                ret.append(sig)
                
        return tuple(ret)
        
        
    def _sppma_ndindex(self, tup, ecut=None):
        """Generates vibrational signatures in SP per mode Aproximation (SPPMA) 
        
        """
        if ecut is not None:
            vec = self.convert_energy_2_internal_u(ecut)

        two = numpy.zeros(len(tup), dtype=numpy.int)
        two[:] = 2
        for i in range(len(two)):
            if two[i] > tup[i] + 1:
                two[i] = tup[i] + 1
                
        shp = tuple(two)
        
        if ecut is not None:
            #with energy_units("int"):    
            for sig in numpy.ndindex(shp):
                en = self.convert_energy_2_internal_u(self.vibenergy(vsig=sig))
                if en <= vec:
                    yield sig
                        
        else:
            return numpy.ndindex(shp)


    def _tpa_ndindex(self, tup, ecut=None):
        """Generates vibrational signature in Two Particle Approximation (TPA)
        
        """
        return self._npa_ndindex(tup, N=2, ecut=ecut)
    

    def _tppma_ndindex(self, tup, ecut=None):
        """Generates vibrational signatures in TP per mode Aproximation (TPPMA) 
        
        """
        return self._nppma_ndindex(tup, N=2, ecut=ecut)


    def _npa_ndindex(self, tup, N=None, ecut=None):
        """Generates vibrational signature in N Particle Approximation (NPA)
        
        """
        
        if N is None:
            raise Exception("Number of states to be used must be defined")
            
        two = numpy.zeros(len(tup), dtype=numpy.int)
        two[:] = N + 1
        shp = tuple(two)
        
        if ecut is not None:
            vec = self.convert_energy_2_internal_u(ecut)
            for tpl in numpy.ndindex(shp):
                en = self.convert_energy_2_internal_u(self.vibenergy(vsig=tpl))
                if (numpy.sum(tpl) <= N) and en <= vec:
                   yield tpl
            
        else:
            
            for tpl in numpy.ndindex(shp):
                if numpy.sum(tpl) <= N:
                   yield tpl


    def _nppma_ndindex(self, tup, N=None, ecut=None):
        """Generates vibrational signatures in NP per mode Aproximation (NPPMA) 
        
        """
        if N is None:
            raise Exception("Number of states to be used must be defined")

        if ecut is not None:
            vec = self.convert_energy_2_internal_u(ecut)

        two = numpy.zeros(len(tup), dtype=numpy.int)
        two[:] = N + 1

        for i in range(len(two)):
            if two[i] > tup[i] + 1:
                two[i] = tup[i] + 1
                
        shp = tuple(two)
        
        if ecut is not None:
            #with energy_units("int"):    
            for sig in numpy.ndindex(shp):
                en = self.convert_energy_2_internal_u(self.vibenergy(vsig=sig))
                if en <= vec:
                    yield sig
                        
        else:
            for sig in numpy.ndindex(shp):
                yield sig        
             
       
    def vsignatures(self, approx=None, N=None, vibenergy_cutoff=None):
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
            'TPA' : two particles approximation. Generates maximum 2 phonon
                    particles per state over all modes
            'TPPMA' : two particles per mode approximation
            'NPA' : N-particles approximation. N is specified
                    as an additional argument
            'NPPMA' : N-particles per mode approximation. N is specified
                    as an additional argument
                        
            The 'per mode' variants respect the maximum number of modes defined
            in the submods. The 'SPA', 'TPA' and 'NPA' approximations do
            not respect the maximum number of states in a mode and always 
            generate the number of states corresponding to the approximation.
            
            
        """
        vibmax = []
        for sm in self.vibmodes:
            vibmax.append(sm.nmax)
            
        if approx is None:    
            return numpy.ndindex(tuple(vibmax))
        elif approx == 'ZPA':
            return numpy.ndindex(tuple([1]*len(vibmax)))
        elif approx == 'SPA':
            return self._spa_ndindex(tuple(vibmax), ecut=vibenergy_cutoff)
        elif approx == 'SPPMA':
            return self._sppma_ndindex(tuple(vibmax), ecut=vibenergy_cutoff)
        elif approx == 'TPA':
            return self._tpa_ndindex(tuple(vibmax), ecut=vibenergy_cutoff)
        elif approx == 'TPPMA':
            return self._tppma_ndindex(tuple(vibmax), ecut=vibenergy_cutoff)
        elif approx == 'NPA':
            return self._npa_ndindex(tuple(vibmax), N=N, ecut=vibenergy_cutoff)
        elif approx == 'NPPMA':
            return self._nppma_ndindex(tuple(vibmax), N=N,
                                       ecut=vibenergy_cutoff)
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
        out  = "\nquantarhei.ElectronicStateobject"
        out += "\n=================================="
        out += "\nenergy = %f" % self.energy()
        out += "\nnumber of substates = %i " % self.number_of_states()
        out += ("\nmember of an aggregate of %i molecules" 
                % self.aggregate.nmono)
        return out
    
        
#class vibronic_state(UnitsManaged):
class VibronicState(UnitsManaged):
    """Represents a vibronic state of an aggregate
    
    Vibronic state of an aggregate is composed of an electronic state
    and its vibronic signature.
    
    """
    
    def __init__(self, elstate, vsig):
        self.elstate = elstate
        self.vsig = vsig
        
        try:
            agg = self.elstate.aggregate
            self.index = agg.vibsigs.index((elstate.elsignature, vsig))
        except:
            self.index = None


    def get_ElectronicState(self):
        """Returns corresponding electronic state
        
        """
        return self.elstate

  
    def get_vibsignature(self):
        """Returns corresponding vibrational signature
        
        """
        return self.vsig


    def energy(self):
        """Returns the energy of the state
        
        """
        return self.elstate.energy(self.vsig)


    def vibenergy(self):
        """ Returns vibrational energy of a state
        
        """
        
        en = 0.0
        
        if self.vsig is not None:
            
            if not (len(self.vsig) == self.elstate.vsiglength):
                raise Exception()  
           
            k = 0
            for nn in self.elstate.vibmodes:
                en += self.vsig[k]*self.convert_energy_2_current_u(nn.omega)
                k += 1
            
        return en        


    def signature(self):
        """Returns a signature of the state
        
        """
        return (self.elstate.elsignature, self.vsig)
    
    
    def __str__(self):
        out  = "\nquantarhei.VibronicState object"
        out += "\n=================================="
        out += "\nenergy = %f" % self.energy()
        out += ("\nmember of an aggregate of %i molecules" 
                % self.elstate.aggregate.nmono)
        return out        
    
    
class Coherence:
    """Class representing quantum mechanical coherence
    
    This class collects methods for analysis of coherence character
    
    
    
    Parameters
    ----------
    
    state1 : {ElectronicState, VibronicState}
        State which participates in the coherence
        
    state2 : {ElectronicState, VibronicState}
        State which participates in the coherence        
    
    
    """
    
    def __init__(self, state1, state2):
        
        self.state1 = state1
        self.state2 = state2
        
        try:
            system1 = state1.system
            system2 = state2.system
        except:
            raise Exception("At least one of the states"+
                            " does not know its system")
            
        if system1 != system2:
            raise Exception("States do not come from the same system")
            
        self.system = system1


    def get_coherence_character(self):
        """Returns a dictionary describing the character of a given coherence
        
        """
        
        return dict(electronic=0.0, vibrational=0.0, mixed=0.0)