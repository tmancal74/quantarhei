# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    systembathinteraction module

"""
import numpy

from ...core.saveable import Saveable
from ...qm.corfunctions.cfmatrix import CorrelationFunctionMatrix
from ...qm.corfunctions.functionstorage import FunctionStorage
from ...qm.corfunctions.correlationfunctions import c2g
from ...core.dfunction import DFunction
from ... import REAL

class SystemBathInteraction(Saveable):
    """Describes interaction of an open quantum system with its environment
    
    Stores the system--bath interaction operator in form of a set of operators
    on the Hilbert space of the system and correlation functions of 
    the operator on the bath Hilbert space,
    
    OR
    
    It stores various relaxation and dephasing rates, to represent e.g. the 
    Linblad form.
    
    
    Parameters
    ----------
    
    sys_operators : list
        List of the system part of the system-bath interaction Hamiltonian
        components
    
    bath_correlation_matrix: CorrelationFunctionMatrix 
        Object of the CorrelationFunctionMatrix type holding all correlation
        functions needed for the description of system bath interaction
    
    rates : list/tuple
        List or tuple of rates. The total number of rates has to be 
        the same as the number of operators in the ``sys_operator`` list
        
    drates : array
        An array of dephasing rates. The dimension of the array must
        correspond to the dimension of the treated system.
        
    dtype : str
        Type of the dehasing defined in `drates`. The types are "Lorentzian"
        which is default, and correponds to a dephasing rate equation with
        the term -\gamma\rho_{ab} on the right hand side. The dephasing
        is exponential. The type "Gaussian" results in a Gaussian dephasing
        and corresponds to the term -\gamma t \rho_{ab} on the right hand 
        side of the rate equation.

    osites : list, array
        List or array of site indices on which the oscillators reside. The 
        indices can repeat, indicating several modes on a single site.
        
    orates : array
        Oscillator decay rates. This rates corresponds to the dephasing
        rate of the oscillations in a harmonic oscillator
        
    system : {Molecule, Aggregate}
        Molecule or Aggregate object in which the system--bath interaction 
        is specified


    """

    def __init__(self, sys_operators=None, bath_correlation_matrix=None,
                 rates=None, drates=None, dtype="Lorentzian", osites=None,
                 orates=None, system=None):

        # information about aggregate is needed when dealing with 
        # multiple excitons
        self.aggregate = None
        self.molecule = None
        self.system = None
        self.rates = None
        self.KK = None
        self.CC = None # correlation function matrix
        self.GG = None # lineshape function storage
        self.TimeAxis = None
        self.drates = None
        self.N = 0
        self.osites = None
        self.orates = None
        self.sbitype = "Linear_Coupling"
        
        self._has_gg_storage = False
        
        #
        # version with bath correlation functions
        #
        if ((sys_operators is not None) 
            and (bath_correlation_matrix is not None)):
            
            self.sbitype = "Linear_Coupling"
            
            # Find the length of the list of operators 
            if isinstance(sys_operators,list):
                self.N = len(sys_operators)
            else:
                raise Exception("sys_operators argument has to a list")
            
            # Second argument has to be a CorrelationFunctionMatrix 
            if not isinstance(bath_correlation_matrix, 
                              CorrelationFunctionMatrix):
                raise Exception("ba_correlation_function argument has to a"+
                                " CorrelationFunctionMatrix")
                
            # Check that sys_operators and bath_correlation matrix has 
            # a compatible number of components
            if bath_correlation_matrix.nob != self.N:
                raise Exception("Incompatile number of bath compoments: " +
                    ("Correlation function matrix - %i vs. operators %i" % 
                    (bath_correlation_matrix.nob,self.N)))
                
            self.TimeAxis = bath_correlation_matrix.timeAxis
            self.set_system(system)

            if self.N > 0:
                
                self._set_operators(sys_operators)
                self.CC = bath_correlation_matrix
             
        
            if rates is not None:
                if len(rates) != self.N:
                    raise Exception("Wrong number of rates specified")
                self.rates = rates
                
        #      
        # version with system-bath operators and rates Lindblad form
        #
        elif ((sys_operators is not None) 
            and (rates is not None)):
            
            self.sbitype = "Lindblad_Form"
            
            # Find the length of the list of operators 
            if isinstance(sys_operators, list):
                self.N = len(sys_operators)
            else:
                raise Exception("First argument has to a list")

            self.set_system(system)
            
            if self.N > 0:
                
                self._set_operators(sys_operators)   
                self.CC = None #bath_correlation_matrix
                
                if len(rates) != self.N:
                    raise Exception("Wrong number of rates specified")
                self.rates = rates

        #
        # version with phenomenological pure dephasing
        #
        elif ((sys_operators is None) and (drates is not None)):
            
            self.sbitype = "Pure_Dephasing"
            
            if len(drates.shape) != 2:
                raise Exception("Pure dephasing rates must"
                                +" be defined by a matrix")
                
            if drates.shape[0] != drates.shape[1]:
                raise Exception("Pure dephasing rates must"
                                +" be defined by a square matrix")
                
            self.N = drates.shape[0]
            self.set_system(system)
            self.CC = None
            
        #
        # version with Lindblad form for vibrational modes
        #
        elif ((sys_operators is None) and (orates is not None)):
            
            self.sbitype = "Vibrational_Lindblad_Form"
            
            if len(orates) != len(osites):
                raise Exception("`orates` and `osites` arguments must"+
                                " have the same lengths")
                
            self.set_system(system)
            self.CC = None
            self.orates = orates
            self.osites = osites
            

    def set_system(self, system):
        """Sets the system attribute
        
        """
        from ...builders.aggregates import Aggregate
        from ...builders.molecules import Molecule  
        from ...builders.opensystem import OpenSystem
        
        if system is not None:
            if isinstance(system, Aggregate):
                self.aggregate = system
                self.molecule = None
                self.system = self.aggregate

            elif isinstance(system, Molecule):
                self.aggregate = None
                self.molecule = system
                self.system = self.molecule
                
            elif isinstance(system, OpenSystem):
                self.aggregate = None
                self.molecule = None
                self.system = system
                
            else:
                raise Exception("Unknown system type")
                

    def _set_operators(self, sys_operators):
        """Sets the system part of the interaction
        
        """
        # First of the operators 
        KK = sys_operators[0]
    
        # Get its dimension 
        dim = KK.data.shape[0]
            
        self.KK = numpy.zeros((self.N, dim, dim), dtype=REAL)
        self.KK[0,:,:] = numpy.real(KK.data)       
        
        # Save other operators and check their dimensions 
        for ii in range(1,self.N):
            
            KK = sys_operators[ii]
            
            if True: #isinstance(KK,Operator):
                if dim == KK.data.shape[0]:
                    self.KK[ii,:,:] = numpy.real(KK.data)
                else:
                    raise Exception("Operators in the list are" 
                    + " not of the same dimension")
            else:
                raise Exception("sys_operators tuple (the first argument)"
                + " has to contain cu.oqs.hilbertspace.Operator")


    def get_time_axis(self):
        """Returns the time axis of the storred correlation functions
        
        """
        return self.TimeAxis
        
    
    def get_correlation_function(self, where):
        """Returns the bath correlation function object defined by a pair of sites (tuple)
         
        """

        return self.CC.get_correlation_function(where[0],where[1])


    def get_coft(self, n, m):
        """Returns bath correlation function corresponding to sites n and m

        
        """
        
        if self.sbitype != "Linear_Coupling":
            raise Exception("Correlation functions only defined for "+
                            "linear microscopic system-bath coupling")
        
        #print("Returning coft" )
        
        if self.system is None:
            #print("Returning coft without the system" )
            return self.CC.get_coft(n,m)
            
        else:
            
            #FIXME: Molecule needs this method
            bn = self.system.which_band[n]
            bm = self.system.which_band[m]
            
            # Ground state
            if ((bn == 0) and (bm == 0)):
                
                #
                # This returns zero correlation function, which is consistent
                #
                return self.CC._cofts[0,:]
                
            # First or higher excited state bands
            elif ((bn >= 1) and (bm >= 1)): 
                #print(bn,bm,"::",n-1,m-1)
                
                #
                # First band starts with n=1, but the correlation functions
                # are stored by site index which starts with 0
                #
                return self.CC.get_coft(n-1,m-1)
                
            # other bands return zero correlation function now
            else:
                
                return self.CC._cofts[0,:]
                
                
                
    def get_coft_elsig(self, n_sig, m_sig):
        """Returns bath correlation based on electronic signatures 


        """ 
        if self.sbitype != "Linear_Coupling":
            raise Exception("Correlation functions only defined for "+
                            "linear microscopic system-bath coupling")
                 
        nb = numpy.sum(n_sig)
        mb = numpy.sum(m_sig)
        
        indices = []
        if mb == nb:
            ni = 0
            for na in n_sig:
                mi = 0
                for ma in m_sig:
                    if ((na == 1) and (ma == 1)):
                        indices.append([ni,mi]) 
                    mi += 1
                ni += 1
            
            ret = numpy.zeros((self.TimeAxis.length),dtype=numpy.complex128)
            for ind in indices:
                #print(nb,":",ind[0],ind[1])
                ret += self.get_coft(ind[0],ind[1]) 
                    
                
            return ret            
            
        else:
            return self.CC._cofts[0,:]


    def get_goft_storage(self, config=None):
        """Returns a lineshape function storage based on correlation functions
        
        The function calculates g(t) fuctions based on the correlation 
        functions specified in this object. This call fails if correlation
        functions are not specified.
        
        Parameters
        ----------
        
        config : dict
            Dictionary of the FunctionStorage configuration. If None submitted
            it wil produced g(t) for a standard 3rd order response calculation.
        
        """
        
        if self._has_gg_storage:
            return self.GG
        
        # number of functions
        Nf = self.CC.nof

        # number of sites
        Nb = self.CC.nob

        # we define storage for lineshape functions with prescribed config
        # FIXME: orinally I had Nb here, but that is wrong. Nf also does not work
        gg = FunctionStorage(Nf, 
                             timeaxis=[self.TimeAxis, self.TimeAxis],
                             config=config)
        
        # create functions to update storage
        fcions = {}
        for kk in range(Nb):
            ii = self.CC.get_index_by_where((kk,kk))
            if ii not in fcions:
                # correlation function in discrete representation
                cf = self.get_coft(kk+1, kk+1) 
                # integration into g(t)
                gf = c2g(self.TimeAxis.data, cf)
                # make it into spline function
                df = DFunction(self.TimeAxis, gf)
                gfunc = df.as_spline_function()
                # save for later
                fcions[ii] = gfunc
                
            # set a function for a give site
            gg.set_goft(kk, func=fcions[ii])    
            
        # make checks
        if Nf != gg.get_number_of_functions():
            raise Exception("Number of functions did not"
                            +" conserve in g(t) calculation")
        
        self.GG = gg
        self._has_gg_storage = True
        
        return gg


    def has_temperature(self):
        """Checkst if the Aggregate has a defined temperature
        
        """
        
        if (self.sbitype == "Lindblad_Form" or
            self.sbitype == "Vibrational_Lindblad_Form"):
            return False
        else:
            try:
                T = self.get_temperature()
                if T >= 0.0:
                    return True
            except:
                return False
                

    def get_temperature(self):
        """Returns temperature associated with the bath
        
        """
        return self.CC.get_temperature()
    
    
    def get_reorganization_energy(self, i, j=None):
        """Returns reorganization energy associated with a given site
        
        If one index `i` is specified, the function returns reorganization
        energy associated with i-th site. If two indices are specified,
        the function returns crosscorrelation function (when i not equal j)
        or a correlation function of the site (when i = j)
        
        """
        if (self.sbitype == "Lindblad_Form" or
            self.sbitype == "Vibrational_Lindblad_Form"):
            return None
        else:
            if j is None:
                j = i
            return self.CC.get_reorganization_energy(i,j)


    def get_correlation_time(self, i, j=None):
        
        if (self.sbitype == "Lindblad_Form" or
            self.sbitype == "Vibrational_Lindblad_Form"):
            return None
        else:
            if j is None:
                j = i
            return self.CC.get_correlation_time(i,j)


    def get_sbitype(self):
        """Returns the type of SystemBathInteraction 
        
        Options are `Lindblad_Form`, `Vibrational_Lindblad_Form`
        and `Linear_Coupling`
         
        """
        return self.sbitype
    
    
