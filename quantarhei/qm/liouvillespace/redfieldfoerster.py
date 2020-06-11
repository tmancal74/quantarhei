# -*- coding: utf-8 -*-
import numpy

from .redfieldtensor import RedfieldRelaxationTensor
from .foerstertensor import FoersterRelaxationTensor
from .rates.foersterrates import _reference_implementation as foerster_rates
from ..corfunctions.correlationfunctions import c2g
from ...core.managers import Manager
from ...core.managers import energy_units


class RedfieldFoersterRelaxationTensor(RedfieldRelaxationTensor,
                                       FoersterRelaxationTensor):
    """Combination of Redfield and Foerster relaxation tensors
    
    Paramaters
    ----------
    ham : cu.oqs.hilbertspace.Hamiltonian
        Hamiltonian of the system
        
    sbi : cu.oqs.liouvillespace.SystemBathInteraction
        Object specifying system bath interaction
        
    initialize : bool
        If True, the tensor is imediately calculated
        
    cutoff_time : float
        Time after which the integration kernel of the Redfield tensor
        is assumed to be zero
        
    coupling_cutoff : float
        The smallest value of coupling which is still considered to cause
        delocalization.
        
    
    """
    def __init__(self, ham, sbi, initialize=True,
                 cutoff_time=None, coupling_cutoff=None):
                     
        # non-zero coupling cut-off requires a calculation of both Redfield
        # and Foerster contributions
        if coupling_cutoff is not None:
            m = Manager()
            self.coupling_cutoff = \
                m.convert_energy_2_internal_u(coupling_cutoff)
        else:
            self.coupling_cutoff = 0.0
            
        super().__init__(ham, sbi, initialize=False, 
                             cutoff_time=cutoff_time)
                
        if initialize: 
            with energy_units("int"):
                self._reference_implementation()

                
    def _reference_implementation(self):
        """ Reference all Python implementation
        
        """
        ham = self.Hamiltonian 
        sbi = self.SystemBathInteraction
        
        tt = sbi.TimeAxis.data
        Nt = sbi.TimeAxis.length
        Na = ham.dim

        if ham._has_remainder_coupling:
            JR = ham.JR
        else:
            JR = numpy.zeros((ham.dim, ham.dim), dtype=numpy.float64)
        
        calcRT = True
        calcFT = True

        # is the remainder coupling different from zero?
        if numpy.allclose(JR, numpy.zeros(JR.shape)):
            calcFT = False

        #
        # calculate Redfield tensor for the strong coupling part
        #
        if calcRT:
            
            if self._has_cutoff_time:
                RT = RedfieldRelaxationTensor(ham, sbi,
                                              cutoff_time=self.cutoff_time)
            else:
                RT = RedfieldRelaxationTensor(ham, sbi)
            
            self.data += RT.data


        #
        # Calculate Foerster for the remainder coupling
        #
        if calcFT:

            hD, SS = numpy.linalg.eigh(ham.data) 

                       
            #
            # identify lineshape functions of excitonic states
            #
            gvals = numpy.zeros((Na,Nt),dtype=numpy.complex128)
            Gt = numpy.zeros((Na,Nt),dtype=numpy.complex128)
            for ii in range(1,Na):
                Gt[ii,:] = c2g(sbi.TimeAxis, sbi.CC.get_coft(ii-1,ii-1))
            for aa in range(Na):
                for bb in range(Na):
                    # Here we assume no correlation between sites 
                    gvals[aa,:] += (SS[bb,aa]**4)*Gt[bb,:]  
                    
            #
            # calculate reorganization energies of exciton states
            #
            lamb = numpy.zeros(Na)
            lamb_sites = numpy.zeros(Na)
            for ii in range(1,Na):
                lamb_sites[ii] = sbi.CC.get_reorganization_energy(ii-1,ii-1)
            for aa in range(1,Na):
                for bb in range(1,Na):
                    # Here we assume no correlation between sites 
                    lamb[aa] += (SS[bb,aa]**4)*lamb_sites[bb]


#            # operators
#            nsbi_op = []
#            for aa in range(1,Na):
#                op = Operator(dim=ham.dim, real=True)
#                op.data[aa,aa] = 1.0
#                nsbi_op.append(op)
#            
#            # correlation function matrix
#            cfm = CorrelationFunctionMatrix(ta, Na-1,Na-1)
#            for ii in range(Na-1):
#                params = dict(ftype="Value-defined",reorg=lamb[ii])
#                fc = CorrelationFunction(ta,params,values=cvals[ii,:])
#                cfm.set_correlation_function(fc,[(ii,ii)],ii+1)
#                
#            nsbi = SystemBathInteraction(nsbi_op, cfm)

            # FIXME: Instead of all the above, we should have transformation
            # of SystemBathInteraction object


                    
            #
            # Hamiltonian matrix
            #
            hj = numpy.dot(numpy.linalg.inv(SS), numpy.dot(ham.JR,SS))
            for i in range(ham.dim):
                hj[i,i] = 0.0
            hh = numpy.diag(hD) + hj

            #
            # Foerster rates
            #
            KF = foerster_rates(Na, hh, tt, gvals, lamb)

#            nham = Hamiltonian(data=hh)
#            with energy_units("1/cm"):
#                print(nham.data)
#            if self._has_cutoff_time:
#                FT = FoersterRelaxationTensor(nham, nsbi,
#                                    cutoff_time=self.cutoff_time)
#            else:
#                FT = FoersterRelaxationTensor(nham, nsbi)                        
            
            
            # 
            # Add the rates to the Redfield
            #
            for b in range(Na):
                gg = 0.0
                for a in range(Na):
                    self.data[a,a,b,b] += KF[a,b]
                    gg += KF[a,b]
                self.data[b,b,b,b] += -gg

            
        self._is_initialized = True
        self._data_initialized = True
        
        
#    def transform(self, SS, inv=None):
#        """Quick fix before I find what is the problem
#        
#        My expectation was that RedfieldFoerster is calculated in 
#        the excitonic basis and needs to undergo a transformation
#        when it is entering code with no basis management context.
#        
#        """
#        super().transform(SS, inv)
#        self._data_initialized = True
        
