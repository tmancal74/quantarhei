# -*- coding: utf-8 -*-
"""

Lindblad form implementation


"""
import numpy

from .redfieldtensor import RedfieldRelaxationTensor
from ...builders.aggregates import Aggregate
from .systembathinteraction import SystemBathInteraction
from ...builders.aggregate_states import VibronicState
from ... import REAL
from ..hilbertspace.operators import ProjectionOperator

class LindbladForm(RedfieldRelaxationTensor):
    """Lindblad form of relaxation tensor

    We use the Redfield tensor class and reimplement its implementation

    """

    def __init__(self, ham, sbi, initialize=True,
                 as_operators=True, name=""):
        super().__init__(ham, sbi, initialize=initialize,
                         as_operators=as_operators, name=name)


    def _implementation(self, ham, sbi):
        """Very simple implementation of Lindblad form

        """
        # dimension of the operators
        Na = ham.dim

        # number of operators
        if sbi is None:
            Nb = 1
        else:
            Nb = sbi.N

        llm = numpy.zeros((Nb, Na, Na), dtype=REAL)
        lld = numpy.zeros((Nb, Na, Na), dtype=REAL)
        
        if sbi is not None:
            for i in range(Nb):
                llm[i, :, :] = sbi.rates[i]*sbi.KK[i, :, :]/2.0
                lld[i, :, :] = numpy.transpose(llm[i, :, :])

        if sbi is None:
            KK = numpy.zeros((1, Na, Na), dtype=REAL)
        else:
            KK = sbi.KK
            
        self._post_implementation(KK, llm, lld)



class ElectronicLindbladForm(LindbladForm):
    """Lindblad for corresponding to purely electronic relaxation


    Purely electronic relaxation defined in the site basis


    Parameters
    ----------
    
    ham : Hamiltonian
        Hamiltonian of the system (aggregate)
        
    sbi : SystemBathInteraction
        Object describing system--bath interaction. For Lindblad form, this
        object holds projection operators describing transitions which
        proceeds with a given rate. Rates are also specified by ``sbi`` object.
        In addition, ``sbi`` has to have an attribute ``system`` which holds
        the Aggregate object of the system
        
    initialize: boolean
        if True, Lindblad form will be immediately initialized
        
    as_operators: boolean
        if True, the Lindblad form will be represented by operators and not
        by a tensor
        
    name : str
        Name of the ElectronicLindbladForm object


    """

    def __init__(self, ham, sbi, initialize=True,
                 as_operators=True, name=""):

        if isinstance(sbi.system, Aggregate):
            agg = sbi.system

            if agg.Nel == agg.Ntot:
                # if the number of electronic states corresponds to the total
                # number of states, the aggregate is purely electronic. There
                # is then no difference between Lindblad form and its
                # electronic version. We call super()
                super().__init__(ham, sbi, initialize=initialize,
                                 as_operators=as_operators, name=name)
            else:
                # check the sbi for rates and (electronic) operators
                
                # sbi operators have to have the same dimension as
                # the single exciton electronic part of the aggregate
                if agg.mult == 1:
                    Nel1 = 1 + agg.nmono
                elif (agg.mult == 2) and (agg.sbi_mult==1):
                    Nel1 = 1 + agg.nmono 
                else:
                    raise Exception("Cannot yet handle"+
                                    "the case of sbi_mult> 1")
                    
                if Nel1 == sbi.KK.shape[1]:
                    
                    # create new interaction operators of higher 
                    # dimensionality
                    Nop = sbi.KK.shape[0]
                    ops = []
                    for k in range(Nop):
                        newkk = numpy.zeros((agg.Ntot, agg.Ntot), 
                                            dtype=numpy.float64)
                        # populate the operator
                        for i_el in range(agg.Nel):
                            for i_vib in agg.vibindices[i_el]:
                                
                                vs_i = agg.vibsigs[i_vib]
                                st_i = agg.get_VibronicState(vs_i[0], vs_i[1])
                                
                                for j_el in range(agg.Nel):
                                    for j_vib in agg.vibindices[j_el]:
                                
                                        vs_j = agg.vibsigs[j_vib]
                                        st_j = agg.get_VibronicState(vs_j[0],
                                                                     vs_j[1])
                                
                                        # electronic transition operator
                                        # dressed in Franck-Condon factors
                                        newkk[i_vib, j_vib] = (
                                        numpy.real(agg.fc_factor(st_i, st_j))*
                                        sbi.KK[k, i_el, j_el])
                    
                        ops.append(newkk)
                    
                    # with the operators constructed, we create Lindblad form
                    newsbi = SystemBathInteraction(sys_operators=ops,
                                                   rates=sbi.rates,
                                                   system=sbi.system)
                    super().__init__(ham, newsbi, initialize=initialize,
                         as_operators=as_operators, name=name)

                
                else:
                    raise Exception("Incompatible dimension of system-bath"+
                                    " interaction operators")

        else:

            raise Exception("SystemBathInteraction has to have `system`"+
                            " attribute set to an Aggregate")


class VibrationalDecayLindbladForm(LindbladForm):
    """Lindblad for corresponding to purely vibrational relaxation


    Purely vibrational relaxation defined in the site basis for an
    oscillator mode residing on a molecule, which is a member of an aggregate


    Parameters
    ----------
    
    ham : Hamiltonian
        Hamiltonian of the system (aggregate)
        
    sbi : SystemBathInteraction
        Object describing system--bath interaction. For Lindblad form, this
        object holds projection operators describing transitions which
        proceeds with a given rate. Rates are also specified by ``sbi`` object.
        In addition, ``sbi`` has to have an attribute ``system`` which holds
        the Aggregate object of the system
        
    initialize: boolean
        if True, Lindblad form will be immediately initialized
        
    as_operators: boolean
        if True, the Lindblad form will be represented by operators and not
        by a tensor
        
    name : str
        Name of the VibrationalDecayLindbladForm object


    """

    def __init__(self, ham, sbi, initialize=True,
                 as_operators=True, name=""):

        if isinstance(sbi.system, Aggregate):
            agg = sbi.system
            
            orates = sbi.orates
            osites = sbi.osites
            Ntot = agg.Ntot
            
            ops = []
            rts = []
            
            zero = numpy.zeros((Ntot,Ntot), dtype=REAL)
            zrs = 0
            
            # we loop over sites in which we want to introduce relaxation
            k = 0
            nops = 0
            for site in osites:
                
                # rates for each site
                rate = orates[k]
                
                # get max number of states for the site 
                nmax = 5
                
                for nk in range(nmax):
                    
                    nn = float(nk+1)
                    pair = [nk, nk+1]
 
                    # projection operator for each transition
                    op = ProjectionOperator(dim=Ntot)
                   
                    # loop over all a matrix of states
                    for a, s1 in agg.all_states:
                    
                        sig1 = s1.signature()
                        el1 = sig1[0]
                        vb1 = sig1[1]
                    
                        for b, s2 in agg.all_states:
                        
                            sig2 = s2.signature()
                            el2 = sig2[0]
                            vb2 = sig2[1]
                        
                            op.data[a,b] = 0.0
                        
                            if ((el1 == el2) and 
                                self._same_regardless_of_one(vb1, vb2, site)):
                            
                                if (vb1[site] == pair[0] 
                                    and vb2[site] == pair[1]):
                                    #print("site:", site, ":", pair[0],"<-",
                                    #      pair[1],"|", sig1, "<--", sig2)
                                
                                    op.data[a,b] = 1.0 #numpy.sqrt(nn)

                    if not numpy.array_equal(op.data, zero):
                        ops.append(op)
                        nops += 1
                        rts.append(rate*numpy.sqrt(nn))
                    else:
                        zrs += 1
                
                k += 1
                
            #print("Number of operators: ", nops)
            #print("Number of rates    : ", len(rts), rts)
            #print("Number of discarted ops", zrs)
    
        else:

            raise Exception("SystemBathInteraction has to have `system`"+
                            " attribute set to an Aggregate")
            
            
        # with the operators constructed, we create Lindblad form
        newsbi = SystemBathInteraction(sys_operators=ops,
                                       rates=rts,
                                       system=sbi.system)
        super().__init__(ham, newsbi, initialize=initialize,
                         as_operators=as_operators, name=name)
        
            
    def _same_regardless_of_one(self, a, b, k):
        
        if (a[:k] == b[:k]) and (a[k+1:] == b[k+1:]):
            return True
        
        
        
        