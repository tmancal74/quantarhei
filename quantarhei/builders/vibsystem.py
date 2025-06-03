# -*- coding: utf-8 -*-
import numpy

from ..core.managers import UnitsManaged
from .opensystem import OpenSystem
from ..core.saveable import Saveable
from .. import REAL
from .. import COMPLEX
from .. import Manager
from ..qm import LindbladForm
from ..qm.oscillators.ho import operator_factory


class VibrationalSystem(UnitsManaged, Saveable, OpenSystem):
    """ Represents a set of coupled oscillators (possibly unharmonic)


    This class forms the basis of the IR spectroscopy treatment
    in Quantarhei

    Parameters
    ----------

    name : str
        Specifies the name of the system

    modes : list or tuple
        List of modes out of which the systems is built

    """

    def __init__(self, modes=None, name=""):
        
        self.modes = modes
        self.Nmodes = len(self.modes)

        self.sbi = None
        self._has_sbi = False

        self.HH = None    
    
        self._mode_couping_init = False

        self._built = False

        ofac = operator_factory()
        self.ad = ofac.creation_operator()
        self.aa = ofac.anihilation_operator()
        self.EE = ofac.unity_operator()   

        #
        #  Relaxation tensor information
        #
        self.RelaxationTensor = None
        self.RelaxationHamiltonian = None
        self._has_relaxation_tensor = False
        self.has_Iterm = False


    def set_mode_coupling(self, N, M, val):
        """Sets bilinear coupling between the modes 
        
    
        >>> import quantarhei as qr
        >>> mod1 = qr.AnharmonicMode(100.0, 1.0/1000.0)
        >>> mod2 = qr.AnharmonicMode(100.0, 1.0/1000.0)
        >>> vsys = qr.VibrationalSystem(modes=[mod1, mod2])
        >>> with qr.energy_units("1/cm"):
        ...     vsys.set_mode_coupling(0,1, 10.0)
        ...     print(vsys.coupling[0,1])
        ...     print(vsys.coupling[1,0])
        0.00188365156731
        0.00188365156731
        >>> print(vsys.coupling[1,0])
        0.00188365156731
        
        """

        if not self._mode_couping_init:
            self.coupling = numpy.zeros((self.Nmodes, self.Nmodes), dtype=REAL)
            self._mode_couping_init = True

        val_int = Manager().convert_energy_2_internal_u(val)
        self.coupling[N,M] = val_int 
        self.coupling[M,N] = val_int


    def get_mode_coupling(self, N, M):
        """Returns the value of the intermode coupling coefficient in current units

        >>> import quantarhei as qr
        >>> mod1 = qr.AnharmonicMode(100.0, 1.0/1000.0)
        >>> mod2 = qr.AnharmonicMode(100.0, 1.0/1000.0)
        >>> vsys = qr.VibrationalSystem(modes=[mod1, mod2])
        >>> with qr.energy_units("1/cm"):
        ...     vsys.set_mode_coupling(0,1, 10.0)
        >>> print(vsys.get_mode_coupling(1,0))
        0.00188365156731
        >>> with qr.energy_units("1/cm"):
        ...     print(vsys.get_mode_coupling(0,1))
        10.0
        """
        return Manager().convert_energy_2_current_u(self.coupling[N,M])
    
        
    def _embed_operator(self, res_tot, member):

        pass

    def _embed_bilinear_interaction(self, i, j, kappa=1.0):
        """
        Embed a bilinear interaction term kappa * K1 ⊗ K2
        between subsystem i and j (0-based indices).

        Parameters:
        - i, j: indices of subsystems
        - dims: list of dimensions of all subsystems
        - kappa: coupling constant

        Returns:
        - Embedded interaction matrix on the full Hilbert space
        """
        N = len(self.dims)
        assert i != j, "Indices i and j must refer to different subsystems"

        # Ensure i < j for simpler indexing
        if i > j:
            i, j = j, i
            #K1, K2 = K2, K1

        # K1 and K2 are both coordinate operators
        qq = (self.ad + self.aa)/numpy.sqrt(2.0)
        K1 = qq[:self.dims[i],:self.dims[i]]
        K2 = qq[:self.dims[j],:self.dims[j]]

        ops = []
        for k in range(N):
            if k == i:
                ops.append(K1)
            elif k == j:
                ops.append(K2)
            else:
                ops.append(numpy.eye(self.dims[k], dtype=complex))

        # Construct the full tensor product
        interaction = ops[0]
        for op in ops[1:]:
            interaction = numpy.kron(interaction, op)

        return kappa * interaction
    

    def build(self):
        """Build the VibrationalSystem based on its components
        
        The VibrationalSystem can be built without system-bath interaction

        >>> import quantarhei as qr
        >>> mod1 = qr.AnharmonicMode(100.0, 1.0/1000.0)
        >>> mod2 = qr.AnharmonicMode(100.0, 1.0/1000.0)
        >>> vsys = qr.VibrationalSystem(modes=[mod1, mod2])
        >>> mod1.set_anharmonicity(0.01)
        >>> mod2.set_anharmonicity(0.01)
        >>> vsys.build()

        In this case, the relaxation tensor should be None

        >>> assert vsys.RelaxationTensor is None

        Vibrational system-bath with interaction not defined for all modes

        >>> mod1 = qr.AnharmonicMode(100.0, 1.0/1000.0, nmax=4)
        >>> mod2 = qr.AnharmonicMode(100.0, 1.0/1000.0, nmax=4)
        >>> mod1.set_anharmonicity(0.01)
        >>> mod2.set_anharmonicity(0.01)
        >>> vsys = qr.VibrationalSystem(modes=[mod1, mod2])
        >>> with qr.energy_units("1/cm"):
        ...     vsys.set_mode_coupling(0,1, 0.0)
        >>> mod1.set_mode_environment(1.0/1000.0, 1.0/500.0)
        >>> mod2.set_mode_environment(1.0/900.0, 1.0/300.0)

        >>> vsys.build()

        >>> assert numpy.amax(numpy.abs(vsys.RelaxationTensor.Ld)) > 0.0
        
        
        """
        from .. import SystemBathInteraction
        from .. import Hamiltonian
        from .. import TransitionDipoleMoment

        Nl = len(self.modes)

        #
        #  Hamiltonian
        #

        # list of Hamiltonian matrices
        H_list = []
        D_list = []
        for md in self.modes:
            #print("Mode:", mcount)
            md.build()
            Ham = md.get_Hamiltonian()
            H_list.append(Ham._data)
            TrDip = md.get_TransitionDipoleMoment()
            D_list.append(TrDip._data)
        
        # list of dimensions
        dims = [H.shape[0] for H in H_list]
        self.dims = dims

        # total Hamiltonian matrix
        H_tot_mx = numpy.zeros((numpy.prod(dims), numpy.prod(dims)), dtype=COMPLEX)

        for i, H in enumerate(H_list):
            ops = []
            for j in range(Nl):
                if j == i:
                    ops.append(H)
                else:
                    ops.append(numpy.eye(self.dims[j], dtype=complex))
            Hi_full = ops[0]
            for op in ops[1:]:
                Hi_full = numpy.kron(Hi_full, op)
            H_tot_mx += Hi_full

        self.HH = H_tot_mx
        self.HamOp = Hamiltonian(data=self.HH)

        # calculate bilinear interaction between modes
        if self._mode_couping_init:
            H_int = numpy.zeros_like(self.HH)
            for ii in range(Nl):
                for jj in range(ii+1, Nl):
                    kap = self.coupling[ii,jj]
                    if numpy.abs(kap) > 0.0:
                        H_int += self._embed_bilinear_interaction(ii, jj, kappa=kap)

            # add the interaction to the mode
            self.HamOp.data += H_int

        #
        # Transition dipole moment
        #
        # total trans dip matrix
        D_tot_mx = numpy.zeros((numpy.prod(dims), numpy.prod(dims), 3), dtype=COMPLEX)

        for i, DD in enumerate(D_list):
            for kk in range(3):
                ops = []
                for j in range(Nl):
                    if j == i:
                        ops.append(DD[:,:,kk])
                    else:
                        ops.append(numpy.eye(self.dims[j], dtype=COMPLEX))
                Di_full = ops[0]
                for op in ops[1:]:
                    Di_full = numpy.kron(Di_full, op)
                D_tot_mx[:,:,kk] += Di_full

        self.DD = D_tot_mx
        self.TrDMOp = TransitionDipoleMoment(data=D_tot_mx)

        #
        # Relaxation
        #

        #print("Dims = ", self.dims)
        ops = []
        rts = []
        mcount = 0
        skip_relaxation = False
        for mm, md in enumerate(self.modes):

            try:
                sbi = md.get_SystemBathInteraction()

            except:
                skip_relaxation = True

            if not skip_relaxation:
                nops = sbi.KK.shape[0]
               
                for ii in range(nops):

                    #print("Current operator is:", ii)
                    oper = sbi.KK[ii,:,:]

                    # here we have to stretch the operator to a new system basis
                    all_ops = []
                    for j in range(Nl):
                        if j == mm:
                            all_ops.append(oper)
                            #print("Shape oper -", j, oper.shape)
                        else:
                            all_ops.append(numpy.eye(self.dims[j], dtype=complex))
                            #print("Dims -", j, self.dims[j])

                    oper_full = all_ops[0]
                    for op in all_ops[1:]:
                        oper_full = numpy.kron(oper_full, op)

                    ops.append(oper_full)
                    rts.append(sbi.rates[ii])
            
            skip_relaxation = False

            mcount += 1

        if len(ops) > 0:

            sbi = SystemBathInteraction(sys_operators=ops, rates=rts)

            self.sbi = sbi
            self._has_sbi = True

        if self._has_sbi:
            self.RelaxationTensor = LindbladForm(self.HamOp, self.sbi, as_operators=True)
        else:
            self.RelaxationTensor = None

        self._built = True
        
