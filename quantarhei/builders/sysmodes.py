
import numpy

from ..utils import Float
from ..utils import Integer

from .opensystem import OpenSystem
from ..qm.hilbertspace.hamiltonian import Hamiltonian 
from ..qm.hilbertspace.dmoment import TransitionDipoleMoment
from ..qm.oscillators.ho import operator_factory
from ..qm import ProjectionOperator
from ..qm import SystemBathInteraction
from ..qm import LindbladForm
from .submodes import SubMode

from .. import REAL


class HarmonicMode(SubMode, OpenSystem):
    """Renaming the SubMode to be used as a standalone Harmonic oscillator mode
    
    
       PES shift is always 0 (will be implemented later)

    """

    def __init__(self, omega=1.0, shift=0.0, nmax=2):
        super().__init__(omega=omega, shift=shift, nmax=nmax)
        self._has_sbi = False
        self._built = False
        self.gamma = 0.0
        self.RT = None


    def build(self, nmax=None, build_relaxation=False):
        """Building all necessary quantities 
        
        
        
        """

        # if provided, we reset nmax
        if nmax is not None:
            self.nmax = nmax

        N = self.nmax

        ofac = operator_factory()
        ad = ofac.creation_operator()
        aa = ofac.anihilation_operator()

        #
        # Hamiltonian
        #
        EE = ofac.unity_operator()
        hh = self.omega*(numpy.dot(ad,aa) + EE*0.5)
        HH = numpy.zeros((N,N), dtype=REAL)
        HH[:,:] = hh[:N,:N]

        HH -= HH[0,0]*EE[:N,:N]
        
        self.HH = HH

        self.HamOp = Hamiltonian(data=HH)
        rind = numpy.array([n for n in range(N)], dtype=int)
        self.HamOp.set_rwa(rwa_indices=rind)

        #
        # Dipole moment
        #
        self._setup_dipole_moment(N, ad, aa)

        #
        # Relaxation 
        #
        if build_relaxation:
            if self._has_sbi:
                self.RT = LindbladForm(self.HamOp, self.sbi, as_operators=True)
            else:
                self.RT = None

        self._built = True
        self._diagonalized = True


    def diagonalize(self):
        """This problem is diagonal from the begining
        
        """
        if self._diagonalized:
            return
        

        self._diagonalized = True


    def _setup_dipole_moment(self, N, ad, aa, ss=None, pfac=1.0):

        DD = numpy.zeros((N,N,3), dtype=REAL)
        
        # FIXME: In what units we define transition dipole moment?
        dip = pfac*(ad + aa)/numpy.sqrt(2.0) 

        # transform if needed
        if ss is not None:
            s1 = numpy.linalg.inv(ss)
            dip = numpy.dot(s1, numpy.dot(dip,ss))

        DD[:,:,0] = numpy.real(dip[:N,:N])
        self.DD = DD

        # FIXME: make this on-demand (if poissible)
        trdata = numpy.zeros((DD.shape[0],DD.shape[1],DD.shape[2]),dtype=REAL)
        trdata[:,:,:] = DD[:,:,:]
        self.TrDMOp = TransitionDipoleMoment(data=trdata)


    # def get_Hamiltonian(self):
    #     """Returns the system Hamiltonian 
        
    #     """
    #     if self._built:
    #         return self.HamOp
    #     else:
    #         raise Exception("The Mode has to be built first.")


    # def get_TransitionDipoleMoment(self):
    #     """Returns the aggregate transition dipole moment operator

    #     """
    #     if self._built:
    #         return self.TrDMOp # TransitionDipoleMoment(data=self.DD)
    #     else:
    #         raise Exception("The Mode has to be built first.")    

    # def get_SystemBathInteraction(self):

    #     if self._built:
    #         return self.sbi
    #     else:
    #         raise Exception("The Mode has to be built first.")         



    def set_mode_environment(self, environ, pdeph=None):
        """
        
        """

        if isinstance(environ, (int, float)):

            self.gamma = environ
            self.pdeph = pdeph

            ops = []
            rts = []

            # relaxation
            if self.gamma > 0.0:
                for kk in range(self.nmax-1):
                    op = ProjectionOperator(to_state=kk, from_state=kk+1, dim=self.nmax)
                    ops.append(op)
                    rts.append(numpy.sqrt(float(kk+1))*self.gamma)

            # pure dephasing
            if (self.pdeph is not None) and (self.pdeph > 0.0):
                for kk in range(self.nmax):
                    op = ProjectionOperator(to_state=kk, from_state=kk, dim=self.nmax)
                    ops.append(op)
                    rts.append(numpy.sqrt(float(kk))*self.pdeph)              


            self.sbi = SystemBathInteraction(sys_operators=ops, rates=rts)

            self._has_sbi = True

        else:

            raise Exception("Environment not implemented")
        



class AnharmonicMode(HarmonicMode):
    """
    
    
    
    """
    def __init__(self, omega=1.0, shift=0.0, nmax=2):
        """
        
        In this case, omega is the difference between levels 1 and 0

        """
        super().__init__(omega=omega, shift=shift, nmax=nmax)

        # unharmonicity
        self.xi = 0.0
        self.om0 = self.omega*(1.0/(1.0-2.0*self.xi))
        self.dom = self.om0 - self.omega

        self._diagonalized = True


    def set_anharmonicity(self, xi):
        """Sets the ocillator anharmonicity
        
        """
        self.xi = xi
        # corresponding harmonic frequency
        self.om0 = self.omega*(1.0/(1.0-2.0*self.xi))
        self.dom = self.om0 - self.omega


    def get_anharmonicity(self, dom=None):
        """Returns the ahnarmonicity set previously, or calculates anharmonicity from frequency difference
        
        """
        if dom is None:
            return self.xi 
        else:
            dom_int = self.convert_energy_2_internal_u(dom)
            return dom_int/(2.0*(self.omega + dom_int))
        

    def build(self, nmax=None, xi=None, build_relaxation=False):
        """Building all necessary quantities """

        # if provided, we reset nmax
        if nmax is not None:
            self.nmax = nmax

        if xi is not None:
            self.set_unharmonicity(xi)

        N = self.nmax

        #
        #  Exact Hamiltonian
        #
        HH = numpy.zeros((N,N), dtype=REAL)

        # energy levels are known exactly
        for nn in range(N):
            HH[nn,nn] = self.om0*(nn + 0.5) - self.om0*self.xi*((nn + 0.5)**2)
            if nn == 0:
                hzer = HH[0,0]
            HH[nn,nn] -= hzer

        # saving the exact Hamiltonian (this line will be removed later)
        self.HamOp_ex = Hamiltonian(data=HH)

        #
        # Morse by diagonalization
        # 
        ofac = operator_factory()
        ad = ofac.creation_operator()
        aa = ofac.anihilation_operator()
        EE = ofac.unity_operator()
        qq = (ad + aa)/numpy.sqrt(2.0)
        pp = 1j*(ad - aa)/numpy.sqrt(2.0)

        hm = (self.omega/2.0)*(numpy.dot(pp,pp) + numpy.dot(qq,qq))

        # potential
        earg = numpy.sqrt(2.0*self.xi)*qq 
        de, ss = numpy.linalg.eigh(earg)

        # e^{-a q}
        e_argd = numpy.dot(ss, numpy.dot(numpy.diag(numpy.exp(-de)), numpy.linalg.inv(ss)))
        # 1 - e^{-a q}
        sqrtV = EE - e_argd
        # De*(1- e^{-a q})
        V_morse = (self.om0/(4.0*self.xi))*numpy.dot(sqrtV,sqrtV)

        hm = (self.om0/2.0)*numpy.dot(pp,pp) + V_morse

        hd, SS = numpy.linalg.eigh(hm)

        # saving transformation coefficients
        self.SS = SS[:N,:N]

        HM = numpy.zeros((N,N), dtype=REAL)
        HM[:,:] = numpy.real(numpy.diag(hd[:N]))

        hzer = HM[0,0]
        for nn in range(N):
            HM[nn,nn] -= hzer

        self.HH = HM

        self.HamOp = Hamiltonian(data=HM)
        rind = numpy.array([n for n in range(N)], dtype=int)
        self.HamOp.set_rwa(rwa_indices=rind)
        
        #
        #  Dipole moment 
        #
        self._setup_dipole_moment(N, ad, aa, ss=SS)


        #
        # Relaxation 
        #
        if build_relaxation:
            if self._has_sbi:
                self.RT = LindbladForm(self.HamOp, self.sbi, as_operators=True)
            else:
                self.RT = None

        self._built = True
        # this system is prepared in its eigenstate basis
        self._diagonalized = True


