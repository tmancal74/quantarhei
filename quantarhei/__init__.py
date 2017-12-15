# -*- coding: utf-8 -*-
"""
*******************************************************************************


    QUANTArhei: Open Quantum System Theory for Molecular Systems 
    ============================================================
    
    (c) 2016 Tomáš Mančal
    
    Charles University
    Faculty of Mathematics and Physics
    Ke Karlovu 5
    CZ-121 16 Prague 2
    Czech Repubic


    For support contact the author at : mancal@karlov.mff.cuni.cz
    
    
*******************************************************************************

    Classes available from quantarhei are devided into several functional
    groups. They can all loaded similarly to the Manager class as
    
    >>> from quantarhei import Manager
    >>> m = qr.Manager()
    >>> print(m.version)
    
    Preferred way of using Quantarhei is to load the package and rename it 
    to something shorter, like below
    
    >>> import quantarhei as qr
    >>> m = qr.Manager()
    >>> print(m.version)
  
    Numeric types
    -------------
    
    REAL ...... real floating point number type, by default numpy.float64
    COMPLEX ... complex floating point number type, by default numpy.complex128 
    
    Builders
    --------
    
    Mode .......... represents a harmonic vibrational mode of a molecule
    Molecule ...... represents a molecule
    Aggregate ..... represents an aggregate of molecules
    PDBFile ....... reader and writter of structures from PDB format
    Disorder ...... class managing static disorder of molecular transition
                    energies
     
    Core classes
    ------------
    
    TimeAxis ......... linear axis of real values representing discrete time
    FrequencyAxis .... linear axis of real values representing discrete
                       frequency axis
    DFunction ........ discrete function
    
    
    Various managers
    ----------------
    
    Manager ............ the main behind-the-scenes manager of the package
    energy_units ....... energy units manager for use with the "with" construct
    frequency_units .... frequency units manager for use with 
                         the "with" construct
    eigenbasis_of ...... manager of the basis transformations to be used with 
                         the "with" construct
    set_current_units .. function to set current units globally
    
    ... to be continued

*******************************************************************************

"""


###############################################################################
#
#
#            Imports of high level classes and functions 
#
#
###############################################################################

#
# Fix used numerical types
#
#import numpy
from .core.managers import Manager
m = Manager()

REAL = m.get_real_type() #numpy.float64
COMPLEX = m.get_complex_type() #numpy.complex128

LOG_URGENT = 0
LOG_REPORT = 3
LOG_INFO = 5
LOG_DETAIL = 7
LOG_QUICK = 9

#
# Builders
#
from .builders.modes import Mode
from .builders.molecules import Molecule
from .builders.aggregates import Aggregate
from .builders.pdb import PDBFile
from .builders.disorder import Disorder

#
# Core classes
#
from .core.time import TimeAxis
from .core.frequency import FrequencyAxis
from .core.dfunction import DFunction
from .core.saveable import Saveable

#
# Various managers
#
from .core.managers import energy_units
from .core.managers import frequency_units
from .core.managers import length_units
from .core.managers import eigenbasis_of
from .core.managers import set_current_units

#
# Parallelization
#
from .core.parallel import distributed_configuration
from .core.parallel import start_parallel_region
from .core.parallel import close_parallel_region
from .core.parallel import parallel_function
from .core.parallel import block_distributed_range

###############################################################################
#                            SPECTROSCOPY
###############################################################################

#
# Linear absorption 
#
from .spectroscopy.abs2 import AbsSpectrum
from .spectroscopy.abs2 import AbsSpectrumContainer
from .spectroscopy.abs2 import AbsSpectrumCalculator

#
# Fourier transform Two-Dimensional Spectra
#
from .spectroscopy.twod2 import TwoDSpectrum
from .spectroscopy.twod2 import TwoDSpectrumContainer
from .spectroscopy.twod2 import TwoDSpectrumCalculator
from .spectroscopy.twod2 import MockTwoDSpectrumCalculator


from .spectroscopy.pathwayanalyzer import LiouvillePathwayAnalyzer

###############################################################################
#                           QUANTUM MECHANICS
###############################################################################


#
# Operators
#
from .qm import StateVector
from .qm import DensityMatrix
from .qm import ReducedDensityMatrix
from .qm import BasisReferenceOperator
from .qm import Hamiltonian
from .qm import TransitionDipoleMoment

#
# Propagators
#
from .qm.propagators.poppropagator import PopulationPropagator
from .qm.propagators.svpropagator import StateVectorPropagator
from .qm import ReducedDensityMatrixPropagator

#
# Evolutions (time-dependent operators)
#
from .qm.propagators.statevectorevolution import StateVectorEvolution
from .qm import DensityMatrixEvolution
from .qm import ReducedDensityMatrixEvolution

#
# System-bath interaction
#
from .qm.corfunctions import CorrelationFunction
from .qm.corfunctions import SpectralDensity



###############################################################################
# Convenience functions
###############################################################################
from .core.saveable import load
from .core.saveable import read_info

from .core.units import convert
from .core.units import in_current_units

from .utils.vectors import normalize2
from .utils.vectors import norm 

from .utils.logging import printlog
from .utils.logging import loglevels2bool
from .utils.logging import log_urgent
from .utils.logging import log_report
from .utils.logging import log_info
from .utils.logging import log_detail
from .utils.logging import log_quick



