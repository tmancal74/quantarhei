"""Quantarhei User Level Classes and Objects
=========================================

In Quantarhei, classes are losely grouped into three categories. First,
there is agroup of classes, which represent basic concepts of quantum
mechanics, provide access to important implementations of spectroscopic
simulations and dynamics of open quantum systems, and classes which allow
basic management of the simulation environment and numerical results.
These classed are called **user level classes**, and they are all
accessible in highest namespace level of the Quantarhei package.
If you import Quantarhei like this:

>>> import quantarhei as qr

you can access user level classes through the qr. prefix, e.g.


>>> manager = qr.Manager()
>>> print(manager.version)
0.0.63

The list of user level classes is provided below. Tue latest and most
uptodate information can be obtained by viewing the source code of the
root `__init__.py` file of the packages. All classes imported there are
considered user level classes.


Other Class Levels
------------------

In this documentation we recognize two more groups (or levels) of classes.
More specialized classes, which normal user does not need as often as the
user level classes are called **advanced level classes**. These use the
second level name space. For instance the class `SystemBathInteraction`
is relatively rarely used directly. It is therefore *hidden* in the name
space `qm` (as quantum mechanics) of the package. This class can be
instantiated e.g. like this

>>> import quantarhei as qr
>>> sbi = qr.qm.SystemBathInteraction()

Advanced level classes are still intendend for relatively frequent use
by the user. However, in order to reduced the *apparent* complexity of
basic usage of Quantarhei, advanced level classes are documented in their
respective sub-packages, one level deeper than user level classes. Complete
documentation of advanced level classes is available in the Advanced Level
Classes section of this documentation.

Everything else in Quantarhei package goes under the banner of
**expert level classes**. This includes all classes and objects used
internally in Quantarhei. We make every effort to document also this part
of the package as completely as possible, but it is the last item on the
list, so to say. The user is welcome to learn and use the expert level
classes, but our aim is to structure Quantarhei in such a way, that this
is not necessary. More on expert level classes in the section in
Quantarhei internals.

User Level Objects and Convenience Functions
============================================

Besides classes, Quantarhei also defines some user level objects and
convenience functions. They are listed here under several categories

Numeric types
-------------

.. toctree::
:maxdepth: 2

functions/numtypes

Convenience Functions
---------------------

.. toctree::
:maxdepth: 2

functions/convenience


Logging Functions and Loglevels
-------------------------------

.. toctree::
:maxdepth: 2

functions/logging


..
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

LOG_URGENT = 1
LOG_REPORT = 3
LOG_INFO = 5
LOG_DETAIL = 7
LOG_QUICK = 9


#
# Non-linear response signals
#
signal_REPH = "rephasing_2D_signal"
signal_NONR = "nonrephasing_2D_signal"
signal_TOTL = "total_2D_signal"
signal_DC = "double_coherence_signal"

TWOD_SIGNALS = dict(signal_REPH = "rephasing_2D_signal",
                    signal_NONR = "nonrephasing_2D_signal",
                    signal_TOTL = "total_2D_signal",
                    signal_DC = "double_coherence_signal")

#
# Parts of the complex data/signal
#
part_REAL = "real_part"
part_IMAGINARY = "imaginary_part"
part_COMPLEX = "complex"
part_ABS = "absolute_value"
part_PHASE = "phase"

SIGNAL_PARTS = dict(part_REAL = "real_part",
                    part_IMAGINARY= "imaginary_part",
                    part_COMPLEX = "complex",
                    part_ABS = "absolute_value",
                    part_PHASE = "phase")

DATA_PARTS = SIGNAL_PARTS

#
# Liouville pathway types
#
ptype_R1g = "pathway_type_R1g"
ptype_R2g = "pathway_type_R2g"
ptype_R3g = "pathway_type_R3g"
ptype_R4g = "pathway_type_R4g"
ptype_R1f = "pathway_type_R1f*"
ptype_R2f = "pathway_type_R2f*"
ptype_R3f = "pathway_type_R3f"
ptype_R4f = "pathway_type_R4f"

PATHWAY_TYPES = dict(ptype_R1g = "pathway_type_R1g",
                     ptype_R2g = "pathway_type_R2g",
                     ptype_R3g = "pathway_type_R3g",
                     ptype_R4g = "pathway_type_R4g",
                     ptype_R1f = "pathway_type_R1f*",
                     ptype_R2f = "pathway_type_R2f*",
                     ptype_R3f = "pathway_type_R3f",
                     ptype_R4f = "pathway_type_R4f")

LIOUVILLE_PATHWAY_TYPES = PATHWAY_TYPES

#
# Builders
#
from .builders.aggregate_test import TestAggregate
from .builders.aggregates import Aggregate
from .builders.disorder import Disorder
from .builders.modes import Mode
from .builders.molecule_test import TestMolecule
from .builders.molecules import Molecule
from .builders.opensystem import OpenSystem
from .builders.pdb import PDBFile
from .builders.sysmodes import AnharmonicMode, HarmonicMode
from .builders.vibsystem import VibrationalSystem
from .core.dfunction import DFunction
from .core.frequency import FrequencyAxis

#
# Various managers
#
from .core.managers import (
    eigenbasis_of,
    energy_units,
    frequency_units,
    length_units,
    set_current_units,
)

#
# Parallelization
#
from .core.parallel import (
    asynchronous_range,
    block_distributed_array,
    block_distributed_list,
    block_distributed_range,
    close_parallel_region,
    collect_block_distributed_data,
    distributed_configuration,
    parallel_function,
    start_parallel_region,
)

###############################################################################
# Convenience functions
###############################################################################
#from .core.saveable import load
#from .core.saveable import read_info
from .core.parcel import Parcel, check_parcel, load_parcel, save_parcel

#from .core.saveable import Saveable
from .core.saveable import Saveable

#
# Core classes
#
from .core.time import TimeAxis
from .core.units import convert, in_current_units
from .core.valueaxis import ValueAxis

###############################################################################
#                           QUANTUM MECHANICS
###############################################################################
#
# State vectors
#
#
# Operators
#
from .qm import (
    BasisReferenceOperator,
    DensityMatrix,
    DensityMatrixEvolution,
    Hamiltonian,
    Liouvillian,
    OQSStateVector,
    OQSStateVectorEvolution,
    OQSStateVectorPropagator,
    ProjectionOperator,
    ReducedDensityMatrix,
    ReducedDensityMatrixEvolution,
    ReducedDensityMatrixPropagator,
    StateVector,
    SystemBathInteraction,
    TransitionDipoleMoment,
    UnityOperator,
)

#
# System-bath interaction
#
#
# LINESHAPE FUNCTIONS
#
from .qm.corfunctions import (
    CorrelationFunction,
    CorrelationFunctionMatrix,
    FastFunctionStorage,
    FunctionStorage,
    LineshapeFunction,
    SpectralDensity,
)
from .qm.corfunctions.correlationfunctions import oscillator_scalled_CorrelationFunction

#
# Evolution operators
#
from .qm.liouvillespace.evolutionsuperoperator import EvolutionSuperOperator
from .qm.liouvillespace.heom import (
    KTHierarchy,
    KTHierarchyPropagator,
    QuTip_KTHierarchyPropagator,
)

#
# Propagators
#
from .qm.propagators.poppropagator import PopulationPropagator

#
# Evolutions (time-dependent operators)
#
from .qm.propagators.statevectorevolution import StateVectorEvolution
from .qm.propagators.svpropagator import StateVectorPropagator

###############################################################################
#                            SPECTROSCOPY
###############################################################################
#
# Linear absorption
#
from .spectroscopy.abs2 import AbsSpectrum
from .spectroscopy.abscalculator import AbsSpectrumCalculator
from .spectroscopy.abscontainer import AbsSpectrumContainer

#
# Circular dichroism
#
from .spectroscopy.circular_dichroism import (
    CircDichSpectrum,
    CircDichSpectrumCalculator,
    CircDichSpectrumContainer,
)
from .spectroscopy.dsfeynman import (
    DSFeynmanDiagram,
    R1f_Diagram,
    R1g_Diagram,
    R1g_R_Diagram,
    R2f_Diagram,
    R2g_Diagram,
    R3g_Diagram,
    R4g_Diagram,
)

#
# Fluorescence
#
from .spectroscopy.fluorescence import (
    FluorSpectrum,
    FluorSpectrumCalculator,
    FluorSpectrumContainer,
)
from .spectroscopy.labsetup import LabField, LabSetup

#
# Linear dichroism
#
from .spectroscopy.linear_dichroism import (
    LinDichSpectrum,
    LinDichSpectrumCalculator,
    LinDichSpectrumContainer,
)
from .spectroscopy.mockabscalculator import MockAbsSpectrumCalculator
from .spectroscopy.mocktwodcalculator import MockTwoDResponseCalculator
from .spectroscopy.pathwayanalyzer import LiouvillePathwayAnalyzer

#
# Pump-probe spectrum
#
from .spectroscopy.pumpprobe import (
    MockPumpProbeSpectrumCalculator,
    PumpProbeSpectrum,
    PumpProbeSpectrumCalculator,
    PumpProbeSpectrumContainer,
)

# 2 deprecated classes
from .spectroscopy.responses import (
    LiouvillePathway,
    NonLinearResponse,
    ResponseFunction,
)
from .spectroscopy.twodcalculator import TwoDResponseCalculator
from .spectroscopy.twodcontainer import TwoDResponseContainer, TwoDSpectrumContainer

#
# Fourier transform Two-Dimensional Spectra
#
from .spectroscopy.twodresponse import TwoDResponse
from .spectroscopy.twodspect import TwoDSpectrum
from .symbolic.cumulant import evaluate_cumulant
from .utils.logging import (
    init_logging,
    log_detail,
    log_info,
    log_quick,
    log_report,
    log_to_file,
    log_urgent,
    loglevels2bool,
    printlog,
    tprint,
)
from .utils.paver import execute_paver
from .utils.timing import done_in, finished_in, timeit, untimeit
from .utils.vectors import norm, normalize2
from .wizard.input.input import Input


def exit(msg=None):
    """Exit to the level above the script with SystemExit exception

    """
    import sys
    if msg is not None:
        printlog("\n(SystemExit) Message: "+msg+"\n", loglevel=0)
    sys.exit()


def stop(msg=None):
    """Stop execution and leave to level above

    """
    exit("Execution stopped")



def show_plot(block=True):
    """Shows current plot

    This function is used to avoid explicit import of matplotlib
    """
    import matplotlib.pyplot as plt
    plt.show(block=block)



def savefig(fname):
    """Saves current plot to a file

    This function is used to avoid explicit import of matplotlib
    """
    import matplotlib.pyplot as plt
    plt.savefig(fname)


def assert_version(check, vno):
    """Throws an exception if the condition is not satisfied

    """
    from packaging import version

    def ext():
        exit("Version requirement not satisfied.")

    if check == ">=":
        if not (version.parse(Manager().version) >= version.parse(vno)):
            ext()

    elif check == "==":
         if not (version.parse(Manager().version) == version.parse(vno)):
            ext()

    elif check == "<=":
        if not (version.parse(Manager().version) <= version.parse(vno)):
            ext()

    elif check == ">":
        if not (version.parse(Manager().version) == version.parse(vno)):
            ext()

    elif check == "<=":
        if not (version.parse(Manager().version) >= version.parse(vno)):
            ext()

    else:
        raise Exception("Unknown comparison operator `"+check+"`")


