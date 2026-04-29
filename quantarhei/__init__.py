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

from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
from importlib.metadata import version as _version

try:
    __version__ = _version("quantarhei")
except _PackageNotFoundError:
    __version__ = "unknown"


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
# import numpy
from .core.managers import Manager

m = Manager()

REAL: type = m.get_real_type()  # numpy.float64
COMPLEX: type = m.get_complex_type()  # numpy.complex128

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

TWOD_SIGNALS: dict[str, str] = dict(
    signal_REPH="rephasing_2D_signal",
    signal_NONR="nonrephasing_2D_signal",
    signal_TOTL="total_2D_signal",
    signal_DC="double_coherence_signal",
)

#
# Parts of the complex data/signal
#
part_REAL = "real_part"
part_IMAGINARY = "imaginary_part"
part_COMPLEX = "complex"
part_ABS = "absolute_value"
part_PHASE = "phase"

SIGNAL_PARTS = dict(
    part_REAL="real_part",
    part_IMAGINARY="imaginary_part",
    part_COMPLEX="complex",
    part_ABS="absolute_value",
    part_PHASE="phase",
)

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

PATHWAY_TYPES = dict(
    ptype_R1g="pathway_type_R1g",
    ptype_R2g="pathway_type_R2g",
    ptype_R3g="pathway_type_R3g",
    ptype_R4g="pathway_type_R4g",
    ptype_R1f="pathway_type_R1f*",
    ptype_R2f="pathway_type_R2f*",
    ptype_R3f="pathway_type_R3f",
    ptype_R4f="pathway_type_R4f",
)

LIOUVILLE_PATHWAY_TYPES = PATHWAY_TYPES

#
# Builders
#
from .builders.aggregate_test import TestAggregate as TestAggregate
from .builders.aggregates import Aggregate as Aggregate
from .builders.disorder import Disorder as Disorder
from .builders.modes import Mode as Mode
from .builders.molecule_test import TestMolecule as TestMolecule
from .builders.molecules import Molecule as Molecule
from .builders.opensystem import OpenSystem as OpenSystem
from .builders.pdb import PDBFile as PDBFile
from .builders.sysmodes import AnharmonicMode as AnharmonicMode
from .builders.sysmodes import HarmonicMode as HarmonicMode
from .builders.vibsystem import VibrationalSystem as VibrationalSystem
from .core.dfunction import DFunction as DFunction
from .core.frequency import FrequencyAxis as FrequencyAxis

#
# Various managers
#
from .core.managers import (
    eigenbasis_of as eigenbasis_of,
)
from .core.managers import (
    energy_units as energy_units,
)
from .core.managers import (
    frequency_units as frequency_units,
)
from .core.managers import (
    length_units as length_units,
)
from .core.managers import (
    set_current_units as set_current_units,
)

#
# Parallelization
#
from .core.parallel import (
    asynchronous_range as asynchronous_range,
)
from .core.parallel import (
    block_distributed_array as block_distributed_array,
)
from .core.parallel import (
    block_distributed_list as block_distributed_list,
)
from .core.parallel import (
    block_distributed_range as block_distributed_range,
)
from .core.parallel import (
    close_parallel_region as close_parallel_region,
)
from .core.parallel import (
    collect_block_distributed_data as collect_block_distributed_data,
)
from .core.parallel import (
    distributed_configuration as distributed_configuration,
)
from .core.parallel import (
    parallel_function as parallel_function,
)
from .core.parallel import (
    start_parallel_region as start_parallel_region,
)

###############################################################################
# Convenience functions
###############################################################################
# from .core.saveable import load
# from .core.saveable import read_info
from .core.parcel import Parcel as Parcel
from .core.parcel import check_parcel as check_parcel
from .core.parcel import load_parcel as load_parcel
from .core.parcel import save_parcel as save_parcel

# from .core.saveable import Saveable
from .core.saveable import Saveable as Saveable

#
# Core classes
#
from .core.time import TimeAxis as TimeAxis
from .core.units import convert as convert
from .core.units import in_current_units as in_current_units
from .core.valueaxis import ValueAxis as ValueAxis

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
    BasisReferenceOperator as BasisReferenceOperator,
)
from .qm import (
    DensityMatrix as DensityMatrix,
)
from .qm import (
    DensityMatrixEvolution as DensityMatrixEvolution,
)
from .qm import (
    Hamiltonian as Hamiltonian,
)
from .qm import (
    Liouvillian as Liouvillian,
)
from .qm import (
    OQSStateVector as OQSStateVector,
)
from .qm import (
    OQSStateVectorEvolution as OQSStateVectorEvolution,
)
from .qm import (
    OQSStateVectorPropagator as OQSStateVectorPropagator,
)
from .qm import (
    ProjectionOperator as ProjectionOperator,
)
from .qm import (
    ReducedDensityMatrix as ReducedDensityMatrix,
)
from .qm import (
    ReducedDensityMatrixEvolution as ReducedDensityMatrixEvolution,
)
from .qm import (
    ReducedDensityMatrixPropagator as ReducedDensityMatrixPropagator,
)
from .qm import (
    StateVector as StateVector,
)
from .qm import (
    SystemBathInteraction as SystemBathInteraction,
)
from .qm import (
    TransitionDipoleMoment as TransitionDipoleMoment,
)
from .qm import (
    UnityOperator as UnityOperator,
)

#
# System-bath interaction
#
#
# LINESHAPE FUNCTIONS
#
from .qm.corfunctions import (
    CorrelationFunction as CorrelationFunction,
)
from .qm.corfunctions import (
    CorrelationFunctionMatrix as CorrelationFunctionMatrix,
)
from .qm.corfunctions import (
    FastFunctionStorage as FastFunctionStorage,
)
from .qm.corfunctions import (
    FunctionStorage as FunctionStorage,
)
from .qm.corfunctions import (
    LineshapeFunction as LineshapeFunction,
)
from .qm.corfunctions import (
    SpectralDensity as SpectralDensity,
)
from .qm.corfunctions.correlationfunctions import (
    oscillator_scalled_CorrelationFunction as oscillator_scalled_CorrelationFunction,
)

#
# Evolution operators
#
from .qm.liouvillespace.evolutionsuperoperator import (
    EvolutionSuperOperator as EvolutionSuperOperator,
)
from .qm.liouvillespace.heom import (
    KTHierarchy as KTHierarchy,
)
from .qm.liouvillespace.heom import (
    KTHierarchyPropagator as KTHierarchyPropagator,
)
from .qm.liouvillespace.heom import (
    QuTip_KTHierarchyPropagator as QuTip_KTHierarchyPropagator,
)

#
# Propagators
#
from .qm.propagators.poppropagator import PopulationPropagator as PopulationPropagator

#
# Evolutions (time-dependent operators)
#
from .qm.propagators.statevectorevolution import (
    StateVectorEvolution as StateVectorEvolution,
)
from .qm.propagators.svpropagator import StateVectorPropagator as StateVectorPropagator

###############################################################################
#                            SPECTROSCOPY
###############################################################################
#
# Linear absorption
#
from .spectroscopy.abs2 import AbsSpectrum as AbsSpectrum
from .spectroscopy.abscalculator import AbsSpectrumCalculator as AbsSpectrumCalculator
from .spectroscopy.abscontainer import AbsSpectrumContainer as AbsSpectrumContainer

#
# Circular dichroism
#
from .spectroscopy.circular_dichroism import (
    CircDichSpectrum as CircDichSpectrum,
)
from .spectroscopy.circular_dichroism import (
    CircDichSpectrumCalculator as CircDichSpectrumCalculator,
)
from .spectroscopy.circular_dichroism import (
    CircDichSpectrumContainer as CircDichSpectrumContainer,
)
from .spectroscopy.dsfeynman import (
    DSFeynmanDiagram as DSFeynmanDiagram,
)
from .spectroscopy.dsfeynman import (
    R1f_Diagram as R1f_Diagram,
)
from .spectroscopy.dsfeynman import (
    R1g_Diagram as R1g_Diagram,
)
from .spectroscopy.dsfeynman import (
    R1g_R_Diagram as R1g_R_Diagram,
)
from .spectroscopy.dsfeynman import (
    R2f_Diagram as R2f_Diagram,
)
from .spectroscopy.dsfeynman import (
    R2g_Diagram as R2g_Diagram,
)
from .spectroscopy.dsfeynman import (
    R3g_Diagram as R3g_Diagram,
)
from .spectroscopy.dsfeynman import (
    R4g_Diagram as R4g_Diagram,
)

#
# Fluorescence
#
from .spectroscopy.fluorescence import (
    FluorSpectrum as FluorSpectrum,
)
from .spectroscopy.fluorescence import (
    FluorSpectrumCalculator as FluorSpectrumCalculator,
)
from .spectroscopy.fluorescence import (
    FluorSpectrumContainer as FluorSpectrumContainer,
)
from .spectroscopy.labsetup import LabField as LabField
from .spectroscopy.labsetup import LabSetup as LabSetup

#
# Linear dichroism
#
from .spectroscopy.linear_dichroism import (
    LinDichSpectrum as LinDichSpectrum,
)
from .spectroscopy.linear_dichroism import (
    LinDichSpectrumCalculator as LinDichSpectrumCalculator,
)
from .spectroscopy.linear_dichroism import (
    LinDichSpectrumContainer as LinDichSpectrumContainer,
)
from .spectroscopy.mockabscalculator import (
    MockAbsSpectrumCalculator as MockAbsSpectrumCalculator,
)
from .spectroscopy.mocktwodcalculator import (
    MockTwoDResponseCalculator as MockTwoDResponseCalculator,
)
from .spectroscopy.pathwayanalyzer import (
    LiouvillePathwayAnalyzer as LiouvillePathwayAnalyzer,
)

#
# Pump-probe spectrum
#
from .spectroscopy.pumpprobe2 import (
    MockPumpProbeSpectrumCalculator as MockPumpProbeSpectrumCalculator,
)
from .spectroscopy.pumpprobe2 import (
    PumpProbeSpectrum as PumpProbeSpectrum,
)
from .spectroscopy.pumpprobe2 import (
    PumpProbeSpectrumCalculator as PumpProbeSpectrumCalculator,
)
from .spectroscopy.pumpprobe2 import (
    PumpProbeSpectrumContainer as PumpProbeSpectrumContainer,
)

# 2 deprecated classes
from .spectroscopy.responses import (
    LiouvillePathway as LiouvillePathway,
)
from .spectroscopy.responses import (
    NonLinearResponse as NonLinearResponse,
)
from .spectroscopy.responses import (
    ResponseFunction as ResponseFunction,
)
from .spectroscopy.twodcalculator import (
    TwoDResponseCalculator as TwoDResponseCalculator,
)
from .spectroscopy.twodcontainer import TwoDResponseContainer as TwoDResponseContainer
from .spectroscopy.twodcontainer import TwoDSpectrumContainer as TwoDSpectrumContainer

#
# Fourier transform Two-Dimensional Spectra
#
from .spectroscopy.twodresponse import TwoDResponse as TwoDResponse
from .spectroscopy.twodspect import TwoDSpectrum as TwoDSpectrum
from .symbolic.cumulant import evaluate_cumulant as evaluate_cumulant
from .utils.logging import (
    init_logging as init_logging,
)
from .utils.logging import (
    log_detail as log_detail,
)
from .utils.logging import (
    log_info as log_info,
)
from .utils.logging import (
    log_quick as log_quick,
)
from .utils.logging import (
    log_report as log_report,
)
from .utils.logging import (
    log_to_file as log_to_file,
)
from .utils.logging import (
    log_urgent as log_urgent,
)
from .utils.logging import (
    loglevels2bool as loglevels2bool,
)
from .utils.logging import (
    printlog as printlog,
)
from .utils.logging import (
    tprint as tprint,
)
from .utils.timing import done_in as done_in
from .utils.timing import finished_in as finished_in
from .utils.timing import timeit as timeit
from .utils.timing import untimeit as untimeit
from .utils.vectors import norm as norm
from .utils.vectors import normalize2 as normalize2
from .wizard.input.input import Input as Input


def exit(msg: str | None = None) -> None:
    """Exit to the level above the script with SystemExit exception"""
    import sys

    if msg is not None:
        printlog("\n(SystemExit) Message: " + msg + "\n", loglevel=0)
    sys.exit()


def stop(msg: str | None = None) -> None:
    """Stop execution and leave to level above"""
    exit("Execution stopped")


def show_plot(block: bool = True) -> None:
    """Shows current plot

    This function is used to avoid explicit import of matplotlib
    """
    import matplotlib.pyplot as plt

    plt.show(block=block)


def savefig(fname: str) -> None:
    """Saves current plot to a file

    This function is used to avoid explicit import of matplotlib
    """
    import matplotlib.pyplot as plt

    plt.savefig(fname)


def assert_version(check: str, vno: str) -> None:
    """Throws an exception if the condition is not satisfied"""
    from packaging import version

    def ext() -> None:
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
        raise Exception("Unknown comparison operator `" + check + "`")
