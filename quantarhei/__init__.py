"""Quantarhei User Level Classes and Objects
=========================================

In Quantarhei, classes are loosely grouped into three categories. First,
there is a group of classes, which represent basic concepts of quantum
mechanics, provide access to important implementations of spectroscopic
simulations and dynamics of open quantum systems, and classes which allow
basic management of the simulation environment and numerical results.
These classes are called **user level classes**, and they are all
accessible in highest namespace level of the Quantarhei package.
If you import Quantarhei like this:

>>> import quantarhei as qr

you can access user level classes through the qr. prefix, e.g.


>>> manager = qr.Manager()
>>> manager.version == qr.__version__
True

The list of user level classes is provided below. The latest and most
up-to-date information can be obtained by viewing the source code of the
root `__init__.py` file of the package, where the imports are grouped
into the same tiers as listed below.


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

Advanced level classes are still intended for relatively frequent use
by the user. However, in order to reduce the *apparent* complexity of
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


Top Level Namespace by Tier
===========================

For historical reasons and for convenience, the top level namespace also
re-exports some classes which belong to the advanced or expert level, as
well as test/mock helpers and deprecated classes. All of them remain
available as ``qr.<name>``, but they are grouped below (and in the source
code) by their intended audience. New scripts should prefer the
user level names.

Tier 1: User level
------------------

Core classes

TimeAxis ......... linear axis of real values representing discrete time
FrequencyAxis .... linear axis of real values representing discrete
frequency axis
ValueAxis ........ linear axis of general real values
DFunction ........ discrete function

Units and basis management

Manager ............ the main behind-the-scenes manager of the package
energy_units ....... energy units manager for use with the "with" construct
frequency_units .... frequency units manager for use with
the "with" construct
length_units ....... length units manager for use with the "with" construct
eigenbasis_of ...... manager of the basis transformations to be used with
the "with" construct
set_current_units .. function to set current units globally
units_state ........ snapshot of the current units settings
convert ............ conversion of values between units
in_current_units ... conversion of values into the current units
EnergyUnit, FrequencyUnit, LengthUnit, TemperatureUnit, TimeUnit
.................... enumerations of supported units

Builders

Mode .......... represents a harmonic vibrational mode of a molecule
HarmonicMode, AnharmonicMode
............... vibrational modes of a VibrationalSystem
Molecule ...... represents a molecule
Aggregate ..... represents an aggregate of molecules
OpenSystem .... common base of open quantum systems
VibrationalSystem
............... represents a system of vibrational modes
PDBFile ....... reader and writer of structures from PDB format
Disorder ...... class managing static disorder of molecular transition
energies

Quantum mechanics

StateVector, DensityMatrix, ReducedDensityMatrix
............... states of the system
Hamiltonian, TransitionDipoleMoment, ProjectionOperator,
BasisReferenceOperator, UnityOperator
............... operators
SystemBathInteraction
............... description of the system-bath coupling
StateVectorPropagator, ReducedDensityMatrixPropagator,
PopulationPropagator, EvolutionSuperOperator
............... propagators and evolution superoperators
StateVectorEvolution, DensityMatrixEvolution,
ReducedDensityMatrixEvolution
............... time evolutions of states

Correlation functions and lineshapes

CorrelationFunction, SpectralDensity, LineshapeFunction,
CorrelationFunctionMatrix, oscillator_scalled_CorrelationFunction

Spectroscopy

AbsSpectrum, AbsSpectrumCalculator, AbsSpectrumContainer
CircDichSpectrum, CircDichSpectrumCalculator, CircDichSpectrumContainer
LinDichSpectrum, LinDichSpectrumCalculator, LinDichSpectrumContainer
FluorSpectrum, FluorSpectrumCalculator, FluorSpectrumContainer
TwoDResponseCalculator, TwoDResponse, TwoDResponseContainer,
TwoDSpectrum, TwoDSpectrumContainer
PumpProbeSpectrum, PumpProbeSpectrumCalculator, PumpProbeSpectrumContainer
LabSetup, LabField

Constants

REAL, COMPLEX ....... numerical types used throughout the package
signal_* ............ types of 2D signals (collected in TWOD_SIGNALS)
part_* .............. parts of complex data (collected in SIGNAL_PARTS,
alias DATA_PARTS)
ptype_* ............. Liouville pathway types (collected in PATHWAY_TYPES,
alias LIOUVILLE_PATHWAY_TYPES)
LOG_* ............... log levels

Saving and loading

save_parcel, load_parcel, check_parcel

Logging, timing and convenience functions

init_logging, printlog, tprint, log_urgent, log_report, log_info,
log_detail, log_quick, log_to_file, loglevels2bool, timeit, untimeit,
finished_in, done_in, norm, normalize2, exit, stop, show_plot, savefig,
assert_version, Input

Exceptions

QuantarheiError, BasisError, BuildError, ConfigurationError,
ImplementationError, UnitsError

Tier 2: Advanced level (re-exported for convenience)
----------------------------------------------------

Liouvillian ......... Liouville superoperator
OQSStateVector, OQSStateVectorPropagator, OQSStateVectorEvolution
..................... open quantum system state vectors
KTHierarchy, KTHierarchyPropagator, QuTip_KTHierarchyPropagator
..................... hierarchical equations of motion (HEOM)
ResponseFunction .... non-linear response function
LiouvillePathwayAnalyzer
..................... analysis of Liouville pathways
DSFeynmanDiagram, R1g_Diagram, R2g_Diagram, R3g_Diagram, R4g_Diagram,
R1f_Diagram, R2f_Diagram, R1g_R_Diagram
..................... double-sided Feynman diagrams (optional dependency)
evaluate_cumulant ... symbolic cumulant evaluation (optional dependency)

Tier 3: Expert level, testing and deprecated (kept for compatibility)
---------------------------------------------------------------------

FunctionStorage, FastFunctionStorage
..................... internal storage of correlation/lineshape functions
Parcel, Saveable, DeserializationWarning
..................... internals of the save/load machinery
TestMolecule, TestAggregate
..................... pre-built test systems
MockAbsSpectrumCalculator, MockTwoDResponseCalculator,
MockPumpProbeSpectrumCalculator
..................... simplified mock calculators
LiouvillePathway, NonLinearResponse
..................... deprecated, use ResponseFunction instead

"""

from __future__ import annotations

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
# The imports below are grouped by API tier (see the module docstring).
# `# isort: split` markers keep Ruff's import sorting within each group.
#
# Module-level constants must be defined before the sub-package imports,
# because several spectroscopy modules import them via ``from .. import``.
#
###############################################################################

import numpy

#
# Managers and exceptions
#
from .core.managers import Manager as Manager
from .exceptions import BasisError as BasisError
from .exceptions import BuildError as BuildError
from .exceptions import ConfigurationError as ConfigurationError
from .exceptions import ImplementationError as ImplementationError
from .exceptions import QuantarheiError as QuantarheiError
from .exceptions import UnitsError as UnitsError

#
# Fix used numerical types
#
REAL: type = numpy.float64
COMPLEX: type = numpy.complex128

#
# Log levels
#
from .utils.logging import (
    LOG_DETAIL,
    LOG_INFO,
    LOG_QUICK,
    LOG_REPORT,
    LOG_URGENT,
)

#
# Non-linear response signals
#
# Each string is defined exactly once, as a standalone constant; the
# collection dicts below are built from these constants.
#
signal_REPH = "rephasing_2D_signal"
signal_NONR = "nonrephasing_2D_signal"
signal_TOTL = "total_2D_signal"
signal_DC = "double_coherence_signal"

TWOD_SIGNALS: dict[str, str] = dict(
    signal_REPH=signal_REPH,
    signal_NONR=signal_NONR,
    signal_TOTL=signal_TOTL,
    signal_DC=signal_DC,
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
    part_REAL=part_REAL,
    part_IMAGINARY=part_IMAGINARY,
    part_COMPLEX=part_COMPLEX,
    part_ABS=part_ABS,
    part_PHASE=part_PHASE,
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
    ptype_R1g=ptype_R1g,
    ptype_R2g=ptype_R2g,
    ptype_R3g=ptype_R3g,
    ptype_R4g=ptype_R4g,
    ptype_R1f=ptype_R1f,
    ptype_R2f=ptype_R2f,
    ptype_R3f=ptype_R3f,
    ptype_R4f=ptype_R4f,
)

LIOUVILLE_PATHWAY_TYPES = PATHWAY_TYPES


###############################################################################
#                       TIER 1: USER LEVEL
###############################################################################

#
# Core classes
#
from .core.dfunction import DFunction as DFunction
from .core.frequency import FrequencyAxis as FrequencyAxis
from .core.time import TimeAxis as TimeAxis
from .core.valueaxis import ValueAxis as ValueAxis

# isort: split

#
# Units and basis management
#
from .core.managers import eigenbasis_of as eigenbasis_of
from .core.managers import energy_units as energy_units
from .core.managers import frequency_units as frequency_units
from .core.managers import length_units as length_units
from .core.managers import set_current_units as set_current_units
from .core.managers import units_state as units_state
from .core.unit_enums import EnergyUnit as EnergyUnit
from .core.unit_enums import FrequencyUnit as FrequencyUnit
from .core.unit_enums import LengthUnit as LengthUnit
from .core.unit_enums import TemperatureUnit as TemperatureUnit
from .core.unit_enums import TimeUnit as TimeUnit
from .core.units import convert as convert
from .core.units import in_current_units as in_current_units

# isort: split

#
# Builders
#
from .builders.aggregates import Aggregate as Aggregate
from .builders.disorder import Disorder as Disorder
from .builders.modes import Mode as Mode
from .builders.molecules import Molecule as Molecule
from .builders.opensystem import OpenSystem as OpenSystem
from .builders.pdb import PDBFile as PDBFile
from .builders.sysmodes import AnharmonicMode as AnharmonicMode
from .builders.sysmodes import HarmonicMode as HarmonicMode
from .builders.vibsystem import VibrationalSystem as VibrationalSystem

# isort: split

#
# Quantum mechanics: states, operators, propagators and evolutions
#
from .qm import BasisReferenceOperator as BasisReferenceOperator
from .qm import DensityMatrix as DensityMatrix
from .qm import DensityMatrixEvolution as DensityMatrixEvolution
from .qm import Hamiltonian as Hamiltonian
from .qm import ProjectionOperator as ProjectionOperator
from .qm import ReducedDensityMatrix as ReducedDensityMatrix
from .qm import ReducedDensityMatrixEvolution as ReducedDensityMatrixEvolution
from .qm import ReducedDensityMatrixPropagator as ReducedDensityMatrixPropagator
from .qm import StateVector as StateVector
from .qm import SystemBathInteraction as SystemBathInteraction
from .qm import TransitionDipoleMoment as TransitionDipoleMoment
from .qm import UnityOperator as UnityOperator
from .qm.liouvillespace.evolutionsuperoperator import (
    EvolutionSuperOperator as EvolutionSuperOperator,
)
from .qm.propagators.poppropagator import PopulationPropagator as PopulationPropagator
from .qm.propagators.statevectorevolution import (
    StateVectorEvolution as StateVectorEvolution,
)
from .qm.propagators.svpropagator import StateVectorPropagator as StateVectorPropagator

# isort: split

#
# Correlation functions and lineshapes
#
from .qm.corfunctions import CorrelationFunction as CorrelationFunction
from .qm.corfunctions import CorrelationFunctionMatrix as CorrelationFunctionMatrix
from .qm.corfunctions import LineshapeFunction as LineshapeFunction
from .qm.corfunctions import SpectralDensity as SpectralDensity
from .qm.corfunctions.correlationfunctions import (
    oscillator_scalled_CorrelationFunction as oscillator_scalled_CorrelationFunction,
)

# isort: split

#
# Spectroscopy: linear absorption, circular and linear dichroism, fluorescence
#
from .spectroscopy.abs2 import AbsSpectrum as AbsSpectrum
from .spectroscopy.abscalculator import AbsSpectrumCalculator as AbsSpectrumCalculator
from .spectroscopy.abscontainer import AbsSpectrumContainer as AbsSpectrumContainer
from .spectroscopy.circular_dichroism import CircDichSpectrum as CircDichSpectrum
from .spectroscopy.circular_dichroism import (
    CircDichSpectrumCalculator as CircDichSpectrumCalculator,
)
from .spectroscopy.circular_dichroism import (
    CircDichSpectrumContainer as CircDichSpectrumContainer,
)
from .spectroscopy.fluorescence import FluorSpectrum as FluorSpectrum
from .spectroscopy.fluorescence import (
    FluorSpectrumCalculator as FluorSpectrumCalculator,
)
from .spectroscopy.fluorescence import (
    FluorSpectrumContainer as FluorSpectrumContainer,
)
from .spectroscopy.linear_dichroism import LinDichSpectrum as LinDichSpectrum
from .spectroscopy.linear_dichroism import (
    LinDichSpectrumCalculator as LinDichSpectrumCalculator,
)
from .spectroscopy.linear_dichroism import (
    LinDichSpectrumContainer as LinDichSpectrumContainer,
)

# isort: split

#
# Spectroscopy: two-dimensional and pump-probe spectra, laboratory setup
#
from .spectroscopy.labsetup import LabField as LabField
from .spectroscopy.labsetup import LabSetup as LabSetup
from .spectroscopy.pumpprobe import PumpProbeSpectrum as PumpProbeSpectrum
from .spectroscopy.pumpprobe import (
    PumpProbeSpectrumCalculator as PumpProbeSpectrumCalculator,
)
from .spectroscopy.pumpprobe import (
    PumpProbeSpectrumContainer as PumpProbeSpectrumContainer,
)
from .spectroscopy.twodcalculator import (
    TwoDResponseCalculator as TwoDResponseCalculator,
)
from .spectroscopy.twodcontainer import TwoDResponseContainer as TwoDResponseContainer
from .spectroscopy.twodcontainer import TwoDSpectrumContainer as TwoDSpectrumContainer
from .spectroscopy.twodresponse import TwoDResponse as TwoDResponse
from .spectroscopy.twodspect import TwoDSpectrum as TwoDSpectrum

# isort: split

#
# Saving and loading
#
from .core.parcel import check_parcel as check_parcel
from .core.parcel import load_parcel as load_parcel
from .core.parcel import save_parcel as save_parcel

# isort: split

#
# Logging, timing, vectors and input
#
from .utils.logging import init_logging as init_logging
from .utils.logging import log_detail as log_detail
from .utils.logging import log_info as log_info
from .utils.logging import log_quick as log_quick
from .utils.logging import log_report as log_report
from .utils.logging import log_to_file as log_to_file
from .utils.logging import log_urgent as log_urgent
from .utils.logging import loglevels2bool as loglevels2bool
from .utils.logging import printlog as printlog
from .utils.logging import tprint as tprint
from .utils.timing import done_in as done_in
from .utils.timing import finished_in as finished_in
from .utils.timing import timeit as timeit
from .utils.timing import untimeit as untimeit
from .utils.vectors import norm as norm
from .utils.vectors import normalize2 as normalize2
from .wizard.input.input import Input as Input

###############################################################################
#          TIER 2: ADVANCED LEVEL (re-exported for convenience)
###############################################################################

# isort: split

#
# Liouville space and open quantum system state vectors
#
from .qm import Liouvillian as Liouvillian
from .qm import OQSStateVector as OQSStateVector
from .qm import OQSStateVectorEvolution as OQSStateVectorEvolution
from .qm import OQSStateVectorPropagator as OQSStateVectorPropagator

# isort: split

#
# Hierarchical equations of motion (HEOM)
#
from .qm.liouvillespace.heom import KTHierarchy as KTHierarchy
from .qm.liouvillespace.heom import KTHierarchyPropagator as KTHierarchyPropagator
from .qm.liouvillespace.heom import (
    QuTip_KTHierarchyPropagator as QuTip_KTHierarchyPropagator,
)

# isort: split

#
# Response functions and Liouville pathway analysis
#
from .spectroscopy.pathwayanalyzer import (
    LiouvillePathwayAnalyzer as LiouvillePathwayAnalyzer,
)
from .spectroscopy.responses import ResponseFunction as ResponseFunction

#
# Double-sided Feynman diagrams (optional dependency)
#
try:
    from .spectroscopy.dsfeynman import DSFeynmanDiagram as DSFeynmanDiagram
    from .spectroscopy.dsfeynman import R1f_Diagram as R1f_Diagram
    from .spectroscopy.dsfeynman import R1g_Diagram as R1g_Diagram
    from .spectroscopy.dsfeynman import R1g_R_Diagram as R1g_R_Diagram
    from .spectroscopy.dsfeynman import R2f_Diagram as R2f_Diagram
    from .spectroscopy.dsfeynman import R2g_Diagram as R2g_Diagram
    from .spectroscopy.dsfeynman import R3g_Diagram as R3g_Diagram
    from .spectroscopy.dsfeynman import R4g_Diagram as R4g_Diagram
except ImportError:
    pass

#
# Symbolic cumulant evaluation (optional dependency)
#
try:
    from .symbolic.cumulant import evaluate_cumulant as evaluate_cumulant
except ImportError:
    pass

###############################################################################
#   TIER 3: EXPERT LEVEL, TESTING AND DEPRECATED (kept for compatibility)
###############################################################################

#
# Internal storage of correlation and lineshape functions
#
from .qm.corfunctions import FastFunctionStorage as FastFunctionStorage
from .qm.corfunctions import FunctionStorage as FunctionStorage

# isort: split

#
# Internals of the save/load machinery
#
from .core.parcel import DeserializationWarning as DeserializationWarning
from .core.parcel import Parcel as Parcel
from .core.saveable import Saveable as Saveable

# isort: split

#
# Pre-built test systems and mock calculators
#
from .builders.aggregate_test import TestAggregate as TestAggregate
from .builders.molecule_test import TestMolecule as TestMolecule
from .spectroscopy.mockabscalculator import (
    MockAbsSpectrumCalculator as MockAbsSpectrumCalculator,
)
from .spectroscopy.mocktwodcalculator import (
    MockTwoDResponseCalculator as MockTwoDResponseCalculator,
)
from .spectroscopy.pumpprobe import (
    MockPumpProbeSpectrumCalculator as MockPumpProbeSpectrumCalculator,
)

# isort: split

#
# Deprecated classes (use ResponseFunction instead)
#
from .spectroscopy.responses import LiouvillePathway as LiouvillePathway
from .spectroscopy.responses import NonLinearResponse as NonLinearResponse


def exit(msg: str | None = None) -> None:
    """Exit to the level above the script by raising :exc:`SystemExit`.

    Parameters
    ----------
    msg : str or None, optional
        Optional message to log before exiting. Default is ``None``.
    """
    import sys

    if msg is not None:
        printlog("\n(SystemExit) Message: " + msg + "\n", loglevel=0)
    sys.exit()


def stop(msg: str | None = None) -> None:
    """Stop execution and exit to the level above.

    Parameters
    ----------
    msg : str or None, optional
        Ignored; present for API compatibility. Default is ``None``.
    """
    exit("Execution stopped")


def show_plot(block: bool = True) -> None:
    """Show the current matplotlib plot.

    Convenience wrapper around ``matplotlib.pyplot.show`` that avoids an
    explicit import of ``matplotlib``.

    Parameters
    ----------
    block : bool, optional
        If ``True``, block until the plot window is closed. Default is
        ``True``.
    """
    import matplotlib.pyplot as plt

    plt.show(block=block)


def savefig(fname: str) -> None:
    """Save the current matplotlib plot to a file.

    Convenience wrapper around ``matplotlib.pyplot.savefig`` that avoids an
    explicit import of ``matplotlib``.

    Parameters
    ----------
    fname : str
        Output file path (format inferred from the extension).
    """
    import matplotlib.pyplot as plt

    plt.savefig(fname)


def assert_version(check: str, vno: str) -> None:
    """Assert that the installed Quantarhei version satisfies a version constraint.

    Parameters
    ----------
    check : str
        Comparison operator string: ``">="``, ``"=="``, or ``"<="``.
    vno : str
        Version string to compare against (e.g. ``"0.0.63"``).

    Raises
    ------
    SystemExit
        If the version constraint is not satisfied.
    """
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
        if not (version.parse(Manager().version) > version.parse(vno)):
            ext()

    elif check == "<":
        if not (version.parse(Manager().version) < version.parse(vno)):
            ext()

    else:
        raise QuantarheiError("Unknown comparison operator `" + check + "`")


#
#  __all__ attribute to define a public API
#
#  Kept as a single sorted list (Ruff RUF022); see the module docstring and
#  the import sections above for the grouping of names by API tier.
#
__all__ = [
    "COMPLEX",
    "DATA_PARTS",
    "LIOUVILLE_PATHWAY_TYPES",
    "LOG_DETAIL",
    "LOG_INFO",
    "LOG_QUICK",
    "LOG_REPORT",
    "LOG_URGENT",
    "PATHWAY_TYPES",
    "REAL",
    "SIGNAL_PARTS",
    "TWOD_SIGNALS",
    "AbsSpectrum",
    "AbsSpectrumCalculator",
    "AbsSpectrumContainer",
    "Aggregate",
    "AnharmonicMode",
    "BasisError",
    "BasisReferenceOperator",
    "BuildError",
    "CircDichSpectrum",
    "CircDichSpectrumCalculator",
    "CircDichSpectrumContainer",
    "ConfigurationError",
    "CorrelationFunction",
    "CorrelationFunctionMatrix",
    "DFunction",
    "DSFeynmanDiagram",
    "DensityMatrix",
    "DensityMatrixEvolution",
    "DeserializationWarning",
    "Disorder",
    "EvolutionSuperOperator",
    "FastFunctionStorage",
    "FluorSpectrum",
    "FluorSpectrumCalculator",
    "FluorSpectrumContainer",
    "FrequencyAxis",
    "FunctionStorage",
    "Hamiltonian",
    "HarmonicMode",
    "ImplementationError",
    "Input",
    "KTHierarchy",
    "KTHierarchyPropagator",
    "LabField",
    "LabSetup",
    "LinDichSpectrum",
    "LinDichSpectrumCalculator",
    "LinDichSpectrumContainer",
    "LineshapeFunction",
    "LiouvillePathway",
    "LiouvillePathwayAnalyzer",
    "Liouvillian",
    "Manager",
    "MockAbsSpectrumCalculator",
    "MockPumpProbeSpectrumCalculator",
    "MockTwoDResponseCalculator",
    "Mode",
    "Molecule",
    "NonLinearResponse",
    "OQSStateVector",
    "OQSStateVectorEvolution",
    "OQSStateVectorPropagator",
    "OpenSystem",
    "PDBFile",
    "Parcel",
    "PopulationPropagator",
    "ProjectionOperator",
    "PumpProbeSpectrum",
    "PumpProbeSpectrumCalculator",
    "PumpProbeSpectrumContainer",
    "QuTip_KTHierarchyPropagator",
    "QuantarheiError",
    "R1f_Diagram",
    "R1g_Diagram",
    "R1g_R_Diagram",
    "R2f_Diagram",
    "R2g_Diagram",
    "R3g_Diagram",
    "R4g_Diagram",
    "ReducedDensityMatrix",
    "ReducedDensityMatrixEvolution",
    "ReducedDensityMatrixPropagator",
    "ResponseFunction",
    "Saveable",
    "SpectralDensity",
    "StateVector",
    "StateVectorEvolution",
    "StateVectorPropagator",
    "SystemBathInteraction",
    "TestAggregate",
    "TestMolecule",
    "TimeAxis",
    "TransitionDipoleMoment",
    "TwoDResponse",
    "TwoDResponseCalculator",
    "TwoDResponseContainer",
    "TwoDSpectrum",
    "TwoDSpectrumContainer",
    "UnitsError",
    "UnityOperator",
    "ValueAxis",
    "VibrationalSystem",
    "assert_version",
    "check_parcel",
    "convert",
    "done_in",
    "eigenbasis_of",
    "energy_units",
    "evaluate_cumulant",
    "exit",
    "finished_in",
    "frequency_units",
    "in_current_units",
    "init_logging",
    "length_units",
    "load_parcel",
    "log_detail",
    "log_info",
    "log_quick",
    "log_report",
    "log_to_file",
    "log_urgent",
    "loglevels2bool",
    "norm",
    "normalize2",
    "oscillator_scalled_CorrelationFunction",
    "part_ABS",
    "part_COMPLEX",
    "part_IMAGINARY",
    "part_PHASE",
    "part_REAL",
    "printlog",
    "ptype_R1f",
    "ptype_R1g",
    "ptype_R2f",
    "ptype_R2g",
    "ptype_R3f",
    "ptype_R3g",
    "ptype_R4f",
    "ptype_R4g",
    "save_parcel",
    "savefig",
    "set_current_units",
    "show_plot",
    "signal_DC",
    "signal_NONR",
    "signal_REPH",
    "signal_TOTL",
    "stop",
    "timeit",
    "tprint",
    "untimeit",
]
