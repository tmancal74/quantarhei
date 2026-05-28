# Quantarhei Class Architecture

## End-to-End Data Flow

How the subsystems connect — the typical pipeline from system construction to observable:

```mermaid
flowchart LR
    subgraph BUILD ["System Construction"]
        Mol[Molecule] --> Agg[Aggregate]
        Mode2[Mode] --> Mol
        Agg -->|build| Ham[Hamiltonian]
        Agg -->|build| TDM[TransitionDipoleMoment]
    end

    subgraph BATH ["Bath & Coupling"]
        CF[CorrelationFunction] --> SBI[SystemBathInteraction]
        SD[SpectralDensity] -->|to CF| CF
        SBI --> RT[RelaxationTensor]
        Ham --> RT
    end

    subgraph DYNAMICS ["Dynamics"]
        RT --> Prop[Propagator]
        Ham --> Prop
        Prop -->|propagate| Evol[Evolution]
        Ham --> ESO[EvolutionSuperOperator]
        RT --> ESO
    end

    subgraph SPECTRA ["Spectroscopy"]
        Agg --> Calc[SpectrumCalculator]
        ESO --> Calc
        Evol --> Calc
        Calc --> Spec[Spectrum]
    end
```

## Typical Usage Pipeline

```mermaid
flowchart TD
    A["1. Create Molecules\n(set energies, dipoles, modes)"] --> B
    B["2. Create Aggregate\n(add molecules, set couplings)"] --> C
    C["3. agg.build()\n(constructs Hamiltonian, state space)"] --> D
    D{"What do you need?"}

    D -->|Dynamics| E["4a. Create CorrelationFunction\n+ SystemBathInteraction"]
    E --> F["5a. get_RelaxationTensor()\nor get_ReducedDensityMatrixPropagator()"]
    F --> G["6a. propagator.propagate(rho0)\n→ DensityMatrixEvolution"]

    D -->|Linear Spectrum| H["4b. AbsSpectrumCalculator(time, agg)"]
    H --> I["5b. calculator.calculate()\n→ AbsSpectrum"]

    D -->|2D Spectrum| J["4c. MockTwoDResponseCalculator(time, agg)"]
    J --> K["5c. calculator.calculate()\n→ TwoDSpectrumContainer"]

    D -->|Exciton Analysis| L["4d. agg.diagonalize()"]
    L --> M["5d. agg.exciton_report()\nagg.get_expansion_squares()"]
```

## Public API (`from quantarhei import ...`)

```mermaid
flowchart TD
    subgraph builders ["Builders"]
        Molecule
        Mode
        Aggregate
        HarmonicMode
        AnharmonicMode
        VibrationalSystem
        OpenSystem
        Disorder
        PDBFile
    end

    subgraph qm_objects ["QM Objects"]
        Hamiltonian
        DensityMatrix
        ReducedDensityMatrix
        TransitionDipoleMoment
        StateVector
        OQSStateVector
        UnityOperator
        ProjectionOperator
        BasisReferenceOperator
    end

    subgraph bath ["Bath & Relaxation"]
        CorrelationFunction
        CorrelationFunctionMatrix
        SpectralDensity
        SystemBathInteraction
        LineshapeFunction
        EvolutionSuperOperator
        Liouvillian
        KTHierarchy
        KTHierarchyPropagator
    end

    subgraph propagators ["Propagators & Evolution"]
        ReducedDensityMatrixPropagator
        ReducedDensityMatrixEvolution
        DensityMatrixEvolution
        StateVectorPropagator
        StateVectorEvolution
        OQSStateVectorPropagator
        OQSStateVectorEvolution
        PopulationPropagator
    end

    subgraph spectra ["Spectroscopy"]
        AbsSpectrum
        AbsSpectrumCalculator
        AbsSpectrumContainer
        MockAbsSpectrumCalculator
        FluorSpectrum
        FluorSpectrumCalculator
        FluorSpectrumContainer
        CircDichSpectrum
        CircDichSpectrumCalculator
        CircDichSpectrumContainer
        LinDichSpectrum
        LinDichSpectrumCalculator
        LinDichSpectrumContainer
        TwoDResponseCalculator
        MockTwoDResponseCalculator
        TwoDResponse
        TwoDSpectrum
        TwoDResponseContainer
        TwoDSpectrumContainer
        PumpProbeSpectrum
        PumpProbeSpectrumCalculator
        PumpProbeSpectrumContainer
        MockPumpProbeSpectrumCalculator
        LiouvillePathway
        LiouvillePathwayAnalyzer
        NonLinearResponse
        ResponseFunction
        LabSetup
        LabField
    end

    subgraph core ["Core"]
        TimeAxis
        FrequencyAxis
        ValueAxis
        DFunction
        Saveable
        energy_units
        frequency_units
        length_units
        eigenbasis_of
    end

    subgraph parallel ["Parallelization"]
        distributed_configuration
        parallel_function
        start_parallel_region
        close_parallel_region
        block_distributed_range
        collect_block_distributed_data
    end

    subgraph wizard_mod ["Wizard"]
        Input
    end
```

## Manager & Context Managers (Global State)

The `Manager` singleton coordinates units and basis across all objects:

```mermaid
flowchart TD
    M[Manager Singleton] -->|tracks| Units["Current energy/frequency/length units"]
    M -->|tracks| Basis["Current basis (site vs exciton)"]

    EU["with energy_units('1/cm'):"] -->|sets/restores| Units
    FU["with frequency_units('THz'):"] -->|sets/restores| Units
    LU["with length_units('nm'):"] -->|sets/restores| Units
    EB["with eigenbasis_of(H):"] -->|sets/restores| Basis

    Units -->|affects| Conv["All EnergyUnitsManaged objects\nautomatically convert on get/set"]
    Basis -->|affects| Trans["All BasisManaged objects\nautomatically transform on get/set"]

    subgraph managed ["Managed Objects"]
        Molecule2[Molecule]
        Aggregate2[Aggregate]
        CF2[CorrelationFunction]
        SD2[SpectralDensity]
        Spec2[Spectra]
    end

    Conv --> managed
    Trans --> managed
```

## build() Lifecycle

What happens when you call `agg.build()`:

```mermaid
flowchart TD
    B["agg.build(mult=1)"] --> S1["Count electronic states\n(ground + N single-exciton + ...)"]
    S1 --> S2["Enumerate vibrational sub-states\nper electronic state (from Mode.nmax)"]
    S2 --> S3["Build full Hamiltonian\n(electronic energies + vib. energies\n+ couplings × FC factors)"]
    S3 --> S4["Build TransitionDipoleMoment\n(dipoles × FC overlaps)"]
    S4 --> S5["Set up state indexing:\nNb[], vibindices[], elinds[], which_band[]"]
    S5 --> S6["Compute Franck-Condon factors\n(FCf matrix)"]
    S6 --> DONE["Aggregate is ready\n(can now diagonalize, propagate, or calculate spectra)"]
```

## Aggregate Inheritance Chain

The core of quantarhei is the `Aggregate` class, built via a deep inheritance chain where each layer adds specific capabilities:

```mermaid
classDiagram
    direction BT

    OpenSystem <|-- AggregateBase : builders/aggregate_base.py
    AggregateBase <|-- AggregateSpectroscopy : builders/aggregate_spectroscopy.py
    AggregateSpectroscopy <|-- AggregateExcitonAnalysis : builders/aggregate_excitonanalysis.py
    AggregateExcitonAnalysis <|-- AggregatePureDephasing : builders/aggregate_pdeph.py
    AggregatePureDephasing <|-- Aggregate : builders/aggregates.py

    class OpenSystem {
        +get_Hamiltonian()
        +get_SystemBathInteraction()
        +get_RelaxationTensor()
        +get_ReducedDensityMatrixPropagator()
        +diagonalize()
    }

    class AggregateBase {
        +build()
        +diagonalize()
        +get_DensityMatrix()
        +convert_to_ground_vibbasis()
        +trace_over_vibrations()
        +set_resonance_coupling()
        +allstates()
        +elsignatures()
        +fc_factor()
        +coupling()
    }

    class AggregateSpectroscopy {
        +liouville_pathways_3()
        +liouville_pathways_3T()
    }

    class AggregateExcitonAnalysis {
        +exciton_report()
        +get_expansion_squares()
        +get_intersite_mixing()
        +get_transition_dipole()
        +get_state_energy()
    }

    class AggregatePureDephasing {
        +get_PureDephasing()
    }

    class Aggregate {
    }
```

## Molecular Building Blocks

```mermaid
classDiagram
    direction BT

    OpenSystem <|-- Molecule
    OpenSystem <|-- Mode
    OpenSystem <|-- HarmonicMode
    OpenSystem <|-- VibrationalSystem
    HarmonicMode <|-- AnharmonicMode

    Molecule --> Mode : has modes
    Mode --> SubMode : has submodes
    Molecule --> ElectronicState : has states
    Aggregate --> Disorder : disorder model

    class Molecule {
        +elenergies
        +dmoments
        +modes
        +position
        +set_dipole()
        +set_energy()
        +add_Mode()
        +get_Hamiltonian()
        +get_TransitionDipoleMoment()
    }

    class Mode {
        +omega
        +set_HR()
        +set_nmax()
        +get_SubMode()
    }

    class SubMode {
        +omega
        +nmax
        +shift
    }

    class ElectronicState {
        +energy
        +transition_dipole
    }

    class Disorder {
        +distribution
        +disorder_type
    }

    class PDBFile {
        +read()
        +get_Molecules()
    }
```

## Quantum Mechanics: Operators & States

```mermaid
classDiagram
    direction BT

    Operator <|-- SelfAdjointOperator
    Operator <|-- ProjectionOperator
    SelfAdjointOperator <|-- Hamiltonian
    SelfAdjointOperator <|-- DensityMatrix
    SelfAdjointOperator <|-- TransitionDipoleMoment
    SelfAdjointOperator <|-- UnityOperator
    SelfAdjointOperator <|-- BasisReferenceOperator
    DensityMatrix <|-- ReducedDensityMatrix

    class Operator {
        +data : ndarray
        +dim : int
        +transform()
    }

    class Hamiltonian {
        +diagonalize()
        +data : (dim x dim)
    }

    class DensityMatrix {
        +data : (dim x dim)
    }

    class ReducedDensityMatrix {
        +data : (dim x dim)
    }

    class TransitionDipoleMoment {
        +data : (dim x dim x 3)
        +dipole_strength()
        +transform()
    }

    class StateVector {
        +data : (dim,)
    }

    class OQSStateVector {
        +data : (dim,)
    }
```

## Quantum Mechanics: Propagators & Evolution

```mermaid
classDiagram
    direction BT

    ReducedDensityMatrixPropagator --> ReducedDensityMatrixEvolution : propagate()
    ReducedDensityMatrixPropagator --> DensityMatrixEvolution : propagate()
    StateVectorPropagator --> StateVectorEvolution : propagate()
    PopulationPropagator --> PopulationEvolution : propagate()
    OQSStateVectorPropagator --> OQSStateVectorEvolution : propagate()

    class ReducedDensityMatrixPropagator {
        +Hamiltonian
        +RelaxationTensor
        +TimeAxis
        +propagate(rho0) → Evolution
    }

    class DensityMatrixEvolution {
        +TimeAxis
        +data : (Nt x dim x dim)
    }

    class ReducedDensityMatrixEvolution {
        +TimeAxis
        +data : (Nt x dim x dim)
    }

    class StateVectorPropagator {
        +Hamiltonian
        +propagate(psi0)
    }

    class StateVectorEvolution {
        +TimeAxis
        +data : (Nt x dim)
    }

    class PopulationPropagator {
        +RateMatrix
        +propagate(pop0)
    }

    class PopulationEvolution {
        +TimeAxis
        +data : (Nt x dim)
    }

    class OQSStateVectorPropagator {
        +propagate(psi0)
    }

    class OQSStateVectorEvolution {
        +TimeAxis
        +data
    }
```

## Quantum Mechanics: Relaxation Tensors

```mermaid
classDiagram
    direction BT

    SuperOperator <|-- RelaxationTensor
    RelaxationTensor <|-- RedfieldRelaxationTensor
    RelaxationTensor <|-- FoersterRelaxationTensor
    RelaxationTensor <|-- ModRedfieldRelaxationTensor
    RelaxationTensor <|-- TDModRedfieldRelaxationTensor
    RedfieldRelaxationTensor <|-- TDRedfieldRelaxationTensor
    RedfieldRelaxationTensor <|-- LindbladForm
    FoersterRelaxationTensor <|-- TDFoersterRelaxationTensor
    FoersterRelaxationTensor <|-- RedfieldFoersterRelaxationTensor
    RedfieldFoersterRelaxationTensor <|-- TDRedfieldFoersterRelaxationTensor
    TDFoersterRelaxationTensor <|-- NEFoersterRelaxationTensor
    LindbladForm <|-- ElectronicLindbladForm
    LindbladForm <|-- VibrationalDecayLindbladForm
    SuperOperator <|-- EvolutionSuperOperator
    SuperOperator <|-- SOpUnity
    SuperOperator <|-- Liouvillian

    class RelaxationTensor {
        +Hamiltonian
        +SystemBathInteraction
        +data : (dim x dim x dim x dim)
    }

    class RedfieldRelaxationTensor {
        Secular approximation
        Weak system-bath coupling
    }

    class FoersterRelaxationTensor {
        Strong coupling limit
        Localized basis
    }

    class ModRedfieldRelaxationTensor {
        Modified Redfield
        Includes reorganization
    }

    class LindbladForm {
        Lindblad master equation
        Completely positive
    }

    class RedfieldFoersterRelaxationTensor {
        Combined theory
        Redfield within clusters
        Foerster between clusters
    }

    class NEFoersterRelaxationTensor {
        Non-equilibrium Foerster
        Time-dependent donor state
    }

    class EvolutionSuperOperator {
        +at(t) → SuperOperator
        +data : (Nt x dim² x dim²)
    }

    class Liouvillian {
        +data
    }
```

## Quantum Mechanics: Rate Matrices

```mermaid
classDiagram
    direction BT

    RateMatrix <|-- RedfieldRateMatrix
    RateMatrix <|-- FoersterRateMatrix
    RateMatrix <|-- ModifiedRedfieldRateMatrix
    RateMatrix <|-- TDRedfieldRateMatrix

    class RateMatrix {
        +data : (dim x dim)
        Population transfer rates
        Used by PopulationPropagator
    }
```

## Quantum Mechanics: HEOM & Advanced Propagators

```mermaid
classDiagram
    direction BT

    KTHierarchy --> KTHierarchyPropagator : used by
    KTHierarchyPropagator <|-- QuTip_KTHierarchyPropagator

    class KTHierarchy {
        Hierarchical Equations of Motion
        Numerically exact for Drude-Lorentz
        +depth
        +SystemBathInteraction
    }

    class KTHierarchyPropagator {
        +propagate()
    }

    class QuTip_KTHierarchyPropagator {
        QuTiP backend
        +propagate()
    }

    class IntegrodiffPropagator {
        Integro-differential equation
        Non-Markovian dynamics
        +propagate()
    }
```

## Correlation Functions & Spectral Densities

```mermaid
classDiagram
    direction BT

    DFunction <|-- CorrelationFunction
    DFunction <|-- SpectralDensity
    DFunction <|-- FTCorrelationFunction
    FTCorrelationFunction <|-- EvenFTCorrelationFunction
    FTCorrelationFunction <|-- OddFTCorrelationFunction
    DFunction <|-- LineshapeFunction

    CorrelationFunction --> CorrelationFunctionMatrix : collected in
    CorrelationFunctionMatrix --> SystemBathInteraction : used by

    class CorrelationFunction {
        +TimeAxis
        +params (type, reorg, cortime, T)
        +temperature
        +get_SpectralDensity()
    }

    class SpectralDensity {
        +FrequencyAxis
        +params
        +reorganization_energy()
        +get_CorrelationFunction()
    }

    class CorrelationFunctionMatrix {
        +data : N x N matrix of CFs
        +nof : number of functions
    }

    class SystemBathInteraction {
        +CC : CorrelationFunctionMatrix
        +system_operators
        +KK : coupling operators
    }

    class LineshapeFunction {
        +data : g(t) lineshape
    }

    class PureDephasing {
        +dephasing_rates
    }

    class ElectronicPureDephasing {
        +dephasing_rates
    }

    class FunctionStorage {
        Cache for lineshape functions
    }
```

## Spectroscopy: Linear Spectra

```mermaid
classDiagram
    direction BT

    DFunction <|-- AbsSpectrumBase
    DFunction <|-- FluorSpectrumBase
    DFunction <|-- CircDichSpectrumBase
    DFunction <|-- LinDichSpectrumBase
    AbsSpectrumBase <|-- AbsSpectrum
    AbsSpectrumBase <|-- LinSpectrum
    FluorSpectrumBase <|-- FluorSpectrum
    CircDichSpectrumBase <|-- CircDichSpectrum
    LinDichSpectrumBase <|-- LinDichSpectrum

    AbsSpectrumCalculator --> AbsSpectrum : produces
    AbsSpectrumCalculator --> AbsSpectrumContainer : collects into
    MockAbsSpectrumCalculator --> AbsSpectrum : produces
    LinSpectrumCalculator --> LinSpectrum : produces
    FluorSpectrumCalculator --> FluorSpectrum : produces
    FluorSpectrumCalculator --> FluorSpectrumContainer : collects into
    CircDichSpectrumCalculator --> CircDichSpectrum : produces
    CircDichSpectrumCalculator --> CircDichSpectrumContainer : collects into
    LinDichSpectrumCalculator --> LinDichSpectrum : produces
    LinDichSpectrumCalculator --> LinDichSpectrumContainer : collects into

    class AbsSpectrumCalculator {
        +TimeAxis
        +system (Aggregate)
        +calculate()
    }

    class MockAbsSpectrumCalculator {
        Liouville pathway based
        +calculate()
    }

    class AbsSpectrum {
        +axis : FrequencyAxis
        +data
        +gaussian_fit()
    }

    class AbsSpectrumContainer {
        +spectra
        +get_spectrum()
        +arange()
    }
```

## Spectroscopy: Nonlinear (2D & Pump-Probe)

```mermaid
classDiagram
    direction BT

    TwoDResponseCalculator <|-- MockTwoDResponseCalculator
    MockTwoDResponseCalculator <|-- MockPumpProbeSpectrumCalculator

    TwoDResponseContainer <|-- TwoDSpectrumContainer
    TwoDSpectrumContainer <|-- PumpProbeSpectrumContainer

    MockTwoDResponseCalculator --> TwoDResponse : produces
    MockTwoDResponseCalculator --> TwoDResponseContainer : collects into
    TwoDResponse --> TwoDSpectrum : FFT
    TwoDSpectrumContainer --> TwoDSpectrum : contains
    PumpProbeSpectrumContainer --> PumpProbeSpectrum : contains

    class TwoDResponseCalculator {
        +TimeAxis (t1, t2, t3)
        +system (Aggregate)
        +calculate()
    }

    class MockTwoDResponseCalculator {
        Liouville pathway based
        Uses EvolutionSuperOperator for t2
        +calculate_pathways()
        +calculate()
    }

    class TwoDResponse {
        +t1axis, t3axis
        +data : complex (Nt1 x Nt3)
        Rephasing + Non-rephasing
    }

    class TwoDSpectrum {
        +xaxis, yaxis : FrequencyAxis
        +data : 2D array
        +get_PumpProbeSpectrum()
    }

    class TwoDSpectrumContainer {
        +t2axis
        +spectra (one per t2)
        +get_spectrum(t2)
    }

    class PumpProbeSpectrum {
        +axis : FrequencyAxis
        +data
    }

    class PumpProbeSpectrumContainer {
        +t2axis
        +spectra
    }
```

## Spectroscopy: Liouville Pathways & Diagrams

```mermaid
classDiagram
    direction BT

    LiouvillePathway --> LiouvillePathwayAnalyzer : analyzed by
    DSFeynmanDiagram <|-- R1g_Diagram
    DSFeynmanDiagram <|-- R2g_Diagram
    DSFeynmanDiagram <|-- R3g_Diagram
    DSFeynmanDiagram <|-- R4g_Diagram
    DSFeynmanDiagram <|-- R1f_Diagram
    DSFeynmanDiagram <|-- R2f_Diagram

    class LiouvillePathway {
        +pathway_type (R1g, R2g, R3g, R4g, R1fs, R2fs)
        +frequency
        +dipole_factor
        +states (initial, coherence, population, final)
    }

    class LiouvillePathwayAnalyzer {
        +select_frequency_window()
        +select_pathways()
        +get_type()
    }

    class NonLinearResponse {
        +pathways
        +calculate()
    }

    class ResponseFunction {
        +TimeAxis
        +data
    }

    class DSFeynmanDiagram {
        Double-sided Feynman diagram
        +draw()
    }

    class LabSetup {
        +pulse_polarizations
        +detection_axis
        +set_polarizations()
    }

    class LabField {
        +polarization
        +frequency
    }
```

## Core Infrastructure

```mermaid
classDiagram
    direction BT

    Managed <|-- UnitsManaged
    Managed <|-- EnergyUnitsManaged
    Managed <|-- LengthUnitsManaged
    Managed <|-- BasisManaged

    Saveable <|-- DFunction
    Saveable <|-- ValueAxis
    ValueAxis <|-- TimeAxis
    ValueAxis <|-- FrequencyAxis
    DFunction <|-- DFunction2

    Manager --> Managed : configures

    class Manager {
        Singleton — global state for all Managed objects
        +current_units : dict
        +current_basis : stack
        +get_current_units()
    }

    class Managed {
        Receives context from Manager
    }

    class UnitsManaged {
        +convert_energy_2_current_u()
        +convert_energy_2_internal_u()
    }

    class EnergyUnitsManaged {
        Energy-specific conversion
    }

    class BasisManaged {
        Automatic basis transformation
    }

    class Saveable {
        +save(filename)
        +load(filename)
    }

    class DFunction {
        +axis : ValueAxis
        +data : ndarray
        +at(x) interpolated value
        +plot()
        +save() / load_data()
    }

    class DFunction2 {
        +xaxis, yaxis
        +data : 2D ndarray
    }

    class TimeAxis {
        +start, step, length
        +data : ndarray
        +atype (complete, upper-half)
    }

    class FrequencyAxis {
        +start, step, length
        +data : ndarray
        +atype
    }
```

## Core: Context Managers

```mermaid
classDiagram
    direction BT

    units_context_manager <|-- energy_units
    energy_units <|-- frequency_units
    units_context_manager <|-- length_units
    basis_context_manager <|-- eigenbasis_of

    class energy_units {
        with energy_units("1/cm"):
        Supported: 1/cm, eV, meV, THz, nm, int
    }

    class frequency_units {
        with frequency_units("THz"):
        Converts frequency axis units
    }

    class length_units {
        with length_units("nm"):
        Converts length/distance units
    }

    class eigenbasis_of {
        with eigenbasis_of(H):
        All operators transform to H eigenbasis
        Restores original basis on exit
    }
```

## Core: Parallelization

```mermaid
classDiagram
    direction BT

    class distributed_configuration {
        Context manager for parallel setup
        +start_parallel_region()
        +close_parallel_region()
    }

    class parallel_function {
        Decorator for parallelizable functions
    }

    class block_distributed_range {
        Splits range across MPI ranks
    }

    class block_distributed_list {
        Splits list across MPI ranks
    }

    class collect_block_distributed_data {
        Gathers results from all ranks
    }

    class Parcel {
        Serializable data package
        for inter-process communication
    }
```

## Models & Databases

```mermaid
classDiagram
    direction BT

    MolecularModel <|-- BacterioChlorophyll
    MolecularModel <|-- BacterioPheophytin
    MolecularModel <|-- ChlorophyllA
    MolecularModel <|-- ChlorophyllB

    SpectralDensityDB --> SpectralDensityDatabaseEntry : contains
    CorrelationFunctionDB --> CorrelationFunctionDatabaseEntry : contains

    class MolecularModel {
        Pre-built molecular parameters
        +get_Molecule()
    }

    class BacterioChlorophyll {
        BChl a parameters
    }

    class ChlorophyllA {
        Chl a parameters
    }

    class SpectralDensityDB {
        Literature spectral densities
        +get_SpectralDensity(name)
    }

    class CorrelationFunctionDB {
        Literature correlation functions
        +get_CorrelationFunction(name)
    }

    class ModelGenerator {
        +generate(model_type)
    }
```

## Wizard (High-Level Interface)

```mermaid
classDiagram
    direction BT

    Simulation <|-- ExcitonDynamics
    Simulation <|-- TwoDPathways

    class Simulation {
        High-level simulation runner
        +input : Input
        +run()
    }

    class Input {
        +parameters : dict
        Configuration for Wizard simulations
    }

    class ExcitonDynamics {
        Population/coherence dynamics
    }

    class TwoDPathways {
        2D spectroscopy via pathways
    }
```

## Module Dependency Map

```
quantarhei/
├── core/              ← Axes, units, managers, DFunction, saving, parallelization
│   ├── managers.py         (Manager singleton, energy_units, eigenbasis_of, ...)
│   ├── time.py             (TimeAxis)
│   ├── frequency.py        (FrequencyAxis)
│   ├── valueaxis.py        (ValueAxis base)
│   ├── dfunction.py        (DFunction - data on axis)
│   ├── saveable.py         (Serialization: save/load)
│   ├── parallel.py         (distributed_configuration, parallel_function)
│   ├── parcel.py           (Parcel - inter-process data)
│   └── units.py            (convert, in_current_units)
│
├── builders/          ← Physical system construction
│   ├── molecules.py          (Molecule)
│   ├── modes.py              (Mode, SubMode)
│   ├── aggregate_base.py     (AggregateBase: build, diagonalize, basis conversion)
│   ├── aggregate_spectroscopy.py  (liouville_pathways_3, _3T)
│   ├── aggregate_excitonanalysis.py (exciton_report, expansion squares)
│   ├── aggregate_pdeph.py    (get_PureDephasing)
│   ├── aggregates.py         (Aggregate - top-level user class)
│   ├── opensystem.py         (OpenSystem - relaxation theory interface)
│   ├── disorder.py           (Disorder models)
│   ├── pdb.py                (PDB file parsing)
│   └── sysmodes.py           (HarmonicMode, AnharmonicMode)
│
├── qm/                ← Quantum mechanical objects
│   ├── hilbertspace/
│   │   ├── operators.py       (Operator, SelfAdjointOperator)
│   │   ├── hamiltonian.py     (Hamiltonian)
│   │   ├── dmoment.py         (TransitionDipoleMoment)
│   │   └── statevector.py     (StateVector, OQSStateVector)
│   ├── liouvillespace/
│   │   ├── superoperator.py   (SuperOperator)
│   │   ├── redfieldtensor.py  (RedfieldRelaxationTensor, TDRedfield)
│   │   ├── foerstertensor.py  (FoersterRelaxationTensor, TDFoerster)
│   │   ├── modrftensr.py      (ModRedfieldRelaxationTensor)
│   │   ├── lindbladform.py    (LindbladForm, Electronic, VibrationalDecay)
│   │   ├── redfieldfoerster.py (RedfieldFoersterRelaxationTensor)
│   │   ├── nefoerstertensor.py (NEFoersterRelaxationTensor)
│   │   ├── evolutionsuperoperator.py (EvolutionSuperOperator)
│   │   ├── puredephasing.py   (PureDephasing, ElectronicPureDephasing)
│   │   ├── rates.py           (RateMatrix variants)
│   │   ├── heom.py            (KTHierarchy, KTHierarchyPropagator)
│   │   ├── systembathinteraction.py (SystemBathInteraction)
│   │   └── supopunity.py      (SOpUnity - identity superoperator)
│   ├── corfunctions/
│   │   ├── correlationfunctions.py  (CorrelationFunction, CF matrix)
│   │   ├── spectraldensities.py     (SpectralDensity)
│   │   ├── lineshapes.py            (LineshapeFunction)
│   │   └── functionstorages.py      (FunctionStorage, FastFunctionStorage)
│   ├── propagators/
│   │   ├── rdmpropagator.py         (ReducedDensityMatrixPropagator)
│   │   ├── svpropagator.py          (StateVectorPropagator)
│   │   ├── poppropagator.py         (PopulationPropagator)
│   │   ├── dmevolution.py           (DensityMatrixEvolution)
│   │   └── statevectorevolution.py  (StateVectorEvolution)
│   └── oscillators/
│       └── ho.py                    (Harmonic oscillator algebra)
│
├── spectroscopy/      ← Spectrum calculation
│   ├── abs2.py                  (AbsSpectrum)
│   ├── abscalculator.py         (AbsSpectrumCalculator)
│   ├── abscontainer.py          (AbsSpectrumContainer)
│   ├── mockabscalculator.py     (MockAbsSpectrumCalculator)
│   ├── fluorescence.py          (FluorSpectrum, Calculator, Container)
│   ├── circular_dichroism.py    (CircDichSpectrum, Calculator, Container)
│   ├── linear_dichroism.py      (LinDichSpectrum, Calculator, Container)
│   ├── twodresponse.py          (TwoDResponse)
│   ├── twodcalculator.py        (TwoDResponseCalculator)
│   ├── mocktwodcalculator.py    (MockTwoDResponseCalculator)
│   ├── twodcontainer.py         (TwoDResponseContainer, TwoDSpectrumContainer)
│   ├── twodspect.py             (TwoDSpectrum)
│   ├── pumpprobe.py             (PumpProbeSpectrum, Calculator, Container)
│   ├── pathwayanalyzer.py       (LiouvillePathwayAnalyzer)
│   ├── responses.py             (LiouvillePathway, NonLinearResponse)
│   ├── diagramatics.py          (DSFeynmanDiagram, R1g/R2g/.../R2f)
│   └── labsetup.py              (LabSetup, LabField)
│
├── models/            ← Pre-built molecular models & spectral density database
│   ├── molecularmodels.py       (BacterioChlorophyll, Chlorophyll, etc.)
│   └── spectral_densities.py   (SpectralDensityDB, CorrelationFunctionDB)
│
├── wizard/            ← High-level simulation interface
│   ├── simulation.py            (Simulation, ExcitonDynamics, TwoDPathways)
│   └── input/input.py           (Input)
│
└── symbolic/          ← Cumulant expansion theory (SymPy-based, advanced)
    └── cumulant.py              (CumulantExpr, symbolic operators dV, Uop, etc.)
```
