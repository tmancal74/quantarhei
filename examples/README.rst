Demo scripts of the package
===========================

This directory contains scripts demonstrating basic usage
of the Quantarhei package.

Examples Index
--------------

Basic Setup
~~~~~~~~~~~

- ``demo_001_Molecule_Hamiltonian.py`` — Building a Molecule and its Hamiltonian
- ``demo_002_Molecule_Aggregate.py`` — Building a molecular Aggregate
- ``demo_005_UnitsManagementHamiltonian.py`` — Units management for Hamiltonians
- ``demo_100_TimeAndFrequencyAxis.py`` — Working with TimeAxis and FrequencyAxis
- ``demo_101_OperatorFactory.py`` — Using the operator factory
- ``demo_102_DFunction_fitting.py`` — Fitting with DFunction

Correlation Functions & Spectral Densities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``demo_003_CorrFcnSpectDens.py`` — Correlation functions and spectral densities
- ``demo_004_SpectDensDatabase.py`` — Spectral density database

Linear Spectroscopy
~~~~~~~~~~~~~~~~~~~

- ``demo_006_Absorption_1.py`` — Absorption spectrum calculation

Relaxation Theory
~~~~~~~~~~~~~~~~~

- ``demo_010_RedfieldTheory_1.py`` — Redfield theory basics
- ``demo_011_LindbladForm_1.py`` — Lindblad form of relaxation
- ``demo_012_Integrodiff.py`` — Integro-differential equation approach
- ``demo_015_RedfieldTheory_2.py`` — Redfield theory with Aggregate object
- ``demo_016_FoersterTheory_1.py`` — Foerster theory with Aggregate object
- ``demo_017_RedfieldTheory_MultiExcitons.py`` — Redfield theory with multi-excitons
- ``demo_018_ModifiedRedfieldTheory_1.py`` — Modified Redfield theory
- ``demo_200_Secular.py`` — Secularization of a relaxation tensor

HEOM
~~~~

- ``demo_013_HEOM.py`` — Hierarchical Equations of Motion (HEOM)
- ``demo_014_HEOM_rates.py`` — HEOM relaxation rates

2D Spectroscopy & Pump-Probe
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``demo_020_EvolutionSuperOperator_1.py`` — Evolution superoperator
- ``demo_030_2DSpectrum1.py`` — 2D spectrum with effective lineshapes
- ``demo_031_PumpProbeSpectrum1.py`` — Pump-probe spectrum
- ``demo_035_2D_Vibronic.py`` — 2D vibronic spectrum
- ``demo_036_2DSpectrumDisorder.py`` — 2D spectrum with disorder
- ``demo_040_2DSpectrumAceto.py`` — 2D spectrum with lineshape functions (Aceto)
- ``demo_854_2DSpectrum_DimerDisorder.py`` — 2D spectrum of a dimer with disorder

PDB & Biological Systems
~~~~~~~~~~~~~~~~~~~~~~~~~

- ``demo_050_PDB_FMO1.py`` — FMO complex from PDB structure

Advanced & Miscellaneous
~~~~~~~~~~~~~~~~~~~~~~~~

- ``demo_300_ParallelIterators.py`` — Parallel iterators
- ``demo_800_DiagProblem.py`` — Diagonalization problem
- ``demo_850_vibrons.py`` — Vibrons
- ``demo_853_RC.py`` — Reaction center simulation
- ``demo_853o_RC.py`` — Reaction center simulation (variant)

Admin
-----

The demo files can be regenerated from within the package by running::

   python admin/make_demos.py

from the ``examples/`` directory.
