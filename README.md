[![Build Status](https://github.com/tmancal74/quantarhei/actions/workflows/python-package.yml/badge.svg)](https://github.com/tmancal74/quantarhei)
[![Documentation Status](https://readthedocs.org/projects/quantarhei/badge/?version=latest)](http://quantarhei.readthedocs.io/en/latest/?badge=latest)
[![Coverage](https://img.shields.io/codecov/c/github/tmancal74/quantarhei.svg)](https://codecov.io/gh/tmancal74/quantarhei)
[![Version](https://img.shields.io/pypi/v/quantarhei.svg)](https://pypi.org/project/quantarhei/)
[![Python](https://img.shields.io/pypi/pyversions/quantarhei.svg)](https://pypi.org/project/quantarhei/)

# QUANTArhei: Open Quantum System Theory for Molecular Systems

Quantarhei is a molecular open quantum systems simulator written in Python. Its name is derived from the aphorism *"Panta rhei"* of Heraclitus of Ephesus — "Everything flows" — fitting for a package centred on quantum dynamics.

It provides tools for building molecular aggregates, computing linear and non-linear spectra, and simulating excitation energy transfer and open quantum system dynamics.

---

## Installation

**From PyPI (recommended):**

```bash
pip install quantarhei
```

**With [uv](https://github.com/astral-sh/uv):**

```bash
uv add quantarhei
```

**From source:**

```bash
git clone https://github.com/tmancal74/quantarhei.git
cd quantarhei
pip install -e .
```

**Requirements:** Python 3.10 or later, NumPy, SciPy.

---

## Quick start

```python
import quantarhei as qr

# Define two molecules with transition energies (in 1/cm)
with qr.energy_units("1/cm"):
    mol1 = qr.Molecule([0.0, 12000.0])
    mol2 = qr.Molecule([0.0, 12200.0])

    mol1.set_dipole(0, 1, [1.0, 0.0, 0.0])
    mol2.set_dipole(0, 1, [0.0, 1.0, 0.0])

# Build a dimer aggregate with resonance coupling
agg = qr.Aggregate([mol1, mol2])
with qr.energy_units("1/cm"):
    agg.set_resonance_coupling(0, 1, 100.0)

agg.build()

# Calculate absorption spectrum
with qr.energy_units("1/cm"):
    time = qr.TimeAxis(0.0, 1000, 1.0)
    calc = qr.AbsSpectrumCalculator(time, system=agg)
    spec = calc.calculate()
    spec.plot()
```

More examples are in the [`docs/examples/`](docs/examples/) directory and on [Read the Docs](http://quantarhei.readthedocs.io).

---

## Features

- **Molecular aggregates** — multi-molecule systems with arbitrary couplings, vibrational modes (Huang-Rhys factors, intra-molecular modes), vibronic coupling, and system-bath interactions
- **Linear spectroscopy** — absorption, circular dichroism, linear dichroism, and fluorescence spectra for monomers and aggregates
- **2D electronic spectroscopy** — non-linear response with Liouville pathway analysis, Feynman diagrams, and pump-probe spectra
- **Energy transfer theories** — Redfield, modified Redfield, and Förster rate theories; time-dependent and non-equilibrium variants; mixed Redfield-Förster
- **Open quantum system dynamics** — Lindblad and Redfield master equations, density matrix propagation, hierarchical equations of motion (HEOM)
- **Bath modeling** — correlation functions, spectral densities, lineshape functions; pre-built models for photosynthetic pigments
- **Laboratory frame** — realistic optical geometries and pulse field configurations via `LabSetup`
- **Extensible** — pure Python reference implementation; performance-critical routines can be replaced with optimised C/Fortran extensions

---

## Documentation

Full API reference and tutorials: [quantarhei.readthedocs.io](http://quantarhei.readthedocs.io)

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full version history.

---

## Acknowledgements

The development of Quantarhei has been supported by:

- [**Neuron Fund for Support of Science**](http://www.nfneuron.cz) — Impuls grant in physics 2014 (2015–2017)
- [**Czech Science Foundation (GACR)**](http://www.gacr.cz) — grants 14-25752S (2014–2016), 17-22160S (2017–2019), 18-18022S (2018–2020)
