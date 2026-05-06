# Quickstart

This page walks through a minimal working simulation in five steps.

## Step 1: Installation

Install Quantarhei from PyPI:

```bash
pip install quantarhei
```

## Step 2: Create a molecule

Define a two-level molecule with a ground state and an excited state at
12 000 cm⁻¹:

```python
import quantarhei as qr

with qr.energy_units("1/cm"):
    m1 = qr.Molecule(name="Mol1", elenergies=[0.0, 12000.0])
    m2 = qr.Molecule(name="Mol2", elenergies=[0.0, 12100.0])
```

`energy_units` sets the unit context so all energy values are interpreted as
inverse centimetres.

## Step 3: Build an aggregate (dimer)

Collect the two molecules into an `Aggregate`, assign transition dipole
moments, set the resonance coupling, and build the single-exciton manifold:

```python
agg = qr.Aggregate(name="Dimer")
agg.add_Molecule(m1)
agg.add_Molecule(m2)

m1.set_dipole(0, 1, [1.0, 0.0, 0.0])
m2.set_dipole(0, 1, [0.0, 1.0, 0.0])

with qr.energy_units("1/cm"):
    agg.set_resonance_coupling(0, 1, 100.0)

agg.build(mult=1)
```

`mult=1` restricts the Hilbert space to at most one excitation at a time.

## Step 4: Add a bath

Create an overdamped Brownian oscillator correlation function and attach it
to each molecule's transition:

```python
ta = qr.TimeAxis(0.0, 1000, 1.0)

params = {
    "ftype": "OverdampedBrownian",
    "reorg": 20.0,
    "cortime": 100.0,
    "T": 300.0,
    "matsubara": 20,
}

with qr.energy_units("1/cm"):
    cf = qr.CorrelationFunction(ta, params)

m1.set_transition_environment((0, 1), cf)
m2.set_transition_environment((0, 1), cf)
```

The reorganisation energy (`reorg`) and correlation time (`cortime`) are given
in cm⁻¹ and femtoseconds respectively; temperature `T` is in Kelvin.

## Step 5: Retrieve the Hamiltonian

Rebuild the aggregate after the bath is attached and inspect the system
Hamiltonian:

```python
agg.build(mult=1)

H = agg.get_Hamiltonian()

with qr.energy_units("1/cm"):
    print(H)
```

The printed matrix shows site energies on the diagonal and resonance
couplings on the off-diagonal, all in cm⁻¹.
