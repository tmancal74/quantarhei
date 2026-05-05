# How to Use `Aggregate` Class

The `Aggregate` class combines multiple `Molecule` objects into a coupled
system. The example below shows how to build a homodimer and diagonalize its
Hamiltonian.

## Creating an aggregate

```python
import quantarhei as qr

en = [0.0, 1.0]
m1 = qr.Molecule(name="Mol1", elenergies=en)
m2 = qr.Molecule(name="Mol2", elenergies=en)

ag = qr.Aggregate(name="Homodimer")
ag.add_Molecule(m1)
ag.add_Molecule(m2)

# set the resonance coupling between molecule 0 and molecule 1
ag.set_resonance_coupling(0, 1, 0.1)

# build the aggregate in the single-excitation manifold
ag.build(mult=1)

H = ag.get_Hamiltonian()
print(H)
```

## Using physical energy units

Energies can be supplied in wavenumbers by wrapping assignments in an
`energy_units` context manager:

```python
import quantarhei as qr

with qr.energy_units("1/cm"):
    m1 = qr.Molecule(name="Mol1", elenergies=[0.0, 10100.0])
    m2 = qr.Molecule(name="Mol2", elenergies=[0.0, 10100.0])

ag = qr.Aggregate(name="Homodimer")
ag.add_Molecule(m1)
ag.add_Molecule(m2)
ag.set_resonance_coupling(0, 1, 0.1)
ag.build(mult=1)

H = ag.get_Hamiltonian()
with qr.energy_units("1/cm"):
    print(H)
```
