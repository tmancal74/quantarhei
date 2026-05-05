# How to Use `Molecule` Class

The `Molecule` class represents a single quantum-mechanical molecule with a set
of electronic energy levels. The examples below show the most common operations.

## Creating a molecule

Provide a list of electronic energies (in internal units, or inside a
`with qr.energy_units(...)` block) and retrieve the Hamiltonian:

```python
import quantarhei as qr

# electronic energies: ground state at 0, excited state at 1 (internal units)
en = [0.0, 1.0]
m = qr.Molecule(elenergies=en)

H = m.get_Hamiltonian()
print(H)
```

## Setting a transition dipole moment

Use `set_dipole` to assign the transition dipole moment between two states:

```python
import quantarhei as qr

en = [0.0, 1.0]
m = qr.Molecule(elenergies=en)

# set the 0->1 transition dipole along the x-axis
m.set_dipole(0, 1, [1.0, 0.0, 0.0])
```

## Adding a vibrational mode

Intramolecular vibrational modes are added with `Mode` and `add_Mode`:

```python
import quantarhei as qr

en = [0.0, 1.0]
m = qr.Molecule(elenergies=en)

# create a mode with frequency 0.01 (internal units)
mod = qr.Mode(frequency=0.01)
m.add_Mode(mod)
```
