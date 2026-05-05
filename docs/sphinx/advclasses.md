(advanced-level-classes)=

# Advanced Level Classes

Advanced level classes are one level below
the user level classes and their usage is for instance as follows

```python
import quantarhei as qr

sbi = qr.qm.SystemBathInteraction
```

## Molecular Environment

```{toctree}
:maxdepth: 2

classes/systembathinteraction
```

## Relaxation Tensors

```{toctree}
:maxdepth: 2

classes/relaxtensors/relaxtensors
classes/relaxtensors/redfieldtensor
classes/relaxtensors/foerstertensor
```

## Open System Propagation

`LindbladForm` is used for Markovian open-system dynamics, providing a
Lindblad master-equation relaxation tensor. `EvolutionSuperOperator`
handles time-evolution of the density matrix and is the primary tool
for propagating open quantum systems over a time grid.

```{toctree}
:maxdepth: 2

classes/evolsupop
```
