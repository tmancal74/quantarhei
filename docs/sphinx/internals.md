(quantarhei-internals)=

# Quantarhei Internals

This section documents the expert-level classes that power Quantarhei's
internal machinery. These components are used throughout the library but are
not typically part of a user's day-to-day workflow. Understanding them is
useful when extending the library, writing custom propagators, or diagnosing
performance-sensitive code.

The main internal subsystems are:

- **`qm.liouvillespace`** — relaxation tensors (Redfield, Förster, HEOM) and
  Liouville-space propagators used to evolve the density matrix.
- **`qm.corfunctions`** — lineshape functions and bath correlation functions
  that underpin energy-gap fluctuation models.
- **`builders.aggregate_base`** — low-level aggregate construction utilities,
  including Hamiltonian assembly and transition-dipole calculations.
- **`implementations`** — back-end implementations (e.g. QuTiP-based HEOM
  solver) that can be swapped depending on the available numerical libraries.

```{note}
This section is under construction. Full API documentation for each
subsystem will be added in future releases.
```
