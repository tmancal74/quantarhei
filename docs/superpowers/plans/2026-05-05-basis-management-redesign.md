# Basis Management Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Redesign the quantarhei basis management system so that tensor constructors and the propagator own their own `eigenbasis_of` contexts, making site-basis the invariant ground truth and eliminating the error-prone `protect_basis`/`unprotect_basis` pattern.

**Architecture:** Add `BasisError`/`BasisContextError` exceptions and two guard functions to `managers.py`. Replace silent auto-transform in `basis_managed_array_property` with a strict raise. Internalize `eigenbasis_of` into `RedfieldRelaxationTensor` and `ModRedfieldRelaxationTensor` constructors, add a `secular=` flag, and record `basis_op`. Update `ReducedDensityMatrixPropagator` to wrap `propagate()` in `eigenbasis_of(basis_op)` automatically. Remove `protect_basis`/`unprotect_basis` from `BasisManaged`. Update `opensystem.py` and affected examples.

**Tech Stack:** Python 3.10+, numpy, pytest, quantarhei internal APIs (`Manager`, `eigenbasis_of`, `BasisManaged`, `RelaxationTensor`, `Secular`)

**Spec:** `docs/superpowers/specs/2026-05-05-basis-management-redesign.md`

---

## File Map

| File | Change |
|---|---|
| `quantarhei/core/managers.py` | Add `BasisError`, `BasisContextError`; add `assert_not_in_eigenbasis_context()`, `assert_in_eigenbasis_of(op)`; remove `protect_basis`/`unprotect_basis` from `BasisManaged` |
| `quantarhei/utils/types.py` | Replace silent auto-transform in `basis_managed_array_property` with raise |
| `quantarhei/qm/liouvillespace/redfieldtensor.py` | Add `secular=False` param to `__init__`; call `self.secularize()` at end of `_implementation` when `secular=True`; set `self.basis_op = ham` after init |
| `quantarhei/qm/liouvillespace/modredfieldtensor.py` | Wrap `initialize()` content in outer `eigenbasis_of(HH)`; set `self.basis_op = ham` |
| `quantarhei/qm/liouvillespace/relaxationtensor.py` | Update `enforce_detailed_balance` — it already uses `eigenbasis_of` internally, verify it still works; add guard at top |
| `quantarhei/qm/propagators/rdmpropagator.py` | Add `basis=None` param to `__init__`; store as `self._basis`; wrap the `has_relaxation` short-exp path in `propagate()` with `eigenbasis_of(basis_op)` |
| `quantarhei/builders/opensystem.py` | Remove all `protect_basis`/`unprotect_basis` calls; remove `eigenbasis_of` wrapping around tensor constructors |
| `quantarhei/wizard/examples/ex_010_RedfieldTheory_1.py` | Update to new API |
| `quantarhei/wizard/examples/ex_018_ModifiedRedfieldTheory_1.py` | Update to new API |
| `quantarhei/wizard/examples/ex_200_Secular.py` | Update to new API |
| `examples/demo_010_RedfieldTheory_1.py` | Update to new API |
| `examples/demo_018_ModifiedRedfieldTheory_1.py` | Update to new API |
| `examples/demo_200_Secular.py` | Update to new API |
| `tests/unit/core/test_basis_error.py` | New: guard and raise-on-mismatch tests |
| `tests/unit/qm/propagators/test_redfield_boltzmann.py` | Remove `xfail`, assert convergence |

---

## Task 1: Add `BasisError`, `BasisContextError`, and guard functions

**Files:**
- Modify: `quantarhei/core/managers.py:910-927` (BasisManaged class)
- Test: `tests/unit/core/test_basis_error.py` (new)

- [ ] **Step 1.1: Write the failing tests**

Create `tests/unit/core/test_basis_error.py`:

```python
import pytest
import numpy
import quantarhei as qr
from quantarhei.core.managers import (
    BasisError,
    BasisContextError,
    assert_not_in_eigenbasis_context,
)


def test_basis_context_error_is_basis_error():
    assert issubclass(BasisContextError, BasisError)


def test_assert_not_in_eigenbasis_context_passes_outside():
    # Should not raise when outside any eigenbasis context
    assert_not_in_eigenbasis_context()


def test_assert_not_in_eigenbasis_context_raises_inside():
    H = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    with pytest.raises(BasisContextError):
        with qr.eigenbasis_of(H):
            assert_not_in_eigenbasis_context()
```

- [ ] **Step 1.2: Run the tests to confirm they fail**

```bash
cd /Users/danherma/projects-personal/quantarhei
pytest tests/unit/core/test_basis_error.py -v
```

Expected: `ImportError` or `AttributeError` — `BasisError`, `BasisContextError`, `assert_not_in_eigenbasis_context` do not exist yet.

- [ ] **Step 1.3: Add exceptions and guard functions to `managers.py`**

In `quantarhei/core/managers.py`, directly after the `BasisManaged` class definition (after line 927), add:

```python
class BasisError(Exception):
    """Raised when an operator is accessed in the wrong basis."""


class BasisContextError(BasisError):
    """Raised when a calculation routine is called in the wrong basis context."""


def assert_not_in_eigenbasis_context() -> None:
    """Raise BasisContextError if currently inside an eigenbasis_of context.

    Call this at the top of any tensor constructor that manages its own
    eigenbasis context internally.
    """
    if Manager()._in_eigenbasis_of_context:
        raise BasisContextError(
            "Tensor constructors must not be called inside an eigenbasis_of "
            "context. Remove the enclosing 'with eigenbasis_of(...)' block — "
            "the constructor manages its own basis context internally."
        )


def assert_in_eigenbasis_of(op: Any) -> None:
    """Raise BasisError if not currently inside eigenbasis_of(op)."""
    mgr = Manager()
    if not mgr._in_eigenbasis_of_context:
        raise BasisError(
            "This operation requires an active eigenbasis_of context."
        )
    if mgr.current_basis_operator is not op:
        raise BasisError(
            "Active eigenbasis context does not match the expected operator."
        )
```

Also add `BasisError`, `BasisContextError`, `assert_not_in_eigenbasis_context`, `assert_in_eigenbasis_of` to the module's public exports (add to `__all__` if present, otherwise just ensure they are importable).

- [ ] **Step 1.4: Run the tests to confirm they pass**

```bash
pytest tests/unit/core/test_basis_error.py -v
```

Expected: all 3 tests pass.

- [ ] **Step 1.5: Commit**

```bash
git add quantarhei/core/managers.py tests/unit/core/test_basis_error.py
git commit -m "#325 feat: add BasisError, BasisContextError and context guard functions"
```

---

## Task 2: Replace silent auto-transform with raise in `.data` property

**Files:**
- Modify: `quantarhei/utils/types.py:87-132`
- Test: `tests/unit/core/test_basis_error.py` (extend)

- [ ] **Step 2.1: Add failing tests**

Append to `tests/unit/core/test_basis_error.py`:

```python
def test_data_access_in_correct_basis_does_not_raise():
    H = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    # Access in site basis — no context, should work
    _ = H.data


def test_data_access_in_wrong_basis_raises():
    H1 = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    H2 = qr.Hamiltonian(data=numpy.array([[0.0, 0.2], [0.2, 2.0]]))
    with qr.eigenbasis_of(H1):
        # H2 was created outside and not part of this context — accessing
        # its .data here should raise because its basis ID doesn't match
        # the current Manager basis.
        # Note: H2 was registered at basis 0 (site basis), Manager is now
        # at a different basis.
        with pytest.raises(BasisError):
            _ = H2.data
```

- [ ] **Step 2.2: Run to confirm the test fails (currently returns silently)**

```bash
pytest tests/unit/core/test_basis_error.py::test_data_access_in_wrong_basis_raises -v
```

Expected: FAIL — currently the property silently transforms instead of raising.

- [ ] **Step 2.3: Update `basis_managed_array_property` in `types.py`**

In `quantarhei/utils/types.py`, in the `prop` getter of `basis_managed_array_property` (lines 93-107), replace the silent transform branch:

```python
# BEFORE (lines 101-106):
if cb == ob:
    pass
else:
    # change basis
    self.manager.transform_to_current_basis(self)
```

with:

```python
if cb != ob:
    from .managers import BasisError  # avoid circular import at module level
    raise BasisError(
        f"Operator '{getattr(self, 'name', type(self).__name__)}' is in "
        f"basis {ob} but the current Manager basis is {cb}. "
        "Access .data only inside the correct eigenbasis_of context, "
        "or outside all contexts (site basis)."
    )
```

Apply the same replacement to the `prop.setter` (lines 109-121) — replace:

```python
if cb == ob:
    pass
else:
    # change basis
    self.manager.transform_to_current_basis(self)
```

with:

```python
if cb != ob:
    from .managers import BasisError
    raise BasisError(
        f"Cannot set .data on operator '{getattr(self, 'name', type(self).__name__)}': "
        f"it is in basis {ob} but the current Manager basis is {cb}."
    )
```

- [ ] **Step 2.4: Run all basis error tests**

```bash
pytest tests/unit/core/test_basis_error.py -v
```

Expected: all tests pass.

- [ ] **Step 2.5: Run the full unit test suite to catch regressions**

```bash
pytest tests/unit/ -x -q 2>&1 | head -60
```

Expected: existing tests that relied on silent transform will now fail — that is intentional. Note which tests fail; they will be fixed in subsequent tasks. If the failure count is overwhelming (>20 tests), stop and reassess before continuing.

- [ ] **Step 2.6: Commit**

```bash
git add quantarhei/utils/types.py tests/unit/core/test_basis_error.py
git commit -m "#325 feat: raise BasisError on .data access in wrong basis context"
```

---

## Task 3: Remove `protect_basis`/`unprotect_basis` from `BasisManaged`

**Files:**
- Modify: `quantarhei/core/managers.py:922-926`
- Test: `tests/unit/core/test_basis_error.py` (extend)

- [ ] **Step 3.1: Add a failing test**

Append to `tests/unit/core/test_basis_error.py`:

```python
def test_protect_basis_no_longer_exists():
    H = qr.Hamiltonian(data=numpy.array([[0.0, 0.1], [0.1, 1.0]]))
    assert not hasattr(H, "protect_basis"), (
        "protect_basis should be removed from BasisManaged"
    )
    assert not hasattr(H, "unprotect_basis"), (
        "unprotect_basis should be removed from BasisManaged"
    )
```

- [ ] **Step 3.2: Run to confirm it fails**

```bash
pytest tests/unit/core/test_basis_error.py::test_protect_basis_no_longer_exists -v
```

Expected: FAIL — `H.protect_basis` currently exists.

- [ ] **Step 3.3: Remove `protect_basis` and `unprotect_basis` from `BasisManaged`**

In `quantarhei/core/managers.py`, delete lines 922-926:

```python
# DELETE these four lines:
    def protect_basis(self) -> None:
        self.is_basis_protected = True

    def unprotect_basis(self) -> None:
        self.is_basis_protected = False
```

Also remove the `is_basis_protected = False` class attribute on line 914, and remove the `is_basis_protected` check in `eigenbasis_of.__exit__` at line 1126:

```python
# DELETE the if-guard (keep the transform call unconditional):
# BEFORE:
if not op.is_basis_protected:
    op.transform(S1, inv=SS)
op.set_current_basis(nb)

# AFTER:
op.transform(S1, inv=SS)
op.set_current_basis(nb)
```

- [ ] **Step 3.4: Run the tests**

```bash
pytest tests/unit/core/test_basis_error.py -v
```

Expected: all tests pass.

- [ ] **Step 3.5: Commit**

```bash
git add quantarhei/core/managers.py tests/unit/core/test_basis_error.py
git commit -m "#325 refactor: remove protect_basis/unprotect_basis from BasisManaged"
```

---

## Task 4: Internalize `eigenbasis_of` in `RedfieldRelaxationTensor`

**Files:**
- Modify: `quantarhei/qm/liouvillespace/redfieldtensor.py:91-148`
- Test: `tests/unit/qm/liouvillespace/test_redfieldtensor.py` (extend or create)

- [ ] **Step 4.1: Write a failing test confirming `basis_op` is set and old pattern raises**

Create or append to `tests/unit/qm/liouvillespace/test_redfieldtensor.py`:

```python
import pytest
import numpy
import quantarhei as qr
from quantarhei.core.managers import BasisContextError


def _make_ham_sbi():
    with qr.energy_units("1/cm"):
        mol1 = qr.Molecule([0.0, 12000.0])
        mol2 = qr.Molecule([0.0, 12200.0])
        mol1.set_transition_environment((0, 1), qr.CorrelationFunction(
            qr.TimeAxis(0.0, 1000, 1.0),
            {"ftype": "OverdampedBrownian", "reorg": 100.0, "cortime": 100.0, "T": 300.0, "matsubara": 20}
        ))
        mol2.set_transition_environment((0, 1), qr.CorrelationFunction(
            qr.TimeAxis(0.0, 1000, 1.0),
            {"ftype": "OverdampedBrownian", "reorg": 100.0, "cortime": 100.0, "T": 300.0, "matsubara": 20}
        ))
        agg = qr.Aggregate([mol1, mol2])
        agg.set_resonance_coupling(0, 1, 100.0)
        agg.build()
        ham = agg.get_Hamiltonian()
        sbi = agg.get_SystemBathInteraction()
    return ham, sbi


def test_redfieldtensor_sets_basis_op():
    ham, sbi = _make_ham_sbi()
    RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi)
    assert RRT.basis_op is ham


def test_redfieldtensor_raises_inside_eigenbasis_context():
    ham, sbi = _make_ham_sbi()
    with pytest.raises(BasisContextError):
        with qr.eigenbasis_of(ham):
            qr.qm.RedfieldRelaxationTensor(ham, sbi)


def test_redfieldtensor_secular_flag():
    ham, sbi = _make_ham_sbi()
    RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi, secular=True)
    assert RRT.is_secular
```

- [ ] **Step 4.2: Run to confirm they fail**

```bash
pytest tests/unit/qm/liouvillespace/test_redfieldtensor.py::test_redfieldtensor_sets_basis_op tests/unit/qm/liouvillespace/test_redfieldtensor.py::test_redfieldtensor_raises_inside_eigenbasis_context tests/unit/qm/liouvillespace/test_redfieldtensor.py::test_redfieldtensor_secular_flag -v
```

Expected: FAIL — `basis_op` does not exist, guard is not in place, `secular=` param does not exist.

- [ ] **Step 4.3: Update `RedfieldRelaxationTensor.__init__` in `redfieldtensor.py`**

At the top of `__init__` (after type checks, before `self.Hamiltonian = ham`), add:

```python
from ...core.managers import assert_not_in_eigenbasis_context
assert_not_in_eigenbasis_context()
```

Add `secular: bool = False` to the `__init__` signature:

```python
def __init__(
    self,
    ham: Hamiltonian,
    sbi: SystemBathInteraction,
    initialize: bool = True,
    cutoff_time: float | None = None,
    as_operators: bool = False,
    name: str = "",
    secular: bool = False,
) -> None:
```

Store the `secular` flag before the `initialize` call:

```python
self._secular_on_init = secular
```

In `_implementation`, after the tensor computation is complete and before the method returns, add:

```python
if self._secular_on_init:
    self.secularize()
```

After the `with energy_units("int"):` block in `__init__` (after `self._implementation(ham, sbi)` is called), set `basis_op`:

```python
self.basis_op = ham
```

The full `__init__` tail should look like:

```python
        if initialize:
            try:
                pass
            except Exception:
                pass
            with energy_units("int"):
                self._implementation(ham, sbi)
            self.basis_op = ham  # <-- add this line

        self.Iterm = None
        self.has_Iterm = False
```

- [ ] **Step 4.4: Run the new tests**

```bash
pytest tests/unit/qm/liouvillespace/test_redfieldtensor.py::test_redfieldtensor_sets_basis_op tests/unit/qm/liouvillespace/test_redfieldtensor.py::test_redfieldtensor_raises_inside_eigenbasis_context tests/unit/qm/liouvillespace/test_redfieldtensor.py::test_redfieldtensor_secular_flag -v
```

Expected: all 3 pass.

- [ ] **Step 4.5: Commit**

```bash
git add quantarhei/qm/liouvillespace/redfieldtensor.py tests/unit/qm/liouvillespace/test_redfieldtensor.py
git commit -m "#325 feat: RedfieldRelaxationTensor owns its eigenbasis context, adds secular= flag"
```

---

## Task 5: Internalize `eigenbasis_of` in `ModRedfieldRelaxationTensor`

**Files:**
- Modify: `quantarhei/qm/liouvillespace/modredfieldtensor.py:68-110`
- Test: `tests/unit/qm/liouvillespace/test_modredfieldtensor.py` (extend or create)

- [ ] **Step 5.1: Write failing tests**

Create or append to `tests/unit/qm/liouvillespace/test_modredfieldtensor.py`:

```python
import pytest
import numpy
import quantarhei as qr
from quantarhei.core.managers import BasisContextError


def _make_ham_sbi():
    with qr.energy_units("1/cm"):
        mol1 = qr.Molecule([0.0, 12000.0])
        mol2 = qr.Molecule([0.0, 12200.0])
        mol1.set_transition_environment((0, 1), qr.CorrelationFunction(
            qr.TimeAxis(0.0, 1000, 1.0),
            {"ftype": "OverdampedBrownian", "reorg": 100.0, "cortime": 100.0, "T": 300.0, "matsubara": 20}
        ))
        mol2.set_transition_environment((0, 1), qr.CorrelationFunction(
            qr.TimeAxis(0.0, 1000, 1.0),
            {"ftype": "OverdampedBrownian", "reorg": 100.0, "cortime": 100.0, "T": 300.0, "matsubara": 20}
        ))
        agg = qr.Aggregate([mol1, mol2])
        agg.set_resonance_coupling(0, 1, 100.0)
        agg.build()
        ham = agg.get_Hamiltonian()
        sbi = agg.get_SystemBathInteraction()
    return ham, sbi


def test_modredfieldtensor_sets_basis_op():
    ham, sbi = _make_ham_sbi()
    MRRT = qr.qm.ModRedfieldRelaxationTensor(ham, sbi)
    assert MRRT.basis_op is ham


def test_modredfieldtensor_raises_inside_eigenbasis_context():
    ham, sbi = _make_ham_sbi()
    with pytest.raises(BasisContextError):
        with qr.eigenbasis_of(ham):
            qr.qm.ModRedfieldRelaxationTensor(ham, sbi)
```

- [ ] **Step 5.2: Run to confirm they fail**

```bash
pytest tests/unit/qm/liouvillespace/test_modredfieldtensor.py::test_modredfieldtensor_sets_basis_op tests/unit/qm/liouvillespace/test_modredfieldtensor.py::test_modredfieldtensor_raises_inside_eigenbasis_context -v
```

Expected: FAIL.

- [ ] **Step 5.3: Update `ModRedfieldRelaxationTensor.__init__` in `modredfieldtensor.py`**

Add guard at the top of `__init__`, after the type checks:

```python
from ...core.managers import assert_not_in_eigenbasis_context
assert_not_in_eigenbasis_context()
```

The `initialize()` method already has an internal `with eigenbasis_of(HH):` at line 95 — this is correct and stays. It is called from `__init__`, not from user code inside an external `eigenbasis_of`.

After the `if initialize:` block in `__init__`, set `basis_op`:

```python
        if initialize:
            self.initialize()
            self._data_initialized = True
            self._is_initialized = True
            self.basis_op = ham   # <-- add this line
        else:
            self._data_initialized = False
```

- [ ] **Step 5.4: Run the new tests**

```bash
pytest tests/unit/qm/liouvillespace/test_modredfieldtensor.py::test_modredfieldtensor_sets_basis_op tests/unit/qm/liouvillespace/test_modredfieldtensor.py::test_modredfieldtensor_raises_inside_eigenbasis_context -v
```

Expected: both pass.

- [ ] **Step 5.5: Commit**

```bash
git add quantarhei/qm/liouvillespace/modredfieldtensor.py tests/unit/qm/liouvillespace/test_modredfieldtensor.py
git commit -m "#325 feat: ModRedfieldRelaxationTensor owns its eigenbasis context"
```

---

## Task 6: Update `ReducedDensityMatrixPropagator` to own its propagation context

**Files:**
- Modify: `quantarhei/qm/propagators/rdmpropagator.py:80-90` (add `basis` param)
- Modify: `quantarhei/qm/propagators/rdmpropagator.py:313-540` (wrap propagation paths)
- Test: `tests/unit/qm/propagators/test_redfield_boltzmann.py`

- [ ] **Step 6.1: Write the regression test (currently xfail)**

Replace the `xfail` marker in `tests/unit/qm/propagators/test_redfield_boltzmann.py`:

```python
import numpy
import pytest
import quantarhei as qr
from quantarhei import TimeAxis, energy_units


def _boltzmann_ratio(dE_cm, T_K=300.0):
    """Expected p2/p1 = exp(-dE / kBT) in wavenumber units."""
    kB_cm = 0.6950356  # cm^-1 / K
    return numpy.exp(-dE_cm / (kB_cm * T_K))


def test_rate_matrix_satisfies_detailed_balance():
    with energy_units("1/cm"):
        mol1 = qr.Molecule([0.0, 12000.0])
        mol2 = qr.Molecule([0.0, 12360.0])
        ta = TimeAxis(0.0, 1000, 1.0)
        for mol in (mol1, mol2):
            mol.set_transition_environment(
                (0, 1),
                qr.CorrelationFunction(
                    ta,
                    {"ftype": "OverdampedBrownian", "reorg": 100.0,
                     "cortime": 100.0, "T": 300.0, "matsubara": 20},
                ),
            )
        agg = qr.Aggregate([mol1, mol2])
        agg.set_resonance_coupling(0, 1, 100.0)
        agg.build()
        ham = agg.get_Hamiltonian()
        sbi = agg.get_SystemBathInteraction()

    RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi, secular=True)
    RRM = qr.qm.RedfieldRateMatrix(ham, sbi)

    with qr.eigenbasis_of(ham):
        k_down = RRM.data[0, 1]
        k_up = RRM.data[1, 0]

    ratio = k_up / k_down
    expected = _boltzmann_ratio(360.0)
    assert abs(ratio - expected) / expected < 0.05, (
        f"Rate ratio {ratio:.4f} deviates >5% from Boltzmann {expected:.4f}"
    )


def test_propagator_converges_to_boltzmann():
    with energy_units("1/cm"):
        mol1 = qr.Molecule([0.0, 12000.0])
        mol2 = qr.Molecule([0.0, 12360.0])
        ta = TimeAxis(0.0, 12000, 1.0)
        for mol in (mol1, mol2):
            mol.set_transition_environment(
                (0, 1),
                qr.CorrelationFunction(
                    ta,
                    {"ftype": "OverdampedBrownian", "reorg": 100.0,
                     "cortime": 100.0, "T": 300.0, "matsubara": 20},
                ),
            )
        agg = qr.Aggregate([mol1, mol2])
        agg.set_resonance_coupling(0, 1, 100.0)
        agg.build()
        ham = agg.get_Hamiltonian()
        sbi = agg.get_SystemBathInteraction()

    RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi, secular=True)
    prop = qr.ReducedDensityMatrixPropagator(ta, ham, RRT)

    rho_i = qr.ReducedDensityMatrix(dim=ham.dim)
    with qr.eigenbasis_of(ham):
        rho_i.data[1, 1] = 1.0   # start in higher exciton state

    rhot = prop.propagate(rho_i)

    # Populations readable directly — no eigenbasis_of context needed
    p1_final = numpy.real(rhot.data[-1, 0, 0])  # lower exciton
    p2_final = numpy.real(rhot.data[-1, 1, 1])  # upper exciton

    ratio = p2_final / p1_final
    expected = _boltzmann_ratio(360.0)
    assert abs(ratio - expected) / expected < 0.05, (
        f"Population ratio {ratio:.4f} deviates >5% from Boltzmann {expected:.4f}"
    )
```

- [ ] **Step 6.2: Run to confirm the second test still fails**

```bash
pytest tests/unit/qm/propagators/test_redfield_boltzmann.py -v
```

Expected: `test_rate_matrix_satisfies_detailed_balance` passes; `test_propagator_converges_to_boltzmann` FAILS.

- [ ] **Step 6.3: Add `basis` parameter to `ReducedDensityMatrixPropagator.__init__`**

In `rdmpropagator.py`, update the `__init__` signature (line 80) to add `basis=None`:

```python
def __init__(
    self,
    timeaxis: Any = None,
    Ham: Any = None,
    RTensor: Any = None,
    Iterm: Any = None,
    Efield: Any = None,
    Trdip: Any = None,
    PDeph: Any = None,
    NonHerm: Any = None,
    basis: Any = None,
) -> None:
```

At the end of the `__init__` body (after all existing attribute assignments), store the basis override:

```python
self._basis = basis
```

- [ ] **Step 6.4: Wrap the `propagate()` short-exp path in `eigenbasis_of`**

In `rdmpropagator.py`, find the `propagate()` method (line 245). The method dispatches to various private `__propagate_*` methods. Wrap the dispatch block in an `eigenbasis_of` context when a basis operator is available.

Locate the section starting at line 313:

```python
if self.has_relaxation:
```

Before this block (but still inside `propagate()`), add basis resolution:

```python
from ..core.managers import eigenbasis_of  # already imported at top of file — confirm
_basis_op = self._basis if hasattr(self, "_basis") and self._basis is not None \
            else getattr(getattr(self, "RelaxationTensor", None), "basis_op", None)
```

Then wrap the entire `if self.has_relaxation:` ... `else:` dispatch block:

```python
if _basis_op is not None:
    _ctx = eigenbasis_of(_basis_op)
    _ctx.__enter__()
    try:
        # --- existing if self.has_relaxation: ... else: block goes here ---
        ...
    finally:
        _ctx.__exit__(None, None, None)
else:
    # --- existing if self.has_relaxation: ... else: block goes here (no context) ---
    ...
```

**Important:** Do not duplicate the dispatch logic. Use a helper or a context-manager `with` statement if the code structure allows. The cleanest approach is:

```python
from contextlib import nullcontext
_ctx = eigenbasis_of(_basis_op) if _basis_op is not None else nullcontext()
with _ctx:
    if self.has_relaxation:
        # ... existing dispatch unchanged ...
    else:
        # ... existing dispatch unchanged ...
```

- [ ] **Step 6.5: Run the Boltzmann test**

```bash
pytest tests/unit/qm/propagators/test_redfield_boltzmann.py -v
```

Expected: both tests pass.

- [ ] **Step 6.6: Run the full propagator unit tests**

```bash
pytest tests/unit/qm/propagators/ -v 2>&1 | tail -30
```

Expected: no regressions.

- [ ] **Step 6.7: Commit**

```bash
git add quantarhei/qm/propagators/rdmpropagator.py tests/unit/qm/propagators/test_redfield_boltzmann.py
git commit -m "#325 feat: propagator wraps propagate() in eigenbasis_of(basis_op) automatically"
```

---

## Task 7: Update `opensystem.py` to remove `protect_basis` and external contexts

**Files:**
- Modify: `quantarhei/builders/opensystem.py`

- [ ] **Step 7.1: Remove all `protect_basis`/`unprotect_basis` calls and wrapping `eigenbasis_of` contexts**

In `opensystem.py`, the pattern repeated ~8 times is:

```python
ham.protect_basis()
with eigenbasis_of(ham):
    SomeTensor = SomeTensorClass(ham, sbi, ...)
    # optional: SomeTensor.secularize()
ham.unprotect_basis()
```

Replace each occurrence with the new API:

```python
SomeTensor = SomeTensorClass(ham, sbi, ..., secular=True)  # if secularize() was called
# or:
SomeTensor = SomeTensorClass(ham, sbi, ...)  # if no secularize()
```

Lines to update: 709-719, 724-734, 747-756, 761-771, 868-878, 884-898, 1043-1046.

Also check lines 926-932 (commented-out protect pattern) — remove the commented-out lines entirely.

- [ ] **Step 7.2: Run the opensystem-related unit tests**

```bash
pytest tests/unit/builders/ -v 2>&1 | tail -30
```

Expected: pass (or reveal further breakage to fix).

- [ ] **Step 7.3: Commit**

```bash
git add quantarhei/builders/opensystem.py
git commit -m "#325 refactor: remove protect_basis pattern from opensystem.py"
```

---

## Task 8: Update wizard examples and `examples/` to use new API

**Files:**
- Modify: `quantarhei/wizard/examples/ex_010_RedfieldTheory_1.py`
- Modify: `quantarhei/wizard/examples/ex_018_ModifiedRedfieldTheory_1.py`
- Modify: `quantarhei/wizard/examples/ex_200_Secular.py`
- Modify: `examples/demo_010_RedfieldTheory_1.py`
- Modify: `examples/demo_018_ModifiedRedfieldTheory_1.py`
- Modify: `examples/demo_200_Secular.py`

Files to update:
- `quantarhei/wizard/examples/ex_010_RedfieldTheory_1.py`
- `quantarhei/wizard/examples/ex_018_ModifiedRedfieldTheory_1.py`
- `quantarhei/wizard/examples/ex_200_Secular.py`
- `examples/demo_010_RedfieldTheory_1.py`
- `examples/demo_015_RedfieldTheory_2.py`
- `examples/demo_017_RedfieldTheory_MultiExcitons.py`
- `examples/demo_018_ModifiedRedfieldTheory_1.py`
- `examples/demo_200_Secular.py`
- `examples/old_demos/demo_redfield_2.py`
- `examples/old_demos/demo_redfield_3.py`

- [ ] **Step 8.1: Update each file**

For each file, replace the old protect pattern:

```python
# OLD — delete this block:
ham.protect_basis()
with qr.eigenbasis_of(ham):
    RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi)
    RRT.secularize()
ham.unprotect_basis()
```

with:

```python
# NEW:
RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi, secular=True)
```

For `eigenbasis_of` blocks that remain only for *reading out results or plotting* (not for tensor construction), leave them in place — those are valid presentation contexts.

- [ ] **Step 8.2: Run the relevant demo scripts to verify no runtime error**

```bash
cd /Users/danherma/projects-personal/quantarhei
python examples/demo_010_RedfieldTheory_1.py 2>&1 | tail -20
python examples/demo_018_ModifiedRedfieldTheory_1.py 2>&1 | tail -20
```

Expected: runs without exception.

- [ ] **Step 8.3: Commit**

```bash
git add quantarhei/wizard/examples/ examples/demo_010_RedfieldTheory_1.py examples/demo_015_RedfieldTheory_2.py examples/demo_017_RedfieldTheory_MultiExcitons.py examples/demo_018_ModifiedRedfieldTheory_1.py examples/demo_200_Secular.py examples/old_demos/demo_redfield_2.py examples/old_demos/demo_redfield_3.py
git commit -m "#325 docs: update examples to new basis management API"
```

---

## Task 8b: Note on `FoersterRateMatrix`, `LindbladForm`, and `enforce_detailed_balance`

**Files:**
- `quantarhei/qm/liouvillespace/foerstertensor.py` — currently does not use `eigenbasis_of` internally; add `basis_op = None` (Förster is computed in site basis, no eigenbasis needed)
- `quantarhei/qm/liouvillespace/lindblad.py` — same: add `basis_op = None` (Lindblad is site-basis)
- `quantarhei/qm/liouvillespace/relaxationtensor.py:152` — `enforce_detailed_balance` already calls `with eigenbasis_of(self.Hamiltonian):` internally; add `assert_not_in_eigenbasis_context()` at entry

- [ ] **Step 8b.1: Add `basis_op = None` to FoersterRateMatrix and LindbladForm**

In `foerstertensor.py` and `lindblad.py`, at the end of their `__init__` methods, add:

```python
self.basis_op = None  # computed in site basis; propagator needs no eigenbasis context
```

- [ ] **Step 8b.2: Guard `enforce_detailed_balance`**

In `relaxationtensor.py`, at the start of `enforce_detailed_balance` (line 152), add:

```python
from ...core.managers import assert_not_in_eigenbasis_context
assert_not_in_eigenbasis_context()
```

- [ ] **Step 8b.3: Run unit tests**

```bash
pytest tests/unit/qm/liouvillespace/ -q 2>&1 | tail -20
```

Expected: pass.

- [ ] **Step 8b.4: Commit**

```bash
git add quantarhei/qm/liouvillespace/foerstertensor.py quantarhei/qm/liouvillespace/lindblad.py quantarhei/qm/liouvillespace/relaxationtensor.py
git commit -m "#325 refactor: set basis_op on Foerster/Lindblad, guard enforce_detailed_balance"
```

---

## Task 9: Final verification

- [ ] **Step 9.1: Run the full unit test suite**

```bash
pytest tests/unit/ -q 2>&1 | tail -30
```

Expected: all tests pass, no regressions beyond those introduced intentionally.

- [ ] **Step 9.2: Run the behave acceptance tests**

```bash
cd /Users/danherma/projects-personal/quantarhei
behave tests/behave/ 2>&1 | tail -30
```

Expected: all scenarios pass.

- [ ] **Step 9.3: Run linting and formatting**

```bash
ruff check quantarhei/ --fix
ruff format quantarhei/
```

Expected: no errors remaining.

- [ ] **Step 9.4: Run mypy**

```bash
mypy quantarhei/core/managers.py quantarhei/utils/types.py quantarhei/qm/liouvillespace/redfieldtensor.py quantarhei/qm/liouvillespace/modredfieldtensor.py quantarhei/qm/propagators/rdmpropagator.py
```

Expected: no new type errors introduced.

- [ ] **Step 9.5: Final commit if any lint/mypy fixes were needed**

```bash
git add -p
git commit -m "#325 chore: lint and type fixes for basis management redesign"
```
