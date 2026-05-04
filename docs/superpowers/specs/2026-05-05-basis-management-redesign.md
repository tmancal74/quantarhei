# Basis Management Redesign — Design Spec

**Date:** 2026-05-05
**Issue:** #325 (root cause), related to tmancal74's Option 1 preference in PR #328
**Status:** Approved, pending implementation plan

---

## Problem

The current basis management system requires users to manually set up `eigenbasis_of` contexts and call `protect_basis()` / `unprotect_basis()` around tensor constructors. This is error-prone:

- Users can forget to protect the Hamiltonian, leaving it in an unexpected state
- Calling `secularize()` outside the correct context produces silently wrong results
- The propagator's commutator step uses `H.data` directly, trusting the caller's context — violating this trust causes PR #325's Boltzmann convergence bug
- The auto-transform-on-access mechanism in `.data` silently masks basis mismatches instead of surfacing them

---

## The New Invariant

> **Outside an `eigenbasis_of` context, every `BasisManaged` operator is always in the site basis. Inside a context, it is in the eigenbasis of the context operator. Accessing `.data` when the operator's recorded basis doesn't match the Manager's current basis raises `BasisError` immediately — no silent transform.**

Site basis is the single ground truth. Eigenbasis is always a scoped, temporary context — entered and exited by library code, never by user code for calculation purposes.

---

## Design

### 1. `BasisError` and Guard Functions

**New exception hierarchy:**
```python
class BasisError(Exception): ...
class BasisContextError(BasisError): ...  # wrong context for this operation
```

Location: `quantarhei/core/managers.py` or new `quantarhei/core/exceptions.py`.

**Two internal guard functions:**

`assert_not_in_eigenbasis_context()` — raises `BasisContextError` if `manager._in_eigenbasis_of_context` is `True`. Called at the top of every tensor constructor. Enforces the hard break: calling a tensor constructor inside a user-managed `eigenbasis_of` is forbidden.

`assert_in_eigenbasis_of(op)` — raises `BasisError` if not currently inside `eigenbasis_of(op)`. Available for internal use where eigenbasis is required.

**`.data` property change:**

`basis_managed_array_property` in `quantarhei/utils/types.py:87-132` — the silent auto-transform branch is replaced with a strict raise:

```python
if cb != ob:
    raise BasisError(
        f"Operator is in basis {ob} but current basis is {cb}. "
        "Access .data only inside the correct eigenbasis_of context."
    )
```

**`protect_basis()` / `unprotect_basis()` removed** from `BasisManaged` (`managers.py:922-926`). No operator ever silently stays diagonal outside a context.

---

### 2. Tensor Constructors Own Their Context

Every tensor constructor that requires eigenbasis computation follows this pattern:

```python
def __init__(self, ham, sbi, ..., secular=False):
    assert_not_in_eigenbasis_context()   # hard break for old protect pattern

    with eigenbasis_of(ham):
        # all existing computation — unchanged numerically
        ...
        if secular:
            self.secularize()            # secularization inside the same context

    # context exited — self.data and ham.data are back in site basis
    self.basis_op = ham                  # tensor remembers which H it was built from
```

**Affected constructors:**
- `RedfieldRelaxationTensor`
- `ModRedfieldRelaxationTensor`
- `FoersterRateMatrix`
- `LindbladForm`
- Any other tensor/rate-matrix constructor currently relying on an external `eigenbasis_of` context

**API change:** `secularize()` as a standalone post-construction call raises `BasisContextError` if called outside an eigenbasis context (which is now always the case for user code). Replaced by `secular=True` constructor flag.

**User-facing before/after:**
```python
# Before (old pattern — now raises BasisContextError)
ham.protect_basis()
with eigenbasis_of(ham):
    RRT = RedfieldRelaxationTensor(ham, sbi)
    RRT.secularize()
ham.unprotect_basis()

# After
RRT = RedfieldRelaxationTensor(ham, sbi, secular=True)
```

---

### 3. Propagator Owns Its Context

`ReducedDensityMatrixPropagator` gains an optional `basis` parameter:

```python
# Default: uses RRT.basis_op
prop = ReducedDensityMatrixPropagator(time, ham, RRT)

# Explicit override
prop = ReducedDensityMatrixPropagator(time, ham, RRT, basis=ham)
```

Internally, `propagate()` resolves the basis operator and wraps all propagation:

```python
def propagate(self, rho_i):
    basis_op = self._basis or getattr(self.RelaxTensor, "basis_op", None)
    if basis_op is not None:
        with eigenbasis_of(basis_op):
            # all existing _INIT_RWA, _COM, short-exp logic — unchanged
            ...
        # result back-transformed to site basis on context exit
    else:
        # no eigenbasis needed (e.g. Lindblad in site basis)
        ...
    return rhot
```

If `basis_op` is `None`, the tensor was built entirely in the site basis (e.g. a pure Lindblad form with no diagonalisation step) and the propagator runs without entering any context.

`_INIT_RWA` requires no `numpy.diag` workaround — inside `eigenbasis_of(ham)`, `H.data` is already diagonal by construction. PR #329's hack is removed.

**Result:** `rhot` is returned in site basis. Users no longer need to extract populations inside an `eigenbasis_of` context — `rhot.data` is directly readable.

---

### 4. Backward Compatibility

This is a **hard break**. The old pattern:

```python
ham.protect_basis()
with eigenbasis_of(ham):
    RRT = RedfieldRelaxationTensor(ham, sbi)
ham.unprotect_basis()
```

raises `BasisContextError` at the `RedfieldRelaxationTensor(...)` line.

The following files use `protect_basis` / `unprotect_basis` or call tensor constructors inside user-managed `eigenbasis_of` contexts and must be updated as part of this change:

**Library source:**
- `quantarhei/builders/opensystem.py` — 10 protect/unprotect call sites
- `quantarhei/wizard/examples/ex_010_RedfieldTheory_1.py`
- `quantarhei/wizard/examples/ex_018_ModifiedRedfieldTheory_1.py`
- `quantarhei/wizard/examples/ex_200_Secular.py`

**Examples:**
- `examples/demo_010_RedfieldTheory_1.py`
- `examples/demo_015_RedfieldTheory_2.py`
- `examples/demo_017_RedfieldTheory_MultiExcitons.py`
- `examples/demo_018_ModifiedRedfieldTheory_1.py`
- `examples/demo_200_Secular.py`
- `examples/old_demos/demo_redfield_2.py`
- `examples/old_demos/demo_redfield_3.py`

Other examples that use `eigenbasis_of` for *presentation* purposes only (reading out results, plotting) are unaffected — those uses remain valid and intentional.

---

## Testing Strategy

### Unit: guard behaviour (`tests/unit/core/test_basis_error.py`)
- `.data` access outside correct context raises `BasisError`
- Tensor constructor called inside `eigenbasis_of` raises `BasisContextError`
- `protect_basis()` call raises `AttributeError` (method removed)

### Unit: tensor constructors (extend `tests/unit/qm/liouvillespace/`)
- Each tensor constructor produces numerically identical results to the old pattern
- `basis_op` is set correctly after construction
- `secular=True` flag produces same result as old external `secularize()` call

### Regression: propagator convergence (`tests/unit/qm/propagators/test_redfield_boltzmann.py`)
- `test_propagator_converges_to_boltzmann` passes (no longer `xfail`)
- Populations readable directly from `rhot.data` without entering a context
- Boltzmann ratio correct within 5% after 12 ps

### Acceptance: existing behave tests (`tests/behave/`)
- All existing scenarios pass after examples are updated to new API

---

## Files Changed

| File | Change |
|---|---|
| `quantarhei/core/managers.py` | Add `BasisError`, `BasisContextError`, guard functions; remove `protect_basis`/`unprotect_basis` from `BasisManaged` |
| `quantarhei/utils/types.py` | Replace silent auto-transform in `.data` property with raise |
| `quantarhei/qm/liouvillespace/redfieldtensor.py` | Internalize `eigenbasis_of`, add `secular=` flag, set `basis_op` |
| `quantarhei/qm/liouvillespace/modredfieldtensor.py` | Same pattern |
| `quantarhei/qm/liouvillespace/foerstertensor.py` | Same pattern |
| `quantarhei/qm/liouvillespace/lindblad.py` | Same pattern |
| `quantarhei/qm/propagators/rdmpropagator.py` | Internalize `eigenbasis_of` in `propagate()`, add `basis=` param, remove `numpy.diag` hack |
| `quantarhei/builders/opensystem.py` | Remove `protect_basis` calls, use new tensor API |
| `examples/demo_010_RedfieldTheory_1.py` | Update to new API |
| (other examples using protect pattern) | Update to new API |
| `tests/unit/core/test_basis_error.py` | New: guard behaviour tests |
| `tests/unit/qm/propagators/test_redfield_boltzmann.py` | Update: remove `xfail`, assert convergence |
