"""Consumers of the Hamiltonian eigensystem must not depend on its lazy basis state.

Inside ``eigenbasis_of(H)`` the data of ``H`` stay in the site basis until
``H.data`` is read, after which ``H._data`` holds the (nearly diagonal)
transformed matrix. Relaxation tensors and rate matrices need the site-basis
eigenvectors of ``H``; they must obtain them independently of that state and,
where the result is labelled with a basis, express it in exactly the basis of
the context (issue #333).

All quantities compared here are computed by identical arithmetic in both
paths up to rounding in the basis transformations, so the tolerances are
~1e-12 relative to the magnitude of the tensors/rates (|R| ~ 1e-3..1e-2 in
internal units), unless stated otherwise.
"""

import numpy
import pytest

import quantarhei as qr
from quantarhei import TimeAxis, energy_units
from quantarhei.qm import (
    ModifiedRedfieldRateMatrix,
    ModRedfieldRelaxationTensor,
    RedfieldFoersterRelaxationTensor,
    RedfieldRateMatrix,
    RedfieldRelaxationTensor,
    TDRedfieldFoersterRelaxationTensor,
    TDRedfieldRateMatrix,
    TDRedfieldRelaxationTensor,
)

ATOL = 1e-12


def _aggregate(ta, energies, couplings):
    """Aggregate of two-level molecules with identical overdamped baths."""
    cpar = {
        "ftype": "OverdampedBrownian",
        "reorg": 20.0,
        "cortime": 100.0,
        "T": 300.0,
        "matsubara": 20,
    }
    with energy_units("1/cm"):
        cfce = qr.CorrelationFunction(ta, cpar)
        mols = [qr.Molecule([0.0, en]) for en in energies]
    for mol in mols:
        mol.set_transition_environment((0, 1), cfce)
    agg = qr.Aggregate(molecules=mols)
    with energy_units("1/cm"):
        for (ii, jj), val in couplings.items():
            agg.set_resonance_coupling(ii, jj, val)
    agg.build()
    return agg


def _ring_trimer(ta):
    """C3-symmetric ring: equal site energies and couplings, so the two
    upper excitons are exactly degenerate."""
    return _aggregate(ta, [12000.0] * 3, {(0, 1): 100.0, (1, 2): 100.0, (0, 2): 100.0})


def _trimer(ta):
    """Non-degenerate trimer with one weak (Foerster-like) coupling."""
    return _aggregate(
        ta, [12000.0, 12150.0, 12400.0], {(0, 1): 100.0, (1, 2): 80.0, (0, 2): 5.0}
    )


def _max_rel_diff(aa, bb):
    return numpy.max(numpy.abs(aa - bb)) / numpy.max(numpy.abs(bb))


# ---------------------------------------------------------------------------
# Relaxation tensors built inside eigenbasis_of(H)
# ---------------------------------------------------------------------------


def _tensor_in_eigenbasis(cls, agg, touch):
    ham = agg.get_Hamiltonian()
    sbi = agg.get_SystemBathInteraction()
    with qr.eigenbasis_of(ham):
        if touch:
            _ = ham.data  # lazily transforms H into its eigenbasis
        RT = cls(ham, sbi)
        return numpy.array(RT.data, copy=True)


@pytest.mark.parametrize("cls", [RedfieldRelaxationTensor, TDRedfieldRelaxationTensor])
def test_tensor_degenerate_ring_trimer_independent_of_lazy_state(cls):
    """Degenerate excitons: eigh of the transformed (noisy diagonal) H would
    rotate the degenerate pair (observed |dR| = 2.4e-4 on |R| = 3e-2). The
    tensor must be expressed in the context's own eigenvectors."""
    ta = TimeAxis(0.0, 500, 2.0)
    RR_ref = _tensor_in_eigenbasis(cls, _ring_trimer(ta), touch=False)
    RR = _tensor_in_eigenbasis(cls, _ring_trimer(ta), touch=True)
    assert numpy.all(numpy.isfinite(RR))
    numpy.testing.assert_allclose(RR, RR_ref, rtol=0, atol=ATOL)


def test_tensor_degenerate_ring_trimer_in_context_basis():
    """The tensor built in eigenbasis_of(H) equals the site-basis tensor
    (built outside and back-transformed) viewed in that same context."""
    ta = TimeAxis(0.0, 500, 2.0)
    agg = _ring_trimer(ta)
    ham = agg.get_Hamiltonian()
    sbi = agg.get_SystemBathInteraction()

    RT_site = RedfieldRelaxationTensor(ham, sbi)
    _, SS = ham.get_site_basis_eigensystem()
    RT_site.transform(numpy.linalg.inv(SS), inv=SS)

    with qr.eigenbasis_of(ham):
        _ = ham.data
        RT = RedfieldRelaxationTensor(ham, sbi)
        RR = numpy.array(RT.data, copy=True)
        RR_ref = numpy.array(RT_site.data, copy=True)

    numpy.testing.assert_allclose(RR, RR_ref, rtol=0, atol=ATOL)


# ---------------------------------------------------------------------------
# Rate matrices (always in the eigenbasis of H, ascending energies)
# ---------------------------------------------------------------------------


def _rates(cls, agg, context, touch):
    ham = agg.get_Hamiltonian()
    sbi = agg.get_SystemBathInteraction()
    if not context:
        return numpy.array(cls(ham, sbi).data, copy=True)
    with qr.eigenbasis_of(ham):
        if touch:
            _ = ham.data
        return numpy.array(cls(ham, sbi).data, copy=True)


@pytest.mark.parametrize(
    "cls, dt, touch",
    [
        # reading H.data first used to give an all-zero rate matrix
        (RedfieldRateMatrix, 2.0, True),
        # the frequency axis must cover the ground-exciton gap: dt = 1 fs
        (ModifiedRedfieldRateMatrix, 1.0, True),
        # eigh(ham.data) was wrong inside any eigenbasis_of(H), even untouched
        (TDRedfieldRateMatrix, 2.0, False),
        (TDRedfieldRateMatrix, 2.0, True),
    ],
)
def test_rate_matrix_in_eigenbasis_matches_outside(cls, dt, touch):
    ta = TimeAxis(0.0, 1000, dt)
    KK_ref = _rates(cls, _trimer(ta), context=False, touch=False)
    KK = _rates(cls, _trimer(ta), context=True, touch=touch)
    assert numpy.max(numpy.abs(KK_ref)) > 0.0
    assert _max_rel_diff(KK, KK_ref) < 1e-10


# ---------------------------------------------------------------------------
# Combined Redfield-Foerster tensors (Foerster part uses site-basis ham.JR)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cls", [RedfieldFoersterRelaxationTensor, TDRedfieldFoersterRelaxationTensor]
)
def test_redfield_foerster_in_eigenbasis_matches_outside(cls):
    """Inside eigenbasis_of(H), eigh(ham.data) returned ~identity, so the
    site-basis remainder coupling JR was not transformed to the excitons."""
    ta = TimeAxis(0.0, 500, 2.0)

    def build(context):
        agg = _trimer(ta)
        ham = agg.get_Hamiltonian()
        sbi = agg.get_SystemBathInteraction()
        with energy_units("1/cm"):
            ham.remove_cutoff_coupling(20.0)  # moves the 5 cm^-1 coupling to JR
        if not context:
            return numpy.array(cls(ham, sbi).data, copy=True)
        with qr.eigenbasis_of(ham):
            return numpy.array(cls(ham, sbi).data, copy=True)

    RR_ref = build(context=False)
    RR = build(context=True)
    numpy.testing.assert_allclose(RR, RR_ref, rtol=0, atol=ATOL)


# ---------------------------------------------------------------------------
# OpenSystem.get_RelaxationTensor
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("make_agg", [_trimer, _ring_trimer])
def test_opensystem_redfield_in_eigenbasis_context(make_agg):
    """Called inside eigenbasis_of(H), the returned tensor must be in that
    context's basis. Building the tensor reads H.data, after which the old
    back-transformation used eigh of the transformed H: ~identity (right by
    accident) without degeneracy, a rotation of the degenerate pair (wrong)
    for the ring trimer."""
    ta = TimeAxis(0.0, 1000, 2.0)
    agg = make_agg(ta)
    ham = agg.get_Hamiltonian()

    RT_out, _ = agg.get_RelaxationTensor(ta, relaxation_theory="stR")
    with qr.eigenbasis_of(ham):
        RT_in, _ = agg.get_RelaxationTensor(ta, relaxation_theory="stR")
        RR_in = numpy.array(RT_in.data, copy=True)
        RR_out = numpy.array(RT_out.data, copy=True)

    numpy.testing.assert_allclose(RR_in, RR_out, rtol=0, atol=ATOL)


def test_opensystem_modified_redfield_is_not_transformed_twice():
    """ModRedfieldRelaxationTensor is already in the site basis when it is
    returned (it fills its data inside its own eigenbasis_of context), so
    get_RelaxationTensor must not back-transform it once more."""
    ta = TimeAxis(0.0, 1000, 1.0)
    agg = _trimer(ta)
    ham = agg.get_Hamiltonian()
    sbi = agg.get_SystemBathInteraction()

    KK = ModifiedRedfieldRateMatrix(ham, sbi).data
    RT_direct = ModRedfieldRelaxationTensor(ham, sbi)
    RT, _ = agg.get_RelaxationTensor(ta, relaxation_theory="mR")

    with qr.eigenbasis_of(ham):
        RR = numpy.array(RT.data, copy=True)
        RR_direct = numpy.array(RT_direct.data, copy=True)

    numpy.testing.assert_allclose(RR, RR_direct, rtol=0, atol=ATOL)
    dim = ham.dim
    for aa in range(dim):
        for bb in range(dim):
            if aa != bb:
                assert abs(RR[aa, aa, bb, bb] - KK[aa, bb]) < ATOL
