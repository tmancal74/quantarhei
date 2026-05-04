"""Regression tests for issue #325.

The secular Redfield propagator used to converge to near-equipartition instead
of the Boltzmann distribution because it used the non-diagonal site-basis
Hamiltonian in the commutator step.

Root cause (fixed)
------------------
``eigenbasis_of(H)`` tells *operator constructors* (RedfieldRelaxationTensor,
etc.) to build their tensors in the eigenbasis but does NOT transform
``H.data`` itself.  The propagator's ``_COM`` step used ``H.data`` directly,
including the off-diagonal coupling J, which drove persistent coherences and
prevented relaxation to the Boltzmann steady state.

Fix
---
``ReducedDensityMatrixPropagator._INIT_RWA`` now detects a secularized
relaxation tensor (``RelaxationTensor.is_secular``) and strips off-diagonal
elements from H before the commutator loop.  The legacy ``secularize()`` path
in ``RelaxationTensor`` was also fixed to set ``is_secular = True`` so the
detection works.

Additional note
---------------
Populations must be extracted from ``rhot.data`` *inside* the
``eigenbasis_of(H)`` context.  Accessing ``rhot.data`` after the context exits
triggers a basis back-transform that mixes populations and coherences, giving
incorrect ratios in the site basis.
"""

import numpy
import scipy.linalg as la

import quantarhei as qr


def _make_heterodimer_redfield(ta):
    cpar = dict(
        ftype="OverdampedBrownian",
        reorg=20.0,
        cortime=100.0,
        T=300.0,
        matsubara=20,
    )
    with qr.energy_units("1/cm"):
        cfce = qr.CorrelationFunction(ta, cpar)
        m1 = qr.Molecule([0.0, 12000.0])
        m2 = qr.Molecule([0.0, 12300.0])
    m1.set_transition_environment((0, 1), cfce)
    m2.set_transition_environment((0, 1), cfce)
    agg = qr.Aggregate(molecules=[m1, m2])
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0, 1, 100.0)
    agg.build()
    H = agg.get_Hamiltonian()
    sbi = agg.get_SystemBathInteraction()
    return H, sbi


class TestRedfieldBoltzmann:
    """#325 — secular Redfield steady-state population tests."""

    def test_rate_matrix_satisfies_detailed_balance(self):
        """Rate matrix satisfies detailed balance: k_up/k_down = exp(-dE/kT)."""
        ta = qr.TimeAxis(0.0, 1001, 2.0)
        H, sbi = _make_heterodimer_redfield(ta)

        H.protect_basis()
        with qr.eigenbasis_of(H):
            RRT = qr.qm.RedfieldRelaxationTensor(H, sbi)
            RRT.secularize()
            k_down = numpy.real(RRT.data[1, 1, 2, 2])
            k_up = numpy.real(RRT.data[2, 2, 1, 1])
        H.unprotect_basis()

        evals = sorted([qr.convert(e, "int", "1/cm") for e in la.eigvalsh(H.data)])
        dE = abs(evals[2] - evals[1])
        kB, T = 0.695, 300.0
        expected = numpy.exp(-dE / (kB * T))

        assert k_down > 0 and k_up > 0
        assert abs(k_up / k_down - expected) / expected < 0.05

    def test_propagator_converges_to_boltzmann(self):
        """Secular Redfield propagation should relax to the Boltzmann distribution."""
        ta = qr.TimeAxis(0.0, 6001, 2.0)
        H, sbi = _make_heterodimer_redfield(ta)

        H.protect_basis()
        with qr.eigenbasis_of(H):
            RRT = qr.qm.RedfieldRelaxationTensor(H, sbi)
            RRT.secularize()
            rho0 = qr.ReducedDensityMatrix(dim=H.dim)
            rho0.data[2, 2] = 1.0
            prop = qr.ReducedDensityMatrixPropagator(ta, H, RRT)
            rhot = prop.propagate(rho0)

            # Extract populations in the eigenbasis while still inside the context.
            # Accessing rhot.data outside the context triggers a basis back-transform
            # that mixes populations and coherences, giving incorrect ratios.
            n = rhot.data.shape[0]
            tail = max(1, n // 10)
            p1 = float(
                numpy.mean([numpy.real(rhot.data[k, 1, 1]) for k in range(n - tail, n)])
            )
            p2 = float(
                numpy.mean([numpy.real(rhot.data[k, 2, 2]) for k in range(n - tail, n)])
            )
        H.unprotect_basis()

        evals = sorted([qr.convert(e, "int", "1/cm") for e in la.eigvalsh(H.data)])
        dE = abs(evals[2] - evals[1])
        kB, T = 0.695, 300.0
        expected = numpy.exp(-dE / (kB * T))

        assert p1 > 0 and p2 > 0
        ratio = p2 / p1
        assert abs(ratio - expected) / expected < 0.05, (
            f"p2/p1={ratio:.4f}, expected {expected:.4f} (dE={dE:.1f} cm⁻¹, T={T} K)"
        )
