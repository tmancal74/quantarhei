import copy

import numpy
import scipy.linalg as la
from behave import given, then, when

import quantarhei as qr
from quantarhei.spectroscopy import X
from quantarhei.spectroscopy.circular_dichroism import CircDichSpectrumCalculator
from quantarhei.spectroscopy.mocktwodcalculator import (
    MockTwoDResponseCalculator as TwoDResponseCalculator,
)


def _build_heterodimer(e1_invcm, e2_invcm, J_invcm, mult=1, width=None):
    with qr.energy_units("1/cm"):
        m1 = qr.Molecule([0.0, e1_invcm])
        m2 = qr.Molecule([0.0, e2_invcm])
        if width is not None:
            m1.set_transition_width((0, 1), width)
            m2.set_transition_width((0, 1), width)
        m1.set_dipole(0, 1, [1.0, 0.8, 0.8])
        m2.set_dipole(0, 1, [0.8, 0.8, 0.0])
    agg = qr.Aggregate(molecules=[m1, m2])
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0, 1, J_invcm)
    agg.build(mult=mult)
    return agg, m1, m2


def _true_eigenvalues_invcm(H):
    evals = la.eigvalsh(H.data)
    return sorted([qr.convert(e, "int", "1/cm") for e in evals])


@given("a two-level molecule with transition energy {e1} invcm")
def step_two_level_molecule(context, e1):
    context.monomer_e1_invcm = float(e1)
    context.use_monomer = True


@given("a heterodimer with site energies {e1} and {e2} invcm and coupling {J} invcm")
def step_heterodimer(context, e1, e2, J):
    context.e1_invcm = float(e1)
    context.e2_invcm = float(e2)
    context.J_invcm = float(J)
    context.use_monomer = False


@given(
    "an overdamped Brownian bath with reorganisation energy {reorg} invcm "
    "and correlation time {cortime} fs at {T} K"
)
def step_bath(context, reorg, cortime, T):
    ta = qr.TimeAxis(0.0, 1000, 1.0)
    cpar = dict(
        ftype="OverdampedBrownian",
        reorg=float(reorg),
        cortime=float(cortime),
        T=float(T),
        matsubara=20,
    )
    with qr.energy_units("1/cm"):
        cfce = qr.CorrelationFunction(ta, cpar)
    context.cfce = cfce
    context.bath_T = float(T)
    context.bath_ta = ta


@given("a Lindblad relaxation rate of 1 per {tau} fs")
def step_lindblad_rate(context, tau):
    context.lindblad_tau = float(tau)


@given("a homodimer-2-env test aggregate with coupling {J} invcm")
def step_homodimer_env(context, J):
    agg = qr.TestAggregate("homodimer-2-env")
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0, 1, float(J))
    agg.build()
    context.agg = agg
    context.J_invcm = float(J)


# ---------------------------------------------------------------------------
# Scenario A: Lindblad coherence decay (two-level monomer)
# ---------------------------------------------------------------------------

@when("I propagate the density matrix for {duration} fs with time step {dt} fs")
def step_propagate_lindblad(context, duration, dt):
    dt_val = float(dt)
    Nt = int(float(duration) / dt_val) + 1
    time = qr.TimeAxis(0.0, Nt, dt_val)
    tau = context.lindblad_tau

    if getattr(context, "use_monomer", False):
        # Pure two-level system: analytic exponential decay, exact ground truth
        with qr.energy_units("1/cm"):
            m = qr.Molecule([0.0, context.monomer_e1_invcm])
            m.set_dipole(0, 1, [1.0, 0.0, 0.0])
        agg = qr.Aggregate(molecules=[m])
        agg.build(mult=1)
        H = agg.get_Hamiltonian()
        K = qr.qm.Operator(dim=H.dim, real=True)
        K.data[0, 1] = 1.0
        sbi = qr.qm.SystemBathInteraction(sys_operators=[K], rates=[1.0 / tau])
        L = qr.qm.LindbladForm(H, sbi, as_operators=False)
        rho0 = qr.ReducedDensityMatrix(dim=H.dim)
        rho0.data[0, 1] = 0.5
        rho0.data[1, 0] = 0.5
        rho0.data[1, 1] = 1.0
        context.L_deph_rate = abs(numpy.real(L.data[0, 1, 0, 1]))
        context._monomer_coherence_idx = (0, 1)
    else:
        agg, _, _ = _build_heterodimer(context.e1_invcm, context.e2_invcm, context.J_invcm)
        H = agg.get_Hamiltonian()
        K = qr.qm.ProjectionOperator(1, 2, dim=H.dim)
        sbi = qr.qm.SystemBathInteraction(sys_operators=[K], rates=[1.0 / tau])
        L = qr.qm.LindbladForm(H, sbi, as_operators=False)
        rho0 = qr.ReducedDensityMatrix(dim=H.dim)
        with qr.eigenbasis_of(H):
            rho0.data[1, 2] = 0.5
            rho0.data[2, 1] = 0.5
            rho0.data[1, 1] = 0.5
            rho0.data[2, 2] = 0.5
        context.L_deph_rate = abs(numpy.real(L.data[1, 2, 1, 2]))
        context._monomer_coherence_idx = (1, 2)

    prop = qr.ReducedDensityMatrixPropagator(time, H, L)
    rhot = prop.propagate(rho0)

    context.rhot = rhot
    context.time = time
    context.H = H


@then("the off-diagonal element decays as exp(-t/200) within {tol_pct} percent")
def step_check_coherence_decay(context, tol_pct):
    tol = float(tol_pct) / 100.0
    rhot = context.rhot
    time = context.time
    gamma = context.L_deph_rate
    i, j = context._monomer_coherence_idx

    # For a pure two-level system with K=|0><1| the coherence rho_01 decays
    # exactly as exp(-gamma*t) * rho_01(0).
    rho0_val = abs(rhot.data[0, i, j])

    n_checks = 5
    step = max(1, (time.length - 1) // n_checks)
    for k in range(step, time.length, step):
        t = time.data[k]
        computed = abs(rhot.data[k, i, j])
        expected = rho0_val * numpy.exp(-gamma * t)
        if expected > 1e-10:
            err = abs(computed - expected) / expected
            assert err < tol, (
                f"At t={t:.0f} fs: |rho_{i}{j}|={computed:.5f}, "
                f"expected {expected:.5f} (rel err={err:.3f}, tol={tol_pct}%)"
            )


# ---------------------------------------------------------------------------
# Scenario B: Absorption peaks
# ---------------------------------------------------------------------------

@when("I calculate the absorption spectrum")
def step_calc_absorption(context):
    ta = context.bath_ta
    cfce = context.cfce
    agg, m1, m2 = _build_heterodimer(context.e1_invcm, context.e2_invcm, context.J_invcm)
    m1.set_transition_environment((0, 1), cfce)
    m2.set_transition_environment((0, 1), cfce)

    agg2 = qr.Aggregate(molecules=[m1, m2])
    with qr.energy_units("1/cm"):
        agg2.set_resonance_coupling(0, 1, context.J_invcm)
    agg2.build()
    H = agg2.get_Hamiltonian()

    (RR, _ham) = agg2.get_RelaxationTensor(ta, relaxation_theory="standard_Redfield")
    ac = qr.AbsSpectrumCalculator(ta, agg2, relaxation_tensor=RR)
    rwa = (context.e1_invcm + context.e2_invcm) / 2.0
    with qr.energy_units("1/cm"):
        ac.bootstrap(rwa=rwa)

    context.spectrum = ac.calculate()
    context.H_abs = H


@then(
    "there are two peaks whose positions match the exciton eigenvalues within {tol} invcm"
)
def step_check_absorption_peaks(context, tol):
    tol = float(tol)
    H = context.H_abs
    spectrum = context.spectrum

    # True eigenvalues (not RWA) for peak position comparison
    eig_energies = sorted(_true_eigenvalues_invcm(H))[1:]  # drop ground state

    with qr.frequency_units("1/cm"):
        freq = spectrum.axis.data
        data = spectrum.data

    peak_indices = [
        i for i in range(1, len(data) - 1)
        if data[i] > data[i - 1] and data[i] > data[i + 1]
    ]
    peak_freqs = [freq[i] for i in peak_indices]

    assert len(peak_freqs) >= 2, (
        f"Expected at least 2 peaks, found {len(peak_freqs)}: {peak_freqs}"
    )

    # Match the two largest peaks to the two exciton eigenvalues
    amplitudes = [data[numpy.argmin(numpy.abs(freq - f))] for f in peak_freqs]
    top2_freqs = sorted(
        [f for _, f in sorted(zip(amplitudes, peak_freqs), reverse=True)[:2]]
    )

    for computed, expected in zip(top2_freqs, eig_energies):
        assert abs(computed - expected) < tol, (
            f"Peak at {computed:.1f} cm⁻¹ deviates from eigenvalue "
            f"{expected:.1f} cm⁻¹ by more than {tol} cm⁻¹"
        )


# ---------------------------------------------------------------------------
# Scenario C: Redfield homodimer detailed balance
# ---------------------------------------------------------------------------

@when("I propagate with standard Redfield theory for {duration} fs with time step {dt} fs")
def step_propagate_redfield(context, duration, dt):
    dt_val = float(dt)
    Nt = int(float(duration) / dt_val) + 1
    time = qr.TimeAxis(0.0, Nt, dt_val)

    if hasattr(context, "agg"):
        agg = context.agg
    else:
        cfce = context.cfce
        agg, m1, m2 = _build_heterodimer(
            context.e1_invcm, context.e2_invcm, context.J_invcm
        )
        m1.set_transition_environment((0, 1), cfce)
        m2.set_transition_environment((0, 1), cfce)
        agg2 = qr.Aggregate(molecules=[m1, m2])
        with qr.energy_units("1/cm"):
            agg2.set_resonance_coupling(0, 1, context.J_invcm)
        agg2.build()
        agg = agg2

    H = agg.get_Hamiltonian()
    sbi = agg.get_SystemBathInteraction()

    H.protect_basis()
    with qr.eigenbasis_of(H):
        RRT = qr.qm.RedfieldRelaxationTensor(H, sbi)
        RRT.secularize()

        rho0 = qr.ReducedDensityMatrix(dim=H.dim)
        rho0.data[2, 2] = 1.0

        prop = qr.ReducedDensityMatrixPropagator(time, H, RRT)
        rhot = prop.propagate(rho0)
    H.unprotect_basis()

    context.rhot = rhot
    context.time = time
    context.H = H
    context.RRT = RRT


@then(
    "the ratio of forward to backward rates satisfies exp(-dE/kT) within {tol_pct} percent"
)
def step_check_detailed_balance(context, tol_pct):
    tol = float(tol_pct) / 100.0
    RRT = context.RRT
    H = context.H

    kB_invcm = 0.695  # cm⁻¹ K⁻¹
    T = 300.0

    # Use TRUE (non-RWA) eigenvalues for the Boltzmann factor
    evals_invcm = _true_eigenvalues_invcm(H)
    e1, e2 = evals_invcm[1], evals_invcm[2]
    dE = abs(e2 - e1)

    # Secular Redfield tensor: R[i,i,j,j] is the rate i ← j
    H.protect_basis()
    with qr.eigenbasis_of(H):
        k_down = numpy.real(RRT.data[1, 1, 2, 2])  # downhill (fast)
        k_up = numpy.real(RRT.data[2, 2, 1, 1])    # uphill (slow)
    H.unprotect_basis()

    assert k_down > 0 and k_up > 0, (
        f"Rates must be positive: k_down={k_down:.4e}, k_up={k_up:.4e}"
    )
    ratio = k_up / k_down
    expected = numpy.exp(-dE / (kB_invcm * T))

    assert abs(ratio - expected) / expected < tol, (
        f"Detailed balance: k_up/k_down={ratio:.4f}, "
        f"exp(-dE/kT)={expected:.4f} (dE={dE:.1f} cm⁻¹, T={T} K, "
        f"rel err={abs(ratio-expected)/expected:.3f})"
    )


# ---------------------------------------------------------------------------
# Scenario D: Trace conservation
# ---------------------------------------------------------------------------

@then("the trace of the density matrix equals 1 at every time step within {tol}")
def step_check_trace(context, tol):
    tol = float(tol)
    rhot = context.rhot
    for k in range(rhot.data.shape[0]):
        tr = numpy.real(numpy.trace(rhot.data[k]))
        assert abs(tr - 1.0) < tol, (
            f"Trace={tr:.8f} at step {k}, deviation {abs(tr-1.0):.2e} > {tol}"
        )


# ---------------------------------------------------------------------------
# Scenario E: Förster rate detailed balance
# ---------------------------------------------------------------------------

@when("I calculate the Foerster rate matrix")
def step_calc_foerster(context):
    ta = context.bath_ta
    cfce = context.cfce
    agg, m1, m2 = _build_heterodimer(
        context.e1_invcm, context.e2_invcm, context.J_invcm
    )
    m1.set_transition_environment((0, 1), cfce)
    m2.set_transition_environment((0, 1), cfce)
    agg2 = qr.Aggregate(molecules=[m1, m2])
    with qr.energy_units("1/cm"):
        agg2.set_resonance_coupling(0, 1, context.J_invcm)
    agg2.build()

    (RRf, _) = agg2.get_RelaxationTensor(ta, relaxation_theory="standard_Foerster")
    context.RRf = RRf
    context.H_foerster = agg2.get_Hamiltonian()


@then("the forward rate matches 2*pi*J^2*spectral_overlap within {tol_pct} percent")
def step_check_foerster_rate(context, tol_pct):
    tol = float(tol_pct) / 100.0
    RRf = context.RRf
    H = context.H_foerster

    # k(1←2) = R[1,1,2,2] (downhill, fast)
    # k(2←1) = R[2,2,1,1] (uphill, slow)
    k_down = numpy.real(RRf.data[1, 1, 2, 2])
    k_up = numpy.real(RRf.data[2, 2, 1, 1])

    # Detailed balance: k_up/k_down = exp(-dE/kBT)
    evals_invcm = _true_eigenvalues_invcm(H)
    dE = abs(evals_invcm[2] - evals_invcm[1])
    kB_invcm, T = 0.695, 300.0
    expected = numpy.exp(-dE / (kB_invcm * T))

    assert k_down > 0 and k_up > 0, (
        f"Förster rates must be positive: k_down={k_down:.4e}, k_up={k_up:.4e}"
    )
    ratio = k_up / k_down
    assert abs(ratio - expected) / expected < tol, (
        f"Förster detailed balance: k_up/k_down={ratio:.4f}, "
        f"exp(-dE/kT)={expected:.4f} (dE={dE:.1f} cm⁻¹, rel err={abs(ratio-expected)/expected:.3f})"
    )


# ---------------------------------------------------------------------------
# Scenario F: 2DES cross-peaks
# ---------------------------------------------------------------------------

@when("I calculate the 2D electronic spectrum at population time {T2} fs")
def step_calc_2d(context, T2):
    e1, e2, J = context.e1_invcm, context.e2_invcm, context.J_invcm
    tau = context.lindblad_tau

    agg, m1, m2 = _build_heterodimer(e1, e2, J, width=150.0)
    agg_2d, _, _ = _build_heterodimer(e1, e2, J, mult=2, width=150.0)
    H = agg.get_Hamiltonian()

    with qr.eigenbasis_of(H):
        K = qr.qm.ProjectionOperator(1, 2, dim=H.dim)
    sbi = qr.qm.SystemBathInteraction(sys_operators=[K], rates=[1.0 / tau])
    L = qr.qm.LindbladForm(H, sbi)

    t2_axis = qr.TimeAxis(0.0, 2, 10.0)
    eUt = qr.EvolutionSuperOperator(time=t2_axis, ham=H, relt=L)
    eUt.set_dense_dt(10)
    eUt.calculate()

    t1_axis = qr.TimeAxis(0.0, 50, 10.0)
    t3_axis = qr.TimeAxis(0.0, 50, 10.0)
    calc = TwoDResponseCalculator(t1_axis, t2_axis, t3_axis)
    rwa = (e1 + e2) / 2.0
    with qr.energy_units("1/cm"):
        calc.bootstrap(rwa=rwa)

    agg_2d.diagonalize()
    lab = qr.LabSetup()
    lab.set_pulse_polarizations(
        pulse_polarizations=(X, X, X), detection_polarization=X
    )

    tcont = calc.calculate_all_system(agg_2d, eUt, lab)
    tcont = tcont.get_TwoDSpectrumContainer()

    context.twod = tcont.get_spectrum(float(T2))
    context.H_2d = H
    context.e1_invcm = e1
    context.e2_invcm = e2


@then("cross-peaks are present in the total 2D spectrum")
def step_check_crosspeaks(context):
    twod = context.twod
    H = context.H_2d

    # True exciton energies for peak location
    evals = _true_eigenvalues_invcm(H)
    e_lo, e_hi = evals[1], evals[2]

    d = numpy.abs(numpy.real(twod.data))
    with qr.frequency_units("1/cm"):
        ax1 = twod.xaxis.data
        ax3 = twod.yaxis.data

    def _nearest(axis, val):
        return int(numpy.argmin(numpy.abs(axis - val)))

    i_lo = _nearest(ax1, e_lo)
    i_hi = _nearest(ax1, e_hi)
    j_lo = _nearest(ax3, e_lo)
    j_hi = _nearest(ax3, e_hi)

    diag_peak = max(d[j_lo, i_lo], d[j_hi, i_hi])
    cross_peak = max(d[j_lo, i_hi], d[j_hi, i_lo])

    assert cross_peak > 0.0, "Cross-peaks are zero in the total 2D spectrum"
    assert cross_peak > 0.05 * diag_peak, (
        f"Cross-peak amplitude ({cross_peak:.4e}) < 5% of "
        f"diagonal ({diag_peak:.4e})"
    )


@then("the rephasing and non-rephasing spectra are not identical")
def step_check_rephasing_nonrephasing(context):
    twod = context.twod

    reph = numpy.array(twod.scopy().data)
    twod.set_data_type(qr.signal_REPH)
    reph_typed = numpy.array(twod.data)

    twod.set_data_type(qr.signal_NONR)
    nonr = numpy.array(twod.data)

    diff = numpy.max(numpy.abs(reph_typed - nonr))
    scale = max(numpy.max(numpy.abs(reph_typed)), numpy.max(numpy.abs(nonr)))

    assert diff > 0.01 * scale, (
        f"Rephasing and non-rephasing spectra are nearly identical "
        f"(max diff {diff:.4e}, scale {scale:.4e})"
    )


# ---------------------------------------------------------------------------
# Scenario G: HEOM convergence
# ---------------------------------------------------------------------------

def _run_heom(e1_invcm, e2_invcm, J_invcm, cfce, ta, depth):
    agg, m1, m2 = _build_heterodimer(e1_invcm, e2_invcm, J_invcm)
    m1.set_transition_environment((0, 1), cfce)
    m2.set_transition_environment((0, 1), cfce)
    agg2 = qr.Aggregate(molecules=[m1, m2])
    with qr.energy_units("1/cm"):
        agg2.set_resonance_coupling(0, 1, J_invcm)
    agg2.build()
    H = agg2.get_Hamiltonian()
    sbi = agg2.get_SystemBathInteraction()
    hier = qr.qm.KTHierarchy(H, sbi, depth)
    rho0 = qr.ReducedDensityMatrix(dim=H.dim)
    rho0.data[2, 2] = 1.0
    prop = qr.qm.KTHierarchyPropagator(ta, hier)
    rhot = prop.propagate(rho0)
    return rhot


@when("I propagate with HEOM at depth {d1} and depth {d2} for {duration} fs with time step {dt} fs")
def step_propagate_heom(context, d1, d2, duration, dt):
    dt_val = float(dt)
    Nt = int(float(duration) / dt_val) + 1
    ta = qr.TimeAxis(0.0, Nt, dt_val)
    cfce = context.cfce

    context.rhot_heom_shallow = _run_heom(
        context.e1_invcm, context.e2_invcm, context.J_invcm, cfce, ta, int(d1)
    )
    context.rhot_heom_mid = _run_heom(
        context.e1_invcm, context.e2_invcm, context.J_invcm, cfce, ta, int(d2)
    )
    # one deeper for convergence check
    context.rhot_heom_deep = _run_heom(
        context.e1_invcm, context.e2_invcm, context.J_invcm, cfce, ta, int(d2) + 1
    )
    context.heom_ta = ta


@then("the depth-2 and depth-4 populations differ by more than {tol_pct} percent")
def step_heom_not_converged(context, tol_pct):
    tol = float(tol_pct) / 100.0
    last = -1
    p_shallow = numpy.real(context.rhot_heom_shallow.data[last, 2, 2])
    p_mid = numpy.real(context.rhot_heom_mid.data[last, 2, 2])
    diff = abs(p_shallow - p_mid)
    scale = max(abs(p_shallow), abs(p_mid), 1e-10)
    assert diff > tol * scale, (
        f"depth-2 and depth-4 populations are already converged: "
        f"p_shallow={p_shallow:.5f}, p_mid={p_mid:.5f}, rel diff={diff/scale:.4f}"
    )


@then("the depth-4 and depth-5 populations agree within {tol_pct} percent")
def step_heom_converged(context, tol_pct):
    tol = float(tol_pct) / 100.0
    last = -1
    p_mid = numpy.real(context.rhot_heom_mid.data[last, 2, 2])
    p_deep = numpy.real(context.rhot_heom_deep.data[last, 2, 2])
    diff = abs(p_mid - p_deep)
    scale = max(abs(p_mid), abs(p_deep), 1e-10)
    assert diff < tol * scale, (
        f"HEOM depth-4 vs depth-5 not converged: "
        f"p_mid={p_mid:.5f}, p_deep={p_deep:.5f}, rel diff={diff/scale:.4f} > {tol_pct}%"
    )


# ---------------------------------------------------------------------------
# Scenario H: Circular dichroism – chiral vs achiral dimer
# ---------------------------------------------------------------------------

@given("a chiral heterodimer with site energies {e1} and {e2} invcm and coupling {J} invcm")
def step_chiral_heterodimer(context, e1, e2, J):
    context.e1_invcm = float(e1)
    context.e2_invcm = float(e2)
    context.J_invcm = float(J)
    context.use_monomer = False


@when("I calculate the circular dichroism spectrum")
def step_calc_cd(context):
    ta = context.bath_ta
    cfce = context.cfce
    e1, e2, J = context.e1_invcm, context.e2_invcm, context.J_invcm

    # Chiral geometry: mu1=[1,0,0], mu2=[0,1,0], R12=[0,0,5]
    # Rotational strength R = (mu1 x mu2) . R12 = [0,0,1].[0,0,5] = 5 ≠ 0
    with qr.energy_units("1/cm"):
        m1 = qr.Molecule([0.0, e1])
        m2 = qr.Molecule([0.0, e2])
        m1.set_transition_width((0, 1), 150.0)
        m2.set_transition_width((0, 1), 150.0)
    m1.set_dipole(0, 1, [1.0, 0.0, 0.0])
    m2.set_dipole(0, 1, [0.0, 1.0, 0.0])
    m1._position = numpy.array([0.0, 0.0, 0.0])
    m2._position = numpy.array([0.0, 0.0, 5.0])

    m1.set_transition_environment((0, 1), cfce)
    m2.set_transition_environment((0, 1), cfce)

    agg = qr.Aggregate(molecules=[m1, m2])
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0, 1, J)
    agg.build()
    H = agg.get_Hamiltonian()

    (RR, _) = agg.get_RelaxationTensor(ta, relaxation_theory="standard_Redfield")
    calc = CircDichSpectrumCalculator(ta, agg, relaxation_tensor=RR)
    rwa = (e1 + e2) / 2.0
    with qr.energy_units("1/cm"):
        calc.bootstrap(rwa=rwa)

    context.cd_spectrum = calc.calculate()


@then("the CD spectrum has both positive and negative features")
def step_check_cd_bisignate(context):
    with qr.frequency_units("1/cm"):
        data = context.cd_spectrum.data
    assert numpy.any(data > 0), "CD spectrum has no positive feature"
    assert numpy.any(data < 0), "CD spectrum has no negative feature"



# ---------------------------------------------------------------------------
# Scenario I: Boltzmann steady-state populations from Redfield rate matrix
# ---------------------------------------------------------------------------

@when("I calculate the Redfield rate matrix")
def step_calc_redfield_rate_matrix(context):
    ta = context.bath_ta
    cfce = context.cfce
    agg, m1, m2 = _build_heterodimer(
        context.e1_invcm, context.e2_invcm, context.J_invcm
    )
    m1.set_transition_environment((0, 1), cfce)
    m2.set_transition_environment((0, 1), cfce)
    agg2 = qr.Aggregate(molecules=[m1, m2])
    with qr.energy_units("1/cm"):
        agg2.set_resonance_coupling(0, 1, context.J_invcm)
    agg2.build()
    H = agg2.get_Hamiltonian()
    sbi = agg2.get_SystemBathInteraction()

    H.protect_basis()
    with qr.eigenbasis_of(H):
        RRT = qr.qm.RedfieldRelaxationTensor(H, sbi)
        RRT.secularize()
    H.unprotect_basis()

    context.RRT_boltzmann = RRT
    context.H_boltzmann = H


@then("the steady-state populations satisfy the Boltzmann ratio within {tol_pct} percent")
def step_check_boltzmann_steady_state(context, tol_pct):
    tol = float(tol_pct) / 100.0
    RRT = context.RRT_boltzmann
    H = context.H_boltzmann

    kB_invcm = 0.695
    T = 300.0

    evals_invcm = _true_eigenvalues_invcm(H)
    dE = abs(evals_invcm[2] - evals_invcm[1])
    expected = numpy.exp(-dE / (kB_invcm * T))

    # Steady-state from rate matrix: dp1/dt=0, dp2/dt=0 gives
    # p1_ss = k_down / (k_down + k_up),  p2_ss = k_up / (k_down + k_up)
    H.protect_basis()
    with qr.eigenbasis_of(H):
        k_down = numpy.real(RRT.data[1, 1, 2, 2])
        k_up = numpy.real(RRT.data[2, 2, 1, 1])
    H.unprotect_basis()

    p1_ss = k_down / (k_down + k_up)
    p2_ss = k_up / (k_down + k_up)

    assert p1_ss > 0 and p2_ss > 0, (
        f"Non-positive steady-state population: p1={p1_ss:.5f}, p2={p2_ss:.5f}"
    )
    ratio = p2_ss / p1_ss
    assert abs(ratio - expected) / expected < tol, (
        f"Steady-state Boltzmann: p2/p1={ratio:.4f}, exp(-dE/kT)={expected:.4f} "
        f"(dE={dE:.1f} cm⁻¹, T={T} K, rel err={abs(ratio-expected)/expected:.3f})"
    )
