import numpy

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
                    {
                        "ftype": "OverdampedBrownian",
                        "reorg": 100.0,
                        "cortime": 100.0,
                        "T": 300.0,
                        "matsubara": 20,
                    },
                ),
            )
        agg = qr.Aggregate([mol1, mol2])
        agg.set_resonance_coupling(0, 1, 100.0)
        agg.build()
        ham = agg.get_Hamiltonian()
        sbi = agg.get_SystemBathInteraction()

    RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi, secular=True)
    RRM = qr.qm.RedfieldRateMatrix(ham, sbi)

    # RRM.data is a plain numpy array already in eigenbasis.
    # Indices: 0=ground, 1=lower exciton, 2=upper exciton
    # k_down: rate from upper (2) -> lower (1), stored at [1, 2]
    # k_up:   rate from lower (1) -> upper (2), stored at [2, 1]
    k_down = RRM.data[1, 2]
    k_up = RRM.data[2, 1]

    ratio = k_up / k_down
    # Compute expected Boltzmann ratio from actual eigenenergies
    SS = ham.get_diagonalization_matrix()
    S1 = numpy.linalg.inv(SS)
    H_eig = S1 @ ham._data @ SS
    dE_int = numpy.real(H_eig[2, 2] - H_eig[1, 1])
    cm_to_int = 2.0 * numpy.pi * 3.0e10 * 1.0e-15
    dE_cm = dE_int / cm_to_int
    expected = _boltzmann_ratio(dE_cm)
    assert abs(ratio - expected) / expected < 0.05, (
        f"Rate ratio {ratio:.4f} deviates >5% from Boltzmann {expected:.4f}"
    )


def test_rate_matrix_populations_sum_to_one():
    """Test that rate matrix elements conserve total population (column sums = 0)."""
    with energy_units("1/cm"):
        mol1 = qr.Molecule([0.0, 12000.0])
        mol2 = qr.Molecule([0.0, 12360.0])
        ta = TimeAxis(0.0, 1000, 1.0)
        for mol in (mol1, mol2):
            mol.set_transition_environment(
                (0, 1),
                qr.CorrelationFunction(
                    ta,
                    {
                        "ftype": "OverdampedBrownian",
                        "reorg": 100.0,
                        "cortime": 100.0,
                        "T": 300.0,
                        "matsubara": 20,
                    },
                ),
            )
        agg = qr.Aggregate([mol1, mol2])
        agg.set_resonance_coupling(0, 1, 100.0)
        agg.build()
        ham = agg.get_Hamiltonian()
        sbi = agg.get_SystemBathInteraction()

    RRM = qr.qm.RedfieldRateMatrix(ham, sbi)
    # Each column should sum to zero (population conservation)
    col_sums = numpy.sum(RRM.data, axis=0)
    numpy.testing.assert_allclose(col_sums, 0.0, atol=1e-10)


def _heterodimer(ta):
    """Dimer with an exciton gap of ~360 cm^-1 (site gap 300, J = 100 cm^-1)."""
    cpar = {
        "ftype": "OverdampedBrownian",
        "reorg": 20.0,
        "cortime": 100.0,
        "T": 300.0,
        "matsubara": 20,
    }
    with energy_units("1/cm"):
        cfce = qr.CorrelationFunction(ta, cpar)
        mol1 = qr.Molecule([0.0, 12000.0])
        mol2 = qr.Molecule([0.0, 12300.0])
    mol1.set_transition_environment((0, 1), cfce)
    mol2.set_transition_environment((0, 1), cfce)
    agg = qr.Aggregate(molecules=[mol1, mol2])
    with energy_units("1/cm"):
        agg.set_resonance_coupling(0, 1, 100.0)
    agg.build()
    return agg.get_Hamiltonian(), agg.get_SystemBathInteraction()


def _redfield_in_eigenbasis(touch_hamiltonian):
    """Build the Redfield tensor inside ``eigenbasis_of(H)``.

    If ``touch_hamiltonian`` is True, ``H.data`` is read before the tensor is
    constructed, which lazily transforms H into its eigenbasis (issue #333).
    Returns H, the tensor, and copies of their data in the eigenbasis.
    """
    ta = TimeAxis(0.0, 1000, 2.0)
    ham, sbi = _heterodimer(ta)
    with qr.eigenbasis_of(ham):
        if touch_hamiltonian:
            _ = ham.data
        RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi)
        HH = numpy.array(ham.data, copy=True)
        RR = numpy.array(RRT.data, copy=True)
    return ham, RRT, HH, RR


def test_redfield_tensor_independent_of_lazy_basis_state():
    """#333: the tensor built in eigenbasis_of(H) must not depend on whether
    H.data was already accessed (i.e. lazily transformed) in the context."""
    _, _, HH_ref, RR_ref = _redfield_in_eigenbasis(touch_hamiltonian=False)
    _, _, HH, RR = _redfield_in_eigenbasis(touch_hamiltonian=True)

    numpy.testing.assert_allclose(HH, HH_ref, rtol=0, atol=1e-12)
    # |R| ~ 4e-2 in internal units; the two paths differ only by rounding in
    # the basis transformations (~1e-17 observed).
    assert numpy.all(numpy.isfinite(RR))
    numpy.testing.assert_allclose(RR, RR_ref, rtol=0, atol=1e-12)


def test_secular_redfield_propagation_in_eigenbasis_reaches_boltzmann():
    """#333: secular Redfield propagation inside eigenbasis_of(H) relaxes to
    the Boltzmann distribution, even when H.data was accessed first."""
    ham, RRT, _, _ = _redfield_in_eigenbasis(touch_hamiltonian=True)

    # Relaxation time is ~2 ps; 20 ps leaves a residual of ~exp(-10) ~ 5e-5.
    ta_long = TimeAxis(0.0, 10001, 2.0)
    with qr.eigenbasis_of(ham):
        RRT.secularize()
        k_down = numpy.real(RRT.data[1, 1, 2, 2])
        k_up = numpy.real(RRT.data[2, 2, 1, 1])
        rho0 = qr.ReducedDensityMatrix(dim=ham.dim)
        rho0.data[2, 2] = 1.0
        prop = qr.ReducedDensityMatrixPropagator(ta_long, ham, RRT)
        rhot = prop.propagate(rho0)
        # copy: the context back-transforms rhot.data in place on exit
        pops = numpy.real(numpy.diag(rhot.data[-1])).copy()
        dE_int = numpy.real(ham.data[2, 2] - ham.data[1, 1])

    assert k_up > 0 and k_down > 0
    ratio = pops[2] / pops[1]
    # For secular Redfield the steady state is exactly k_up / k_down.
    assert abs(ratio - k_up / k_down) / (k_up / k_down) < 1e-3

    cm_to_int = 2.0 * numpy.pi * 3.0e10 * 1.0e-15
    expected = _boltzmann_ratio(dE_int / cm_to_int)
    # Redfield rates obey detailed balance up to the time discretization of
    # the rate integrals (~2% for dt = 2 fs here, shrinking with dt).
    assert abs(ratio - expected) / expected < 0.05, (
        f"p2/p1 = {ratio:.4f}, Boltzmann {expected:.4f}"
    )
