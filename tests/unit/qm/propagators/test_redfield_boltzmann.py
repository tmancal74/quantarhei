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
