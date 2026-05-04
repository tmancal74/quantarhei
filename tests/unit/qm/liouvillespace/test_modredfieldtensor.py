import pytest

import quantarhei as qr
from quantarhei.core.managers import BasisContextError


def _make_ham_sbi():
    with qr.energy_units("1/cm"):
        mol1 = qr.Molecule([0.0, 12000.0])
        mol2 = qr.Molecule([0.0, 12200.0])
        mol1.set_transition_environment(
            (0, 1),
            qr.CorrelationFunction(
                qr.TimeAxis(0.0, 1000, 1.0),
                {
                    "ftype": "OverdampedBrownian",
                    "reorg": 100.0,
                    "cortime": 100.0,
                    "T": 300.0,
                    "matsubara": 20,
                },
            ),
        )
        mol2.set_transition_environment(
            (0, 1),
            qr.CorrelationFunction(
                qr.TimeAxis(0.0, 1000, 1.0),
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
