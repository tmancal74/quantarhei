import quantarhei as qr


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


def test_redfieldtensor_secular_flag():
    ham, sbi = _make_ham_sbi()
    RRT = qr.qm.RedfieldRelaxationTensor(ham, sbi, secular=True)
    assert RRT.is_secular
