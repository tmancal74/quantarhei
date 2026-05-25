import unittest
from pathlib import Path

import numpy
import numpy.testing as npt

import quantarhei as qr
from quantarhei.spectroscopy.pumpprobe import PumpProbeSpectrum
from quantarhei.symbolic.cumulant import GFInitiator
from quantarhei.utils.vectors import X

# import matplotlib.pyplot as plt


"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.responses package


*******************************************************************************
"""

TEST_DIR = Path(__file__).parent


def efce_R1g_vec(t2, t1s, t3s, aa):
    """Our R1g response: vectorial version"""
    import numpy as np

    aa.set_t2(t2)

    return np.exp(
        -aa.gg_t1["ab"]
        + aa.gg_t1_t2["ab"]
        - aa.gg_t1_t2_t3["aa"]
        - np.conj(aa.gg_t2["bb"])
        + np.conj(aa.gg_t2_t3["ab"])
        - np.conj(aa.gg_t3["ab"])
    )


def efce_R2g_vec(t2, t1s, t3s, aa):
    """Our R2g response: vectorial version"""
    import numpy as np

    aa.set_t2(t2)

    return np.exp(
        aa.gg_t2["ab"]
        - aa.gg_t2_t3["bb"]
        - np.conj(aa.gg_t1["ab"])
        - np.conj(aa.gg_t1_t2["aa"])
        + np.conj(aa.gg_t1_t2_t3["ab"])
        - np.conj(aa.gg_t3["ab"])
    )


def efce_R3g_vec(t2, t1s, t3s, aa):
    """Our R23g response: vectorial version"""
    import numpy as np

    aa.set_t2(t2)

    return np.exp(
        -aa.gg_t3["bb"]
        - np.conj(aa.gg_t1["aa"])
        - np.conj(aa.gg_t1_t2["ab"])
        + np.conj(aa.gg_t1_t2_t3["ab"])
        + np.conj(aa.gg_t2["ab"])
        - np.conj(aa.gg_t2_t3["ab"])
    )


def efce_R4g_vec(t2, t1s, t3s, aa):
    """Our R4g response: vectorial version"""
    import numpy as np

    aa.set_t2(t2)

    return np.exp(
        -aa.gg_t1["aa"]
        + aa.gg_t1_t2["ab"]
        - aa.gg_t1_t2_t3["ab"]
        - aa.gg_t2["ab"]
        + aa.gg_t2_t3["ab"]
        - aa.gg_t3["bb"]
    )


def _reference_pump_probe_aggregate():
    """Small dimer with enough bath data for the pump-probe reference test."""
    env_axis = qr.TimeAxis(0.0, 400, 5.0)
    params = {
        "ftype": "OverdampedBrownian",
        "reorg": 40.0,
        "cortime": 100.0,
        "T": 100.0,
        "matsubara": 20,
    }

    with qr.energy_units("1/cm"):
        mol1 = qr.Molecule([0.0, 12000.0])
        mol2 = qr.Molecule([0.0, 12300.0])
        cf1 = qr.CorrelationFunction(env_axis, params)
        cf2 = qr.CorrelationFunction(env_axis, params)
        mol1.set_transition_environment((0, 1), cf1)
        mol2.set_transition_environment((0, 1), cf2)

    mol1.set_dipole(0, 1, [1.0, 0.8, 0.8])
    mol2.set_dipole(0, 1, [0.8, 0.8, 0.0])

    agg = qr.Aggregate(molecules=[mol1, mol2])
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0, 1, 100.0)

    agg.build(mult=2, sbi_for_higher_ex=True)
    agg.diagonalize()
    return agg


def _reference_calculate_pathways_rdm(calc, rdm0, rdm, tau, lab, ptol=1.0e-6):
    """Copied pump-probe RDM formula used as a regression reference."""
    onepp = PumpProbeSpectrum()
    onepp.set_axis(calc.oa3)

    ss = calc.system.SS.copy()

    if calc.goft_matrix is not None:
        gt3s = calc.goft_matrix
    else:
        gt3s = calc._SE_excitonic_gofts(
            ss, calc.system, tau=0.0, _diag_double_only=True
        )
        calc.goft_matrix = gt3s
    gt3tau = calc._SE_excitonic_gofts(ss, calc.system, tau=tau, _diag_double_only=True)

    ppspec = numpy.zeros(calc.t3axis.length, dtype=numpy.complex128)
    dim = calc.system.Nb[1] + calc.system.Nb[0]

    for jj in range(1, dim):
        pref_gsb = lab.F4eM4[0]
        pref_gsb *= 2
        pref_gsb *= numpy.sum(numpy.diag(rdm0))
        pref_gsb *= numpy.dot(calc.system.DD[jj, 0, :], calc.system.DD[jj, 0, :])

        om = calc.system.HH[jj, jj] - calc.system.HH[0, 0] - calc.rwa
        ft = -1j * om * calc.t3axis.data
        ft -= gt3s[jj, jj]
        ppspec += pref_gsb * numpy.exp(ft)

    for ii in range(1, dim):
        om = calc.system.HH[ii, ii] - calc.system.HH[0, 0] - calc.rwa

        for jj in range(1, dim):
            if rdm[ii, jj] < ptol:
                continue

            omtau = calc.system.HH[ii, ii] - calc.system.HH[jj, jj]
            pref_se = lab.F4eM4[0]
            pref_se *= 2
            pref_se *= rdm[ii, jj]
            pref_se *= numpy.dot(calc.system.DD[ii, 0, :], calc.system.DD[jj, 0, :])

            ft = -1j * om * calc.t3axis.data - 1j * omtau * tau
            gt1 = gt3tau[ii, ii] - numpy.conj(gt3tau[ii, jj])
            gt2 = numpy.conj(gt3tau[jj, jj, 0]) - gt3tau[ii, jj, 0]
            gt3 = numpy.conj(gt3s[jj, ii])
            ft -= gt1 + gt2 + gt3

            ppspec += pref_se * numpy.exp(ft)

    for ii in range(1, dim):
        for jj in range(1, dim):
            if numpy.abs(rdm[ii, jj]) < ptol:
                continue

            for ll in range(dim, calc.system.Ntot):
                om = calc.system.HH[ll, ll] - calc.system.HH[jj, jj] - calc.rwa
                omtau = calc.system.HH[ii, ii] - calc.system.HH[jj, jj]

                pref_esa = lab.F4eM4[0]
                pref_esa *= 2
                pref_esa *= rdm[ii, jj]
                pref_esa *= numpy.dot(
                    calc.system.DD[ll, ii, :], calc.system.DD[jj, ll, :]
                )

                fl = ll
                ek = jj
                ej = ii

                ft = -1j * om * calc.t3axis.data - 1j * omtau * tau
                gt1 = gt3s[fl, fl] - gt3s[fl, ek] + gt3s[ej, ek] - gt3s[ej, fl]
                gt2 = (
                    gt3tau[ej, ej, 0]
                    - numpy.conj(gt3tau[ej, ek, 0])
                    - gt3tau[fl, ej, 0]
                    + numpy.conj(gt3tau[fl, ek, 0])
                )
                gt3 = (
                    numpy.conj(gt3tau[ek, ek])
                    - gt3tau[ek, ej]
                    + gt3tau[fl, ej]
                    - numpy.conj(gt3tau[fl, ek])
                )
                ft -= gt1 + gt2 + gt3

                ppspec -= pref_esa * numpy.exp(ft)

    ppspec = -ppspec
    ft = numpy.fft.hfft(ppspec) * calc.t3axis.step
    ft = numpy.fft.fftshift(ft)
    ft = numpy.flipud(ft)
    nt = calc.t3axis.length

    data = numpy.real(ft[nt // 2 : nt + nt // 2])
    data = calc.oa3.data * data

    onepp._add_data(data)
    onepp.set_t2(tau)

    return onepp


class TestPumpProbe(unittest.TestCase):
    """Tests for the response package"""

    def setUp(self, verbose=False):
        pass

    def test_PumpProbe_1(self):
        """Testing basic functions of the TwoDSpectrumBase class"""
        #
        # Lineshape function
        #
        Nt1 = 50
        dt1 = 10
        Nt3 = 50
        dt3 = 10

        Nt2 = 2
        dt2 = 10.0

        t1_axis = qr.TimeAxis(0.0, Nt1, dt1)
        t3_axis = qr.TimeAxis(0.0, Nt3, dt3)
        t2_axis = qr.TimeAxis(0.0, Nt2, dt2)

        temperature = 300.0
        # parameters of the correlation function
        params = {
            "ftype": "OverdampedBrownian",
            "reorg": 80.0,
            "cortime": 150.0,
            "T": temperature,
            "matsubara": 20,
        }

        deg = [1.0, 0.0, 0.0]
        dfe = [1.0, 0.0, 0.0]

        d1 = deg
        d2 = deg
        d3 = deg
        d4 = deg
        d3f = dfe
        d4f = dfe

        omeg = 12100.0
        omfe = 12100.0

        #
        # we get the g-functions as spline interpolated functions
        #
        tt = qr.TimeAxis(0.0, 3 * t1_axis.length, t1_axis.step)
        with qr.energy_units("1/cm"):
            cf = qr.LineshapeFunction(tt, params)

        sfce = cf.as_spline_function()

        responses = []

        response_types = ["R1g", "R2g", "R3g", "R4g"]
        response_functions = [efce_R1g_vec, efce_R2g_vec, efce_R3g_vec, efce_R4g_vec]

        gdict = dict(aa=sfce, bb=sfce, ab=sfce, ba=sfce)

        gf = GFInitiator(t1_axis.data, t3_axis.data, gdict)

        for rftype, rffce in zip(response_types, response_functions):
            rfce = qr.ResponseFunction(rftype)
            rfce.set_evaluation_function(rffce)
            responses.append(rfce)

        #
        # This part is the same for all two-level system responses
        #
        for rfce in responses:
            rfce.set_dipoles(d1, d2, d3, d4)
            args = (gf,)
            rfce.set_auxliary_arguments(args)

        #
        # Set transition frequencies for all responses (beware of units)
        #
        with qr.energy_units("1/cm"):
            for rsp in responses:
                rsp.set_frequencies(omeg, omeg)

        #
        # Calculation of the 2D spectra from response functions
        #

        rwa = 12000.0
        Npad = 0

        # laboratory settings
        lab = qr.LabSetup()
        lab.set_pulse_polarizations(
            pulse_polarizations=(X, X, X), detection_polarization=X
        )

        calc = qr.TwoDResponseCalculator(t1_axis, t2_axis, t3_axis, responses=responses)
        with qr.energy_units("1/cm"):
            calc.bootstrap(rwa=rwa, pad=Npad, lab=lab, verbose=True)

        one = False
        T2 = 0.0

        if one:
            print("Calculating one spectrum")
            twod = calc.calculate_one(0)

        else:
            print("Calculating", Nt2, "spectra")
            tcont = calc.calculate()
            tcont = tcont.get_TwoDSpectrumContainer()
            twod = tcont.get_spectrum(T2)

        file_path_1 = TEST_DIR / "pumpprobe_test_data_0.dat"
        save_data = False
        if save_data:
            twod.save_data(file_path_1)

        #
        # Here we load spectrum for comparison
        #
        twod0 = qr.TwoDSpectrum()
        twod0.load_data(file_path_1)
        twod0.set_axis_1(t1_axis)
        twod0.set_axis_3(t3_axis)

        npt.assert_allclose(twod0.data, twod.data)

    def test_PumpProbe_2(self):
        """Testing basic functions of the TwoDSpectrumCalculator class"""
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)

        t2 = qr.TimeAxis(30, 10, 10.0)

        twod_calc = qr.TwoDResponseCalculator(t1, t2, t3)

    def test_pump_probe_rdm_reference_formula(self):
        """Compare current pump-probe RDM calculation with a copied reference."""
        t2_axis = qr.TimeAxis(0.0, 3, 10.0)
        t3_axis = qr.TimeAxis(0.0, 32, 5.0)
        agg = _reference_pump_probe_aggregate()

        calc = qr.PumpProbeSpectrumCalculator(t2_axis, t3_axis, system=agg)
        lab = qr.LabSetup()
        lab.set_pulse_polarizations(
            pulse_polarizations=(X, X, X), detection_polarization=X
        )
        with qr.energy_units("1/cm"):
            calc.bootstrap(rwa=12100.0, lab=lab)

        rdm0 = numpy.zeros((calc.system.Ntot, calc.system.Ntot), dtype=numpy.float64)
        rdm0[1, 1] = 0.65
        rdm0[2, 2] = 0.35

        rdm = numpy.zeros((calc.system.Ntot, calc.system.Ntot), dtype=numpy.float64)
        rdm[1, 1] = 0.55
        rdm[2, 2] = 0.25
        rdm[1, 2] = 0.05
        rdm[2, 1] = 0.05

        tau = 20.0
        calc.goft_matrix = None
        current = calc.calculate_pathways_rdm(rdm0, rdm, tau, lab)

        calc.goft_matrix = None
        reference = _reference_calculate_pathways_rdm(calc, rdm0, rdm, tau, lab)

        npt.assert_allclose(current.axis.data, reference.axis.data)
        npt.assert_allclose(current.data, reference.data, rtol=1.0e-12, atol=1.0e-12)
        self.assertEqual(current.get_t2(), reference.get_t2())
