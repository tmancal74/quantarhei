import unittest
from pathlib import Path

import numpy
import numpy.testing as npt

import quantarhei as qr
from quantarhei.spectroscopy.response_implementations import get_implementation
from quantarhei.spectroscopy.responses import RESPONSE_DIAGRAMS, NonLinearResponse
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


class TestResponses(unittest.TestCase):
    """Tests for the response package"""

    def setUp(self, verbose=False):
        pass

    def test_response_diagram_mapping(self):
        """Testing explicit response diagram metadata"""
        self.assertEqual(RESPONSE_DIAGRAMS["R1g"]["process"], "SE")
        self.assertEqual(RESPONSE_DIAGRAMS["R1g"]["signal"], qr.signal_NONR)
        self.assertEqual(RESPONSE_DIAGRAMS["R2g"]["process"], "SE")
        self.assertEqual(RESPONSE_DIAGRAMS["R2g"]["signal"], qr.signal_REPH)

        self.assertEqual(RESPONSE_DIAGRAMS["R3g"]["process"], "GSB")
        self.assertEqual(RESPONSE_DIAGRAMS["R3g"]["signal"], qr.signal_REPH)
        self.assertEqual(RESPONSE_DIAGRAMS["R4g"]["process"], "GSB")
        self.assertEqual(RESPONSE_DIAGRAMS["R4g"]["signal"], qr.signal_NONR)

        self.assertEqual(RESPONSE_DIAGRAMS["R1f"]["process"], "ESA")
        self.assertEqual(RESPONSE_DIAGRAMS["R1f"]["signal"], qr.signal_NONR)
        self.assertEqual(RESPONSE_DIAGRAMS["R2f"]["process"], "ESA")
        self.assertEqual(RESPONSE_DIAGRAMS["R2f"]["signal"], qr.signal_REPH)

        self.assertTrue(RESPONSE_DIAGRAMS["R1g_scM0g"]["transfer"])
        self.assertEqual(RESPONSE_DIAGRAMS["R1g_scM0g"]["process"], "SE")
        self.assertTrue(RESPONSE_DIAGRAMS["R2f_scM0e"]["transfer"])
        self.assertEqual(RESPONSE_DIAGRAMS["R2f_scM0e"]["process"], "ESA")
        self.assertTrue(RESPONSE_DIAGRAMS["R1g_scM1g"]["transfer"])
        self.assertEqual(RESPONSE_DIAGRAMS["R1g_scM1g"]["process"], "SE")
        self.assertTrue(RESPONSE_DIAGRAMS["R2f_scM1e"]["transfer"])
        self.assertEqual(RESPONSE_DIAGRAMS["R2f_scM1e"]["process"], "ESA")

    def test_Responses(self):
        """Testing basic functions of the TwoDSpectrumBase class"""
        #
        # Lineshape function
        #
        Nt1 = 100
        dt1 = 11
        Nt3 = 100
        dt3 = 12

        Nt2 = 50
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

        file_path_1 = TEST_DIR / "responses_test_twod_0.dat"
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

    def test_LouvillePathway(self):
        """Testing basic functions of the TwoDSpectrumCalculator class"""
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)

        t2 = qr.TimeAxis(30, 10, 10.0)

        twod_calc = qr.TwoDResponseCalculator(t1, t2, t3)

    def test_time_dependent_rate_matrix_setup(self):
        """Testing response setup with time-dependent population rates"""

        class System:
            Nb = [0, 2, 0]
            Ntot = 2
            mult = 1

        t1 = qr.TimeAxis(0.0, 10, 1.0)
        t2 = qr.TimeAxis(0.0, 11, 1.0)
        t3 = qr.TimeAxis(0.0, 10, 1.0)

        response = NonLinearResponse(None, System(), "R1g", t1, t2, t3)

        KK = numpy.zeros((t2.length, 2, 2), dtype=numpy.float64)
        KK[:, 0, 1] = 0.001
        KK[:, 1, 0] = 0.002
        KK[:, 0, 0] = -KK[:, 1, 0]
        KK[:, 1, 1] = -KK[:, 0, 1]

        response.set_rate_matrix(KK)

        self.assertEqual(response.Uee.shape, (2, 2, t2.length))
        self.assertEqual(response.U1_t2.shape, (2, 2, t2.length))
        npt.assert_allclose(response.U1_t2[:, :, 0], numpy.eye(2))
        npt.assert_allclose(response.Uremainder_t2, response.Uee - response.U1_t2)
        npt.assert_allclose(response.Utransfer_t2, response.Uremainder_t2)

    def test_zero_jump_dephasing_factors(self):
        """Testing half-rate damping factors on all response time axes"""

        class System:
            Nb = [0, 2, 0]
            Ntot = 2
            mult = 1

        t1 = qr.TimeAxis(0.0, 10, 1.0)
        t2 = qr.TimeAxis(0.0, 6, 2.0)
        t3 = qr.TimeAxis(0.0, 10, 1.0)

        response = NonLinearResponse(None, System(), "R1g", t1, t2, t3)

        KK = numpy.zeros((2, 2), dtype=numpy.float64)
        KK[0, 0] = -0.002
        KK[1, 1] = -0.004

        response.set_rate_matrix(KK)

        npt.assert_allclose(response.U0_t1[0], numpy.exp(0.5 * KK[0, 0] * t1.data))
        npt.assert_allclose(response.U0_t2[1], numpy.exp(0.5 * KK[1, 1] * t2.data))
        npt.assert_allclose(response.U0_t3[0], numpy.exp(0.5 * KK[0, 0] * t3.data))
        npt.assert_allclose(response.U1_t2[0, 0], numpy.exp(KK[0, 0] * t2.data))
        npt.assert_allclose(response.Uremainder_t2, response.Uee - response.U1_t2)
        npt.assert_allclose(response.Utransfer_t2, response.Uremainder_t2)

    def test_single_jump_transfer_channel(self):
        """Testing one-jump transfer channel and remainder bookkeeping"""

        class LineshapeStorage:
            time_mapping = {"t1": 0, "t2": 1, "s2": 2}

        class System:
            Nb = [0, 2, 0]
            Ntot = 2
            mult = 1

            def get_lineshape_functions(self):
                return LineshapeStorage()

        t1 = qr.TimeAxis(0.0, 10, 1.0)
        t2 = qr.TimeAxis(0.0, 6, 2.0)
        t3 = qr.TimeAxis(0.0, 10, 1.0)

        response0 = NonLinearResponse(None, System(), "R1g_scM0g", t1, t2, t3)
        response1 = NonLinearResponse(None, System(), "R1g_scM1g", t1, t2, t3)

        KK = numpy.array([[-0.002, 0.003], [0.002, -0.003]], dtype=numpy.float64)
        response0.set_rate_matrix(KK)
        response1.set_rate_matrix(KK)

        self.assertEqual(len(response0.Ujump_t2), 2)
        npt.assert_allclose(
            response0.Uremainder_t2,
            response0.Uee - response0.U1_t2 - response0.Usingle_t2,
        )
        npt.assert_allclose(
            response0._get_transfer_matrix(1), response0.Uremainder_t2[:, :, 1]
        )
        npt.assert_allclose(
            response1._get_transfer_matrix(1), response1.Usingle_t2[:, :, 1]
        )

    def test_single_jump_response_requires_storage_labels(self):
        """Testing storage validation for single-jump response functions"""

        class LineshapeStorage:
            time_mapping = {"t1": 0, "t2": 1, "t3": 2}

        class System:
            def get_lineshape_functions(self):
                return LineshapeStorage()

        response = get_implementation("R1g_scM1g")
        with self.assertRaisesRegex(Exception, "FunctionStorage\\(config=1\\)"):
            response(0.0, numpy.zeros(1), numpy.zeros(1), None, System(), None, None)

    def test_response_rate_matrix_relaxation_theory(self):
        """Testing response setup from OpenSystem rate matrix theory names"""

        class RateMatrix:
            data = numpy.array(
                [
                    [0.0, 0.0, 0.0],
                    [0.0, -0.002, 0.001],
                    [0.0, 0.002, -0.001],
                ],
                dtype=numpy.float64,
            )

        class System:
            Nb = [1, 2, 0]
            Ntot = 3
            mult = 1
            requested_theory = None

            def get_band(self, band):
                if band == 1:
                    return (1, 2)
                return (0,)

            def get_RateMatrix(
                self,
                relaxation_theory=None,
                time_dependent=False,
                relaxation_cutoff_time=None,
            ):
                self.requested_theory = relaxation_theory
                self.requested_time_dependent = time_dependent
                self.requested_cutoff_time = relaxation_cutoff_time
                return RateMatrix()

        t1 = qr.TimeAxis(0.0, 10, 1.0)
        t2 = qr.TimeAxis(0.0, 11, 1.0)
        t3 = qr.TimeAxis(0.0, 10, 1.0)
        system = System()

        response = NonLinearResponse(
            None,
            system,
            "R1g",
            t1,
            t2,
            t3,
            relaxation_theory="standard_Foerster",
            relaxation_cutoff_time=100.0,
        )

        self.assertEqual(system.requested_theory, "standard_Foerster")
        self.assertFalse(system.requested_time_dependent)
        self.assertEqual(system.requested_cutoff_time, 100.0)
        npt.assert_allclose(response.KK, RateMatrix.data[1:3, 1:3])
