import unittest
from pathlib import Path

import numpy
import numpy.testing as npt

import quantarhei as qr
from quantarhei import (
    REAL,
    CorrelationFunctionMatrix,
    SystemBathInteraction,
    TransitionDipoleMoment,
)
from quantarhei.qm import Operator

_show_spectra_ = False
_show_loaded_ = False
_save_data_ = False

if _show_spectra_ or _show_loaded_:
    import matplotlib.pyplot as plt


"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.responses package


*******************************************************************************
"""

TEST_DIR = Path(__file__).parent


class TestTwoDSpectrum(unittest.TestCase):
    """Tests for the response package"""

    def setUp(self, verbose=False):
        #
        #  Chlorophyll parameters
        #

        E1 = 16807.0
        dE = 100.0

        d_ge = 1.0
        d_ef = 1.0

        E2 = E1 + E1 + dE  # 23600.0

        d_0_1 = [d_ge, 0.0, 0.0]
        # d_0_Qx = [0.2, 0.0, 0.0]
        d_1_2 = [d_ef, 0.0, 0.0]
        # d_Qx_S2 = [1.0, 0.0, 0.0]

        #
        # Other parameters
        #

        # temperature in K
        temperature = 300.0

        # specifying the pulse shapes
        pulse_fwhm = 359.0
        pulse_amplitude = 1.0
        pulse_frequency = 16807.0

        #
        # Time axis parameters
        #
        Nt = 50
        dt = 5.0

        Nt2 = 2
        dt2 = 100.0

        #
        #  Time axes
        #

        t1_axis = qr.TimeAxis(0.0, Nt, dt)
        t3_axis = qr.TimeAxis(0.0, Nt, dt)

        t2_axis = qr.TimeAxis(0.0, Nt2, dt2)

        #
        # Bath correlation functions for the molecular transitions
        #

        cf01_params = {
            "ftype": "OverdampedBrownian",
            "reorg": 140.0,
            "cortime": 60.0,
            "T": temperature,
            "matsubara": 20,
        }

        # md01_params = {"ftype":  "UnderdampedBrownian",
        #               "reorg": 70.0,
        #               "freq":  500.0,
        #               "gamma": 100.0,
        #               "T": temperature}

        cf02_params = {
            "ftype": "OverdampedBrownian",
            "reorg": 140.0,
            "cortime": 60.0,
            "T": temperature,
            "matsubara": 20,
        }

        # md02_params = {"ftype":  "UnderdampedBrownian",
        #               "reorg": 60.0,
        #               "freq":  500.0,
        #               "gamma": 100.0,
        #               "T": temperature}

        # Build the correlation function
        with qr.energy_units("1/cm"):
            cfce1 = qr.CorrelationFunction(t1_axis, cf01_params)
            # c1_under = qr.CorrelationFunction(t1_axis, md01_params)
            cfce2 = qr.CorrelationFunction(t1_axis, cf02_params)
            # c2_under = qr.CorrelationFunction(t1_axis, md02_params)

        # cfce1 += c1_under
        # cfce2 += c2_under
        cfce3 = cfce1

        from quantarhei.core.managers import UnitsManaged

        #
        # custom three level molecule class
        #

        class Chlorophyll(UnitsManaged, qr.OpenSystem):
            def build(self):
                """Building internal information of the system"""
                self.el_rwa_indices = numpy.array(
                    [ii for ii in range(self.Nel)], dtype=numpy.int32
                )

                # Hamiltonian
                ham = qr.Hamiltonian(data=numpy.diag(self.elenergies))

                self.HH = numpy.zeros((self.Nel, self.Nel), dtype=REAL)
                self.HH[:, :] = ham.data[:, :]
                self.HamOp = ham

                # transition dipole moment operator

                DD = self.dmoments
                self.DD = DD
                trdata = numpy.zeros(
                    (DD.shape[0], DD.shape[1], DD.shape[2]), dtype=REAL
                )
                trdata[:, :, :] = DD[:, :, :]
                self.TrDMOp = TransitionDipoleMoment(data=trdata)

                #
                # system-bath interaction
                #

                #
                # Correlation fuctions
                #
                # Assuming that transition environments were already set
                Ncf = (
                    self.Nel - 1
                )  # all transitions from the ground state will be specified
                egcf1 = self.get_transition_environment(
                    (0, 1)
                )  # we get the time axis like this
                timea = egcf1.axis

                self.egcf_matrix = CorrelationFunctionMatrix(timea, Ncf)

                # only have two transitions so far
                fc1 = self.get_transition_environment((0, 1))
                if self.Nel == 3:
                    fc2 = self.get_transition_environment((0, 2))
                # fc3 = self.get_transition_environment((0,3))

                # !!! in the egcf matrix the excited states are labeled from 0 (it is like numbering the monomers in an aggregate)
                self.egcf_matrix.set_correlation_function(fc1, [(0, 0)])
                if self.Nel == 3:
                    self.egcf_matrix.set_correlation_function(fc2, [(1, 1)])
                # self.egcf_matrix.set_correlation_function(fc3,[(2,2)])

                Nop = Ncf

                iops = []
                for i in range(1, Nop + 1):
                    op1 = Operator(dim=self.HH.shape[0], real=True)
                    op1.data[i, i] = 1.0
                    iops.append(op1)

                self.sbi = SystemBathInteraction(iops, self.egcf_matrix, system=self)

                self._built = True

            def get_weighted_participation(self):
                MM = super().get_weighted_participation()
                MM[1, 0, 0] = 0.0
                MM[0, 1, 0] = 0.0
                return MM

        #
        # Definition of the Chlorophyll molecule
        #
        with qr.energy_units("1/cm"):
            chl = Chlorophyll([0.0, E1, E2])

        chl.set_dipole(0, 1, d_0_1)
        chl.set_dipole(1, 2, d_1_2)

        chl.set_transition_environment((0, 1), cfce1)
        chl.set_transition_environment((0, 2), cfce2)

        chl.build()

        #
        # Verification
        #

        sys = chl

        Ham = sys.get_Hamiltonian()

        print("\nHamiltonian")
        print(Ham)

        Dip = sys.get_TransitionDipoleMoment()
        print("\nDipole moment x\n")
        print(Dip.data[:, :, 0])

        rwa = sys.get_RWA_suggestion()

        print(rwa)

        #
        # Lab settings (Pulse polarization for the 2D spectrum calculation)
        #
        from quantarhei.spectroscopy import X

        # laboratory settings class
        lab = qr.LabSetup()
        lab.set_pulse_polarizations(
            pulse_polarizations=(X, X, X), detection_polarization=X
        )

        self.sys = sys
        self.lab = lab
        self.Nt = Nt
        self.dt = dt
        self.Nt2 = Nt2
        self.dt2 = dt2

    def test_twod_1(self):

        Nt1 = self.Nt
        dt1 = self.dt
        Nt3 = self.Nt
        dt3 = self.dt

        Nt2 = self.Nt2
        dt2 = self.dt2
        t1_axis = qr.TimeAxis(0.0, Nt1, dt1)
        t3_axis = qr.TimeAxis(0.0, Nt3, dt3)
        t2_axis = qr.TimeAxis(0.0, Nt2, dt2)

        sys = self.sys
        lab = self.lab

        #
        # Setting up 2D spectrum calculation
        #
        calc = qr.TwoDResponseCalculator(t1_axis, t2_axis, t3_axis, system=sys)
        with qr.energy_units("1/cm"):
            rwa = sys.get_RWA_suggestion()
            print("RWA = ", rwa, "1/cm")
            calc.bootstrap(rwa=rwa, pad=0, lab=lab, verbose=True)

        print("RWA = ", calc.rwa)

        #
        # calculate one or many 2D spectra
        #

        # T2 time for verification ploting
        T2 = 100.0  # fs

        calc.reset_t2_time()
        # print("Calculating", Nt2,"spectra")
        # t1 = time.time()
        tcont = calc.calculate()
        # t2 = time.time()
        # print("...done in", t2-t1,"sec")

        #
        # Plot one spectrum for verification
        #

        scont = tcont.get_TwoDSpectrumContainer()

        twod1 = scont.get_spectrum(0.0)
        twod2 = scont.get_spectrum(T2)

        if _show_spectra_:
            # plot_window = [11000,13000,11000,13000]
            plot_window = [3.07, 3.27, 3.07, 3.27]
            with qr.energy_units("1/fs"):
                twod1.plot(Npos_contours=10, window=plot_window)

            # plot_window = [11000,13000,11000,13000]
            plot_window = [3.07, 3.27, 3.07, 3.27]
            with qr.energy_units("1/fs"):
                twod2.plot(Npos_contours=10, window=plot_window)

            plt.show()

        file_path_1 = TEST_DIR / "twodspectrum_test_data_0.dat"
        file_path_2 = TEST_DIR / "twodspectrum_test_data_100.dat"

        if _save_data_:
            twod1.save_data(file_path_1)
            twod2.save_data(file_path_2)

        #
        # Here we load spectrum for comparison
        #
        twod01 = qr.TwoDSpectrum()
        twod01.load_data(file_path_1)
        twod01.set_axis_1(t1_axis)
        twod01.set_axis_3(t3_axis)

        twod02 = qr.TwoDSpectrum()
        twod02.load_data(file_path_2)
        twod02.set_axis_1(t1_axis)
        twod02.set_axis_3(t3_axis)

        # if _show_loaded_:
        #    #plot_window = [11000,13000,11000,13000]
        #    plot_window = [3.07,3.27,3.07,3.27]
        #    with qr.energy_units("1/fs"):
        #        twod01.plot(Npos_contours=10, window=plot_window)

        #    #plot_window = [11000,13000,11000,13000]
        #    plot_window = [3.07,3.27,3.07,3.27]
        #    with qr.energy_units("1/fs"):
        #        twod02.plot(Npos_contours=10, window=plot_window)
        #
        #    plt.show()

        npt.assert_allclose(twod01.data, twod1.data)
        npt.assert_allclose(twod02.data, twod2.data)

    def test_twod_2(self):
        """Testing basic functions of the TwoDSpectrumCalculator class"""
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)

        t2 = qr.TimeAxis(30, 10, 10.0)

        twod_calc = qr.TwoDResponseCalculator(t1, t2, t3)
