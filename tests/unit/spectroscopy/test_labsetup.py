import unittest

import matplotlib.pyplot as plt
import numpy
import numpy.testing as npt

from quantarhei import (
    Aggregate,
    CorrelationFunction,
    DFunction,
    LabSetup,
    Molecule,
    ReducedDensityMatrixPropagator,
    TimeAxis,
    convert,
    eigenbasis_of,
    energy_units,
)
from quantarhei.exceptions import QuantarheiError

# import quantarhei as qr
from quantarhei.utils.vectors import X, Y, Z


class TestLabSetup(unittest.TestCase):
    """Test of the laboratory setup"""

    def setUp(self):

        self._plot_ = False

        ################################################################################
        #
        # Set up laboratory/experiment configuration
        #
        ################################################################################

        # three pulses will be used
        lab = LabSetup(nopulses=3)

        # on a time axis starting at a specified time, with a certain number of steps
        # and a step size
        Nfr = 10
        time = TimeAxis(-500.0, Nfr * 1500, 1.0 / Nfr, atype="complete")

        # pulse shapes are specified below
        pulse2 = dict(ptype="Gaussian", FWHM=20, amplitude=0.1)
        params = (pulse2, pulse2, pulse2)

        # with self.assertRaises(Exception) as context:
        #
        #    lab.set_pulse_shapes(time, params)
        #
        # self.assertTrue("Pulse arrival times have to specified"
        #                in str(context.exception))

        # each pulse has a defined frequency
        lab.set_pulse_polarizations(
            pulse_polarizations=(X, Y, Z), detection_polarization=X
        )

        # time of arrival
        lab.set_pulse_arrival_times([0.0, 0.0, 100.0])

        Om = 0.0
        # and polarization
        ome = convert(10200.0 - Om, "1/cm", "int")
        lab.set_pulse_frequencies([ome, ome, ome])

        # additional phases can be also controlled
        lab.set_pulse_phases([0.0, 1.0, 0.0])

        lab.set_pulse_shapes(time, params)

        self.lab = lab
        self.time = time

        with energy_units("1/cm"):
            m1 = Molecule([0.0, 10000.0 - Om])
            m1.set_dipole((0, 1), [1.0, 0.0, 0.0])
            m2 = Molecule([0.0, 10000.0 - Om])
            m2.set_dipole((0, 1), [0.0, 0.0, 0.0])

            agg = Aggregate(molecules=[m1, m2])

            agg.set_resonance_coupling(0, 1, 200.0)

        agg.build()

        time_r = TimeAxis(0.0, 1000, 1.0)
        with energy_units("1/cm"):
            m1 = Molecule([0.0, 10000.0 - Om])
            m1.set_dipole((0, 1), [1.0, 0.0, 0.0])
            m2 = Molecule([0.0, 10000.0 - Om])
            m2.set_dipole((0, 1), [0.0, 0.0, 0.0])

            agg2 = Aggregate(molecules=[m1, m2])

            params = dict(
                ftype="OverdampedBrownian",
                T=300.0,
                reorg=30.0,
                cortime=30.0,
                matsubara=30,
            )
            cf = CorrelationFunction(time_r, params)

            agg2.set_resonance_coupling(0, 1, 200.0)

        m1.set_transition_environment((0, 1), cf)
        m2.set_transition_environment((0, 1), cf)

        agg2.build()

        self.agg = agg
        self.aggB = agg2  # aggregate with bath
        self.time_r = time_r

    def test_lab_pulse_setters(self):
        """(LabSetup) Testing LabSetup pulse properties setters"""
        lab = self.lab
        time = self.time

        fields = lab.get_labfields()

        fld = lab.get_labfield(1)

        self.assertTrue(fld.om == fields[1].om)
        self.assertTrue(fld.phi == fields[1].phi)

        _plot_ = self._plot_

        if _plot_:
            fields[2].set_rwa(0.9)
            fields[2].phi = numpy.pi

            fld = fields[2].get_field()
            plt.plot(time.data, numpy.real(fld))
            plt.show()

            fields[2].tc = 300.0

            fld = fields[2].get_field()
            plt.plot(time.data, numpy.real(fld))
            plt.show()

    def test_gaussian_uses_peak_amplitude_and_intensity_fwhm_by_default(self):
        """Finite Gaussian pulses use optical pulse conventions by default."""
        time = TimeAxis(-100.0, 20001, 0.01, atype="complete")
        pulse = dict(ptype="Gaussian", FWHM=20.0, amplitude=0.3)
        lab = LabSetup(nopulses=1)

        lab.set_pulse_shapes(time, (pulse,))
        envelope = lab.pulse_t[0].data

        center = time.nearest(0.0)
        half_width = time.nearest(10.0)
        self.assertAlmostEqual(envelope[center], 0.3)
        self.assertAlmostEqual(abs(envelope[half_width]) ** 2, 0.3**2 / 2.0)

    def test_gaussian_legacy_area_and_amplitude_fwhm_are_explicit(self):
        """Explicit options reproduce the historical Gaussian definition."""
        time = TimeAxis(-100.0, 20001, 0.01, atype="complete")
        pulse = dict(
            ptype="Gaussian",
            FWHM=20.0,
            amplitude=0.3,
            amplitude_type="area",
            FWHM_type="amplitude",
        )
        lab = LabSetup(nopulses=1)

        lab.set_pulse_shapes(time, (pulse,))
        envelope = lab.pulse_t[0].data
        expected = (
            (2.0 / pulse["FWHM"])
            * numpy.sqrt(numpy.log(2.0) / numpy.pi)
            * pulse["amplitude"]
            * numpy.exp(-4.0 * numpy.log(2.0) * (time.data / pulse["FWHM"]) ** 2)
        )

        npt.assert_allclose(envelope, expected)
        self.assertAlmostEqual(numpy.sum(envelope) * time.step, 0.3)

    def test_delta_uses_area_and_accepts_deprecated_amplitude(self):
        """Delta-pulse strength is its sampled integral."""
        time = TimeAxis(-10.0, 201, 0.1, atype="complete")
        lab = LabSetup(nopulses=1)
        lab.set_pulse_shapes(time, ({"ptype": "delta", "area": 0.4},))

        self.assertAlmostEqual(numpy.sum(lab.pulse_t[0].data) * time.step, 0.4)

        with self.assertWarns(DeprecationWarning):
            lab.set_pulse_shapes(
                time,
                ({"ptype": "delta", "amplitude": 0.2},),
            )
        self.assertAlmostEqual(numpy.sum(lab.pulse_t[0].data) * time.step, 0.2)

    def test_field_phase_is_defined_at_pulse_center(self):
        """Translation preserves the configured carrier phase at the peak."""
        time = TimeAxis(-100.0, 2001, 0.1, atype="complete")
        pulse = dict(ptype="Gaussian", FWHM=20.0, amplitude=0.3)
        lab = LabSetup(nopulses=1)
        lab.set_pulse_arrival_times([12.0])
        lab.set_pulse_frequencies([0.25])
        lab.set_pulse_phases([0.4])
        lab.set_pulse_shapes(time, (pulse,))
        field = lab.get_labfield(0)

        expected = 0.3 * numpy.exp(1j * 0.4)
        npt.assert_allclose(field.field_p_at(12.0), expected)

        field.set_center(-17.0)
        npt.assert_allclose(field.field_p_at(-17.0), expected)

    def test_envelope_at_accepts_scalar_and_array_times(self):
        """Envelope evaluation preserves scalar and array input shape."""
        field = self.lab.get_labfield(2)
        times = numpy.array([[95.0, 100.0], [105.0, 110.0]])

        scalar = field.envelope_at(100.0)
        values = field.envelope_at(times)

        self.assertTrue(numpy.isscalar(scalar))
        self.assertEqual(values.shape, times.shape)
        npt.assert_allclose(values.ravel(), self.lab.pulse_t[2].at(times.ravel()))
        npt.assert_allclose(field.envelope_at(), self.lab.pulse_t[2].data)

    def test_numeric_envelope_is_complex_and_zero_outside_support(self):
        """Sampled complex envelopes have finite, zero-padded support."""
        time = TimeAxis(-2.0, 5, 1.0, atype="complete")
        data = numpy.array([0.0, 1.0 + 2.0j, 2.0 - 1.0j, 1.0j, 0.0])
        pulse = dict(ptype="numeric", function=DFunction(time, data))
        lab = LabSetup(nopulses=1)
        lab.set_pulse_shapes(time, (pulse,))
        field = lab.get_labfield(0)

        times = numpy.array([-3.0, -2.0, -0.5, 2.0, 3.0])
        expected = numpy.array([0.0, 0.0, 1.5 + 0.5j, 0.0, 0.0])

        npt.assert_allclose(field.envelope_at(times), expected)
        self.assertEqual(field.envelope_at(times).dtype, data.dtype)
        self.assertEqual(field.envelope_at(-3.0), 0.0j)
        self.assertEqual(field.envelope_at(3.0), 0.0j)

    def test_empty_chirp_is_compatible_and_nonempty_chirp_is_rejected(self):
        """A chirp must not be silently ignored by a finite-pulse setup."""
        time = TimeAxis(-2.0, 5, 1.0, atype="complete")
        lab = LabSetup(nopulses=1)
        pulse = {"ptype": "Gaussian", "FWHM": 1.0, "amplitude": 1.0, "chirp": []}

        lab.set_pulse_shapes(time, (pulse,))

        with self.assertRaisesRegex(
            QuantarheiError, "Chirped pulses are not implemented"
        ):
            lab.set_pulse_shapes(
                time,
                (
                    {
                        "ptype": "Gaussian",
                        "FWHM": 1.0,
                        "amplitude": 1.0,
                        "chirp": [0.1],
                    },
                ),
            )

    def test_get_field_at_time_wraps_envelope_at(self):
        """The legacy time-evaluation call delegates to envelope_at()."""
        field = self.lab.get_labfield(2)
        times = numpy.array([95.0, 100.0, 105.0])

        with self.assertWarns(DeprecationWarning):
            legacy = field.get_field(times)

        npt.assert_allclose(legacy, field.envelope_at(times))

    def test_field_components_follow_the_analytic_signal_convention(self):
        """Negative-frequency and real fields derive from the analytic field."""
        field = self.lab.get_labfield(2)
        times = numpy.array([95.0, 100.0, 105.0])

        field_p = field.field_p_at(times)
        field_m = field.field_m_at(times)
        real_field = field.real_field_at(times)

        npt.assert_allclose(field_m, numpy.conj(field_p))
        npt.assert_allclose(real_field, (field_p + field_m) / 2.0)
        self.assertTrue(numpy.isrealobj(real_field))

    def test_rwa_field_evaluation_is_local_and_non_mutating(self):
        """RWA evaluation changes only the local carrier detuning."""
        time = TimeAxis(-100.0, 2001, 0.1, atype="complete")
        pulse = dict(ptype="Gaussian", FWHM=20.0, amplitude=0.3)
        lab = LabSetup(nopulses=1)
        lab.set_pulse_arrival_times([12.0])
        lab.set_pulse_frequencies([0.25])
        lab.set_pulse_phases([0.4])
        lab.set_pulse_shapes(time, (pulse,))
        field = lab.get_labfield(0)
        times = numpy.array([10.0, 12.0, 14.0])
        omega_before = lab.omega.copy()

        actual = field.field_p_at(times, rwa_frequency=0.1)
        envelope = lab.pulse_t[0].at(times)
        expected = envelope * numpy.exp(-1j * (0.25 - 0.1) * (times - 12.0) + 1j * 0.4)

        npt.assert_allclose(actual, expected)
        npt.assert_allclose(lab.omega, omega_before)
        npt.assert_allclose(
            lab.get_field(0, rwa_frequency=0.1),
            field.field_p_at(rwa_frequency=0.1),
        )

    def test_carrier_frequency_accessors_use_active_energy_units(self):
        """LabSetup and LabField expose the same unit-safe carrier values."""
        lab = LabSetup(nopulses=1)
        with energy_units("1/cm"):
            lab.set_pulse_frequencies([12000.0])
            field = lab.get_labfield(0)

            self.assertEqual(lab.get_pulse_frequency(0), 12000.0)
            self.assertEqual(field.get_frequency(), 12000.0)
            self.assertEqual(field.om, 12000.0)

            field.set_frequency(12500.0)
            self.assertEqual(lab.get_pulse_frequency(0), 12500.0)
            self.assertEqual(field.om, 12500.0)

            field.om = 13000.0
            self.assertEqual(lab.get_pulse_frequency(0), 13000.0)

        with energy_units("int"):
            expected = convert(13000.0, "1/cm", "int")
            self.assertAlmostEqual(lab.get_pulse_frequency(0), expected)
            self.assertAlmostEqual(field.get_frequency(), expected)

    def test_delay_phase_storage_does_not_affect_the_field(self):
        """Legacy delay-phase storage is not part of field evaluation."""
        field = self.lab.get_labfield(2)
        original = field.field_p_at()

        self.lab.delay_phases[2] += 10.0

        npt.assert_allclose(field.field_p_at(), original)

    def test_mutating_rwa_methods_are_deprecated(self):
        """Legacy mutating RWA methods remain available during migration."""
        lab = self.lab
        omega_before = lab.omega.copy()

        with self.assertWarns(DeprecationWarning):
            lab.set_rwa(0.1)
        npt.assert_allclose(lab.omega, omega_before - 0.1)

        with self.assertWarns(DeprecationWarning):
            lab.restore_rwa()
        npt.assert_allclose(lab.omega, omega_before)

    def test_dm_propagation_with_fields(self):
        """(LabSetup) Time evolution with explicit electric field"""
        from quantarhei.qm import LindbladForm, Operator, SystemBathInteraction

        lab = self.lab
        time = self.time
        agg = self.agg

        HH = agg.get_Hamiltonian()
        DD = agg.get_TransitionDipoleMoment()
        # print(DD.data.shape)

        ops = []
        KK = Operator(dim=HH.dim)
        with eigenbasis_of(HH):
            KK.data[1, 2] = 1.0

        ops.append(KK)
        rates = []
        rates.append(1.0 / 1000.0)

        SBI = SystemBathInteraction(sys_operators=ops, rates=rates, system=agg)

        LT = LindbladForm(HH, SBI, as_operators=False)

        ef = lab.get_labfield(0)

        rhoi = agg.get_thermal_ReducedDensityMatrix()

        #######################################################################
        #
        # Relaxation time-independent + LabField
        #
        #######################################################################
        prop = ReducedDensityMatrixPropagator(
            timeaxis=time, Ham=HH, Efield=ef, Trdip=DD, RTensor=LT
        )

        #
        #
        # propagation has to be reimplemented with LabFields
        #
        rhot = prop.propagate(rhoi)
        self.assertTrue(rhot.is_in_rwa)
        # rhot.convert_from_RWA(HH)

        ef1 = ef.field

        #######################################################################
        #
        # Relaxation time-independent + field as an array
        #
        #######################################################################

        prop2 = ReducedDensityMatrixPropagator(
            timeaxis=time, Ham=HH, Efield=ef1, Trdip=DD, RTensor=LT
        )

        #
        # propagation has to be reimplemented with LabFields
        #
        rhot2 = prop2.propagate(rhoi)
        self.assertFalse(rhot2.is_in_rwa)

        _plot_ = self._plot_

        if _plot_:
            om = HH.rwa_energies[HH.rwa_indices[1]]

            with eigenbasis_of(HH):
                plt.plot(time.data, numpy.real(rhot.data[:, 1, 1]), "-b")
                plt.plot(time.data, numpy.real(rhot2.data[:, 1, 1]), "--g")
                plt.plot(time.data, numpy.real(rhot.data[:, 2, 2]), "-r")
                plt.plot(time.data, numpy.real(rhot2.data[:, 2, 2]), "--k")
                plt.plot(time.data, numpy.real(rhot.data[:, 1, 2]), "-b")
                plt.plot(time.data, numpy.real(rhot2.data[:, 1, 2]), "--g")
                # ef.set_rwa(om)
                # plt.plot(time.data, ef.field_p, "-m")
                # ef.restore_rwa()

            plt.show()

            with eigenbasis_of(HH):
                plt.plot(time.data, numpy.real(rhot.data[:, 0, 2]), "-k")
                plt.plot(
                    time.data,
                    numpy.real(rhot2.data[:, 0, 2] * numpy.exp(-1j * om * time.data)),
                    "--g",
                )

            plt.show()

        aggB = self.aggB
        time_r = self.time_r
        HH = aggB.get_Hamiltonian()
        DD = aggB.get_TransitionDipoleMoment()
        (RT, ham) = aggB.get_RelaxationTensor(
            time_r, relaxation_theory="stR", time_dependent=False
        )

        propB = ReducedDensityMatrixPropagator(
            timeaxis=time, Ham=HH, Efield=ef1, Trdip=DD, RTensor=RT
        )

        rhotB = propB.propagate(rhoi)
        self.assertFalse(rhot2.is_in_rwa)

        _plot_ = False
        if _plot_:
            with eigenbasis_of(HH):
                plt.plot(time.data, numpy.real(rhotB.data[:, 1, 1]), "-b")
                plt.plot(time.data, numpy.real(rhotB.data[:, 2, 2]), "-r")

            plt.show()

    def test_phase_setting_and_time_shifts(self):
        """(Labsetup) Phase and time-shift setting"""
        lab = LabSetup(nopulses=4)

        Nfr = 10
        time = TimeAxis(-500.0, Nfr * 1500, 1.0 / Nfr, atype="complete")

        # pulse shapes are specified below
        pulse2 = dict(ptype="Gaussian", FWHM=20, amplitude=0.1)
        params = (pulse2, pulse2, pulse2, pulse2)

        # pulse envelops
        lab.set_pulse_shapes(time, params)

        # pulse arrival times
        lab.set_pulse_arrival_times([0.0, 100.0, -80.0, 10.0])

        # pulse polarizations
        X = numpy.zeros(3, dtype=float)
        X[0] = 1
        lab.set_pulse_polarizations(
            pulse_polarizations=(X, X, X, X), detection_polarization=X
        )

        # pulse frequencies
        ome = convert(10200.0, "1/cm", "int")
        lab.set_pulse_frequencies([ome, ome, ome, ome])

        # additional phases can be also controlled
        lab.set_pulse_phases([0.0, 1.0, 0.0, 2.0])

        #
        # Test delay settings
        #
        fields = lab.get_labfields()
        cntr = fields[0].get_center()
        self.assertTrue(cntr == 0.0)

        fields[0].set_center(10.0)
        cntr = fields[0].get_center()
        self.assertTrue(cntr == 10.0)

        #
        # Test phases
        #
        phs = fields[0].get_phase()
        self.assertTrue(phs == 0.0)

        fields[0].set_phase(1.53)
        phs = fields[0].get_phase()
        self.assertTrue(phs == 1.53)

        lab.set_pulse_phases([0.6, 1.2, 0.1, 2.3])
        phs = fields[0].get_phase()
        self.assertTrue(phs == 0.6)
        phs = fields[1].get_phase()
        self.assertTrue(phs == 1.2)

        fields[3].set_phase(1.8)

        phses = lab.get_pulse_phases()
        self.assertAlmostEqual(phses, [0.6, 1.2, 0.1, 1.8])

        for ii in range(4):
            self.assertTrue(phses[ii] == fields[ii].get_phase())

        #
        #   delays again
        #
        cntrs = lab.get_pulse_arrival_times()
        for ii in range(4):
            cntr = fields[ii].get_center()
            self.assertTrue(cntr == cntrs[ii])

        #
        # Check indiviual field values
        #
        for kk in range(4):
            fld_f = fields[kk].get_field()
            fld_c = lab.get_field(kk)
            npt.assert_allclose(fld_f, fld_c)

        #
        # Check field values when resetting phase and time delay
        #
        setthis = [0.67, 1.53, 2.14, 3.0]
        for kk in range(4):
            fld = fields[kk].get_field()
            phs = fields[kk].get_phase()

            nphs = setthis[kk]
            phs_diff = nphs - phs
            fields[kk].set_phase(nphs)
            nfld = fields[kk].get_field()

            npt.assert_allclose(
                fld, nfld * numpy.exp(-1j * phs_diff), rtol=0, atol=1.0e-10
            )

        setthis = [20.0, 100.0, 80.0, 180.0]
        for kk in range(4):
            fld = fields[kk].get_field()
            cnt = fields[kk].get_center()
            dph = fields[kk].get_delay_phase()

            ncnt = setthis[kk]

            cnt_diff = ncnt - cnt

            fields[kk].set_center(ncnt)
            ndph = fields[kk].get_delay_phase()

            dph_diff = dph - ndph
            om = fields[kk].get_frequency()

            exp_diff = -cnt_diff * om

            # print(om, cnt_diff, dph_diff, exp_diff)

            npt.assert_allclose(exp_diff, dph_diff)

    def test_from_scratch_vs_by_objects(self):
        """(Labsetup) Phase and time-shift setting"""
        Nfr = 10
        time = TimeAxis(-500.0, Nfr * 1500, 1.0 / Nfr, atype="complete")
        # pulse shapes are specified below
        pulse2 = dict(ptype="Gaussian", FWHM=20, amplitude=0.1)
        params = (pulse2, pulse2, pulse2, pulse2)

        ome = convert(10200.0, "1/cm", "int")

        lab = LabSetup(nopulses=4)

        # pulse envelops
        lab.set_pulse_shapes(time, params)

        # pulse arrival times
        arr_times = [0.0, 100.0, -80.0, 10.0]
        lab.set_pulse_arrival_times(arr_times)

        # pulse frequencies
        lab.set_pulse_frequencies([ome, ome, ome, ome])

        # additional phases can be also controlled
        # lab.set_pulse_phases([0.0, 1.0, 0.0, 2.0])

        for ii in range(4):
            self.assertEqual(0.0, lab.get_pulse_phase(ii))
            fld = lab.get_labfield(ii)
            self.assertEqual(0.0, fld.get_phase())
            fld.set_phase(10.0)
            self.assertEqual(10.0, lab.get_pulse_phase(ii))

        npt.assert_allclose(arr_times, lab.get_pulse_arrival_times())

        for ii in range(4):
            tp = lab.get_pulse_arrival_time(ii)
            self.assertEqual(arr_times[ii], tp)

        tset = [-30.0, 0.0, 100.0, 200.0]
        lab.set_pulse_arrival_times(tset)

        fld1 = lab.get_field()

        _plot = False
        if _plot:
            plt.plot(time.data, numpy.real(fld1), "-r")
            plt.show()

        zeros = [0.0, 0.0, 0.0, 0.0]
        lab.set_pulse_arrival_times(zeros)

        fields = lab.get_labfields()
        for ii in range(4):
            self.assertEqual(0.0, fields[ii].get_center())
            self.assertEqual(0.0, fields[ii].get_delay_phase())

        for ii in range(4):
            fields[ii].set_center(tset[ii])

        istset = lab.get_pulse_arrival_times()
        npt.assert_allclose(istset, tset)

        fld2 = lab.get_field()

        if _plot:
            plt.plot(time.data, numpy.real(fld1), "-r")
            plt.plot(time.data, numpy.real(fld2), "-b")
            plt.show()

        npt.assert_allclose(fld1, fld2, rtol=0, atol=1.0e-7)

    def test_comparing_fields(self):
        """(Labsetup) Phase and time-shift setting"""
        Nfr = 10
        time = TimeAxis(-500.0, Nfr * 1500, 1.0 / Nfr, atype="complete")
        # pulse shapes are specified below
        pulse2 = dict(ptype="Gaussian", FWHM=20, amplitude=0.1)
        params = (pulse2, pulse2, pulse2, pulse2)

        ome = convert(10200.0, "1/cm", "int")

        lab = LabSetup(nopulses=4)

        # pulse envelops
        lab.set_pulse_shapes(time, params)

        # pulse arrival times
        arr_times = [0.0, 100.0, -80.0, 10.0]
        lab.set_pulse_arrival_times(arr_times)

        # pulse frequencies
        lab.set_pulse_frequencies([ome, ome, ome, ome])

        lab.set_pulse_phases([0.0, 1.0, 0.0, 2.0])

        fields = lab.get_labfields()

        setthis = [10.0, 20.0, 15.0, -10.0]
        for ii in range(4):
            fld = fields[ii]

            fwhm = fld.get_fwhm()
            om = fld.get_frequency()

            npt.assert_almost_equal(20.0, fwhm)
            t1 = fld.get_center()

            ef1 = fld.get_field()

            t2 = t1 + setthis[ii]
            fld.set_center(t2)

            lfc = 2.0 * numpy.log(2.0)
            kappa = numpy.exp(
                -lfc * (2.0 * (t1 - t2) * time.data - (t1**2 - t2**2)) / (fwhm**2)
                - 1j * om * (t1 - t2)
            )

            ef2 = fld.get_field()

            npt.assert_allclose(ef2, ef1 * kappa, rtol=0, atol=1.0e-7)


if __name__ == "__main__":
    unittest.main()
