import unittest

import numpy

from quantarhei import SpectralDensity, TimeAxis, energy_units


class TestSpectralDensityRoundTrip(unittest.TestCase):
    """Regression test for #393: SD -> CF -> SD round-trip accuracy."""

    def _roundtrip_error(self, params, time, temperature=300):
        with energy_units("1/cm"):
            sd_orig = SpectralDensity(time, params)
            cf = sd_orig.get_CorrelationFunction(temperature=temperature)
            sd_back = cf.get_SpectralDensity()

        with energy_units("int"):
            freq = sd_orig.axis.data
            positive = freq > 0
            orig_data = sd_orig.data[positive]
            back_data = sd_back.data[positive]

        max_val = numpy.max(numpy.abs(orig_data))
        if max_val == 0:
            return 0.0
        return numpy.max(numpy.abs(orig_data - back_data)) / max_val

    def test_overdamped_brownian_roundtrip(self):
        time = TimeAxis(0.0, 10000, 1.0)
        params = dict(ftype="OverdampedBrownian", reorg=200.0, cortime=100.0, T=300.0)

        rdiff = self._roundtrip_error(params, time)
        self.assertLess(rdiff, 0.05)

    def test_underdamped_brownian_roundtrip(self):
        time = TimeAxis(0.0, 10000, 1.0)
        params = dict(
            ftype="UnderdampedBrownian",
            reorg=50.0,
            freq=500.0,
            gamma=1.0 / 200.0,
            T=300.0,
        )

        rdiff = self._roundtrip_error(params, time)
        self.assertLess(rdiff, 0.05)

    def test_higher_resolution_improves_accuracy(self):
        params = dict(ftype="OverdampedBrownian", reorg=200.0, cortime=100.0, T=300.0)

        time_coarse = TimeAxis(0.0, 2000, 2.0)
        time_fine = TimeAxis(0.0, 10000, 1.0)

        err_coarse = self._roundtrip_error(params, time_coarse)
        err_fine = self._roundtrip_error(params, time_fine)

        self.assertLess(err_fine, err_coarse)


if __name__ == "__main__":
    unittest.main()
