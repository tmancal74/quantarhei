import os
import unittest
from importlib.resources import files as resource_files

import numpy

from quantarhei import Input


class TestInput(unittest.TestCase):
    """Tests for the units package


    """

    def test_input(self):
        """Testing Input class
        """
        files = ["input.json", "input.yaml"]

        for fl in files:

            extsplt = os.path.splitext(fl)

            filejson = str(resource_files(__package__).joinpath(fl))

            inp = Input(filejson)

            ecalc = inp.ECalc

            if extsplt[1] in [".json"]:
                wd0 = inp.width_dist["0"]
                wd1 = inp.width_dist["1"]
            elif extsplt[1] in [".yaml", ".yml"]:
                wd0 = inp.width_dist[0]
                wd1 = inp.width_dist[1]
            make_movie = inp.make_movie

            numpy.testing.assert_equal(wd0, 50.0)
            numpy.testing.assert_equal(wd1, 50.0)
            numpy.testing.assert_equal(ecalc, 10000.0)
            numpy.testing.assert_equal(make_movie, False)
