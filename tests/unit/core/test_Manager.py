import unittest
from unittest.mock import patch

"""
*******************************************************************************


    Tests of the quantarhei.Manager class


*******************************************************************************
"""

import quantarhei
from quantarhei import Manager


class TestManager(unittest.TestCase):
    """Tests for the Manager class"""

    def setUp(self):
        pass

    def test_that_Manager_is_a_singleton(self):
        """Testing that Manager object is a singleton"""
        m = Manager()
        n = Manager()

        if m is not n:
            raise Exception()

    def test_assert_version_gt_passes_when_current_is_greater(self):
        """assert_version('>') must not exit when current version is strictly greater"""
        with patch.object(Manager, "version", new="1.0.0"):
            try:
                quantarhei.assert_version(">", "0.9.0")
            except SystemExit:
                self.fail(
                    "assert_version('>') called exit() even though current > required"
                )

    def test_assert_version_gt_exits_when_current_is_equal(self):
        """assert_version('>') must exit when current version equals required"""
        with patch.object(Manager, "version", new="1.0.0"):
            with self.assertRaises(SystemExit):
                quantarhei.assert_version(">", "1.0.0")

    def test_assert_version_gt_exits_when_current_is_less(self):
        """assert_version('>') must exit when current version is less than required"""
        with patch.object(Manager, "version", new="0.9.0"):
            with self.assertRaises(SystemExit):
                quantarhei.assert_version(">", "1.0.0")
