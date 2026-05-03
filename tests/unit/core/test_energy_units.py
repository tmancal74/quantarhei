import unittest

"""
*******************************************************************************


    Tests of the quantarhei.energy_units class


*******************************************************************************
"""

from quantarhei import Manager, energy_units, set_current_units
from quantarhei.core.units import cm2int, conversion_facs_energy

# let us reuse a class from previous test
from .test_UnitsManaged import UnitsManagedObject


class TestEnergyUnits(unittest.TestCase):
    def setUp(self):
        # reusining class from previous test
        self.u = UnitsManagedObject()

    def test_using_different_units(self):
        """Testing that 'energy_units' context manager works"""
        # set value in 1/cm
        with energy_units("1/cm"):
            val = 100.0
            val2 = 50.0
            self.u.set_value(val)

            w = UnitsManagedObject()
            w.set_value(val2)

        # compare it in internal units
        self.assertEqual(self.u.value, val * cm2int)
        self.assertEqual(w.value, val2 * cm2int)

    def test_using_global_units_testing(self):
        """Testing that 'set_current_units' works"""
        a = UnitsManagedObject()
        a.set_value(1.0)

        #
        # Test of energy units
        #
        units = dict(energy="1/cm")

        set_current_units(units)

        self.assertEqual(a.get_value(), 1.0 / cm2int)

        val1 = 100
        val2 = 50.0
        w = UnitsManagedObject()
        w.set_value(val1)
        self.u.set_value(val2)

        self.assertEqual(self.u.get_value(), val2)
        self.assertEqual(w.get_value(), val1)

        set_current_units()

        self.assertEqual(self.u.value, val2 * cm2int)
        self.assertEqual(w.value, val1 * cm2int)

        self.assertEqual(self.u.get_value(), val2 * cm2int)
        self.assertEqual(w.get_value(), val1 * cm2int)

        #
        # Test of frequency units
        #
        units = dict(frequency="1/cm")

        set_current_units(units)

        self.assertEqual(a.get_value(), 1.0 / cm2int)

        val1 = 100
        val2 = 50.0
        w = UnitsManagedObject()
        w.set_value(val1)
        self.u.set_value(val2)

        self.assertEqual(self.u.get_value(), val2)
        self.assertEqual(w.get_value(), val1)

        set_current_units()

        self.assertEqual(self.u.value, val2 * cm2int)
        self.assertEqual(w.value, val1 * cm2int)

        self.assertEqual(self.u.get_value(), val2 * cm2int)
        self.assertEqual(w.get_value(), val1 * cm2int)

    def test_energy_units_state_restored_after_exception(self):
        """energy_units.__exit__ must restore Manager state even if body raises"""
        manager = Manager()

        try:
            with energy_units("1/cm"):
                raise RuntimeError("body error")
        except RuntimeError:
            pass

        self.assertEqual(manager._in_eu_count, 0)
        self.assertFalse(manager._in_energy_units_context)

    def test_energy_units_nested_state_restored_after_exception(self):
        """Nested energy_units: inner exception must not corrupt outer counter"""
        manager = Manager()

        with energy_units("1/cm"):
            try:
                with energy_units("eV"):
                    raise RuntimeError("inner error")
            except RuntimeError:
                pass

            self.assertEqual(manager._in_eu_count, 1)
            self.assertTrue(manager._in_energy_units_context)

        self.assertEqual(manager._in_eu_count, 0)
        self.assertFalse(manager._in_energy_units_context)
    def test_cu_energy_roundtrip(self):
        """cu_energy() must convert input units -> current units correctly.

        100 1/cm expressed in eV, then cu_energy(..., units="eV") with
        current_units="1/cm" must return 100.0.  Before the fix this
        returned the original eV value unchanged (identity bug).
        """
        m = Manager()
        cm_fac = conversion_facs_energy["1/cm"]
        eV_fac = conversion_facs_energy["eV"]

        # 100 1/cm in internal units, then expressed as eV
        val_eV = (100.0 * cm_fac) / eV_fac

        set_current_units({"energy": "1/cm"})
        try:
            result = m.cu_energy(val_eV, units="eV")
            self.assertAlmostEqual(result, 100.0, places=8)
        finally:
            set_current_units()
