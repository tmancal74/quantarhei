import numpy
import pytest

from quantarhei.builders.interactions import HARTREE_TO_INVCM, dipole_dipole

# CODATA 2018 hartree energy expressed in wavenumbers
HARTREE_IN_INVCM = 219474.6313632


def test_hartree_to_invcm_matches_codata():
    # units.py uses the CODATA 2014 hartree (27.21138602 eV), 1e-8 relative
    # below the CODATA 2018 value
    assert HARTREE_TO_INVCM == pytest.approx(HARTREE_IN_INVCM, rel=1e-7)


@pytest.mark.parametrize(
    ("dipole", "expected_hartree"),
    [
        ([0.0, 0.0, 1.0], 1.0e-3),
        ([1.0, 0.0, 0.0], -2.0e-3),
    ],
    ids=["side_by_side", "head_to_tail"],
)
def test_dipole_dipole_units(dipole, expected_hartree):
    center1 = numpy.zeros(3)
    center2 = numpy.array([10.0, 0.0, 0.0])
    dipole = numpy.array(dipole)

    args = (center1, dipole, center2, dipole)
    assert dipole_dipole(*args, "Hartree") == pytest.approx(expected_hartree)
    assert dipole_dipole(*args, "cm-1") == pytest.approx(
        expected_hartree * HARTREE_IN_INVCM, rel=1e-7
    )
    assert dipole_dipole(*args) == dipole_dipole(*args, "cm-1")
