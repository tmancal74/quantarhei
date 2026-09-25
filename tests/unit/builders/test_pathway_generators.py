"""Regression tests for the third-order Liouville pathway generators.

The ten ``generate_*`` functions in
:mod:`quantarhei.builders.aggregate_spectroscopy` (R1g, R2g, R3g, R4g, R1f*,
R2f*, R1gE, R2gE, R1f*E, R2f*E) were rewritten from hand-written nested loops
into descriptor-driven generators (issue #366). These tests pin the *full*
content of every generated pathway (type, name, transitions, sides, states,
relaxations, evolution factor, sign, bands, dipoles, energies, frequencies,
widths, dephasings, orientational-averaging vector) against the output of the
original, unrefactored implementation.

Test systems
------------
``excitonic``
    Trimer of two-level molecules built with ``mult=2``, i.e. with a
    two-exciton band (exercises R1f*/R2f* ESA detection in band 2).
``vibronic``
    Dimer with one underdamped mode per molecule (two vibrational levels in
    each electronic state, different Huang-Rhys factors), built with
    ``mult=2``. Several vibrational ground states are populated in ``rho0``,
    so the ground-state bleach/relaxation variants (R3g, R4g, R*E) iterate
    over more than one ground state. One ground state is deliberately kept
    below the population tolerance.

Both systems are used in the site (non-diagonalized) basis so that the
reference data do not depend on the eigenvector phase convention of the
local LAPACK build. The evolution superoperator ``eUt2`` is a sparse complex
tensor defined by a deterministic integer formula (no RNG), with magnitudes
straddling ``EVF_TOL`` so the evolution-amplitude filter is exercised.

Reference fixture
-----------------
``test_pathway_generators.npz`` was produced from upstream master commit
0b7d6d9 (the last commit before the refactoring) with::

    git worktree add --detach /tmp/wt401-master 0b7d6d9
    PYTHONPATH=/tmp/wt401-master python \\
        tests/unit/builders/test_pathway_generators.py --write-fixture
    git worktree remove /tmp/wt401-master

The script checks that it imports quantarhei from the worktree. Do not
regenerate the fixture from refactored code: that would make the test
compare the implementation with itself.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy
import pytest

FIXTURE = Path(__file__).with_name("test_pathway_generators.npz")

POP_TOL = 1.0e-3
DIP_TOL = 0.05
EVF_TOL = 0.03

# (fixture key, generator function name, uses evolution-amplitude tolerance)
GENERATORS = (
    ("R1g", "generate_R1g", True),
    ("R2g", "generate_R2g", True),
    ("R3g", "generate_R3g", False),
    ("R4g", "generate_R4g", False),
    ("R1f*", "generate_R1f", True),
    ("R2f*", "generate_R2f", True),
    ("R1gE", "generate_R1gE", True),
    ("R2gE", "generate_R2gE", True),
    ("R1f*E", "generate_R1fE", True),
    ("R2f*E", "generate_R2fE", True),
)

SYSTEMS = ("excitonic", "vibronic")

INT_FIELDS = ("sinit", "transitions", "sides", "states", "relaxations")
FLOAT_FIELDS = (
    "sign",
    "popt_band",
    "relax_order",
    "dmoments",
    "energy",
    "frequency",
    "widths",
    "dephs",
    "F4n",
)
STR_FIELDS = ("pathway_type", "pathway_name", "event")


def _build_excitonic() -> Any:
    from quantarhei import Aggregate, Molecule

    energies = (1.00, 1.05, 0.97)
    dipoles = ([1.0, 0.0, 0.0], [0.6, 0.8, 0.0], [0.2, -0.3, 0.9])
    widths = (0.010, 0.012, 0.015)
    dephs = (0.020, 0.018, 0.025)
    mols = []
    for k in range(3):
        m = Molecule(name=f"M{k}", elenergies=[0.0, energies[k]])
        m.set_dipole(0, 1, dipoles[k])
        m.position = [float(k), 0.0, 0.0]
        m.set_transition_width((0, 1), widths[k])
        m.set_transition_dephasing((0, 1), dephs[k])
        mols.append(m)
    agg = Aggregate(molecules=mols)
    agg.set_resonance_coupling(0, 1, 0.010)
    agg.set_resonance_coupling(1, 2, 0.007)
    agg.set_resonance_coupling(0, 2, -0.003)
    agg.build(mult=2)
    rho0 = numpy.zeros(agg.HH.shape, dtype=complex)
    rho0[0, 0] = 1.0
    agg.rho0 = rho0
    return agg


def _build_vibronic() -> Any:
    from quantarhei import Aggregate, Mode, Molecule

    energies = (1.00, 1.04)
    dipoles = ([1.0, 0.0, 0.0], [0.5, 0.7, 0.0])
    freqs = (0.010, 0.013)
    hrs = (0.3, 0.6)
    widths = (0.011, 0.014)
    dephs = (0.021, 0.017)
    mols = []
    for k in range(2):
        m = Molecule(name=f"V{k}", elenergies=[0.0, energies[k]])
        m.set_dipole(0, 1, dipoles[k])
        m.position = [float(k), 0.0, 0.0]
        m.set_transition_width((0, 1), widths[k])
        m.set_transition_dephasing((0, 1), dephs[k])
        mode = Mode(freqs[k])
        m.add_Mode(mode)
        mode.set_nmax(0, 2)
        mode.set_nmax(1, 2)
        mode.set_HR(1, hrs[k])
        mols.append(m)
    agg = Aggregate(molecules=mols)
    agg.set_resonance_coupling(0, 1, 0.008)
    agg.build(mult=2)
    ngs = len(agg.get_electronic_groundstate())
    assert ngs == 4
    rho0 = numpy.zeros(agg.HH.shape, dtype=complex)
    # three populated ground states, one below POP_TOL
    for idx, pop in enumerate((0.55, 0.30, 0.15, 1.0e-4)):
        rho0[idx, idx] = pop
    agg.rho0 = rho0
    return agg


def _sparse_eUt2(dim: int) -> numpy.ndarray:
    """Deterministic sparse complex tensor with magnitudes around EVF_TOL."""
    idx = numpy.indices((dim, dim, dim, dim))
    a, b, c, d = idx
    key = 7 * a + 11 * b + 13 * c + 19 * d
    mask = (key % 3) == 0
    mag = 0.01 * (1 + (a + 2 * b + 3 * c + 4 * d) % 8)  # 0.01 ... 0.08
    phase = 0.37 * (a - b + 2 * c - 3 * d)
    out = numpy.zeros((dim,) * 4, dtype=complex)
    out[mask] = (mag * numpy.exp(1j * phase))[mask]
    return out


def _system(name: str) -> Any:
    agg = _build_excitonic() if name == "excitonic" else _build_vibronic()
    return agg, _sparse_eUt2(agg.HH.shape[0])


def _generate(agg: Any, eUt2: numpy.ndarray, func_name: str, evf: bool) -> list:
    from quantarhei.builders import aggregate_spectroscopy as asp

    func = getattr(asp, func_name)
    pathways: list = []
    if evf:
        func(agg, pathways, eUt2, POP_TOL, DIP_TOL, EVF_TOL, 0)
    else:
        func(agg, pathways, eUt2, POP_TOL, DIP_TOL, 0)
    return pathways


def _serialize(pathways: list) -> dict[str, numpy.ndarray]:
    """Stack the relevant attributes of all pathways into arrays."""
    fields: dict[str, list] = {
        f: [] for f in INT_FIELDS + FLOAT_FIELDS + STR_FIELDS + ("evolfac",)
    }
    for p in pathways:
        fields["sinit"].append(numpy.asarray(p.sinit))
        fields["transitions"].append(numpy.asarray(p.transitions))
        fields["sides"].append(numpy.asarray(p.sides))
        fields["states"].append(numpy.asarray(p.states))
        fields["relaxations"].append(
            numpy.asarray([list(r[0]) + list(r[1]) for r in p.relaxations])
        )
        fields["evolfac"].append(complex(p.evolfac))
        fields["sign"].append(float(p.sign))
        fields["popt_band"].append(float(p.popt_band))
        fields["relax_order"].append(float(p.relax_order))
        fields["dmoments"].append(numpy.asarray(p.dmoments))
        fields["energy"].append(numpy.asarray(p.energy))
        fields["frequency"].append(numpy.asarray(p.frequency))
        fields["widths"].append(numpy.asarray(p.widths, dtype=float))
        fields["dephs"].append(numpy.asarray(p.dephs, dtype=float))
        fields["F4n"].append(numpy.asarray(p.F4n))
        fields["pathway_type"].append(str(p.pathway_type))
        fields["pathway_name"].append(str(p.pathway_name))
        fields["event"].append("".join(str(e) for e in p.event))
    out: dict[str, numpy.ndarray] = {}
    for f, vals in fields.items():
        if f in STR_FIELDS:
            out[f] = numpy.asarray(vals, dtype=str)
        elif f in INT_FIELDS:
            out[f] = numpy.asarray(vals, dtype=numpy.int64)
        elif f == "evolfac":
            out[f] = numpy.asarray(vals, dtype=complex)
        else:
            out[f] = numpy.asarray(vals, dtype=float)
    return out


@pytest.fixture(scope="module")
def reference() -> Any:
    with numpy.load(FIXTURE) as data:
        return {k: data[k] for k in data.files}


@pytest.fixture(scope="module")
def systems() -> dict[str, Any]:
    return {name: _system(name) for name in SYSTEMS}


@pytest.mark.parametrize("sysname", SYSTEMS)
@pytest.mark.parametrize(
    ("ptype", "func_name", "evf"), GENERATORS, ids=[g[0] for g in GENERATORS]
)
def test_pathways_match_master(
    reference: dict, systems: dict, sysname: str, ptype: str, func_name: str, evf: bool
) -> None:
    agg, eUt2 = systems[sysname]
    got = _serialize(_generate(agg, eUt2, func_name, evf))
    prefix = f"{sysname}/{ptype}/"
    count = int(reference[prefix + "count"])
    assert len(got["sign"]) == count
    for field in INT_FIELDS + STR_FIELDS:
        numpy.testing.assert_array_equal(
            got[field], reference[prefix + field], err_msg=field
        )
    # Values are copied from the same Hamiltonian/dipole/width arrays and the
    # same eUt2 entries, so only floating-point round-off is tolerated.
    for field in FLOAT_FIELDS + ("evolfac",):
        numpy.testing.assert_allclose(
            got[field],
            reference[prefix + field],
            rtol=1e-12,
            atol=1e-14,
            err_msg=field,
        )


def test_fixture_is_nontrivial(reference: dict) -> None:
    """Each generator must produce pathways in at least one test system,
    and the vibronic system must start pathways from several ground states."""
    for ptype, _, _ in GENERATORS:
        total = sum(int(reference[f"{s}/{ptype}/count"]) for s in SYSTEMS)
        assert total > 0, ptype
    for ptype in ("R1g", "R3g", "R4g", "R1gE"):
        sinit = reference[f"vibronic/{ptype}/sinit"]
        assert len(numpy.unique(sinit[:, 0])) >= 3, ptype


def test_band2_missing_raises_quantarhei_error() -> None:
    from quantarhei.builders import aggregate_spectroscopy as asp
    from quantarhei.exceptions import QuantarheiError

    agg = _build_excitonic()

    def no_band(band: int = 1) -> tuple:
        if band == 2:
            raise IndexError("no band 2")
        return tuple(range(1, 4))

    agg.get_excitonic_band = no_band  # type: ignore[method-assign]
    eUt2 = _sparse_eUt2(agg.HH.shape[0])
    with pytest.raises(QuantarheiError, match="R1f\\* pathway generation"):
        asp.generate_R1f(agg, [], eUt2, POP_TOL, DIP_TOL, EVF_TOL, 0)
    with pytest.raises(QuantarheiError, match="R2f\\* pathway generation"):
        asp.generate_R2f(agg, [], eUt2, POP_TOL, DIP_TOL, EVF_TOL, 0)


# Behaviour of the original generators when building a single pathway raises:
# a QuantarheiError with the given message, or (None) silently skipping it.
FAILURE_MESSAGES = {
    "R1g": "Pathway generation failed",
    "R2g": "",
    "R3g": "Generation of pathway failed",
    "R4g": None,
    "R1f*": "Constructionrelaxation pathway failed",
    "R2f*": None,
    "R1gE": "",
    "R2gE": "",
    "R1f*E": "Constructionrelaxation pathway failed",
    "R2f*E": None,
}


@pytest.mark.parametrize(
    ("ptype", "func_name", "evf"), GENERATORS, ids=[g[0] for g in GENERATORS]
)
def test_pathway_construction_failure(
    monkeypatch: pytest.MonkeyPatch, ptype: str, func_name: str, evf: bool
) -> None:
    from quantarhei.exceptions import QuantarheiError
    from quantarhei.spectroscopy import diagramatics

    def broken(*args: Any, **kwargs: Any) -> Any:
        raise ValueError("boom")

    monkeypatch.setattr(diagramatics, "liouville_pathway", broken)
    agg, eUt2 = _system("excitonic")
    message = FAILURE_MESSAGES[ptype]
    if message is None:
        assert _generate(agg, eUt2, func_name, evf) == []
    else:
        with pytest.raises(QuantarheiError) as info:
            _generate(agg, eUt2, func_name, evf)
        assert str(info.value) == message


def _write_fixture() -> None:  # pragma: no cover - maintenance helper
    import quantarhei

    print("Using quantarhei from", quantarhei.__file__)
    data: dict[str, numpy.ndarray] = {}
    for sysname in SYSTEMS:
        agg, eUt2 = _system(sysname)
        for ptype, func_name, evf in GENERATORS:
            ser = _serialize(_generate(agg, eUt2, func_name, evf))
            prefix = f"{sysname}/{ptype}/"
            data[prefix + "count"] = numpy.asarray(len(ser["sign"]))
            for k, v in ser.items():
                data[prefix + k] = v
            print(f"{sysname:10s} {ptype:6s} {len(ser['sign']):5d} pathways")
    numpy.savez_compressed(FIXTURE, **data)
    print("Wrote", FIXTURE)


if __name__ == "__main__":  # pragma: no cover
    if "--write-fixture" in sys.argv:
        _write_fixture()
