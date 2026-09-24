"""Tests for exciton analysis reports on aggregates."""

import io

import quantarhei as qr


def _dimer():
    agg = qr.TestAggregate("dimer-2")
    agg.set_coupling_by_dipole_dipole()
    agg.build()
    agg.diagonalize()
    return agg


def test_report_on_expansion_accepts_state_positionally(capsys):
    """``report_on_expansion(state)`` reports on ``state`` (issue #335)."""
    agg = _dimer()
    agg.report_on_expansion(1)
    positional = capsys.readouterr().out
    agg.report_on_expansion(state=1)
    keyword = capsys.readouterr().out
    assert positional == keyword
    assert "0.90998848" in positional


def test_report_on_expansion_writes_to_file():
    agg = _dimer()
    buf = io.StringIO()
    agg.report_on_expansion(1, N=2, file=buf)
    lines = buf.getvalue().splitlines()
    # header rule, header, rule, N rows, closing rule
    assert len(lines) == 6
    assert "0.90998848" in lines[3]
