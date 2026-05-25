"""Benchmark 2D response calculation with jump-order approximations.

The script compares the current zero-jump transfer approximation
(``jump_order=0``) with the one-jump-ready approximation (``jump_order=1``).
For one-jump runs, the calculator selects ``FunctionStorage(config=1)`` and
splits transfer pathways into explicit single-jump and remainder contributions.

Run from the repository root, for example:

    python examples/benchmark_twod_storage.py --nt 32 --nt2 12 --dt 10.0
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path
from typing import Any

os.environ.setdefault(
    "MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "quantarhei-mpl-cache")
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import quantarhei as qr
from quantarhei.spectroscopy import X


def _rss_mb() -> float:
    """Return peak resident set size in MB on macOS/Linux."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024.0 * 1024.0)
    return rss / 1024.0


def _build_system(nt: int, nt2: int, dt: float) -> Any:
    env_axis = qr.TimeAxis(0.0, nt, dt)

    with qr.energy_units("1/cm"):
        m1 = qr.Molecule([0.0, 12000.0])
        m2 = qr.Molecule([0.0, 12300.0])

        params = dict(
            ftype="OverdampedBrownian",
            reorg=40.0,
            cortime=100.0,
            T=100,
            matsubara=20,
        )
        m1.set_transition_environment((0, 1), qr.CorrelationFunction(env_axis, params))
        m2.set_transition_environment((0, 1), qr.CorrelationFunction(env_axis, params))

    m1.set_dipole(0, 1, [1.0, 0.8, 0.8])
    m2.set_dipole(0, 1, [0.8, 0.8, 0.0])

    agg = qr.Aggregate(molecules=[m1, m2])
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0, 1, 100.0)

    agg.build(mult=2)
    agg.diagonalize()
    return agg


def _s2_point_count(t2_axis: Any, jump_time_graining: int) -> int:
    """Return the number of s2 grid points used for the largest t2."""
    if t2_axis.length == 0:
        return 0
    t2_index = t2_axis.length - 1
    indices = list(range(0, t2_index + 1, jump_time_graining))
    if not indices or indices[-1] != t2_index:
        indices.append(t2_index)
    return len(indices)


def _run_single(
    jump_order: int,
    nt: int,
    nt2: int,
    dt: float,
    jump_time_graining: int,
    jump_kernel_cutoff: float,
    jump_kernel_zero_cutoff: float,
) -> dict[str, Any]:
    t1_axis = qr.TimeAxis(0.0, nt, dt)
    t2_axis = qr.TimeAxis(0.0, nt2, dt)
    t3_axis = qr.TimeAxis(0.0, nt, dt)

    system = _build_system(nt, nt2, dt)

    lab = qr.LabSetup()
    lab.set_pulse_polarizations(pulse_polarizations=(X, X, X), detection_polarization=X)

    tracemalloc.start()
    start = time.perf_counter()

    calc = qr.TwoDResponseCalculator(
        t1_axis,
        t2_axis,
        t3_axis,
        system=system,
        jump_order=jump_order,
        jump_time_graining=jump_time_graining,
        jump_kernel_cutoff=jump_kernel_cutoff,
        jump_kernel_zero_cutoff=jump_kernel_zero_cutoff,
    )
    with qr.energy_units("1/cm"):
        calc.bootstrap(rwa=12100.0, pad=0, lab=lab)

    gg = system.get_lineshape_functions()
    storage_size_mb = float(gg.data_size)
    storage_stride = int(gg.data_stride)
    storage_arguments = len(gg.time_mapping)

    response = calc.calculate()
    container = response.get_TwoDSpectrumContainer()
    spectrum = container.get_spectrum(t2_axis.data[0])
    checksum = float(abs(spectrum.data).sum())
    jump_diagnostics = [
        diag.get("jump", {})
        for diag in calc.response_diagnostics
        if str(diag.get("transfer_channel", "")).startswith("scM1")
    ]
    skipped_jumps = sum(1 for diag in jump_diagnostics if diag.get("skipped", False))
    used_s2_points = sum(int(diag.get("used_points", 0)) for diag in jump_diagnostics)
    remainder_changes = [
        diag["remainder_relative_change"]
        for diag in calc.response_diagnostics
        if diag.get("remainder_relative_change") is not None
    ]
    max_remainder_change = max(remainder_changes, default=0.0)

    elapsed = time.perf_counter() - start
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "jump_order": jump_order,
        "jump_time_graining": jump_time_graining,
        "jump_kernel_cutoff": jump_kernel_cutoff,
        "jump_kernel_zero_cutoff": jump_kernel_zero_cutoff,
        "s2_point_count": _s2_point_count(t2_axis, jump_time_graining),
        "used_s2_points": used_s2_points,
        "skipped_jumps": skipped_jumps,
        "max_remainder_change": max_remainder_change,
        "nt": nt,
        "nt2": nt2,
        "dt": dt,
        "elapsed_s": elapsed,
        "tracemalloc_current_mb": current / (1024.0 * 1024.0),
        "tracemalloc_peak_mb": peak / (1024.0 * 1024.0),
        "max_rss_mb": _rss_mb(),
        "storage_data_size_mb": storage_size_mb,
        "storage_stride": storage_stride,
        "storage_arguments": storage_arguments,
        "checksum": checksum,
    }


def _print_table(results: list[dict[str, Any]]) -> None:
    print(
        "jump  grain  cutoff  zero  s2pts  used  skip  rem_dR  args  stride  storage_MB  elapsed_s  trace_peak_MB  max_RSS_MB  checksum"
    )
    for result in results:
        print(
            f"{result['jump_order']:>4}  "
            f"{result['jump_time_graining']:>5}  "
            f"{result['jump_kernel_cutoff']:>6.1e}  "
            f"{result['jump_kernel_zero_cutoff']:>4.1e}  "
            f"{result['s2_point_count']:>5}  "
            f"{result['used_s2_points']:>4}  "
            f"{result['skipped_jumps']:>4}  "
            f"{result['max_remainder_change']:>6.1e}  "
            f"{result['storage_arguments']:>4}  "
            f"{result['storage_stride']:>6}  "
            f"{result['storage_data_size_mb']:>10.3f}  "
            f"{result['elapsed_s']:>9.3f}  "
            f"{result['tracemalloc_peak_mb']:>13.3f}  "
            f"{result['max_rss_mb']:>10.3f}  "
            f"{result['checksum']:.6e}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nt", type=int, default=24)
    parser.add_argument("--nt2", type=int, default=8)
    parser.add_argument("--dt", type=float, default=10.0)
    parser.add_argument("--jump-time-graining", type=int, default=1)
    parser.add_argument("--jump-kernel-cutoff", type=float, default=0.0)
    parser.add_argument("--jump-kernel-zero-cutoff", type=float, default=0.0)
    parser.add_argument("--single-jump-order", type=int, choices=(0, 1))
    args = parser.parse_args()

    if args.single_jump_order is not None:
        print(
            json.dumps(
                _run_single(
                    args.single_jump_order,
                    args.nt,
                    args.nt2,
                    args.dt,
                    args.jump_time_graining,
                    args.jump_kernel_cutoff,
                    args.jump_kernel_zero_cutoff,
                )
            )
        )
        return

    results = []
    for jump_order in (0, 1):
        command = [
            sys.executable,
            __file__,
            "--single-jump-order",
            str(jump_order),
            "--nt",
            str(args.nt),
            "--nt2",
            str(args.nt2),
            "--dt",
            str(args.dt),
            "--jump-time-graining",
            str(args.jump_time_graining),
            "--jump-kernel-cutoff",
            str(args.jump_kernel_cutoff),
            "--jump-kernel-zero-cutoff",
            str(args.jump_kernel_zero_cutoff),
        ]
        completed = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
        results.append(json.loads(completed.stdout.strip().splitlines()[-1]))

    _print_table(results)


if __name__ == "__main__":
    main()
