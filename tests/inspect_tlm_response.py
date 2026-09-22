"""Inspect TLM responses: uv run --no-sync python -m tests.inspect_tlm_response."""

import argparse

import matplotlib.pyplot as plt
import numpy as np

import quantarhei as qr
from tests.unit.spectroscopy.twod_tlm_test import make_tlm_calculator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--component",
        choices=["rephasing", "nonrephasing", "R2g", "R3g", "R1g", "R4g"],
        default="rephasing",
    )
    parser.add_argument("--save", help="Save the plot instead of opening a window")
    args = parser.parse_args()
    calc = make_tlm_calculator(underdamped=True)
    responses = calc.calculate()
    dtype = {"rephasing": qr.signal_REPH, "nonrephasing": qr.signal_NONR}.get(
        args.component, args.component
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    frequency = qr.convert(calc.oa1.data, "int", "1/cm")
    for index, tau in enumerate((0.0, 100.0)):
        response = responses.get_response(tau)
        for name in (qr.signal_REPH, qr.signal_NONR, "R2g", "R3g", "R1g", "R4g"):
            response.set_data_flag(name)
            print(
                f"t2={tau:g} fs, {name}: max |spectrum| = {np.max(np.abs(response.data)):.6g}"
            )
        for name in ("rTot", "nTot", "rSE", "rGSB"):
            data = calc.responses[index][name]
            print(f"  time-domain {name}: max |response| = {np.max(np.abs(data)):.6g}")
        response.set_data_flag(dtype)
        data = response.data.real
        vmax = np.max(np.abs(data))
        ax = axes[index]
        im = ax.pcolormesh(
            frequency,
            frequency,
            data,
            shading="auto",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
        )
        ax.set(
            title=f"{args.component}, t₂={tau:g} fs (real)",
            xlabel="Excitation (cm⁻¹)",
            ylabel="Detection (cm⁻¹)",
        )
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    main()
