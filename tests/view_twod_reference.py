"""View stored 2D test spectra; run from the repository root with uv run."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

import quantarhei as qr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--system", choices=["chlorophyll", "tlm"], default="chlorophyll"
    )
    parser.add_argument(
        "--bath", choices=["overdamped", "underdamped"], default="underdamped"
    )
    parser.add_argument("--t2", type=int, choices=[0, 100], default=0)
    parser.add_argument("--save", help="Save a figure instead of opening a window")
    args = parser.parse_args()

    prefix = (
        "twodspectrum_underdamped"
        if args.bath == "underdamped"
        else "twodspectrum_test"
    )
    if args.system == "tlm":
        prefix = f"twod_tlm_{args.bath}"
    path = (
        Path(__file__).parent / "unit" / "spectroscopy" / f"{prefix}_data_{args.t2}.dat"
    )
    spectrum = qr.TwoDSpectrum()
    spectrum.load_data(path)
    # The .dat files store complex data only. Reconstruct the test's FFT axes.
    time = qr.TimeAxis(0.0, 50, 5.0, atype="complete")
    axis = time.get_FrequencyAxis()
    rwa = qr.convert(16807.0, "1/cm", "int")
    axis.data += rwa
    axis.start += rwa
    spectrum.set_axis_1(axis)
    spectrum.set_axis_3(axis)
    spectrum.set_t2(args.t2)
    with qr.energy_units("1/cm"):
        spectrum.plot(window=[15000, 18500, 15000, 18500], show=False)
    plt.title(f"{args.system}: {args.bath}, t₂ = {args.t2} fs (real total signal)")
    if args.save:
        plt.savefig(args.save, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()
