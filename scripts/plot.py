"""Plot results CSVs.

Usage: python3 scripts/plot.py results/example_spectrum.csv [--band 10 500] [-o out.png]

Milestone 3.2.1 covers the spectrum plot only (one window of P[k] in dB).
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE = "#fcfcfb"
SERIES = "#2a78d6"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"


def read_spectrum(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return [float(r["freq_hz"]) for r in rows], [float(r["power"]) for r in rows]


def plot_spectrum(path, band, out):
    import math

    freq, power = read_spectrum(path)
    power_db = [10.0 * math.log10(p) if p > 0 else float("nan") for p in power]

    fig, ax = plt.subplots(figsize=(9, 4), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axvspan(band[0], band[1], color=GRID, alpha=0.5, lw=0, zorder=0)
    ax.text(band[0], 0.97, f" search band {band[0]:g}-{band[1]:g} Hz", color=INK_2,
            fontsize=9, transform=ax.get_xaxis_transform(), va="top")
    ax.plot(freq, power_db, color=SERIES, lw=1.0, zorder=2)

    k = max(range(len(power)), key=power.__getitem__)
    ax.annotate(f"peak {freq[k]:.2f} Hz", (freq[k], power_db[k]), xytext=(8, -4),
                textcoords="offset points", color=INK, fontsize=9)

    ax.set_xlim(0, min(freq[-1], band[1] * 1.2))
    ax.set_xlabel("Frequency (Hz)", color=INK_2)
    ax.set_ylabel("P[k] (dB re sigma^2)", color=INK_2)
    ax.set_title("Hann periodogram, one window", color=INK, loc="left")
    ax.grid(True, color=GRID, lw=0.6)
    ax.tick_params(colors=MUTED)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)

    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("csv", type=Path)
    ap.add_argument("--band", nargs=2, type=float, default=(10.0, 1000.0),
                    metavar=("F_LO", "F_HI"))
    ap.add_argument("-o", "--out", type=Path)
    args = ap.parse_args()
    plot_spectrum(args.csv, args.band, args.out or args.csv.with_suffix(".png"))


if __name__ == "__main__":
    main()
