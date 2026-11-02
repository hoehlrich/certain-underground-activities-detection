"""Plot results CSVs.

Usage: python3 scripts/plot.py results/example_spectrum.csv [--band 10 1000] [-o out.png]

Draws one window's P[k] in dB with the CA-CFAR threshold and detections.
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE = "#fcfcfb"
SERIES = "#2a78d6"     # categorical slot 1
THRESHOLD = "#eb6834"  # categorical slot 2
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"


def to_db(v):
    return 10.0 * math.log10(v) if v > 0 else float("nan")


def read_spectrum(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    freq = [float(r["freq_hz"]) for r in rows]
    power = [float(r["power"]) for r in rows]
    threshold = [float(r["threshold"]) if r.get("threshold") else float("nan") for r in rows]
    detect = [r.get("detect") == "1" for r in rows]
    return freq, power, threshold, detect


def plot_spectrum(path, band, out):
    freq, power, threshold, detect = read_spectrum(path)
    power_db = [to_db(p) for p in power]
    thresh_db = [to_db(t) if not math.isnan(t) else float("nan") for t in threshold]

    fig, ax = plt.subplots(figsize=(9, 4), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axvspan(band[0], band[1], color=GRID, alpha=0.5, lw=0, zorder=0)
    ax.text(band[0], 0.97, f" search band {band[0]:g}-{band[1]:g} Hz", color=INK_2,
            fontsize=9, transform=ax.get_xaxis_transform(), va="top")

    ax.plot(freq, power_db, color=SERIES, lw=1.0, zorder=2, label="P[k]")
    if any(not math.isnan(t) for t in thresh_db):
        ax.plot(freq, thresh_db, color=THRESHOLD, lw=2.0, zorder=3, label="CA-CFAR threshold")

    hits = [k for k, d in enumerate(detect) if d]
    if hits:
        ax.scatter([freq[k] for k in hits], [power_db[k] for k in hits], s=64, color=SERIES,
                   edgecolors=SURFACE, linewidths=2, zorder=4, label="Detections")
        k = max(hits, key=power.__getitem__)
        ax.annotate(f"{len(hits)} bins over threshold, peak {freq[k]:.2f} Hz",
                    (freq[k], power_db[k]), xytext=(10, -4), textcoords="offset points",
                    color=INK, fontsize=9)

    ax.set_xlim(0, min(freq[-1], band[1] * 1.2))
    ax.set_xlabel("Frequency (Hz)", color=INK_2)
    ax.set_ylabel("P[k] (dB re sigma^2)", color=INK_2)
    ax.set_title("Hann periodogram with CA-CFAR threshold, one window", color=INK, loc="left")
    ax.grid(True, color=GRID, lw=0.6)
    ax.tick_params(colors=MUTED)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    legend = ax.legend(loc="lower right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(INK_2)

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
