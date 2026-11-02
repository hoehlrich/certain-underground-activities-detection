"""Plot results CSVs.

Usage: python3 scripts/plot.py results/example_spectrum.csv [--band 10 1000] [-o out.png]
       python3 scripts/plot.py results/pd_vs_snr.csv [--pd-target 0.9] [-o out.png]
       python3 scripts/plot.py results/roc.csv [-o out.png]

The CSV's columns pick the plot: one window's P[k] with the CA-CFAR threshold
and detections, Pd vs C/N0 with the known-noise bound, or ROC curves at fixed SNRs.
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
SLOT_3 = "#1baf7a"     # categorical slot 3
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
    style_axes(ax)
    legend = ax.legend(loc="lower right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(INK_2)

    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    print(f"wrote {out}")


def crossing(x, y, target):
    for i in range(1, len(y)):
        if y[i - 1] < target <= y[i]:
            return x[i - 1] + (target - y[i - 1]) / (y[i] - y[i - 1]) * (x[i] - x[i - 1])
    return float("nan")


def style_axes(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.tick_params(colors=MUTED)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)


def plot_pd(path, pd_target, out):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    cn0 = [float(r["cn0_dbhz"]) for r in rows]
    pd = [float(r["pd"]) for r in rows]
    lo = [float(r["ci95_lo"]) for r in rows]
    hi = [float(r["ci95_hi"]) for r in rows]
    ideal = [float(r["pd_ideal"]) for r in rows]
    offset = float(rows[0]["snr_bin_db"]) - cn0[0]
    trials = rows[0]["trials"]

    fig, ax = plt.subplots(figsize=(9, 4.5), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axhline(pd_target, color=MUTED, lw=1.0, zorder=1)
    ax.text(cn0[0], pd_target, f" Pd = {pd_target:g}", color=INK_2, fontsize=9, va="bottom")

    ax.plot(cn0, ideal, color=THRESHOLD, lw=2.0, zorder=2, label="Known noise (ideal bound)")
    ax.fill_between(cn0, lo, hi, color=SERIES, alpha=0.18, lw=0, zorder=2)
    ax.plot(cn0, pd, color=SERIES, lw=2.0, zorder=3, label=f"CA-CFAR, {trials} trials per point")
    ax.scatter(cn0, pd, s=64, color=SERIES, edgecolors=SURFACE, linewidths=2, zorder=4)

    # Exact crossings from pd_sweep's summary when present; else interpolate the grid.
    x_ca = crossing(cn0, pd, pd_target)
    x_id = crossing(cn0, ideal, pd_target)
    summary = path.with_name(path.stem + "_summary.csv")
    if summary.exists():
        with open(summary, newline="") as f:
            row = next(csv.DictReader(f))
        if float(row["pd_target"]) == pd_target:
            x_ca = float(row["cn0_ca_dbhz"])
            x_id = float(row["cn0_ideal_dbhz"])
    for x, color, text, dy, ha, va in ((x_id, THRESHOLD, "ideal", 6, "right", "bottom"),
                                       (x_ca, SERIES, "CA-CFAR", -6, "left", "top")):
        if math.isnan(x):
            continue
        ax.scatter([x], [pd_target], s=64, color=color, edgecolors=SURFACE, linewidths=2, zorder=5)
        ax.annotate(f"{text} {x:.1f} dB-Hz", (x, pd_target), xytext=(-dy, dy),
                    textcoords="offset points", ha=ha, va=va, color=INK, fontsize=9)
    if not (math.isnan(x_ca) or math.isnan(x_id)):
        ax.text(0.99, 0.04, f"CFAR loss at Pd = {pd_target:g}: {x_ca - x_id:.1f} dB",
                transform=ax.transAxes, ha="right", color=INK_2, fontsize=9)

    ax.set_xlim(cn0[0], cn0[-1])
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(f"C/N0 (dB-Hz)    per-bin SNR = C/N0 + {offset:.2f} dB", color=INK_2)
    ax.set_ylabel("Probability of detection", color=INK_2)
    ax.set_title("Pd vs C/N0, tone uniform over the band, white noise", color=INK, loc="left")
    style_axes(ax)
    legend = ax.legend(loc="center right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(INK_2)

    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    print(f"wrote {out}")


def plot_roc(path, out):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    curves = {}
    for r in rows:
        curves.setdefault(float(r["snr_bin_db"]), []).append(r)
    if len(curves) > 3:
        raise SystemExit("plot_roc: at most 3 SNRs per plot")
    trials = rows[0]["trials"]

    fig, ax = plt.subplots(figsize=(9, 4.5), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axvline(0.01, color=MUTED, lw=1.0, zorder=1)
    ax.text(0.01, 0.02, " 1% per window\n (draft target)", color=INK_2, fontsize=9,
            transform=ax.get_xaxis_transform(), va="bottom")

    for color, (snr, rs) in zip((SERIES, THRESHOLD, SLOT_3), sorted(curves.items())):
        x = [float(r["pfa_window"]) for r in rs]
        pd = [float(r["pd"]) for r in rs]
        ideal = [float(r["pd_ideal"]) for r in rs]
        cn0 = float(rs[0]["cn0_dbhz"])
        ax.fill_between(x, [float(r["ci95_lo"]) for r in rs], [float(r["ci95_hi"]) for r in rs],
                        color=color, alpha=0.18, lw=0, zorder=2)
        ax.plot(x, ideal, color=color, lw=1.5, ls=(0, (4, 3)), zorder=2)
        ax.plot(x, pd, color=color, lw=2.0, zorder=3,
                label=f"C/N0 {cn0:.1f} dB-Hz (SNR_bin {snr:g} dB)")
        ax.scatter(x, pd, s=64, color=color, edgecolors=SURFACE, linewidths=2, zorder=4)

    ax.plot([], [], color=INK_2, lw=2.0, label=f"CA-CFAR, {trials} trials per point")
    ax.plot([], [], color=INK_2, lw=1.5, ls=(0, (4, 3)), label="Known noise (ideal bound)")
    ax.set_xscale("log")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("False-alarm probability per window (design)", color=INK_2)
    ax.set_ylabel("Probability of detection", color=INK_2)
    ax.set_title("ROC at fixed SNR, tone uniform over the band, white noise", color=INK,
                 loc="left")
    style_axes(ax)
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
    ap.add_argument("--pd-target", type=float, default=0.9)
    ap.add_argument("-o", "--out", type=Path)
    args = ap.parse_args()
    out = args.out or args.csv.with_suffix(".png")
    with open(args.csv, newline="") as f:
        header = next(csv.reader(f))
    if "pfa_window" in header:
        plot_roc(args.csv, out)
    elif "pd" in header:
        plot_pd(args.csv, args.pd_target, out)
    else:
        plot_spectrum(args.csv, args.band, out)


if __name__ == "__main__":
    main()
