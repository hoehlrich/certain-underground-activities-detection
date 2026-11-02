#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cfar/config.hpp"
#include "cfar/sim/config.hpp"
#include "cfar/sim/pfa.hpp"
#include "cfar/sim/theory.hpp"

namespace cfar::sim {

// Wilson 95% interval for hits out of n trials; holds up near 0 and 1, where
// the normal approximation does not.
Interval wilson_interval_95(std::uint64_t hits, std::uint64_t n);

// 10 log10(T) - 10 log10(ENBW): add to C/N0 [dB-Hz] to get the bin-centred
// per-bin SNR [dB]. 1.35 dB at the baseline (T = 2.05 s, Hann).
double snr_bin_offset_db(const ProcessingConfig& proc);

struct PdPoint {
    double snr_bin_db;  // bin-centred per-bin SNR
    double cn0_dbhz;
    double pfa_bin;     // design false-alarm probability per bin
    double pfa_window;  // 1 - (1 - pfa_bin)^bins_in_band
    std::uint64_t trials;
    std::uint64_t hits;        // windows with a detection within one bin of f0
    std::uint64_t fa_windows;  // windows with a detection more than one bin from f0
    double ideal_pd;           // known-noise bound, same f0 distribution
    double pd() const { return static_cast<double>(hits) / static_cast<double>(trials); }
    Interval interval() const { return wilson_interval_95(hits, trials); }
};

// Known-noise bound for cfg at pfa_bin (default: the design Pfa per bin), with
// tone offsets uniform over a bin, or the single offset of a fixed f0.
IdealBound ideal_bound(const SimConfig& cfg, double pfa_bin);
IdealBound ideal_bound(const SimConfig& cfg);

// Bin-centred SNR grid from cfg's [pd_sweep] section.
std::vector<double> pd_snr_grid(const SimConfig& cfg);

// Tone in white noise through the full detector (cfar::Detector) at every
// (SNR, per-bin Pfa) pair; points are ordered SNR-major. Trial t uses
// make_trial_rng(cfg.seed, t) for f0, phase and noise at every point, so points
// differ only in tone amplitude and threshold, and curves are smooth.
std::vector<PdPoint> run_pd_grid(const SimConfig& cfg, std::span<const double> snr_bin_db,
                                 std::span<const double> pfa_bins, std::size_t trials);

// Pd vs SNR at the design Pfa per bin.
std::vector<PdPoint> run_pd_sweep(const SimConfig& cfg, std::span<const double> snr_bin_db,
                                  std::size_t trials);

// SNR [dB] where Pd first reaches target, by linear interpolation between grid
// points. NaN if the sweep never reaches it.
double crossing_snr_bin_db(std::span<const PdPoint> points, double target);

}  // namespace cfar::sim
