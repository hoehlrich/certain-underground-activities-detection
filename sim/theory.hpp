#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace cfar::sim {

// Marcum Q-function of order 1, Q1(a, b) = P(|a + n| > b) for n complex
// Gaussian with unit variance per component.
double marcum_q1(double a, double b);

// Tone power gain of a window for a tone offset_bins from a bin centre,
// relative to a bin-centred tone: |W(offset)|^2 / |W(0)|^2. Hann at 0.5 bins
// gives 1.42 dB of loss.
double scallop_gain(std::span<const double> window, double offset_bins);

// Pd of the known-noise detector on one bin: threshold -ln(pfa) sigma^2,
// per-bin SNR snr_bin (linear).
double ideal_pd(double snr_bin, double pfa);

// Known-noise bound for a tone tested in its nearest bin, averaged over tone
// offsets from the bin centre. snr_bin_db is the bin-centred per-bin SNR
// (C/N0 + 10 log10(T) - 10 log10(ENBW)); scalloping is applied per offset.
class IdealBound {
public:
    // offsets_bins: tone offsets to average over (one value for a fixed f0).
    IdealBound(std::span<const double> window, double pfa_bin,
               std::span<const double> offsets_bins);

    double pd(double snr_bin_db) const;
    // Smallest bin-centred SNR [dB] with pd >= target, by bisection.
    double snr_bin_db_for(double pd_target) const;

private:
    double pfa_;
    std::vector<double> gains_;
};

// n offsets at the midpoints of [0, 0.5) bins: a tone uniform in frequency,
// for a window whose response is symmetric in offset.
std::vector<double> uniform_offsets(std::size_t n = 64);

}  // namespace cfar::sim
