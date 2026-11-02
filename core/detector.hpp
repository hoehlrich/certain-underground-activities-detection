#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "core/cfar.hpp"
#include "core/config.hpp"
#include "core/spectrum.hpp"

namespace cfar {

// One detection in a window: a run of adjacent bins over threshold, reported at
// its peak bin with its detection statistic and estimated SNR.
struct Detection {
    std::size_t first_bin;  // cluster extent
    std::size_t last_bin;
    std::size_t peak_bin;   // bin with the largest P[k] in the cluster
    double freq_hz;         // peak refined by parabolic interpolation on log power
    double statistic;       // P[peak] / noise[peak]
    double snr_est;         // (P[peak] - noise[peak]) / noise[peak], linear
};

// Groups adjacent in-band bins with power > threshold into detections. noise is
// the per-bin noise level the threshold was built from. NaN thresholds (outside
// the band) never detect.
std::vector<Detection> find_detections(std::span<const double> power,
                                       std::span<const double> threshold,
                                       std::span<const double> noise, BinRange band,
                                       double bin_hz);

// Peak offset in bins, in [-0.5, 0.5], from a parabola through the log power
// of bins k-1, k, k+1.
double parabolic_offset(double p_minus, double p_peak, double p_plus);

// The per-window pipeline: window, periodogram, CA-CFAR over the band, clustering.
// Built from processing settings only, so it runs the same on simulated and
// recorded data.
class Detector {
public:
    explicit Detector(const ProcessingConfig& cfg);

    // samples.size() must equal cfg.nfft.
    std::vector<Detection> process(std::span<const double> samples) const;
    // power is a periodogram of nfft / 2 + 1 bins made with window().
    std::vector<Detection> detect(std::span<const double> power) const;

    const std::vector<double>& window() const { return window_; }
    BinRange band() const { return band_; }
    double bin_hz() const { return bin_hz_; }
    const CaCfar& cfar() const { return ca_; }

private:
    std::vector<double> window_;
    BinRange band_;
    double bin_hz_;
    CaCfar ca_;
};

}  // namespace cfar
