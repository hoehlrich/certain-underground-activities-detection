#include "cfar/detector.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfar {

double parabolic_offset(double p_minus, double p_peak, double p_plus) {
    if (!(p_minus > 0.0 && p_peak > 0.0 && p_plus > 0.0)) {
        return 0.0;
    }
    const double lm = std::log(p_minus);
    const double l0 = std::log(p_peak);
    const double lp = std::log(p_plus);
    const double curvature = lm - 2.0 * l0 + lp;
    if (!(curvature < 0.0)) {
        return 0.0;  // not a local maximum
    }
    return std::clamp(0.5 * (lm - lp) / curvature, -0.5, 0.5);
}

std::vector<Detection> find_detections(std::span<const double> power,
                                       std::span<const double> threshold,
                                       std::span<const double> noise, BinRange band,
                                       double bin_hz) {
    if (threshold.size() != power.size() || noise.size() != power.size() ||
        band.last >= power.size()) {
        throw std::invalid_argument("find_detections: spectrum sizes do not match the band");
    }
    std::vector<Detection> out;
    std::size_t k = band.first;
    while (k <= band.last) {
        if (!(power[k] > threshold[k])) {
            ++k;
            continue;
        }
        Detection d{};
        d.first_bin = k;
        d.peak_bin = k;
        while (k <= band.last && power[k] > threshold[k]) {
            if (power[k] > power[d.peak_bin]) {
                d.peak_bin = k;
            }
            d.last_bin = k++;
        }
        const std::size_t p = d.peak_bin;
        const double offset = p > 0 && p + 1 < power.size()
                                  ? parabolic_offset(power[p - 1], power[p], power[p + 1])
                                  : 0.0;
        d.freq_hz = (static_cast<double>(p) + offset) * bin_hz;
        d.statistic = power[p] / noise[p];
        d.snr_est = d.statistic - 1.0;
        out.push_back(d);
    }
    return out;
}

Detector::Detector(const ProcessingConfig& cfg)
    : window_(make_window(cfg.window, cfg.nfft)),
      band_(band_bins(cfg)),
      bin_hz_(cfar::bin_hz(cfg.nfft, cfg.fs_hz)),
      ca_(ca_params(cfg), window_) {}

std::vector<Detection> Detector::process(std::span<const double> samples) const {
    if (samples.size() != window_.size()) {
        throw std::invalid_argument("Detector::process: expected nfft samples");
    }
    return detect(periodogram(samples, window_));
}

std::vector<Detection> Detector::detect(std::span<const double> power) const {
    const auto& p = ca_.params();
    const NoiseEstimate noise = ca_noise(power, band_, p.guard_per_side, p.train_per_side);
    const auto threshold = ca_.threshold(noise);
    return find_detections(power, threshold, noise.level, band_, bin_hz_);
}

}  // namespace cfar
