#include "cfar/sim/pd.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>

#include "cfar/detector.hpp"
#include "cfar/sim/noise.hpp"
#include "cfar/sim/signal.hpp"
#include "cfar/spectrum.hpp"

namespace cfar::sim {

Interval wilson_interval_95(std::uint64_t hits, std::uint64_t n) {
    constexpr double z = 1.96;
    const double nn = static_cast<double>(n);
    const double p = static_cast<double>(hits) / nn;
    const double denom = 1.0 + z * z / nn;
    const double centre = (p + z * z / (2.0 * nn)) / denom;
    const double half = z * std::sqrt(p * (1.0 - p) / nn + z * z / (4.0 * nn * nn)) / denom;
    return {centre - half, centre + half};
}

double snr_bin_offset_db(const ProcessingConfig& proc) {
    const double t = static_cast<double>(proc.nfft) / proc.fs_hz;
    const double enbw = enbw_bins(make_window(proc.window, proc.nfft));
    return 10.0 * std::log10(t) - 10.0 * std::log10(enbw);
}

IdealBound ideal_bound(const SimConfig& cfg) {
    const auto& proc = cfg.proc;
    std::vector<double> offsets = uniform_offsets();
    if (cfg.f0_hz) {
        const double bins = *cfg.f0_hz / bin_hz(proc.nfft, proc.fs_hz);
        offsets = {bins - std::round(bins)};
    }
    return {make_window(proc.window, proc.nfft), design_pfa_bin(proc), offsets};
}

std::vector<double> pd_snr_grid(const SimConfig& cfg) {
    std::vector<double> grid;
    const double eps = 1e-9 * cfg.pd_snr_step_db;
    for (double s = cfg.pd_snr_min_db; s <= cfg.pd_snr_max_db + eps; s += cfg.pd_snr_step_db) {
        grid.push_back(s);
    }
    return grid;
}

std::vector<PdPoint> run_pd_sweep(const SimConfig& cfg, std::span<const double> snr_bin_db,
                                  std::size_t trials) {
    const auto& proc = cfg.proc;
    const Detector detector(proc);
    const double df = detector.bin_hz();
    const double offset_db = snr_bin_offset_db(proc);

    const IdealBound ideal = ideal_bound(cfg);

    std::vector<PdPoint> points;
    std::vector<double> amplitudes;
    for (double s : snr_bin_db) {
        const double cn0 = s - offset_db;
        points.push_back({s, cn0, trials, 0, 0, ideal.pd(s)});
        amplitudes.push_back(amplitude_for_cn0(cn0, cfg.noise_sigma, proc.fs_hz));
    }

    std::vector<double> noise(proc.nfft);
    std::vector<double> unit_tone(proc.nfft);
    std::vector<double> x(proc.nfft);
    for (std::size_t t = 0; t < trials; ++t) {
        auto rng = make_trial_rng(cfg.seed, t);
        ToneParams tone;
        tone.freq_hz = cfg.f0_hz.value_or(
            std::uniform_real_distribution<double>(proc.f_lo_hz, proc.f_hi_hz)(rng));
        tone.phase_rad = std::uniform_real_distribution<double>(0.0, 2.0 * std::numbers::pi)(rng);
        std::fill(noise.begin(), noise.end(), 0.0);
        add_white_noise(noise, cfg.noise_sigma, rng);
        std::fill(unit_tone.begin(), unit_tone.end(), 0.0);
        add_tone(unit_tone, proc.fs_hz, tone);

        for (std::size_t i = 0; i < points.size(); ++i) {
            for (std::size_t k = 0; k < x.size(); ++k) {
                x[k] = noise[k] + amplitudes[i] * unit_tone[k];
            }
            bool hit = false;
            bool fa = false;
            for (const auto& d : detector.process(x)) {
                (std::abs(d.freq_hz - tone.freq_hz) <= df ? hit : fa) = true;
            }
            points[i].hits += hit ? 1 : 0;
            points[i].fa_windows += fa ? 1 : 0;
        }
    }
    return points;
}

double crossing_snr_bin_db(std::span<const PdPoint> points, double target) {
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (points[i].pd() < target) {
            continue;
        }
        if (i == 0) {
            return points[0].snr_bin_db;
        }
        const auto& a = points[i - 1];
        const auto& b = points[i];
        const double frac = (target - a.pd()) / (b.pd() - a.pd());
        return a.snr_bin_db + frac * (b.snr_bin_db - a.snr_bin_db);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace cfar::sim
