#include "sim/pd.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>

#include "core/detector.hpp"
#include "sim/noise.hpp"
#include "sim/signal.hpp"
#include "core/spectrum.hpp"

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

IdealBound ideal_bound(const SimConfig& cfg, double pfa_bin) {
    const auto& proc = cfg.proc;
    std::vector<double> offsets = uniform_offsets();
    if (cfg.f0_hz) {
        const double bins = *cfg.f0_hz / bin_hz(proc.nfft, proc.fs_hz);
        offsets = {bins - std::round(bins)};
    }
    return {make_window(proc.window, proc.nfft), pfa_bin, offsets};
}

IdealBound ideal_bound(const SimConfig& cfg) { return ideal_bound(cfg, design_pfa_bin(cfg.proc)); }

std::vector<double> pd_snr_grid(const SimConfig& cfg) {
    std::vector<double> grid;
    const double eps = 1e-9 * cfg.pd_snr_step_db;
    for (double s = cfg.pd_snr_min_db; s <= cfg.pd_snr_max_db + eps; s += cfg.pd_snr_step_db) {
        grid.push_back(s);
    }
    return grid;
}

std::vector<PdPoint> run_pd_grid(const SimConfig& cfg, std::span<const double> snr_bin_db,
                                 std::span<const double> pfa_bins, std::size_t trials) {
    const auto& proc = cfg.proc;
    const double df = bin_hz(proc.nfft, proc.fs_hz);
    const double offset_db = snr_bin_offset_db(proc);
    const std::size_t n_bins = band_bins(proc).size();

    // One detector and one bound per false-alarm setting.
    std::vector<Detector> detectors;
    std::vector<IdealBound> bounds;
    for (double pfa : pfa_bins) {
        ProcessingConfig p = proc;
        p.pfa_bin = pfa;
        detectors.emplace_back(p);
        bounds.push_back(ideal_bound(cfg, pfa));
    }
    const auto& window = detectors.front().window();

    std::vector<PdPoint> points;
    std::vector<double> amplitudes;
    for (double s : snr_bin_db) {
        const double cn0 = s - offset_db;
        amplitudes.push_back(amplitude_for_cn0(cn0, cfg.noise_sigma, proc.fs_hz));
        for (std::size_t j = 0; j < pfa_bins.size(); ++j) {
            const double pfa_window =
                -std::expm1(static_cast<double>(n_bins) * std::log1p(-pfa_bins[j]));
            points.push_back({s, cn0, pfa_bins[j], pfa_window, trials, 0, 0, bounds[j].pd(s)});
        }
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

        for (std::size_t i = 0; i < amplitudes.size(); ++i) {
            for (std::size_t k = 0; k < x.size(); ++k) {
                x[k] = noise[k] + amplitudes[i] * unit_tone[k];
            }
            const auto power = periodogram(x, window);
            for (std::size_t j = 0; j < detectors.size(); ++j) {
                bool hit = false;
                bool fa = false;
                for (const auto& d : detectors[j].detect(power)) {
                    (std::abs(d.freq_hz - tone.freq_hz) <= df ? hit : fa) = true;
                }
                auto& pt = points[i * detectors.size() + j];
                pt.hits += hit ? 1 : 0;
                pt.fa_windows += fa ? 1 : 0;
            }
        }
    }
    return points;
}

std::vector<PdPoint> run_pd_sweep(const SimConfig& cfg, std::span<const double> snr_bin_db,
                                  std::size_t trials) {
    const double pfa = design_pfa_bin(cfg.proc);
    return run_pd_grid(cfg, snr_bin_db, std::span(&pfa, 1), trials);
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
