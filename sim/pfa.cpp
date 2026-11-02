#include "sim/pfa.hpp"

#include <algorithm>
#include <cmath>

#include "core/cfar.hpp"
#include "sim/noise.hpp"
#include "core/spectrum.hpp"

namespace cfar::sim {

Interval binomial_interval_95(double p, std::uint64_t n) {
    const double half = 1.96 * std::sqrt(p * (1.0 - p) / static_cast<double>(n));
    return {p - half, p + half};
}

std::vector<PfaResult> run_pfa_check(const SimConfig& cfg, std::span<const double> design_pfas,
                                     std::size_t windows) {
    const auto& proc = cfg.proc;
    const auto window = make_window(proc.window, proc.nfft);
    const BinRange band = band_bins(proc);

    std::vector<CaCfar> detectors;
    std::vector<PfaResult> results;
    for (double pfa : design_pfas) {
        detectors.emplace_back(CaParams{proc.guard_per_side, proc.train_per_side, pfa}, window);
        results.push_back({pfa, 0, 0, 0, windows});
    }

    std::vector<double> x(proc.nfft);
    for (std::size_t w = 0; w < windows; ++w) {
        auto rng = make_trial_rng(cfg.seed, w);
        std::fill(x.begin(), x.end(), 0.0);
        add_white_noise(x, cfg.noise_sigma, rng);
        const auto power = periodogram(x, window);
        const auto noise = ca_noise(power, band, proc.guard_per_side, proc.train_per_side);

        for (std::size_t d = 0; d < detectors.size(); ++d) {
            std::uint64_t hits = 0;
            for (std::size_t k = band.first; k <= band.last; ++k) {
                hits += power[k] > detectors[d].alpha(noise.sides[k]) * noise.level[k] ? 1 : 0;
            }
            results[d].hits += hits;
            results[d].window_hits += hits > 0 ? 1 : 0;
            results[d].cells += band.size();
        }
    }
    return results;
}

}  // namespace cfar::sim
