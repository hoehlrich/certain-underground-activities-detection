#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "sim/config.hpp"

namespace cfar::sim {

struct Interval {
    double lo;
    double hi;
    bool contains(double x) const { return lo <= x && x <= hi; }
};

// 95% interval for the fraction of hits in n Bernoulli(p) trials (normal
// approximation; fine once n p is a few tens).
Interval binomial_interval_95(double p, std::uint64_t n);

struct PfaResult {
    double design_pfa;
    std::uint64_t hits;         // bins over threshold
    std::uint64_t window_hits;  // windows with at least one bin over threshold
    std::uint64_t cells;        // bins tested
    std::uint64_t windows;
    double measured() const { return static_cast<double>(hits) / static_cast<double>(cells); }
    Interval interval() const { return binomial_interval_95(design_pfa, cells); }
    bool pass() const { return interval().contains(measured()); }
};

// Noise-only Monte Carlo of CA-CFAR over the search band in cfg.proc: one CA
// threshold per design Pfa, all applied to the same windows. Window w uses
// make_trial_rng(cfg.seed, w).
std::vector<PfaResult> run_pfa_check(const SimConfig& cfg, std::span<const double> design_pfas,
                                     std::size_t windows);

}  // namespace cfar::sim
