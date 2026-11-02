#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "core/config.hpp"

namespace cfar::sim {

// Simulation config: the processing sections plus [noise], [source], [run] and
// [pfa_check], [pd_sweep] and [roc]. None of the defaults are confirmed by the sponsor.
struct SimConfig {
    ProcessingConfig proc;
    // [noise]
    double noise_sigma = 1.0;
    // [source]
    std::string source_kind = "tone";
    std::optional<double> f0_hz;  // empty means "uniform": drawn over the band each trial
    double cn0_dbhz = 30.0;
    // [run]
    std::uint64_t seed = 1;
    // [pfa_check] design values checked besides the derived baseline
    std::size_t pfa_check_windows = 10000;
    std::vector<double> pfa_check_design = {1e-3, 1e-4, 1e-5};
    // [pd_sweep] grid of bin-centred per-bin SNR, trials per point, Pd to report
    double pd_snr_min_db = 0.0;
    double pd_snr_max_db = 25.0;
    double pd_snr_step_db = 1.0;
    std::size_t pd_trials = 2000;
    double pd_target = 0.9;
    // [roc] fixed SNRs; false-alarm probabilities per window (per-bin values are derived)
    std::vector<double> roc_snr_bin_db = {12.0, 14.0, 16.0};
    std::vector<double> roc_pfa_window = {1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 0.01,
                                          0.02, 0.05, 0.1,  0.2,  0.5};
    std::size_t roc_trials = 2000;
};

// Throws std::runtime_error on a parse error, a bad value or a SUB-03.2 violation.
SimConfig load_config(const std::filesystem::path& path);

// Checks the values in cfg, including cfg.proc; load_config calls this.
void validate(const SimConfig& cfg);

}  // namespace cfar::sim
