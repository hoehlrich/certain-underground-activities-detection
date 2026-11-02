#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "cfar/config.hpp"

namespace cfar::sim {

// Simulation config: the processing sections plus [noise], [source], [run] and
// [pfa_check] and [pd_sweep]. Defaults match docs/cfar-sim-plan.md; none are confirmed by the sponsor.
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
};

// Throws std::runtime_error on a parse error, a bad value or a SUB-03.2 violation.
SimConfig load_config(const std::filesystem::path& path);

// Checks the values in cfg, including cfg.proc; load_config calls this.
void validate(const SimConfig& cfg);

}  // namespace cfar::sim
