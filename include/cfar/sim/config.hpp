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
// [pfa_check]. Defaults match docs/cfar-sim-plan.md; none are confirmed by the sponsor.
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
};

// Throws std::runtime_error on a parse error, a bad value or a SUB-03.2 violation.
SimConfig load_config(const std::filesystem::path& path);

// Checks the values in cfg, including cfg.proc; load_config calls this.
void validate(const SimConfig& cfg);

}  // namespace cfar::sim
