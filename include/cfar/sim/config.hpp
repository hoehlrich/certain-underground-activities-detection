#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>

namespace cfar::sim {

// Simulation parameters, loaded from TOML (configs/baseline.toml). Defaults
// match docs/cfar-sim-plan.md; none are confirmed by the sponsor.
struct SimConfig {
    // [band] SYS-01 search band, TBD
    double f_lo_hz = 10.0;
    double f_hi_hz = 1000.0;
    // [acq]
    double fs_hz = 4000.0;  // SUB-03.2: at least 2.5 x f_hi_hz
    std::size_t nfft = 8192;
    std::string window = "hann";
    // [noise]
    double noise_sigma = 1.0;
    // [source]
    std::string source_kind = "tone";
    std::optional<double> f0_hz;  // empty means "uniform": drawn over the band each trial
    double cn0_dbhz = 30.0;
    // [run]
    std::uint64_t seed = 1;
};

// Throws std::runtime_error on a parse error, a bad value or a SUB-03.2 violation.
SimConfig load_config(const std::filesystem::path& path);

// Checks the values in cfg; load_config calls this.
void validate(const SimConfig& cfg);

}  // namespace cfar::sim
