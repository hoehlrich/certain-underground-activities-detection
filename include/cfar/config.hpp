#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>

#include "cfar/cfar.hpp"
#include "cfar/spectrum.hpp"

namespace cfar {

// Processing parameters shared by simulation and real data: the [band], [acq],
// [false_alarm] and [cfar] sections of a config. Defaults match
// docs/cfar-sim-plan.md; none are confirmed by the sponsor.
struct ProcessingConfig {
    // [band] SYS-01 search band, TBD
    double f_lo_hz = 10.0;
    double f_hi_hz = 1000.0;
    // [acq]
    double fs_hz = 4000.0;  // SUB-03.2: at least 2.5 x f_hi_hz
    std::size_t nfft = 8192;
    std::string window = "hann";
    // [false_alarm]
    double pfa_window = 0.01;     // SUB-05.1 draft
    std::optional<double> pfa_bin;  // empty means "derived" from pfa_window and the band
    // [cfar]
    std::size_t guard_per_side = 2;
    std::size_t train_per_side = 16;
};

// Reads the processing sections of a TOML config; other sections are ignored.
// Throws std::runtime_error on a parse error or a bad value.
ProcessingConfig load_processing_config(const std::filesystem::path& path);

// Checks the values in cfg, including SUB-03.2; load_processing_config calls this.
void validate(const ProcessingConfig& cfg);

BinRange band_bins(const ProcessingConfig& cfg);

// pfa_bin if set, otherwise derived so the band meets pfa_window.
double design_pfa_bin(const ProcessingConfig& cfg);

CaParams ca_params(const ProcessingConfig& cfg);

}  // namespace cfar
