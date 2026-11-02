#pragma once

#include <cstddef>
#include <optional>
#include <string>

#include "core/cfar.hpp"
#include "core/spectrum.hpp"

namespace cfar {

// Processing parameters shared by simulation and real data: the [band], [acq],
// [false_alarm] and [cfar] sections of a config (loaded by io/config_loader).
// None of the defaults are confirmed by the sponsor.
struct ProcessingConfig {
    // [band] search band, TBD with the sponsor
    double f_lo_hz = 10.0;
    double f_hi_hz = 1000.0;
    // [acq]
    double fs_hz = 4000.0;  // at least 2.5 x f_hi_hz, leaving room for anti-alias roll-off
    std::size_t nfft = 8192;
    std::string window = "hann";
    // [false_alarm]
    double pfa_window = 0.01;     // draft target, not agreed with the sponsor
    std::optional<double> pfa_bin;  // empty means "derived" from pfa_window and the band
    // [cfar]
    std::size_t guard_per_side = 2;
    std::size_t train_per_side = 16;
};

// Checks the values in cfg, including fs_hz >= 2.5 x f_hi_hz. Throws std::runtime_error.
void validate(const ProcessingConfig& cfg);

BinRange band_bins(const ProcessingConfig& cfg);

// pfa_bin if set, otherwise derived so the band meets pfa_window.
double design_pfa_bin(const ProcessingConfig& cfg);

CaParams ca_params(const ProcessingConfig& cfg);

}  // namespace cfar
