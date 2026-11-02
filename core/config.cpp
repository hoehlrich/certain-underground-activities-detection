#include "core/config.hpp"

#include <stdexcept>

namespace cfar {

void validate(const ProcessingConfig& cfg) {
    if (cfg.fs_hz < 2.5 * cfg.f_hi_hz) {
        throw std::runtime_error("config: fs_hz must be at least 2.5 x f_hi_hz");
    }
    if (!(0.0 < cfg.f_lo_hz && cfg.f_lo_hz < cfg.f_hi_hz)) {
        throw std::runtime_error("config: need 0 < f_lo_hz < f_hi_hz");
    }
    if (cfg.nfft < 2 || cfg.nfft % 2 != 0) {
        throw std::runtime_error("config: nfft must be even and at least 2");
    }
    if (cfg.window != "hann") {
        throw std::runtime_error("config: only window = \"hann\" is supported");
    }
    if (!(0.0 < cfg.pfa_window && cfg.pfa_window < 1.0)) {
        throw std::runtime_error("config: pfa_window must be in (0, 1)");
    }
    if (cfg.pfa_bin && !(0.0 < *cfg.pfa_bin && *cfg.pfa_bin < 1.0)) {
        throw std::runtime_error("config: pfa_bin must be in (0, 1)");
    }
    if (cfg.train_per_side == 0) {
        throw std::runtime_error("config: train_per_side must be at least 1");
    }
}

BinRange band_bins(const ProcessingConfig& cfg) {
    return band_bins(cfg.f_lo_hz, cfg.f_hi_hz, cfg.nfft, cfg.fs_hz);
}

double design_pfa_bin(const ProcessingConfig& cfg) {
    return cfg.pfa_bin.value_or(pfa_bin_for_window(cfg.pfa_window, band_bins(cfg).size()));
}

CaParams ca_params(const ProcessingConfig& cfg) {
    return {cfg.guard_per_side, cfg.train_per_side, design_pfa_bin(cfg)};
}

}  // namespace cfar
