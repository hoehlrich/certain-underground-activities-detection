#include "io/config_loader.hpp"

#include <stdexcept>
#include <string>

#include <toml++/toml.hpp>

namespace cfar::io {

ProcessingConfig load_processing_config(const std::filesystem::path& path) {
    toml::table tbl;
    try {
        tbl = toml::parse_file(path.string());
    } catch (const toml::parse_error& e) {
        throw std::runtime_error("config " + path.string() + ": " + std::string(e.description()));
    }

    ProcessingConfig cfg;
    cfg.f_lo_hz = tbl["band"]["f_lo_hz"].value_or(cfg.f_lo_hz);
    cfg.f_hi_hz = tbl["band"]["f_hi_hz"].value_or(cfg.f_hi_hz);
    cfg.fs_hz = tbl["acq"]["fs_hz"].value_or(cfg.fs_hz);
    cfg.nfft = tbl["acq"]["nfft"].value_or(cfg.nfft);
    cfg.window = tbl["acq"]["window"].value_or(cfg.window);
    cfg.pfa_window = tbl["false_alarm"]["pfa_window"].value_or(cfg.pfa_window);
    cfg.guard_per_side = tbl["cfar"]["guard_per_side"].value_or(cfg.guard_per_side);
    cfg.train_per_side = tbl["cfar"]["train_per_side"].value_or(cfg.train_per_side);

    // pfa_bin is either "derived" or a fixed value.
    const auto pfa_bin = tbl["false_alarm"]["pfa_bin"];
    if (auto p = pfa_bin.value<double>()) {
        cfg.pfa_bin = *p;
    } else if (pfa_bin && pfa_bin.value_or(std::string{}) != "derived") {
        throw std::runtime_error("config: false_alarm.pfa_bin must be \"derived\" or a number");
    }

    validate(cfg);
    return cfg;
}

}  // namespace cfar::io
