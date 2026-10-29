#include "cfar/sim/config.hpp"

#include <stdexcept>

#include <toml++/toml.hpp>

namespace cfar::sim {

SimConfig load_config(const std::filesystem::path& path) {
    toml::table tbl;
    try {
        tbl = toml::parse_file(path.string());
    } catch (const toml::parse_error& e) {
        throw std::runtime_error("config " + path.string() + ": " + std::string(e.description()));
    }

    SimConfig cfg;
    cfg.f_lo_hz = tbl["band"]["f_lo_hz"].value_or(cfg.f_lo_hz);
    cfg.f_hi_hz = tbl["band"]["f_hi_hz"].value_or(cfg.f_hi_hz);
    cfg.fs_hz = tbl["acq"]["fs_hz"].value_or(cfg.fs_hz);
    cfg.nfft = tbl["acq"]["nfft"].value_or(cfg.nfft);
    cfg.window = tbl["acq"]["window"].value_or(cfg.window);
    cfg.noise_sigma = tbl["noise"]["sigma"].value_or(cfg.noise_sigma);
    cfg.source_kind = tbl["source"]["kind"].value_or(cfg.source_kind);
    cfg.cn0_dbhz = tbl["source"]["cn0_dbhz"].value_or(cfg.cn0_dbhz);
    cfg.seed = tbl["run"]["seed"].value_or(cfg.seed);

    // f0_hz is either "uniform" or a fixed frequency.
    const auto f0 = tbl["source"]["f0_hz"];
    if (auto hz = f0.value<double>()) {
        cfg.f0_hz = *hz;
    } else if (f0 && f0.value_or(std::string{}) != "uniform") {
        throw std::runtime_error("config: source.f0_hz must be \"uniform\" or a number");
    }

    validate(cfg);
    return cfg;
}

void validate(const SimConfig& cfg) {
    if (cfg.fs_hz < 2.5 * cfg.f_hi_hz) {
        throw std::runtime_error("config: fs_hz must be at least 2.5 x f_hi_hz (SUB-03.2)");
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
    if (!(cfg.noise_sigma > 0.0)) {
        throw std::runtime_error("config: noise sigma must be positive");
    }
    if (cfg.source_kind != "tone") {
        throw std::runtime_error("config: only source.kind = \"tone\" is supported");
    }
    if (cfg.f0_hz && !(cfg.f_lo_hz <= *cfg.f0_hz && *cfg.f0_hz <= cfg.f_hi_hz)) {
        throw std::runtime_error("config: source.f0_hz must lie inside the band");
    }
}

}  // namespace cfar::sim
