#include "cfar/sim/config.hpp"

#include <stdexcept>

#include <toml++/toml.hpp>

namespace cfar::sim {

SimConfig load_config(const std::filesystem::path& path) {
    SimConfig cfg;
    cfg.proc = load_processing_config(path);

    toml::table tbl;
    try {
        tbl = toml::parse_file(path.string());
    } catch (const toml::parse_error& e) {
        throw std::runtime_error("config " + path.string() + ": " + std::string(e.description()));
    }

    cfg.noise_sigma = tbl["noise"]["sigma"].value_or(cfg.noise_sigma);
    cfg.source_kind = tbl["source"]["kind"].value_or(cfg.source_kind);
    cfg.cn0_dbhz = tbl["source"]["cn0_dbhz"].value_or(cfg.cn0_dbhz);
    cfg.seed = tbl["run"]["seed"].value_or(cfg.seed);
    cfg.pfa_check_windows = tbl["pfa_check"]["windows"].value_or(cfg.pfa_check_windows);

    // f0_hz is either "uniform" or a fixed frequency.
    const auto f0 = tbl["source"]["f0_hz"];
    if (auto hz = f0.value<double>()) {
        cfg.f0_hz = *hz;
    } else if (f0 && f0.value_or(std::string{}) != "uniform") {
        throw std::runtime_error("config: source.f0_hz must be \"uniform\" or a number");
    }

    if (const auto* design = tbl["pfa_check"]["design_pfa"].as_array()) {
        cfg.pfa_check_design.clear();
        for (const auto& v : *design) {
            const auto p = v.value<double>();
            if (!p) {
                throw std::runtime_error("config: pfa_check.design_pfa must hold numbers");
            }
            cfg.pfa_check_design.push_back(*p);
        }
    }

    validate(cfg);
    return cfg;
}

void validate(const SimConfig& cfg) {
    validate(cfg.proc);
    if (!(cfg.noise_sigma > 0.0)) {
        throw std::runtime_error("config: noise sigma must be positive");
    }
    if (cfg.source_kind != "tone") {
        throw std::runtime_error("config: only source.kind = \"tone\" is supported");
    }
    if (cfg.f0_hz && !(cfg.proc.f_lo_hz <= *cfg.f0_hz && *cfg.f0_hz <= cfg.proc.f_hi_hz)) {
        throw std::runtime_error("config: source.f0_hz must lie inside the band");
    }
    if (cfg.pfa_check_windows == 0) {
        throw std::runtime_error("config: pfa_check.windows must be at least 1");
    }
    for (double p : cfg.pfa_check_design) {
        if (!(0.0 < p && p < 1.0)) {
            throw std::runtime_error("config: pfa_check.design_pfa values must be in (0, 1)");
        }
    }
}

}  // namespace cfar::sim
