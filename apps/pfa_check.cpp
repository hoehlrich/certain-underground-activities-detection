// Milestone 3.2.2 Pfa check: noise-only CA-CFAR, measured vs design Pfa per bin.
//
// Usage: pfa_check [config.toml] [output.csv]

#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "core/config.hpp"
#include "sim/config.hpp"
#include "sim/pfa.hpp"

int main(int argc, char** argv) {
    const std::filesystem::path config_path = argc > 1 ? argv[1] : "configs/baseline.toml";
    const std::filesystem::path out_path = argc > 2 ? argv[2] : "results/pfa_check.csv";

    try {
        const cfar::sim::SimConfig cfg = cfar::sim::load_config(config_path);
        const cfar::ProcessingConfig& proc = cfg.proc;

        std::vector<double> design = cfg.pfa_check_design;
        design.push_back(cfar::design_pfa_bin(proc));
        const auto results = cfar::sim::run_pfa_check(cfg, design, cfg.pfa_check_windows);

        if (out_path.has_parent_path()) {
            std::filesystem::create_directories(out_path.parent_path());
        }
        std::ofstream out(out_path);
        out.precision(6);
        out << "nfft,guard_per_side,train_per_side,design_pfa,cells,hits,measured_pfa,"
               "ci95_lo,ci95_hi,pass,windows,window_fa_rate\n";

        std::printf("%zu windows, %zu bins in band, guard %zu, train %zu per side\n",
                    cfg.pfa_check_windows, cfar::band_bins(proc).size(), proc.guard_per_side,
                    proc.train_per_side);
        std::printf("%-10s %-11s %-25s %-6s %s\n", "design", "measured", "95% interval", "pass",
                    "window FA rate");
        bool all_pass = true;
        for (const auto& r : results) {
            const auto ci = r.interval();
            const double window_rate =
                static_cast<double>(r.window_hits) / static_cast<double>(r.windows);
            all_pass = all_pass && r.pass();
            out << proc.nfft << ',' << proc.guard_per_side << ',' << proc.train_per_side << ','
                << r.design_pfa << ',' << r.cells << ',' << r.hits << ',' << r.measured() << ','
                << ci.lo << ',' << ci.hi << ',' << (r.pass() ? 1 : 0) << ',' << r.windows << ','
                << window_rate << '\n';
            std::printf("%-10.3g %-11.4g [%-10.4g, %-10.4g]  %-6s %.4g\n", r.design_pfa,
                        r.measured(), ci.lo, ci.hi, r.pass() ? "yes" : "NO", window_rate);
        }
        std::cout << "wrote " << out_path.string() << '\n';
        return all_pass ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "pfa_check: " << e.what() << '\n';
        return 1;
    }
}
