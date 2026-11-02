// ROC for CA-CFAR on a tone in white noise: Pd against the per-window false-alarm
// probability at fixed SNRs, with the known-noise bound.
//
// Usage: roc [config.toml] [output.csv]

#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "cfar/cfar.hpp"
#include "cfar/config.hpp"
#include "cfar/sim/config.hpp"
#include "cfar/sim/pd.hpp"

int main(int argc, char** argv) {
    const std::filesystem::path config_path = argc > 1 ? argv[1] : "configs/baseline.toml";
    const std::filesystem::path out_path = argc > 2 ? argv[2] : "results/roc.csv";

    try {
        const auto cfg = cfar::sim::load_config(config_path);
        const auto& proc = cfg.proc;
        const std::size_t n_bins = cfar::band_bins(proc).size();
        std::vector<double> pfa_bins;
        for (double p : cfg.roc_pfa_window) {
            pfa_bins.push_back(cfar::pfa_bin_for_window(p, n_bins));
        }

        std::printf("%zu trials per point, %zu SNRs x %zu Pfa values, %zu bins in band\n",
                    cfg.roc_trials, cfg.roc_snr_bin_db.size(), pfa_bins.size(), n_bins);
        const auto points =
            cfar::sim::run_pd_grid(cfg, cfg.roc_snr_bin_db, pfa_bins, cfg.roc_trials);

        if (out_path.has_parent_path()) {
            std::filesystem::create_directories(out_path.parent_path());
        }
        std::ofstream out(out_path);
        out << "scenario,snr_bin_db,cn0_dbhz,pfa_window,pfa_bin,trials,hits,pd,ci95_lo,ci95_hi,"
               "pd_ideal,fa_windows,nfft,guard,train\n";
        out.precision(8);
        std::printf("SNR_bin  C/N0    Pfa/window  Pfa/bin    Pd      95%% interval      ideal   "
                    "FA windows\n");
        for (const auto& p : points) {
            const auto ci = p.interval();
            out << "white," << p.snr_bin_db << ',' << p.cn0_dbhz << ',' << p.pfa_window << ','
                << p.pfa_bin << ',' << p.trials << ',' << p.hits << ',' << p.pd() << ',' << ci.lo
                << ',' << ci.hi << ',' << p.ideal_pd << ',' << p.fa_windows << ',' << proc.nfft
                << ',' << proc.guard_per_side << ',' << proc.train_per_side << '\n';
            std::printf("%5.1f   %6.2f   %-10.3g  %-9.3g  %.4f  [%.4f, %.4f]  %.4f  %llu\n",
                        p.snr_bin_db, p.cn0_dbhz, p.pfa_window, p.pfa_bin, p.pd(), ci.lo, ci.hi,
                        p.ideal_pd, static_cast<unsigned long long>(p.fa_windows));
        }
        std::cout << "wrote " << out_path.string() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "roc: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
