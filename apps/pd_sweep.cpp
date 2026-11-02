// Pd vs SNR for CA-CFAR on a tone in white noise, with the known-noise bound.
// Reports the C/N0 needed to reach pd_target (0.9 by default).
//
// Usage: pd_sweep [config.toml] [output.csv]
// Also writes <output>_summary.csv with the C/N0 at pd_target.

#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "core/config.hpp"
#include "sim/config.hpp"
#include "sim/pd.hpp"

int main(int argc, char** argv) {
    const std::filesystem::path config_path = argc > 1 ? argv[1] : "configs/baseline.toml";
    const std::filesystem::path out_path = argc > 2 ? argv[2] : "results/pd_vs_snr.csv";

    try {
        const auto cfg = cfar::sim::load_config(config_path);
        const auto& proc = cfg.proc;
        const auto grid = cfar::sim::pd_snr_grid(cfg);
        const double offset_db = cfar::sim::snr_bin_offset_db(proc);
        const double pfa_bin = cfar::design_pfa_bin(proc);

        std::printf("%zu trials per point, %zu points, design Pfa per bin %.3g, f0 %s\n",
                    cfg.pd_trials, grid.size(), pfa_bin, cfg.f0_hz ? "fixed" : "uniform");
        const auto points = cfar::sim::run_pd_sweep(cfg, grid, cfg.pd_trials);

        if (out_path.has_parent_path()) {
            std::filesystem::create_directories(out_path.parent_path());
        }
        std::ofstream out(out_path);
        out << "scenario,snr_bin_db,cn0_dbhz,trials,hits,pd,ci95_lo,ci95_hi,pd_ideal,"
               "fa_windows,nfft,guard,train,pfa_bin\n";
        out.precision(8);
        std::printf("SNR_bin  C/N0     Pd      95%% interval        ideal   FA windows\n");
        for (const auto& p : points) {
            const auto ci = p.interval();
            out << "white," << p.snr_bin_db << ',' << p.cn0_dbhz << ',' << p.trials << ','
                << p.hits << ',' << p.pd() << ',' << ci.lo << ',' << ci.hi << ',' << p.ideal_pd
                << ',' << p.fa_windows << ',' << proc.nfft << ',' << proc.guard_per_side << ','
                << proc.train_per_side << ',' << pfa_bin << '\n';
            std::printf("%5.1f   %6.2f   %.4f  [%.4f, %.4f]  %.4f  %llu\n", p.snr_bin_db,
                        p.cn0_dbhz, p.pd(), ci.lo, ci.hi, p.ideal_pd,
                        static_cast<unsigned long long>(p.fa_windows));
        }

        const double ca_snr = cfar::sim::crossing_snr_bin_db(points, cfg.pd_target);
        const double ideal_snr = cfar::sim::ideal_bound(cfg).snr_bin_db_for(cfg.pd_target);

        std::printf("\nPd = %.2f at Pfa per bin %.3g (%.3g per window):\n", cfg.pd_target,
                    pfa_bin, proc.pfa_window);
        if (std::isnan(ca_snr)) {
            std::printf("  CA-CFAR  not reached in the sweep\n");
        } else {
            std::printf("  CA-CFAR  C/N0 %.2f dB-Hz (SNR_bin %.2f dB)\n", ca_snr - offset_db,
                        ca_snr);
        }
        std::printf("  ideal    C/N0 %.2f dB-Hz (SNR_bin %.2f dB)\n", ideal_snr - offset_db,
                    ideal_snr);
        if (!std::isnan(ca_snr)) {
            std::printf("  CFAR loss %.2f dB\n", ca_snr - ideal_snr);
        }

        auto summary_path = out_path;
        summary_path.replace_filename(out_path.stem().string() + "_summary.csv");
        std::ofstream summary(summary_path);
        summary.precision(8);
        summary << "pd_target,pfa_bin,pfa_window,snr_offset_db,cn0_ca_dbhz,cn0_ideal_dbhz,"
                   "cfar_loss_db\n"
                << cfg.pd_target << ',' << pfa_bin << ',' << proc.pfa_window << ',' << offset_db
                << ',' << ca_snr - offset_db << ',' << ideal_snr - offset_db << ','
                << ca_snr - ideal_snr << '\n';
        std::cout << "wrote " << out_path.string() << " and " << summary_path.string() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "pd_sweep: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
