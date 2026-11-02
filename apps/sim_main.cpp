// One window of tone plus white noise through the detector, written as a spectrum CSV.
//
// Usage: sim_main [config.toml] [output.csv]

#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <random>
#include <vector>

#include "core/cfar.hpp"
#include "core/config.hpp"
#include "core/detector.hpp"
#include "sim/config.hpp"
#include "sim/noise.hpp"
#include "sim/signal.hpp"
#include "core/spectrum.hpp"

int main(int argc, char** argv) {
    const std::filesystem::path config_path = argc > 1 ? argv[1] : "configs/baseline.toml";
    const std::filesystem::path out_path = argc > 2 ? argv[2] : "results/example_spectrum.csv";

    try {
        const cfar::sim::SimConfig cfg = cfar::sim::load_config(config_path);
        const cfar::ProcessingConfig& proc = cfg.proc;
        auto rng = cfar::sim::make_trial_rng(cfg.seed, 0);

        cfar::sim::ToneParams tone;
        tone.freq_hz = cfg.f0_hz.value_or(
            std::uniform_real_distribution<double>(proc.f_lo_hz, proc.f_hi_hz)(rng));
        tone.phase_rad = std::uniform_real_distribution<double>(0.0, 2.0 * std::numbers::pi)(rng);
        tone.amplitude = cfar::sim::amplitude_for_cn0(cfg.cn0_dbhz, cfg.noise_sigma, proc.fs_hz);

        std::vector<double> x(proc.nfft, 0.0);
        cfar::sim::add_tone(x, proc.fs_hz, tone);
        cfar::sim::add_white_noise(x, cfg.noise_sigma, rng);

        const cfar::Detector det(proc);
        const auto power = cfar::periodogram(x, det.window());
        const double df = det.bin_hz();
        const cfar::CaCfar& ca = det.cfar();
        const auto noise = cfar::ca_noise(power, det.band(), proc.guard_per_side,
                                          proc.train_per_side);
        const auto threshold = ca.threshold(noise);
        const auto detections =
            cfar::find_detections(power, threshold, noise.level, det.band(), df);

        if (out_path.has_parent_path()) {
            std::filesystem::create_directories(out_path.parent_path());
        }
        std::ofstream out(out_path);
        out << "freq_hz,power,threshold,detect\n";
        out.precision(10);
        std::size_t hits = 0;
        for (std::size_t k = 0; k < power.size(); ++k) {
            const bool detect = power[k] > threshold[k];  // false outside the band (NaN)
            hits += detect ? 1 : 0;
            out << static_cast<double>(k) * df << ',' << power[k] << ',';
            if (std::isnan(threshold[k])) {
                out << ",0\n";
            } else {
                out << threshold[k] << ',' << (detect ? 1 : 0) << '\n';
            }
        }

        const auto tone_bin = static_cast<std::size_t>(std::lround(tone.freq_hz / df));
        std::cout << "f0 = " << tone.freq_hz << " Hz, expected peak bin " << tone_bin << " ("
                  << df << " Hz bins)\n"
                  << "design Pfa per bin " << ca.params().pfa_bin << ", alpha "
                  << ca.alpha(cfar::Sides::both) << " (two-sided), "
                  << ca.alpha(cfar::Sides::lead_only) << " (one-sided)\n"
                  << hits << " bins over threshold in " << detections.size()
                  << " detection(s):\n";
        for (const auto& d : detections) {
            std::cout << "  " << d.freq_hz << " Hz (bins " << d.first_bin << '-' << d.last_bin
                      << "), statistic " << 10.0 * std::log10(d.statistic) << " dB, SNR est "
                      << 10.0 * std::log10(d.snr_est) << " dB"
                      << (std::abs(d.freq_hz - tone.freq_hz) <= df ? ", the tone" : ", false alarm")
                      << '\n';
        }
        std::cout << "wrote " << out_path.string() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "sim_main: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
