// Milestone 3.2.1: one window of tone plus white noise, written as a spectrum CSV.
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

#include "cfar/sim/config.hpp"
#include "cfar/sim/noise.hpp"
#include "cfar/sim/signal.hpp"
#include "cfar/spectrum.hpp"

int main(int argc, char** argv) {
    const std::filesystem::path config_path = argc > 1 ? argv[1] : "configs/baseline.toml";
    const std::filesystem::path out_path = argc > 2 ? argv[2] : "results/example_spectrum.csv";

    try {
        const cfar::sim::SimConfig cfg = cfar::sim::load_config(config_path);
        auto rng = cfar::sim::make_trial_rng(cfg.seed, 0);

        cfar::sim::ToneParams tone;
        tone.freq_hz = cfg.f0_hz.value_or(
            std::uniform_real_distribution<double>(cfg.f_lo_hz, cfg.f_hi_hz)(rng));
        tone.phase_rad = std::uniform_real_distribution<double>(0.0, 2.0 * std::numbers::pi)(rng);
        tone.amplitude = cfar::sim::amplitude_for_cn0(cfg.cn0_dbhz, cfg.noise_sigma, cfg.fs_hz);

        std::vector<double> x(cfg.nfft, 0.0);
        cfar::sim::add_tone(x, cfg.fs_hz, tone);
        cfar::sim::add_white_noise(x, cfg.noise_sigma, rng);

        const auto window = cfar::hann_window(cfg.nfft);
        const auto power = cfar::periodogram(x, window);
        const double df = cfar::bin_hz(cfg.nfft, cfg.fs_hz);

        if (out_path.has_parent_path()) {
            std::filesystem::create_directories(out_path.parent_path());
        }
        std::ofstream out(out_path);
        out << "freq_hz,power\n";
        out.precision(10);
        for (std::size_t k = 0; k < power.size(); ++k) {
            out << static_cast<double>(k) * df << ',' << power[k] << '\n';
        }

        std::cout << "f0 = " << tone.freq_hz << " Hz, expected peak bin "
                  << std::lround(tone.freq_hz / df) << " (" << df << " Hz bins)\n"
                  << "wrote " << out_path.string() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "sim_main: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
