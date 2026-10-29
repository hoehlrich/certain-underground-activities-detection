// Milestone 3.2.1: tone placement and ENBW of the Hann periodogram.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "cfar/sim/config.hpp"
#include "cfar/sim/noise.hpp"
#include "cfar/sim/signal.hpp"
#include "cfar/spectrum.hpp"

namespace {

constexpr std::size_t kN = 8192;
constexpr double kFs = 4000.0;

std::vector<double> tone_spectrum(double freq_hz, double amplitude) {
    std::vector<double> x(kN, 0.0);
    cfar::sim::add_tone(x, kFs, {.amplitude = amplitude, .freq_hz = freq_hz});
    return cfar::periodogram(x, cfar::hann_window(kN));
}

std::size_t argmax(const std::vector<double>& v) {
    return static_cast<std::size_t>(std::distance(v.begin(), std::ranges::max_element(v)));
}

// Mean and variance of noise-only P[k] over interior bins (no DC, no Nyquist).
struct NoiseStats {
    double mean;
    double var;
};

NoiseStats noise_stats(double sigma, std::size_t windows) {
    const auto w = cfar::hann_window(kN);
    double sum = 0.0;
    double sum_sq = 0.0;
    std::size_t count = 0;
    for (std::size_t trial = 0; trial < windows; ++trial) {
        auto rng = cfar::sim::make_trial_rng(42, trial);
        std::vector<double> x(kN, 0.0);
        cfar::sim::add_white_noise(x, sigma, rng);
        const auto p = cfar::periodogram(x, w);
        for (std::size_t k = 1; k + 1 < p.size(); ++k) {
            sum += p[k];
            sum_sq += p[k] * p[k];
            ++count;
        }
    }
    const double mean = sum / static_cast<double>(count);
    return {mean, sum_sq / static_cast<double>(count) - mean * mean};
}

}  // namespace

TEST(HannWindow, PeriodicAndEnbwIsOnePointFive) {
    const auto w = cfar::hann_window(kN);
    ASSERT_EQ(w.size(), kN);
    EXPECT_DOUBLE_EQ(w[0], 0.0);
    EXPECT_NEAR(w[kN / 2], 1.0, 1e-15);
    for (std::size_t k = 1; k < kN; ++k) {
        EXPECT_NEAR(w[k], w[kN - k], 1e-15);
    }
    EXPECT_NEAR(cfar::enbw_bins(w), 1.5, 1e-12);
}

TEST(Periodogram, OnBinToneLandsInItsBin) {
    const std::size_t k0 = 1000;
    const double a = 2.0;
    const auto p = tone_spectrum(static_cast<double>(k0) * cfar::bin_hz(kN, kFs), a);

    ASSERT_EQ(p.size(), kN / 2 + 1);
    EXPECT_EQ(argmax(p), k0);
    const double expected_peak = a * a / 4.0 * static_cast<double>(kN) / 1.5;
    EXPECT_NEAR(p[k0] / expected_peak, 1.0, 1e-9);
    // Hann main lobe: neighbours are half the amplitude, a quarter of the power.
    EXPECT_NEAR(p[k0 - 1] / p[k0], 0.25, 1e-9);
    EXPECT_NEAR(p[k0 + 1] / p[k0], 0.25, 1e-9);
}

TEST(Periodogram, ToneInHzLandsInNearestBin) {
    const double f0 = 123.4;
    const auto p = tone_spectrum(f0, 1.0);
    EXPECT_EQ(argmax(p), static_cast<std::size_t>(std::lround(f0 / cfar::bin_hz(kN, kFs))));
}

TEST(Periodogram, HalfBinScallopingLoss) {
    const double df = cfar::bin_hz(kN, kFs);
    const std::size_t k0 = 1000;
    const auto centred = tone_spectrum(static_cast<double>(k0) * df, 1.0);
    const auto offset = tone_spectrum((static_cast<double>(k0) + 0.5) * df, 1.0);

    const double peak = *std::ranges::max_element(offset);
    const double loss_db = 10.0 * std::log10(centred[k0] / peak);
    EXPECT_NEAR(loss_db, 1.42, 0.01);
}

TEST(Periodogram, WhiteNoiseBinsAreExponentialWithMeanSigmaSquared) {
    const double sigma = 1.7;
    const auto s = noise_stats(sigma, 200);
    EXPECT_NEAR(s.mean / (sigma * sigma), 1.0, 0.01);
    // Exponential distribution: variance equals mean squared.
    EXPECT_NEAR(s.var / (s.mean * s.mean), 1.0, 0.03);
}

// 3.2.1 done criterion: ENBW measured from the spectrum is 1.5 bins.
TEST(Periodogram, MeasuredEnbwIsOnePointFive) {
    const double sigma = 1.0;
    const std::size_t k0 = 1000;
    const double a = 1.0;
    const auto p = tone_spectrum(static_cast<double>(k0) * cfar::bin_hz(kN, kFs), a);
    const auto s = noise_stats(sigma, 200);

    // Tone peak gain relative to the noise power each bin collects.
    const double enbw_meas =
        a * a / 4.0 * static_cast<double>(kN) / p[k0] * (s.mean / (sigma * sigma));
    EXPECT_NEAR(enbw_meas, 1.5, 0.015);
}

TEST(Periodogram, Cn0MapsToPerBinSnr) {
    const double sigma = 1.0;
    const double cn0_dbhz = 30.0;
    const std::size_t k0 = 1000;
    const double a = cfar::sim::amplitude_for_cn0(cn0_dbhz, sigma, kFs);
    const auto p = tone_spectrum(static_cast<double>(k0) * cfar::bin_hz(kN, kFs), a);

    // SNR_bin [dB] = C/N0 + 10 log10(T) - 10 log10(ENBW), on-bin so no scalloping.
    const double t = static_cast<double>(kN) / kFs;
    const double snr_bin_db = cn0_dbhz + 10.0 * std::log10(t) - 10.0 * std::log10(1.5);
    EXPECT_NEAR(10.0 * std::log10(p[k0] / (sigma * sigma)), snr_bin_db, 1e-9);
}

TEST(Noise, SameSeedAndTrialReproduces) {
    auto draw = [](std::uint64_t seed, std::uint64_t trial) {
        auto rng = cfar::sim::make_trial_rng(seed, trial);
        std::vector<double> x(256, 0.0);
        cfar::sim::add_white_noise(x, 1.0, rng);
        return x;
    };
    EXPECT_EQ(draw(7, 3), draw(7, 3));
    EXPECT_NE(draw(7, 3), draw(7, 4));
    EXPECT_NE(draw(7, 3), draw(8, 3));
}

TEST(Config, BaselineLoads) {
    const auto cfg = cfar::sim::load_config(CFAR_SOURCE_DIR "/configs/baseline.toml");
    EXPECT_DOUBLE_EQ(cfg.f_lo_hz, 10.0);
    EXPECT_DOUBLE_EQ(cfg.f_hi_hz, 1000.0);
    EXPECT_DOUBLE_EQ(cfg.fs_hz, 4000.0);
    EXPECT_EQ(cfg.nfft, 8192u);
    EXPECT_EQ(cfg.window, "hann");
    EXPECT_EQ(cfg.source_kind, "tone");
    EXPECT_FALSE(cfg.f0_hz.has_value());  // "uniform"
}

TEST(Config, FixedToneOutsideBandThrows) {
    cfar::sim::SimConfig cfg;
    cfg.f0_hz = 1500.0;
    EXPECT_THROW(cfar::sim::validate(cfg), std::runtime_error);
    cfg.f0_hz = 60.0;
    EXPECT_NO_THROW(cfar::sim::validate(cfg));
}

TEST(Config, SampleRateBelowSub032Throws) {
    cfar::sim::SimConfig cfg;
    cfg.fs_hz = 1000.0;
    cfg.f_hi_hz = 500.0;
    EXPECT_THROW(cfar::sim::validate(cfg), std::runtime_error);
}
