// Milestone 3.2.3: detector logic, the ideal bound, and Pd vs SNR.

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "cfar/config.hpp"
#include "cfar/detector.hpp"
#include "cfar/sim/config.hpp"
#include "cfar/sim/noise.hpp"
#include "cfar/sim/pd.hpp"
#include "cfar/sim/signal.hpp"
#include "cfar/sim/theory.hpp"
#include "cfar/spectrum.hpp"

namespace {

constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

}  // namespace

TEST(MarcumQ1, MatchesReferenceValues) {
    // Reference values from scipy.stats.ncx2.sf(b^2, 2, a^2).
    EXPECT_NEAR(cfar::sim::marcum_q1(1.0, 1.0), 0.7328798037968203, 1e-12);
    EXPECT_NEAR(cfar::sim::marcum_q1(2.0, 3.0), 0.21436208816264943, 1e-12);
    EXPECT_NEAR(cfar::sim::marcum_q1(5.0, 4.94), 0.5640714892980341, 1e-12);
    EXPECT_NEAR(cfar::sim::marcum_q1(0.5, 4.94) / 1.5405022371008636e-05, 1.0, 1e-9);
    EXPECT_NEAR(cfar::sim::marcum_q1(10.0, 4.94), 0.9999998550828115, 1e-12);
    // Edge cases: no signal gives the Rayleigh tail; zero threshold always detects.
    EXPECT_NEAR(cfar::sim::marcum_q1(0.0, 3.0), std::exp(-4.5), 1e-15);
    EXPECT_EQ(cfar::sim::marcum_q1(3.0, 0.0), 1.0);
}

TEST(IdealBound, HannScallopingAndOnBinPd) {
    const auto w = cfar::hann_window(8192);
    EXPECT_NEAR(cfar::sim::scallop_gain(w, 0.0), 1.0, 1e-12);
    EXPECT_NEAR(-10.0 * std::log10(cfar::sim::scallop_gain(w, 0.5)), 1.42, 0.01);

    // On-bin bound at 5e-6 reaches Pd = 0.9 at 12.74 dB (scipy brentq on ncx2).
    const std::vector<double> on_bin{0.0};
    EXPECT_NEAR(cfar::sim::IdealBound(w, 4.956e-6, on_bin).snr_bin_db_for(0.9), 12.743, 0.005);
    // Averaging over offsets adds the mean scalloping penalty: 13.28 dB.
    EXPECT_NEAR(cfar::sim::IdealBound(w, 4.956e-6, cfar::sim::uniform_offsets())
                    .snr_bin_db_for(0.9),
                13.285, 0.005);
}

// The bound only means something if the sim's tone and noise scaling matches it:
// an on-bin tone with a known-noise threshold must cross at the Q1 rate.
TEST(IdealBound, MatchesKnownNoiseMonteCarlo) {
    const std::size_t n = 8192;
    const double fs = 4000.0;
    const std::size_t k0 = 1000;
    const double pfa = 1e-3;
    const double snr_bin_db = 10.0;
    const auto w = cfar::hann_window(n);
    const double t = static_cast<double>(n) / fs;
    const double cn0 = snr_bin_db - 10.0 * std::log10(t) + 10.0 * std::log10(1.5);
    const double amp = cfar::sim::amplitude_for_cn0(cn0, 1.0, fs);
    const double threshold = -std::log(pfa);  // sigma = 1

    const std::size_t trials = 4000;
    std::uint64_t hits = 0;
    for (std::size_t i = 0; i < trials; ++i) {
        auto rng = cfar::sim::make_trial_rng(11, i);
        std::vector<double> x(n, 0.0);
        const double phase = std::uniform_real_distribution<double>(0.0, 6.283185307179586)(rng);
        cfar::sim::add_tone(x, fs, {.amplitude = amp,
                                    .freq_hz = static_cast<double>(k0) * fs / static_cast<double>(n),
                                    .phase_rad = phase});
        cfar::sim::add_white_noise(x, 1.0, rng);
        hits += cfar::periodogram(x, w)[k0] > threshold ? 1 : 0;
    }
    const double expected = cfar::sim::ideal_pd(std::pow(10.0, snr_bin_db / 10.0), pfa);
    const auto ci = cfar::sim::wilson_interval_95(hits, trials);
    EXPECT_TRUE(ci.contains(expected)) << "Q1 " << expected << " outside [" << ci.lo << ", "
                                       << ci.hi << "]";
}

TEST(WilsonInterval, StaysInsideZeroOne) {
    const auto none = cfar::sim::wilson_interval_95(0, 2000);
    EXPECT_NEAR(none.lo, 0.0, 1e-15);
    EXPECT_GT(none.hi, 0.0);
    const auto all = cfar::sim::wilson_interval_95(2000, 2000);
    EXPECT_NEAR(all.hi, 1.0, 1e-15);
    EXPECT_LT(all.lo, 1.0);
}

TEST(Detections, ClustersAdjacentBins) {
    std::vector<double> p(64, 1.0);
    std::vector<double> noise(64, 1.0);
    std::vector<double> thr(64, 5.0);
    thr[0] = kNan;
    p[10] = 8.0;
    p[11] = 20.0;
    p[12] = 6.0;
    p[30] = 9.0;
    p[0] = 100.0;  // NaN threshold never detects
    const auto d = cfar::find_detections(p, thr, noise, {1, 62}, 0.5);
    ASSERT_EQ(d.size(), 2u);
    EXPECT_EQ(d[0].first_bin, 10u);
    EXPECT_EQ(d[0].last_bin, 12u);
    EXPECT_EQ(d[0].peak_bin, 11u);
    EXPECT_DOUBLE_EQ(d[0].statistic, 20.0);
    EXPECT_DOUBLE_EQ(d[0].snr_est, 19.0);
    EXPECT_EQ(d[1].first_bin, 30u);
    EXPECT_EQ(d[1].last_bin, 30u);
    EXPECT_DOUBLE_EQ(d[1].freq_hz, 15.0);  // equal neighbours: no offset
}

TEST(Detections, ParabolicInterpolationRefinesHannTone) {
    const std::size_t n = 8192;
    const double fs = 4000.0;
    const double df = cfar::bin_hz(n, fs);
    const auto w = cfar::hann_window(n);
    for (double off : {0.0, 0.1, 0.25, 0.4, -0.3}) {
        std::vector<double> x(n, 0.0);
        const double f0 = (1000.0 + off) * df;
        cfar::sim::add_tone(x, fs, {.amplitude = 1.0, .freq_hz = f0});
        const auto p = cfar::periodogram(x, w);
        const double est = 1000.0 + cfar::parabolic_offset(p[999], p[1000], p[1001]);
        EXPECT_NEAR(est, 1000.0 + off, 0.05) << "offset " << off;
    }
}

TEST(Detector, StrongToneGivesOneDetectionNearF0) {
    cfar::ProcessingConfig proc;
    const cfar::Detector det(proc);
    const double f0 = 423.9;
    const double snr_bin_db = 30.0;
    const double cn0 = snr_bin_db - cfar::sim::snr_bin_offset_db(proc);
    std::vector<double> x(proc.nfft, 0.0);
    cfar::sim::add_tone(x, proc.fs_hz,
                        {.amplitude = cfar::sim::amplitude_for_cn0(cn0, 1.0, proc.fs_hz),
                         .freq_hz = f0});
    auto rng = cfar::sim::make_trial_rng(3, 0);
    cfar::sim::add_white_noise(x, 1.0, rng);

    const auto d = det.process(x);
    ASSERT_EQ(d.size(), 1u);
    EXPECT_NEAR(d[0].freq_hz, f0, 0.1 * det.bin_hz());
    const double offset = f0 / det.bin_hz() - std::round(f0 / det.bin_hz());
    const double expected_db = snr_bin_db + 10.0 * std::log10(cfar::sim::scallop_gain(det.window(), offset));
    EXPECT_NEAR(10.0 * std::log10(d[0].snr_est), expected_db, 2.0);
}

// 3.2.3 sanity: CA-CFAR stays below the known-noise bound and within 2 dB of it
// around Pd = 0.9 at the baseline config.
TEST(PdSweep, CaSitsJustBelowIdealBound) {
    auto cfg = cfar::sim::load_config(CFAR_SOURCE_DIR "/configs/baseline.toml");
    const std::vector<double> grid{13.0, 15.0, 17.0};
    const auto points = cfar::sim::run_pd_sweep(cfg, grid, 1000);
    const auto ideal = cfar::sim::ideal_bound(cfg);
    for (const auto& p : points) {
        const auto ci = p.interval();
        EXPECT_LT(ci.lo, p.ideal_pd) << "SNR " << p.snr_bin_db;
        EXPECT_GT(ci.hi, ideal.pd(p.snr_bin_db - 2.0)) << "SNR " << p.snr_bin_db;
        EXPECT_LT(p.fa_windows, 30u) << "SNR " << p.snr_bin_db;  // about 1% of 1000 expected
    }
}
