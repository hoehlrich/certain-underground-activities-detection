// Milestone 3.2.2: CA-CFAR noise estimate, alpha, and the Pfa check.

#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "cfar/cfar.hpp"
#include "cfar/config.hpp"
#include "cfar/sim/config.hpp"
#include "cfar/sim/pfa.hpp"
#include "cfar/spectrum.hpp"

namespace {

constexpr std::size_t kN = 8192;
constexpr double kFs = 4000.0;

}  // namespace

TEST(Band, BaselineBandHas2028Bins) {
    const auto band = cfar::band_bins(10.0, 1000.0, kN, kFs);
    EXPECT_EQ(band.first, 21u);    // 10 Hz / 0.488 Hz = 20.48, rounded up
    EXPECT_EQ(band.last, 2048u);   // exactly 1000 Hz
    EXPECT_EQ(band.size(), 2028u);
    EXPECT_EQ(cfar::band_bins(10.0, 1000.0, 4096, kFs).size(), 1014u);
    EXPECT_THROW(cfar::band_bins(10.0, 2500.0, kN, kFs), std::invalid_argument);  // past Nyquist
}

TEST(FalseAlarmBudget, DerivedPfaPerBin) {
    const double p = cfar::pfa_bin_for_window(0.01, 2028);
    EXPECT_NEAR(p, 4.956e-6, 0.001e-6);
    EXPECT_NEAR(1.0 - std::pow(1.0 - p, 2028), 0.01, 1e-12);

    cfar::ProcessingConfig cfg;
    EXPECT_DOUBLE_EQ(cfar::design_pfa_bin(cfg), p);
    cfg.pfa_bin = 1e-4;
    EXPECT_DOUBLE_EQ(cfar::design_pfa_bin(cfg), 1e-4);
}

TEST(CaNoise, MatchesBruteForceIncludingEdges) {
    std::mt19937_64 rng(5);
    std::exponential_distribution<double> exp1(1.0);
    std::vector<double> p(400);
    for (double& v : p) {
        v = exp1(rng);
    }
    const cfar::BinRange band{50, 300};
    const std::size_t g = 2;
    const std::size_t t = 16;
    const auto est = cfar::ca_noise(p, band, g, t);

    for (std::size_t k = 0; k < p.size(); ++k) {
        if (k < band.first || k > band.last) {
            EXPECT_TRUE(std::isnan(est.level[k]));
            EXPECT_EQ(est.sides[k], cfar::Sides::none);
            continue;
        }
        const bool lead = k >= band.first + g + t;
        const bool lag = k + g + t <= band.last;
        double sum = 0.0;
        std::size_t cells = 0;
        for (std::size_t i = 1; i <= t; ++i) {
            if (lead) {
                sum += p[k - g - i];
                ++cells;
            }
            if (lag) {
                sum += p[k + g + i];
                ++cells;
            }
        }
        ASSERT_GT(cells, 0u);
        EXPECT_NEAR(est.level[k], sum / static_cast<double>(cells), 1e-12) << "bin " << k;
        EXPECT_EQ(est.sides[k], lead && lag ? cfar::Sides::both
                                            : (lead ? cfar::Sides::lead_only : cfar::Sides::lag_only));
    }
}

TEST(CaNoise, BandTooNarrowThrows) {
    std::vector<double> p(100, 1.0);
    EXPECT_THROW(cfar::ca_noise(p, {40, 60}, 2, 16), std::invalid_argument);
}

TEST(CaAlpha, IndependentBinsGiveTextbookFormula) {
    const std::vector<double> independent{1.0};
    for (double pfa : {1e-3, 1e-5, 5e-6}) {
        for (std::size_t n : {16u, 32u}) {
            const double textbook =
                static_cast<double>(n) * (std::pow(pfa, -1.0 / static_cast<double>(n)) - 1.0);
            EXPECT_NEAR(cfar::ca_alpha(pfa, n / 2, n - n / 2, independent) / textbook, 1.0, 1e-9);
            EXPECT_NEAR(cfar::ca_pfa(textbook, n / 2, n - n / 2, independent) / pfa, 1.0, 1e-9);
        }
    }
}

TEST(CaAlpha, HannBinCorrelation) {
    const auto rho = cfar::bin_correlation(cfar::hann_window(kN), 4);
    EXPECT_NEAR(rho[0], 1.0, 1e-12);
    EXPECT_NEAR(rho[1], -2.0 / 3.0, 1e-12);
    EXPECT_NEAR(rho[2], 1.0 / 6.0, 1e-12);
    EXPECT_EQ(rho[3], 0.0);
    EXPECT_EQ(rho[4], 0.0);
}

TEST(CaAlpha, HannCorrectionRaisesAlpha) {
    // Correlated training cells average less noise away, so alpha must rise to
    // hold the same Pfa. Reference values from an independent numpy eigenvalue
    // computation that matched Monte Carlo.
    const auto rho = cfar::bin_correlation(cfar::hann_window(kN), 16);
    EXPECT_NEAR(cfar::ca_alpha(1e-3, 16, 16, rho), 8.46, 0.01);
    EXPECT_NEAR(cfar::ca_alpha(5e-6, 16, 16, rho), 17.44, 0.01);
    // The textbook alpha (14.86) gives about 4x the design Pfa at 5e-6.
    EXPECT_NEAR(cfar::ca_pfa(14.86, 16, 16, rho) / 5e-6, 4.0, 0.1);
}

TEST(CaCfar, GuardTooSmallForHannThrows) {
    const auto w = cfar::hann_window(kN);
    EXPECT_THROW(cfar::CaCfar({1, 16, 1e-3}, w), std::invalid_argument);
    EXPECT_NO_THROW(cfar::CaCfar({2, 16, 1e-3}, w));
}

TEST(CaCfar, ThresholdIsAlphaTimesNoise) {
    const auto w = cfar::hann_window(kN);
    const cfar::CaCfar ca({2, 16, 1e-4}, w);
    std::vector<double> p(kN / 2 + 1, 2.0);
    const cfar::BinRange band{21, 2048};
    const auto t = ca.threshold(p, band);
    EXPECT_TRUE(std::isnan(t[20]));
    EXPECT_NEAR(t[21], 2.0 * ca.alpha(cfar::Sides::lag_only), 1e-12);
    EXPECT_NEAR(t[1000], 2.0 * ca.alpha(cfar::Sides::both), 1e-12);
    EXPECT_NEAR(t[2048], 2.0 * ca.alpha(cfar::Sides::lead_only), 1e-12);
    EXPECT_GT(ca.alpha(cfar::Sides::lead_only), ca.alpha(cfar::Sides::both));
}

// The 3.2.2 done criterion: in white noise, measured Pfa per bin falls inside the
// 95% binomial interval of the design value for 1e-3, 1e-4, 1e-5 and the derived
// baseline (about 5e-6). Uses the baseline config.
TEST(PfaCheck, MeasuredMatchesDesignInWhiteNoise) {
    const auto cfg = cfar::sim::load_config(CFAR_SOURCE_DIR "/configs/baseline.toml");
    std::vector<double> design = cfg.pfa_check_design;
    design.push_back(cfar::design_pfa_bin(cfg.proc));

    const auto results = cfar::sim::run_pfa_check(cfg, design, cfg.pfa_check_windows);
    for (const auto& r : results) {
        const auto ci = r.interval();
        EXPECT_TRUE(r.pass()) << "design " << r.design_pfa << ": measured " << r.measured()
                              << " outside [" << ci.lo << ", " << ci.hi << "] (" << r.hits
                              << " hits in " << r.cells << " cells)";
    }
}
