#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cfar/spectrum.hpp"

namespace cfar {

// Which training sides a bin used. Near the band edges only one side fits.
enum class Sides : std::uint8_t { none, both, lead_only, lag_only };

// Per-bin noise estimate over the search band. level[k] and sides[k] cover the
// whole spectrum; bins outside the band get NaN and Sides::none.
//
// This is the interface other estimators (OS, GO, SO) will share later: a
// function from P[k] and the band to a noise level per bin.
struct NoiseEstimate {
    std::vector<double> level;
    std::vector<Sides> sides;
};

// Cell-averaging estimate: mean of `train` cells on each side of bin k, skipping
// `guard` cells next to it. Training cells stay inside the band, so a bin within
// guard + train of a band edge uses the other side only. O(1) per bin via prefix
// sums. Throws std::invalid_argument if some bin in the band fits neither side.
NoiseEstimate ca_noise(std::span<const double> power, BinRange band, std::size_t guard,
                       std::size_t train);

// Complex correlation between windowed DFT bins m apart, for white input:
// rho[m] = sum(w^2 cos(2 pi m n / N)) / sum(w^2), m = 0..max_lag.
// Hann gives 1, -2/3, 1/6, then 0.
std::vector<double> bin_correlation(std::span<const double> window, std::size_t max_lag);

// Exact CA false-alarm probability in white Gaussian noise for threshold
// alpha * (mean of the n_lead + n_lag training cells), with training cells
// correlated as rho (the test cell must be independent of them, which a guard of
// at least rho.size() - 1 bins ensures):
//
//   Pfa = 1 / (det(I + (alpha/N) C_lead) * det(I + (alpha/N) C_lag)),  N = n_lead + n_lag
//
// C_side is the Toeplitz correlation matrix of one side's cells. With rho = {1}
// (independent bins) this is the textbook (1 + alpha/N)^-N.
double ca_pfa(double alpha, std::size_t n_lead, std::size_t n_lag, std::span<const double> rho);

// Inverts ca_pfa for alpha by bisection.
double ca_alpha(double pfa, std::size_t n_lead, std::size_t n_lag, std::span<const double> rho);

// Per-bin Pfa giving pfa_window over n_bins independent tests:
// 1 - (1 - pfa_window)^(1 / n_bins).
double pfa_bin_for_window(double pfa_window, std::size_t n_bins);

struct CaParams {
    std::size_t guard_per_side = 2;
    std::size_t train_per_side = 16;
    double pfa_bin = 5e-6;
};

// CA-CFAR with alpha corrected for the window's bin correlation. Alpha is
// computed once, for two-sided and one-sided training.
class CaCfar {
public:
    CaCfar(const CaParams& params, std::span<const double> window);

    // Noise estimate times alpha for each bin in the band; NaN outside it.
    std::vector<double> threshold(std::span<const double> power, BinRange band) const;
    std::vector<double> threshold(const NoiseEstimate& noise) const;

    double alpha(Sides sides) const;
    const CaParams& params() const { return params_; }

private:
    CaParams params_;
    double alpha_both_;
    double alpha_one_;
};

}  // namespace cfar
