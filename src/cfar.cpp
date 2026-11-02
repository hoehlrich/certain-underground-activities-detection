#include "cfar/cfar.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace cfar {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// log det(I + t C) for the n x n Toeplitz correlation matrix C built from rho,
// by Cholesky. I + t C is positive definite for t >= 0 because C is a covariance.
double log_det_i_plus_tc(std::size_t n, double t, std::span<const double> rho) {
    if (n == 0) {
        return 0.0;
    }
    auto c = [&](std::size_t i, std::size_t j) {
        const std::size_t lag = i > j ? i - j : j - i;
        const double r = lag < rho.size() ? rho[lag] : 0.0;
        return (i == j ? 1.0 : 0.0) + t * r;
    };
    std::vector<double> l(n * n, 0.0);
    double log_det = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        double d = c(j, j);
        for (std::size_t k = 0; k < j; ++k) {
            d -= l[j * n + k] * l[j * n + k];
        }
        if (!(d > 0.0)) {
            throw std::domain_error("ca_pfa: correlation matrix is not positive definite");
        }
        const double ljj = std::sqrt(d);
        l[j * n + j] = ljj;
        log_det += 2.0 * std::log(ljj);
        for (std::size_t i = j + 1; i < n; ++i) {
            double s = c(i, j);
            for (std::size_t k = 0; k < j; ++k) {
                s -= l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = s / ljj;
        }
    }
    return log_det;
}

}  // namespace

NoiseEstimate ca_noise(std::span<const double> power, BinRange band, std::size_t guard,
                       std::size_t train) {
    if (band.last >= power.size() || band.first > band.last) {
        throw std::invalid_argument("ca_noise: band outside the spectrum");
    }
    if (train == 0) {
        throw std::invalid_argument("ca_noise: need at least one training cell per side");
    }

    // prefix[i] = sum of power[0 .. i-1], so sum over [a, b] = prefix[b + 1] - prefix[a].
    std::vector<double> prefix(power.size() + 1, 0.0);
    for (std::size_t i = 0; i < power.size(); ++i) {
        prefix[i + 1] = prefix[i] + power[i];
    }

    NoiseEstimate est{std::vector<double>(power.size(), kNaN),
                      std::vector<Sides>(power.size(), Sides::none)};
    const std::size_t reach = guard + train;
    for (std::size_t k = band.first; k <= band.last; ++k) {
        const bool lead = k >= band.first + reach;
        const bool lag = k + reach <= band.last;
        double sum = 0.0;
        std::size_t cells = 0;
        if (lead) {
            sum += prefix[k - guard] - prefix[k - reach];
            cells += train;
        }
        if (lag) {
            sum += prefix[k + reach + 1] - prefix[k + guard + 1];
            cells += train;
        }
        if (cells == 0) {
            throw std::invalid_argument("ca_noise: band too narrow for guard + training cells");
        }
        est.level[k] = sum / static_cast<double>(cells);
        est.sides[k] = lead && lag ? Sides::both : (lead ? Sides::lead_only : Sides::lag_only);
    }
    return est;
}

std::vector<double> bin_correlation(std::span<const double> window, std::size_t max_lag) {
    const double n = static_cast<double>(window.size());
    double sum_sq = 0.0;
    for (double w : window) {
        sum_sq += w * w;
    }
    std::vector<double> rho(max_lag + 1);
    for (std::size_t m = 0; m <= max_lag; ++m) {
        double acc = 0.0;
        for (std::size_t i = 0; i < window.size(); ++i) {
            const double phase = 2.0 * std::numbers::pi * static_cast<double>(m) *
                                 static_cast<double>(i) / n;
            acc += window[i] * window[i] * std::cos(phase);
        }
        rho[m] = acc / sum_sq;
        // Clean up rounding so lags a window does not couple are exactly zero.
        if (std::abs(rho[m]) < 1e-12) {
            rho[m] = 0.0;
        }
    }
    return rho;
}

double ca_pfa(double alpha, std::size_t n_lead, std::size_t n_lag, std::span<const double> rho) {
    const std::size_t n = n_lead + n_lag;
    if (n == 0) {
        throw std::invalid_argument("ca_pfa: no training cells");
    }
    const double t = alpha / static_cast<double>(n);
    return std::exp(-(log_det_i_plus_tc(n_lead, t, rho) + log_det_i_plus_tc(n_lag, t, rho)));
}

double ca_alpha(double pfa, std::size_t n_lead, std::size_t n_lag, std::span<const double> rho) {
    if (!(pfa > 0.0 && pfa < 1.0)) {
        throw std::invalid_argument("ca_alpha: pfa must be in (0, 1)");
    }
    // ca_pfa falls monotonically from 1 at alpha = 0. Bracket, then bisect.
    double lo = 0.0;
    double hi = 1.0;
    while (ca_pfa(hi, n_lead, n_lag, rho) > pfa) {
        lo = hi;
        hi *= 2.0;
    }
    for (int i = 0; i < 200 && hi - lo > 1e-12 * hi; ++i) {
        const double mid = 0.5 * (lo + hi);
        (ca_pfa(mid, n_lead, n_lag, rho) > pfa ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

double pfa_bin_for_window(double pfa_window, std::size_t n_bins) {
    if (!(pfa_window > 0.0 && pfa_window < 1.0) || n_bins == 0) {
        throw std::invalid_argument("pfa_bin_for_window: need 0 < pfa_window < 1 and n_bins > 0");
    }
    return -std::expm1(std::log1p(-pfa_window) / static_cast<double>(n_bins));
}

CaCfar::CaCfar(const CaParams& params, std::span<const double> window) : params_(params) {
    const std::size_t reach = params.guard_per_side + params.train_per_side;
    const auto rho = bin_correlation(window, reach);
    // The closed form needs the test cell independent of every training cell.
    for (std::size_t m = params.guard_per_side + 1; m <= reach; ++m) {
        if (rho[m] != 0.0) {
            throw std::invalid_argument("CaCfar: guard too small for this window's bin correlation");
        }
    }
    alpha_both_ = ca_alpha(params.pfa_bin, params.train_per_side, params.train_per_side, rho);
    alpha_one_ = ca_alpha(params.pfa_bin, params.train_per_side, 0, rho);
}

double CaCfar::alpha(Sides sides) const {
    switch (sides) {
        case Sides::both:
            return alpha_both_;
        case Sides::lead_only:
        case Sides::lag_only:
            return alpha_one_;
        case Sides::none:
            break;
    }
    return kNaN;
}

std::vector<double> CaCfar::threshold(const NoiseEstimate& noise) const {
    std::vector<double> t(noise.level.size(), kNaN);
    for (std::size_t k = 0; k < t.size(); ++k) {
        if (noise.sides[k] != Sides::none) {
            t[k] = alpha(noise.sides[k]) * noise.level[k];
        }
    }
    return t;
}

std::vector<double> CaCfar::threshold(std::span<const double> power, BinRange band) const {
    return threshold(ca_noise(power, band, params_.guard_per_side, params_.train_per_side));
}

}  // namespace cfar
