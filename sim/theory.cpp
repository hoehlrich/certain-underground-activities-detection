#include "sim/theory.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace cfar::sim {

double marcum_q1(double a, double b) {
    if (b <= 0.0) {
        return 1.0;
    }
    // |a + n|^2 is noncentral chi-square with 2 degrees of freedom, a Poisson
    // (mean a^2/2) mixture of chi-square(2 + 2j). Its tail above b^2 is
    //   sum_j Pois(j; a^2/2) * P(Pois(b^2/2) <= j),
    // a sum of positive terms that stays accurate for small and large Q1.
    const double lam = 0.5 * a * a;
    const double y = 0.5 * b * b;
    if (lam == 0.0) {
        return std::exp(-y);
    }
    const auto j_max = static_cast<std::size_t>(lam + 12.0 * std::sqrt(lam) + 40.0);
    const double log_lam = std::log(lam);
    const double log_y = std::log(y);
    double cdf_y = 0.0;        // P(Pois(y) <= j)
    double log_term_y = -y;    // log P(Pois(y) = j)
    double sum = 0.0;
    for (std::size_t j = 0; j <= j_max; ++j) {
        const double dj = static_cast<double>(j);
        cdf_y = std::min(1.0, cdf_y + std::exp(log_term_y));
        sum += std::exp(-lam + dj * log_lam - std::lgamma(dj + 1.0)) * cdf_y;
        log_term_y += log_y - std::log(dj + 1.0);
    }
    return std::min(sum, 1.0);
}

double scallop_gain(std::span<const double> window, double offset_bins) {
    const double n = static_cast<double>(window.size());
    std::complex<double> w_off{0.0, 0.0};
    double w_sum = 0.0;
    for (std::size_t i = 0; i < window.size(); ++i) {
        const double phase = -2.0 * std::numbers::pi * offset_bins * static_cast<double>(i) / n;
        w_off += window[i] * std::polar(1.0, phase);
        w_sum += window[i];
    }
    return std::norm(w_off) / (w_sum * w_sum);
}

double ideal_pd(double snr_bin, double pfa) {
    return marcum_q1(std::sqrt(2.0 * snr_bin), std::sqrt(-2.0 * std::log(pfa)));
}

IdealBound::IdealBound(std::span<const double> window, double pfa_bin,
                       std::span<const double> offsets_bins)
    : pfa_(pfa_bin) {
    if (offsets_bins.empty()) {
        throw std::invalid_argument("IdealBound: no tone offsets");
    }
    for (double off : offsets_bins) {
        gains_.push_back(scallop_gain(window, off));
    }
}

double IdealBound::pd(double snr_bin_db) const {
    const double snr = std::pow(10.0, snr_bin_db / 10.0);
    double sum = 0.0;
    for (double g : gains_) {
        sum += ideal_pd(snr * g, pfa_);
    }
    return sum / static_cast<double>(gains_.size());
}

double IdealBound::snr_bin_db_for(double pd_target) const {
    double lo = -20.0;
    double hi = 60.0;
    for (int i = 0; i < 60; ++i) {
        const double mid = 0.5 * (lo + hi);
        (pd(mid) >= pd_target ? hi : lo) = mid;
    }
    return hi;
}

std::vector<double> uniform_offsets(std::size_t n) {
    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = 0.5 * (static_cast<double>(i) + 0.5) / static_cast<double>(n);
    }
    return out;
}

}  // namespace cfar::sim
