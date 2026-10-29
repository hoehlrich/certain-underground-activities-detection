#include "cfar/spectrum.hpp"

#include <cmath>
#include <complex>
#include <numbers>
#include <numeric>
#include <stdexcept>

#include "pocketfft_hdronly.h"

namespace cfar {

std::vector<double> hann_window(std::size_t n) {
    std::vector<double> w(n);
    for (std::size_t k = 0; k < n; ++k) {
        w[k] = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * static_cast<double>(k) /
                                    static_cast<double>(n));
    }
    return w;
}

double enbw_bins(std::span<const double> window) {
    const double sum = std::accumulate(window.begin(), window.end(), 0.0);
    const double sum_sq = std::inner_product(window.begin(), window.end(), window.begin(), 0.0);
    return static_cast<double>(window.size()) * sum_sq / (sum * sum);
}

double bin_hz(std::size_t n, double fs) { return fs / static_cast<double>(n); }

std::vector<double> periodogram(std::span<const double> x, std::span<const double> window) {
    if (x.size() != window.size()) {
        throw std::invalid_argument("periodogram: signal and window lengths differ");
    }
    const std::size_t n = x.size();

    std::vector<double> xw(n);
    double sum_sq = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        xw[k] = x[k] * window[k];
        sum_sq += window[k] * window[k];
    }

    std::vector<std::complex<double>> spec(n / 2 + 1);
    pocketfft::r2c<double>({n}, {sizeof(double)}, {sizeof(std::complex<double>)}, 0, true,
                           xw.data(), spec.data(), 1.0);

    std::vector<double> power(spec.size());
    for (std::size_t k = 0; k < spec.size(); ++k) {
        power[k] = std::norm(spec[k]) / sum_sq;
    }
    return power;
}

}  // namespace cfar
