#pragma once

#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

namespace cfar {

// Inclusive range of FFT bins [first, last].
struct BinRange {
    std::size_t first = 0;
    std::size_t last = 0;
    std::size_t size() const { return last - first + 1; }
};

// Periodic (DFT-even) Hann window of length n: w[k] = 0.5 - 0.5 cos(2 pi k / n).
// ENBW is exactly 1.5 bins in this form.
std::vector<double> hann_window(std::size_t n);

// Window by config name. Only "hann" for now; throws std::invalid_argument otherwise.
std::vector<double> make_window(std::string_view name, std::size_t n);

// Bins whose centre frequency lies inside [f_lo_hz, f_hi_hz] for an n-point FFT
// at fs. Throws std::invalid_argument if the band holds no bins or passes Nyquist.
BinRange band_bins(double f_lo_hz, double f_hi_hz, std::size_t n, double fs);

// Equivalent noise bandwidth in bins: n * sum(w^2) / sum(w)^2.
double enbw_bins(std::span<const double> window);

// Frequency spacing of an n-point FFT at sample rate fs.
double bin_hz(std::size_t n, double fs);

// One-sided windowed periodogram, n/2 + 1 bins (DC to Nyquist).
//
// P[k] = |X[k]|^2 / sum(w^2), where X is the DFT of x * w. With this scaling,
// noise-only bins of white noise with variance sigma^2 have mean sigma^2, and a
// bin-centred tone A cos(...) peaks at (A^2 / 4) * n / ENBW.
std::vector<double> periodogram(std::span<const double> x, std::span<const double> window);

}  // namespace cfar
