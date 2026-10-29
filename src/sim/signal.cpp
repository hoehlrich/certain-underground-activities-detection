#include "cfar/sim/signal.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>

namespace cfar::sim {

void add_tone(std::span<double> x, double fs, const ToneParams& tone) {
    for (std::size_t k = 0; k < x.size(); ++k) {
        const double t = static_cast<double>(k) / fs;
        const double cycles = tone.freq_hz * t + 0.5 * tone.drift_hz_per_s * t * t;
        x[k] += tone.amplitude * std::cos(2.0 * std::numbers::pi * cycles + tone.phase_rad);
    }
}

double amplitude_for_cn0(double cn0_dbhz, double sigma, double fs) {
    const double n0 = 2.0 * sigma * sigma / fs;
    const double carrier = n0 * std::pow(10.0, cn0_dbhz / 10.0);
    return std::sqrt(2.0 * carrier);
}

}  // namespace cfar::sim
