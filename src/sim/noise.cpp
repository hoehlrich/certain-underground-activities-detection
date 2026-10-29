#include "cfar/sim/noise.hpp"

namespace cfar::sim {

std::mt19937_64 make_trial_rng(std::uint64_t seed, std::uint64_t trial) {
    std::seed_seq seq{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32),
                      static_cast<std::uint32_t>(trial), static_cast<std::uint32_t>(trial >> 32)};
    return std::mt19937_64(seq);
}

void add_white_noise(std::span<double> x, double sigma, std::mt19937_64& rng) {
    std::normal_distribution<double> dist(0.0, sigma);
    for (double& v : x) {
        v += dist(rng);
    }
}

}  // namespace cfar::sim
