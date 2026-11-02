#pragma once

#include <cstdint>
#include <random>
#include <span>

namespace cfar::sim {

// Per-trial generator (SYS-09): the same (seed, trial) always gives the same
// stream, independent of how trials are ordered or split across threads.
// std::normal_distribution is implementation-defined, so runs reproduce exactly
// on the same toolchain, not necessarily across compilers.
std::mt19937_64 make_trial_rng(std::uint64_t seed, std::uint64_t trial);

// Adds zero-mean white Gaussian noise of standard deviation sigma to x in place.
void add_white_noise(std::span<double> x, double sigma, std::mt19937_64& rng);

}  // namespace cfar::sim
