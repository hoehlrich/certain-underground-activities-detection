#pragma once

#include <span>

namespace cfar::sim {

// s(t) = A cos(2 pi (f0 t + 0.5 drift t^2) + phi), with t = k / fs.
struct ToneParams {
    double amplitude = 1.0;
    double freq_hz = 0.0;
    double phase_rad = 0.0;
    double drift_hz_per_s = 0.0;
};

// Adds the tone to x in place.
void add_tone(std::span<double> x, double fs, const ToneParams& tone);

// Amplitude giving the requested C/N0 [dB-Hz] in real white noise of standard
// deviation sigma sampled at fs. Uses C = A^2 / 2 and one-sided N0 = 2 sigma^2 / fs.
double amplitude_for_cn0(double cn0_dbhz, double sigma, double fs);

}  // namespace cfar::sim
