# CFAR detection sim

CA-CFAR detector for a low-frequency magnetic tone in HFR-1 rod recordings from an Analog Discovery 3, with a Monte Carlo simulation to validate it. Part of the 8.9 Quantum Sensors senior design capstone (Colorado School of Mines, NNSS sponsor).

The detector library processes FFT power spectra and has no dependency on the simulation code. Until rod recordings are available, it is tested on simulated tones in noise. The tool for running it on recorded files is not written yet.

## Build and test

Needs CMake 3.24+, a C++20 compiler and GoogleTest. pocketfft and toml++ are downloaded at the first configure.

```
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Run

```
./build/sim_main configs/baseline.toml          # writes results/example_spectrum.csv
python3 scripts/plot.py results/example_spectrum.csv   # writes results/example_spectrum.png
./build/pfa_check configs/baseline.toml         # writes results/pfa_check.csv
```

`sim_main` runs one window of a tone in white noise through CA-CFAR and reports which bins crossed the threshold. `pfa_check` runs noise-only windows and compares the measured false-alarm rate with the design value; it exits non-zero if any design value fails. Parameters are in `configs/baseline.toml`; none are confirmed by the sponsor yet.

## Layout

| Path | Contents |
| --- | --- |
| `include/cfar/`, `src/` | Detector library (target `cfar`): spectrum, CA-CFAR, processing config. No dependency on the simulation code. |
| `include/cfar/sim/`, `src/sim/` | Signal, noise, simulation config, Pfa Monte Carlo (target `cfar_sim`). |
| `apps/` | `sim_main`, `pfa_check` |
| `tests/` | GoogleTest suite |
| `configs/` | TOML configs |
| `scripts/` | Plotting |
| `docs/cfar-sim-plan.md` | Signal model, parameters, CFAR math, validation plan |

## Status

| Milestone | State |
| --- | --- |
| 3.2.1 Tone, white noise, Hann periodogram | Done |
| 3.2.2 CA-CFAR and Pfa check | Done |
| 3.2.3 Pd vs SNR | Next |
