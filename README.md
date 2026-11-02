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
./build/pd_sweep configs/baseline.toml          # writes results/pd_vs_snr.csv (about 4 s)
python3 scripts/plot.py results/pd_vs_snr.csv   # writes results/pd_vs_snr.png
./build/roc configs/baseline.toml               # writes results/roc.csv
python3 scripts/plot.py results/roc.csv         # writes results/roc.png
```

`sim_main` runs one window of a tone in white noise through the detector and lists the detections. `pfa_check` runs noise-only windows and compares the measured false-alarm rate with the design value; it exits non-zero if any design value fails. `pd_sweep` measures probability of detection against SNR, compares it with the known-noise bound, and prints the C/N0 needed for Pd = 0.9. `roc` holds SNR fixed and sweeps the false-alarm target. Parameters are in `configs/baseline.toml`; none are confirmed by the sponsor yet.

## Layout

| Path | Contents |
| --- | --- |
| `include/cfar/`, `src/` | Detector library (target `cfar`): spectrum, CA-CFAR, per-window detector, processing config. No dependency on the simulation code. |
| `include/cfar/sim/`, `src/sim/` | Signal, noise, simulation config, Pfa and Pd Monte Carlo, detection theory (target `cfar_sim`). |
| `apps/` | `sim_main`, `pfa_check`, `pd_sweep`, `roc` |
| `tests/` | GoogleTest suite |
| `configs/` | TOML configs |
| `scripts/` | Plotting |
| `docs/cfar-sim-plan.md` | Signal model, parameters, CFAR math, validation plan |
| `docs/pd-results.md` | Pd vs SNR and ROC results, and the C/N0 needed for Pd = 0.9 |

## Status

| Milestone | State |
| --- | --- |
| 3.2.1 Tone, white noise, Hann periodogram | Done |
| 3.2.2 CA-CFAR and Pfa check | Done |
| 3.2.3 Pd vs SNR | Done: Pd = 0.9 at C/N0 = 13.5 dB-Hz in white noise ([results](docs/pd-results.md)) |
