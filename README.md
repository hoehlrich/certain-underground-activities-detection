# CFAR detection sim

CA-CFAR detector for a low-frequency magnetic tone in HFR-1 rod recordings from
an Analog Discovery 3. Part of the 8.9 Quantum Sensors senior design capstone
(Colorado School of Mines, NNSS sponsor).

The detection code in `core/` has no dependency on `io/` or `sim/`; a CTest
check enforces it.

## Build and test

Needs CMake 3.24+, a C++20 compiler and GoogleTest. pocketfft and toml++ are
downloaded at the first configure.

```
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Run

```
./build/sim_main configs/baseline.toml                  # writes results/example_spectrum.csv
python3 scripts/plot.py results/example_spectrum.csv    # writes results/example_spectrum.png

./build/pfa_check configs/baseline.toml                 # writes results/pfa_check.csv
./build/pd_sweep configs/baseline.toml                  # writes results/pd_vs_snr.csv (about 4 s)
python3 scripts/plot.py results/pd_vs_snr.csv           # writes results/pd_vs_snr.png

./build/roc configs/baseline.toml                       # writes results/roc.csv
python3 scripts/plot.py results/roc.csv                 # writes results/roc.png
```

`sim_main` runs one window of a tone in white noise through the detector and
lists the detections.

`pfa_check` runs noise-only windows and compares the measured false-alarm rate
with the design value.

`pd_sweep` measures probability of detection against SNR, compares it with the
known-noise bound, and prints the C/N0 needed for Pd = 0.9.

`roc` holds SNR fixed and sweeps the false-alarm target. Parameters are in
`configs/baseline.toml`

## Layout

| Path | Contents |
| --- | --- |
| `core/` | Detection core (target `core`): spectrum, CA-CFAR, detector, config structs and validation. No I/O. |
| `io/` | Config loading (target `io`). Depends on `core` only. |
| `sim/` | Tone and noise generators, simulation config, Marcum Q bound, Pfa, Pd and ROC Monte Carlo (target `sim`). |
| `apps/` | `sim_main`, `pfa_check`, `pd_sweep`, `roc` |
| `tests/` | GoogleTest suite |
| `cmake/` | Include-boundary check run by CTest |
| `configs/` | TOML configs |
| `scripts/` | Plotting |
| `docs/img/` | Pd vs SNR and ROC figures |
