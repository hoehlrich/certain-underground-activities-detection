# CFAR detection sim plan

Technical reference for the CA-CFAR detector and its validation simulation. Scope for now is CA-CFAR in white Gaussian noise; other variants and noise scenarios are listed under "Deferred" below. The team-facing version is the "CFAR Detection Sim Plan" doc; this file tracks it for what the repo builds.

## Scope update (sponsor slides, 2026-10-29)

- Quantum magnetometers are out for now (export control). HFR-1 rods are the receive sensors; NNSS proposes lending two with amplifiers.
- The project goal shifts to building and characterizing a tunable source, comparing propagation models with measurements on traditional sensors, and working out the sensitivity a sensor would need. Testing runs lab, then outdoors on campus, then the mine.
- The source options include a VFD reaching 1000 Hz, and the MagPi's 500 Hz ceiling went with the quantum sensors, so the default search band widens to 10 to 1000 Hz.
- For this sim, the C/N0 needed for Pd = 0.9 becomes a headline output, since it converts directly into a required sensor noise floor.

## Detection approach (from the team's Sprint 3 Concept Definition)

- FFT over a window shorter than the signal's coherence time, then 1D CFAR across frequency bins.
- Tx and Rx share no frequency reference (free-running crystals on both ends), so the detector searches the whole band rather than one known bin.
- SNR is peak bin power over mean non-signal power.

## Default parameters

None of these are confirmed by the sponsor. All are config values.

| Parameter | Default | Reason | Confirm with |
| --- | --- | --- | --- |
| Search band | 10 to 1000 Hz | VFD source reaches 1000 Hz; stays above DC and the worst 1/f | NNSS (SYS-01) |
| Sample rate | 4 kHz after decimation | At least 2.5 x 1000 Hz (SUB-03.2), leaving 1 to 2 kHz for anti-alias roll-off | AD3 config |
| FFT length | 8192 samples (2.05 s, 0.49 Hz bins, 2028 bins in band) | Later swept over 1, 2, 4, 8 s windows | Sim results |
| Window | Hann | Keeps 60 Hz harmonic leakage out of nearby bins; ENBW 1.5 bins, scalloping up to 1.42 dB | Sim results |
| Pfa per window | 0.01 | SUB-05.1 draft | NNSS |
| Pfa per bin | 5e-6, derived | 2028 bins searched: 1 - (1 - 5e-6)^2028 is about 0.01 (was 1e-5 over 1004 bins for the 500 Hz band) | Derived |
| Tone frequency | Unknown, uniform over the band; drift under 0.05 Hz per window | Free-running oscillators | Bench data |
| Rod response | Flat after antenna-factor correction | HFR-1 datasheet not yet in hand; gain and band of the loaner amplifiers unknown | NNSS (HFR-1 and amp specs) |
| Guard cells per side | 2 (sweep 1 to 4 later) | Hann main lobe spans +/-2 bins | Sim results |
| Training cells per side | 16 (sweep 8, 16, 32 later) | At 0.49 Hz bins, 16 per side covers +/-8 Hz | Sim results |

## False-alarm budget

Pfa per bin is derived from the per-window target and the number of bins searched, not fixed, so a new f2 or window length is a config change. With B = f2 - f1 and window T, N_bins is about B x T:

```
Pfa_bin = 1 - (1 - Pfa_win)^(1 / N_bins)  ~  Pfa_win / N_bins
```

| Window T (FFT length) | Bin width | Bins in 10 to 1000 Hz | Pfa per bin | Ideal threshold over noise mean |
| --- | --- | --- | --- | --- |
| 1.02 s (4096) | 0.98 Hz | 1014 | 9.9e-6 | 10.6 dB |
| 2.05 s (8192, baseline) | 0.49 Hz | 2028 | 5.0e-6 | 10.9 dB |
| 4.10 s (16384) | 0.24 Hz | 4056 | 2.5e-6 | 11.1 dB |
| 8.19 s (32768) | 0.12 Hz | 8111 | 1.2e-6 | 11.3 dB |

The threshold column is -ln(Pfa_bin) in dB, which CA approaches as training grows. Hann bins overlap, so the independent-bin formula slightly overstates window Pfa; the Pfa check measures the real value.

## Signal model

```
s(t) = A cos(2 pi (f0 t + 0.5 drift t^2) + phi)
```

- f0 drawn uniformly over the band each trial, so off-bin scalloping loss is included. A fixed f0 can be set in the config.
- phi uniform on [0, 2 pi).
- Optional linear frequency drift.
- A is set from the target C/N0, so path loss stays in the team's propagation model and the sim doesn't depend on it.

Report SNR as C/N0 in dB-Hz (tone power over noise density) so results compare across window lengths. With C = A^2 / 2 and one-sided N0 = 2 sigma^2 / fs for real white noise of variance sigma^2, convert to per-bin SNR for theory overlays:

```
SNR_bin [dB] = C/N0 [dB-Hz] + 10 log10(T) - 10 log10(ENBW_bins) - L_scallop
```

where T is the window length in seconds and L_scallop is the scalloping loss for that trial's f0.

## Test statistic

Per bin, the one-sided Hann periodogram (periodic Hann, n/2 + 1 bins):

```
P[k] = |X[k]|^2 / sum(w^2)
```

where X is the DFT of x * w. With this scaling, noise-only bins have mean sigma^2 and a bin-centred tone peaks at (A^2 / 4) n / ENBW. With white Gaussian noise and no Welch averaging, P[k] in noise-only bins is exponentially distributed, which gives CA-CFAR a closed-form threshold.

## CA-CFAR

For bin k, take N training cells total (N/2 leading, N/2 lagging), skipping G guard cells on each side of k.

```
noise_k   = (1/N) * sum over training cells of P[i]
alpha_CA  = N * (Pfa^(-1/N) - 1)
T_k       = alpha_CA * noise_k
detect    = P[k] > T_k
```

- Compute the training sums with a running sum so each bin costs O(1).
- Compute alpha once per config.
- At band edges, use one-sided training (the available side only) and rescale alpha for the reduced N.

Detector logic per window:

1. Threshold every bin in the search band.
2. Cluster adjacent hits and keep each cluster's peak bin.
3. Refine frequency by parabolic interpolation on log power around the peak.
4. Estimate SNR as (peak - noise_k) / noise_k.
5. Output per window (SYS-05): decision, detection statistic P[k] / noise_k, estimated SNR, refined frequency.

## Ideal bound

The known-noise detector for a steady tone, which every CFAR curve is compared against:

```
Pd = Q1( sqrt(2 * SNR_bin), sqrt(-2 ln Pfa) )
```

Q1 is the Marcum Q-function of order 1 (SNR_bin linear, not dB). The horizontal gap between this curve and the CA curve at Pd = 0.9 is the CFAR loss.

## Validation

1. **Pfa check:** noise only, at least 1e7 bin trials per config (about 1e4 windows). Pass when measured Pfa falls inside the 95% binomial interval of the design value, for design Pfa of 1e-3, 1e-4 and 1e-5, plus the derived baseline 5e-6.
2. **Pd vs SNR:** per-bin SNR from 0 to 25 dB in 1 dB steps, 2000 trials each, ideal curve overlaid. Report the C/N0 needed for Pd = 0.9.
3. **ROC (optional):** at three SNRs near the knee of the Pd curve, sweep Pfa from 1e-6 to 1e-1.

## Status

- **3.2.1 framework: done 2026-10-29.** Tone and white noise, Hann periodogram, TOML config, example spectrum CSV and plot. `tests/test_spectrum.cpp` checks tone placement, measured ENBW of 1.5 bins, 1.42 dB half-bin scalloping, exponential noise bins with mean sigma^2, the C/N0 to SNR_bin conversion, per-trial reproducibility and config validation.
- **3.2.2 CA-CFAR:** next.
- **3.2.3 Pd vs SNR:** after that.

## Repo layout

```
CMakeLists.txt     targets: cfar (detector library), cfar_sim (simulation), sim_main, test_spectrum
include/cfar/      spectrum.hpp; later cfar.hpp detector.hpp metrics.hpp
include/cfar/sim/  signal.hpp noise.hpp config.hpp
src/, src/sim/     matching .cpp files
apps/              sim_main.cpp (one window to CSV now, sweeps later), run_file.cpp (AD3 captures, later)
tests/             test_spectrum.cpp; test_cfar.cpp in 3.2.2
configs/           baseline.toml
scripts/           plot.py (spectrum now; Pd vs SNR and ROC later)
docs/              this file
results/           CSV and PNG output, git-ignored
```

`cfar` must not depend on `cfar_sim`; it becomes the delivered pipeline for real AD3 data (D-09). `cfar.hpp` will expose a function that fills a threshold array from P[k] and a params struct, written so other noise estimators can be added later behind the same signature.

Dependencies: pocketfft and toml++ come in through CMake FetchContent (pinned); GoogleTest is the system package.

## Config

`configs/baseline.toml` follows the section names in the shared plan doc. Built so far:

| Section | Keys |
| --- | --- |
| `[band]` | `f_lo_hz`, `f_hi_hz` |
| `[acq]` | `fs_hz`, `nfft`, `window` (`"hann"` only) |
| `[noise]` | `sigma` |
| `[source]` | `kind` (`"tone"` only), `f0_hz` (`"uniform"` or a value in Hz), `cn0_dbhz` |
| `[run]` | `seed` |

To add in 3.2.2: `[false_alarm]` with `pfa_window = 0.01` and `pfa_bin = "derived"`, and `[cfar]` with `guard_per_side` and `train_per_side`. The shared doc's `variant = "os"`, `os_rank_frac`, `welch_segments`, `[interference]` and `source.kind = "vfd"` wait for the deferred work below.

| Output file | Contents |
| --- | --- |
| example_spectrum.csv | One window's P[k] (threshold and detections added in 3.2.2), for a slide |
| pfa_check.csv | Config, design Pfa, measured Pfa, 95% interval |
| pd_vs_snr.csv | Scenario, SNR, Pd, trials |
| roc.csv | SNR, Pfa, Pd |

## Deferred until real data (do not build now)

- **OS-CFAR:** noise estimate is the k-th smallest training cell (k about 3/4 of N). Robust to 60 Hz harmonics and nearby tones in the training cells. Alpha solved by bisection from `Pfa = prod_{i=0}^{k-1} (N - i) / (N - i + alpha)`. The shared doc names OS as the likely fielded choice.
- **GO-CFAR:** larger of the leading and lagging means. Holds Pfa on sloped noise floors (1/f) and edges.
- **SO-CFAR:** smaller of the two means. Resists a one-sided interferer; Pfa rises at edges. Comparison only.
- **Noise scenarios:** colored 1/f^a background; 60 Hz harmonics to 960 Hz (16 lines in the band) at 20 to 40 dB over the floor with slight wander; impulsive Gaussian mixture (about 1% of samples at 10 to 30 x variance); all combined.
- **VFD source (Option B):** harmonics and slip sidebands; the fundamental can sit anywhere up to 1000 Hz, and harmonics above the band are dropped.
- **Known-tone gate (Option A):** an optional +/-Delta Hz gate around the expected frequency, which shrinks the search and the per-window Pfa.
- **Welch averaging:** M = 1, 4, 16 segments at 50% overlap.
- **Real data:** noise-only rod recordings through the AD3 with synthetic tones injected (TP-10), via `apps/run_file.cpp`.

## Open questions for the sponsor

- Which band and tone frequencies to plan around (SYS-01 f1 and f2)? 10 to 1000 Hz is the working default.
- What Pd, per-window Pfa, window length and SNR define "detectable" (SUB-05.1)?
- HFR-1 usable band, antenna factor, and the gain and band of the loaner amplifiers?
- Any noise-only or source-on recordings from their amp, loop and AD3 setup to inject into?
