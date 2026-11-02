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
| Guard cells per side | 2 (sweep 2 to 4 later) | Hann main lobe spans +/-2 bins, and Hann bins up to 2 apart are correlated, so 2 is the minimum for the exact alpha | Sim results |
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
T_k       = alpha * noise_k
detect    = P[k] > T_k
```

The textbook alpha, `N * (Pfa^(-1/N) - 1)`, assumes independent training cells. Hann-windowed bins are not independent: the complex correlation between bins m apart is

```
rho[m] = sum(w^2 cos(2 pi m n / N)) / sum(w^2)      Hann: 1, -2/3, 1/6, 0, 0, ...
```

Correlated training cells average away less noise, so the textbook alpha lets too many false alarms through. With 16 cells per side, Monte Carlo gave 1.3x the design Pfa at 1e-2 and 1.7x at 1e-3, and the exact formula below predicts 4x at the 5e-6 baseline. The exact Pfa, for white Gaussian noise and a test cell independent of its training cells, is

```
Pfa = 1 / ( det(I + (alpha/N) C_lead) * det(I + (alpha/N) C_lag) )
```

where C_side is the Toeplitz matrix of rho over one side's cells (the two sides are far enough apart to be independent). With rho = {1} this reduces to the textbook (1 + alpha/N)^-N. Alpha is solved from it by bisection, once per config, using a Cholesky factorization for the determinants. At the baseline (16 + 16 cells, 4.96e-6), alpha is 17.46 instead of the textbook 14.86, a threshold 0.7 dB higher.

- The test cell must be independent of every training cell, so G must be at least the largest lag with nonzero rho (2 for Hann). `CaCfar` throws if it isn't.
- Training sums come from prefix sums, so each bin costs O(1).
- Training cells stay inside the search band. A bin within G + N/2 of a band edge uses the side that fits, with alpha recomputed for N/2 cells (24.85 at baseline).
- The formula assumes circular complex Gaussian bins, which holds away from DC and Nyquist. The band stays clear of both.

Detector logic per window (`cfar::Detector` in `detector.hpp`):

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

Q1 is the Marcum Q-function of order 1 (SNR_bin linear, not dB). With f0 uniform over the band, the sim averages this over tone offsets from the bin centre, applying each offset's scalloping loss, so the bound and the CA trials see the same signals. The horizontal gap between this curve and the CA curve at Pd = 0.9 is the CFAR loss.

A CA trial counts as detected when a reported peak lies within one bin of f0, while the bound tests only the nearest bin. This favours CA slightly, so the measured loss is if anything a little low.

## Validation

1. **Pfa check:** noise only, at least 1e7 bin trials per config (about 1e4 windows). Pass when measured Pfa falls inside the 95% binomial interval of the design value, for design Pfa of 1e-3, 1e-4 and 1e-5, plus the derived baseline 5e-6. The interval uses the normal approximation, and the expected count is at least 100 hits at the lowest design value. Adjacent test cells are correlated, but a false alarm rarely spills into the next bin, so the binomial interval is close.
2. **Pd vs SNR:** per-bin SNR (bin-centred) from 0 to 25 dB in 1 dB steps, 2000 trials each, ideal curve overlaid. Report the C/N0 needed for Pd = 0.9. Each trial reuses its noise, f0 and phase at every SNR, so the curve is smooth.
3. **ROC:** at three SNRs near the knee of the Pd curve (12, 14, 16 dB per bin), sweep the per-window Pfa from 1e-4 to 0.5 (5e-8 to 3.4e-4 per bin), using the same trials as the Pd sweep.

## Status

- **3.2.1 framework: done 2026-10-29.** Tone and white noise, Hann periodogram, TOML config, example spectrum CSV and plot. `tests/test_spectrum.cpp` checks tone placement, measured ENBW of 1.5 bins, 1.42 dB half-bin scalloping, exponential noise bins with mean sigma^2, the C/N0 to SNR_bin conversion, per-trial reproducibility and config validation.
- **3.2.2 CA-CFAR: done 2026-10-29.** Prefix-sum noise estimate, correlation-corrected alpha, derived per-bin Pfa, one-sided training at band edges. `tests/test_cfar.cpp` checks the noise estimate against brute force, the textbook limit of alpha, the Hann correlation values, and the Pfa check below. Baseline config, 10,000 noise-only windows, 2.03e7 bin trials:

  | Design Pfa | Measured | 95% interval | Windows with a false alarm |
  | --- | --- | --- | --- |
  | 1e-3 | 9.94e-4 | [9.86e-4, 1.01e-3] | 82% |
  | 1e-4 | 9.78e-5 | [9.57e-5, 1.04e-4] | 17% |
  | 1e-5 | 9.62e-6 | [8.62e-6, 1.14e-5] | 1.8% |
  | 4.96e-6 (derived) | 4.88e-6 | [3.99e-6, 5.93e-6] | 0.96% |

  The derived value holds the per-window target of 1%.
- **3.2.3 Pd vs SNR: done 2026-11-02.** Per-window detector (clustering, parabolic frequency refinement, SNR estimate), Marcum Q1 bound averaged over scalloping, and the sweep. `tests/test_pd.cpp` checks Q1 against scipy, the bound against a known-noise Monte Carlo, clustering and interpolation, and that CA sits within 2 dB below the bound. At baseline, Pd = 0.9 needs **C/N0 = 13.5 dB-Hz** for CA-CFAR and 11.9 dB-Hz for the known-noise bound, a CFAR loss of 1.6 dB. Window false alarms in the sweep are 1.2% (24 of 2000). ROC at 12, 14 and 16 dB per-bin SNR: near Pd = 0.9, each decade of per-window Pfa costs about 1 dB of C/N0. Write-up in `docs/pd-results.md`.

## Repo layout

```
CMakeLists.txt     targets: cfar (detector library), cfar_sim (simulation), apps, tests
include/cfar/      spectrum.hpp cfar.hpp config.hpp detector.hpp
include/cfar/sim/  signal.hpp noise.hpp config.hpp pfa.hpp theory.hpp pd.hpp
src/, src/sim/     matching .cpp files
apps/              sim_main.cpp (one window to CSV), pfa_check.cpp (Pfa check to CSV),
                   pd_sweep.cpp (Pd vs SNR to CSV), roc.cpp (ROC to CSV),
                   run_file.cpp (AD3 captures, later)
tests/             test_spectrum.cpp test_cfar.cpp test_pd.cpp
configs/           baseline.toml
scripts/           plot.py (spectrum, Pd vs SNR, ROC)
docs/              this file, pd-results.md (3.2.3 write-up), img/
results/           CSV and PNG output, git-ignored
```

`cfar` must not depend on `cfar_sim`; it becomes the delivered pipeline for real AD3 data (D-09). `cfar::Detector` is built from `ProcessingConfig` alone and turns nfft samples into a list of detections (SYS-05 outputs). In `cfar.hpp`, `ca_noise` maps P[k] and the band to a `NoiseEstimate` (level and sides used per bin). OS, GO and SO estimators can be added later behind the same signature, each with its own alpha.

Dependencies: pocketfft and toml++ come in through CMake FetchContent (pinned); GoogleTest is the system package.

## Config

`configs/baseline.toml` follows the section names in the shared plan doc. Processing sections are read by the detector library (`cfar::load_processing_config`), so a real-data app needs no simulation code. The simulation sections are read by `cfar::sim::load_config`.

| Section | Read by | Keys |
| --- | --- | --- |
| `[band]` | library | `f_lo_hz`, `f_hi_hz` |
| `[acq]` | library | `fs_hz`, `nfft`, `window` (`"hann"` only) |
| `[false_alarm]` | library | `pfa_window`, `pfa_bin` (`"derived"` or a value) |
| `[cfar]` | library | `guard_per_side`, `train_per_side` |
| `[noise]` | sim | `sigma` |
| `[source]` | sim | `kind` (`"tone"` only), `f0_hz` (`"uniform"` or a value in Hz), `cn0_dbhz` |
| `[run]` | sim | `seed` |
| `[pfa_check]` | sim | `windows`, `design_pfa` (list; the derived value is always added) |
| `[pd_sweep]` | sim | `snr_bin_db_min`, `snr_bin_db_max`, `snr_bin_db_step` (bin-centred per-bin SNR), `trials`, `pd_target` |
| `[roc]` | sim | `snr_bin_db` (list), `pfa_window` (list; per-bin values are derived), `trials` |

The shared doc's `variant = "os"`, `os_rank_frac`, `welch_segments`, `[interference]` and `source.kind = "vfd"` wait for the deferred work below.

| Output file | Contents |
| --- | --- |
| example_spectrum.csv | One window's P[k], CA threshold and per-bin detections, for a slide |
| pfa_check.csv | Config, design Pfa, cells, hits, measured Pfa, 95% interval, pass, per-window false-alarm rate |
| pd_vs_snr.csv | Scenario, per-bin SNR, C/N0, trials, hits, Pd, 95% Wilson interval, ideal Pd, windows with a false alarm, config |
| pd_vs_snr_summary.csv | Pd target, C/N0 for CA and for the bound at that Pd, CFAR loss |
| roc.csv | Scenario, per-bin SNR, C/N0, Pfa per window and per bin, trials, hits, Pd, 95% Wilson interval, ideal Pd, windows with a false alarm, config |

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
