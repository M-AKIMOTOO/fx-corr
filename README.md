# fx-corr

FX correlator toolkit for VLBI-style dual-station processing.

This repository builds four binaries from the same codebase:

- `yi-acf`: auto-correlation only (A1xA1, A2xA2)
- `yi-xcf`: cross-correlation only (A1xA2)
- `yi-corr`: auto + cross correlation (A1xA1, A1xA2, A2xA2)
- `yi-phasedarray`: phased-array synthesis (`.raw`, phased `.cor`, plots)

## Features

- XML schedule driven operation (`--sc <schedule.xml>`)
- Multi-process scan support (`<process>` entries in XML)
- Per-antenna bit-depth/level/shuffle/sideband/rotation handling
- Geometric delay + clock delay/rate + residual delay/rate model
- Band alignment for different station center frequencies
- ACF/XCF/Phased product generation in one executable
- Throughput-oriented implementation:
  - rayon parallel compute
  - chunked I/O pipeline
  - process-window prefetch thread
  - FFT scratch reuse
  - Linux `posix_fadvise(SEQUENTIAL|WILLNEED)`

## Build

```bash
cd /home/akimoto/program/rust/fx-corr
cargo build --release
```

`fx-corr` now performs configure-like host probing at build time (`build.rs`):

- logical CPU count
- L3 cache size (Linux `/sys/devices/system/cpu/cpu0/cache/...`)

These build-time values are used as defaults for automatic `yi-corr` tuning
(`--cpu` omitted, `--chunk-frames` omitted). Manual CLI values still override.

Binaries are created under `target/release/`:

- `target/release/yi-acf`
- `target/release/yi-xcf`
- `target/release/yi-corr`
- `target/release/yi-phasedarray`

## Quick Start

### 1) Create a template XML

```bash
target/release/yi-corr --mkxml
```

This writes `example.xml` in the current directory.

### 2) Run correlation

```bash
target/release/yi-corr \
  --sc r24131a_cor.xml \
  --raw raw \
  --cor test/
```

### 3) Run per mode

```bash
# ACF only
target/release/yi-acf --sc r24131a_cor.xml --raw raw --cor test/

# XCF only
target/release/yi-xcf --sc r24131a_cor.xml --raw raw --cor test/

# ACF + XCF
target/release/yi-corr --sc r24131a_cor.xml --raw raw --cor test/

# Phased-array (writes outputs to current directory)
target/release/yi-phasedarray --sc r24131a_cor.xml --raw raw
```

## How input files are resolved

By default, input raw files are auto-resolved from:

- `--raw <DIR>`
- station names from XML (`<station><name>...`)
- epoch tag (`YYYYDDDhhmmss`)

Expected pattern:

- `<DIR>/<ANT1_NAME>_<YYYYDDDhhmmss>.raw`
- `<DIR>/<ANT2_NAME>_<YYYYDDDhhmmss>.raw`

Fallback station names are `YAMAGU32` and `YAMAGU34`.

You can override auto-resolution with explicit files:

```bash
--ant1 /path/to/ant1.raw --ant2 /path/to/ant2.raw
```

## Schedule XML (supported structure)

The parser reads the first `<process>` as default, and can run all `<process>` entries sequentially.

Important nodes:

- `<station key="A">`
  - `<name>`, `<pos-x>`, `<pos-y>`, `<pos-z>`, `<terminal>`
- `<clock key="A">`
  - `<delay>`, `<rate>` (seconds, seconds/second)
  - optional `<epoch>` clock-reference UTC (`delay` is interpreted at this epoch, then propagated to each process epoch using `rate`)
- `<terminal name="...">`
  - `<speed>` (Hz), `<bit>`, `<level>`
- `<shuffle key="A">`
  - 32-entry permutation
- `<source name="...">`
  - `<ra>`, `<dec>`
- `<stream>`
  - `<frequency>` (Hz), `<fft>` (optional), `<label>` (optional)
  - `<special key="A"><rotation>...</rotation><sideband>LSB|USB</sideband></special>`
- `<process>`
  - `<epoch>YYYY/DDD HH:MM:SS</epoch>`
  - `<skip>` (optional, s)
  - `<length>` (optional, s)
  - `<object>` (optional)
  - `<stations>` (optional station key string, e.g. `AB`)

## Multi-process behavior

If XML contains multiple `<process>` entries and you do not pass overriding scan controls, the tool runs all scans sequentially.

Automatic all-process mode is enabled only when all of these are unset/default:

- `--process-index` (hidden advanced option)
- `--epoch`
- `--ra`
- `--dec`
- `--length`
- `--skip` must be `0`

For process 2+ logs:

- compact logs are used
- `Antenna Parameters` are still printed for each process
- `Global Parameters` are omitted in compact mode

## Output files

### yi-acf

- `<ANT1>_<ANT1>_<TAG>_<LABEL>.cor`
- `<ANT2>_<ANT2>_<TAG>_<LABEL>.cor`

### yi-xcf

- `<ANT1>_<ANT2>_<TAG>_<LABEL>.cor`

### yi-corr

- `<ANT1>_<ANT1>_<TAG>_<LABEL>.cor`
- `<ANT1>_<ANT2>_<TAG>_<LABEL>.cor`
- `<ANT2>_<ANT2>_<TAG>_<LABEL>.cor`

### yi-phasedarray

- `YAMAGU66_<TAG>.raw` (phased time-series)
- `YAMAGU66_YAMAGU66_<TAG>_phasedarray.cor`
- `YAMAGU66_<TAG>_phased_auto_spectrum.png`
- `YAMAGU66_<TAG>_phased_spectrum_amplitude.png`
- `YAMAGU66_<TAG>_phased_autocorrelation.png`

Where:

- `TAG`: `YYYYDDDhhmmss` from process epoch (+ effective skip)
- `LABEL`: stream label from XML (`<stream><label>`), sanitized for file safety

Additionally, a run summary log is appended to:

- `<schedule.xml>.stdout.txt`

For `yi-corr`, full runtime stdout is also streamed when `--stdout` is specified:

- `./stdout/stdout_<yyyydddhhmmss>.log` (created at command start time)

## Key CLI options

Core:

- `--sc, --schedule <XML>`
- `--raw, --raw-directory <DIR>`
- `--cor, --cor-directory <DIR>` (required for `yi-acf/yi-xcf/yi-corr`)

Signal/processing:

- `--fft <N>`
- `--sampling <MSPS>`
- `--obsfreq <MHz>`
- `--bin <BITS>`
- `--level ...`
- `--shuffle ...`
- `--sideband LSB|USB`
- `--coarse <s>` (default: `0.0`)
- `--delay <samples>`
- `--rate <Hz>`
- `--resdelay <samples>`
- `--resrate <Hz>`
- `--rotation ...`
- `--skip <s>`
- `--length <s>`

### Delay Compensation Model

`yi-corr` applies delay correction in two parts for each FFT frame:

1. Integer-sample shift in time domain
2. Fractional-sample phase rotation in frequency domain

Per frame, total delay is evaluated at frame midpoint (`t_mid`) from geometric/clock/residual terms:

- `tau(t) = d + r*t + 0.5*a*t^2` (plus configured clock/residual terms)

Then it is split as:

- `n = round(tau * fs)` (integer samples)
- `eps = tau - n/fs` (fractional seconds)

This `round` split is intentional because it minimizes fractional residual:

- `eps` is always within about `[-0.5, +0.5]` sample
- Example: `12.983` samples is handled as `n=13`, `eps=-0.017` sample
- Example: `12.345` samples is handled as `n=12`, `eps=+0.345` sample

Fractional correction is applied as per-bin complex phase rotation:

- `X(f) <- X(f) * exp(-j*2*pi*f*eps)`

For packed quantized input (for example 2-bit/sample), integer sample shift corresponds to a bit offset conceptually:

- `bit_shift = n * bits_per_sample`
- For 2-bit sampling: `n=12` samples means `24` bits (`3` bytes)

Implementation note:

- The correlator decodes packed bits to samples first, then applies integer sample shift on decoded samples.
- This is equivalent to sample-aligned packed-bit shifting for delay compensation, while keeping the code path simple and robust.

### Geometric Delay / Rate / Accel and Doppler Tracking

Geometric terms are computed from:

- antenna ECEF coordinates (`<station><pos-x/pos-y/pos-z>`)
- source RA/Dec
- epoch (MJD)
- Earth rotation (sidereal-time rotation via ERFA)

In other words, geometric delay/rate/accel include Earth-rotation-driven fringe-rate (geometric Doppler tracking term).

Conceptually:

- `tau_geom = -((b2 - b1) · s) / c`
- `rate_geom` and `accel_geom` are time derivatives of `tau_geom`

Important behavior:

- If both antennas are set to the same coordinates (including both `[0, 0, 0]`), baseline vector is zero.
- Then geometric delay/rate/accel become approximately zero, so geometric Doppler tracking is effectively disabled.
- Clock/rate/rotation/user terms are still applied independently; only the geometric component is removed.

Note:

- If station coordinates are omitted in XML, built-in defaults are used (not zero), so geometric tracking remains enabled.

### Cor Normalization (Quantization-Loss Corrected)

`yi-corr` writes `.cor` with a fixed, physically motivated normalization (GICO3-compatible), including quantization-loss correction.

- Core formula:
  - `inv = 1 / (0.5 * P * nf * fft^2)`
  - `nf`: number of integrated FFT frames in the sector
  - `fft`: FFT length
  - `P`: quantization power term from level map
- Quantization power terms:
  - `P11 = mean(level1^2)` for `A1xA1`
  - `P22 = mean(level2^2)` for `A2xA2`
  - `P12 = sqrt(P11 * P22)` for `A1xA2`
- Example (typical 2-bit levels `[-1.5, 0.5, -0.5, 1.5]`):
  - `mean(level^2) = 1.25` (quantization-loss correction term)

About the `0.5` factor:

- This `0.5` is for **one-sided FFT storage** (`0..Nyquist` bins only in `.cor`).
- It is **not** the same concept as receiver **LSB/USB sideband** setting.

Performance:

- `--cpu <N>`
- `--chunk-frames <N>`
- `--pipeline-depth <N>`

Diagnostics:

- `--stdout` (yi-corr only; write runtime stdout log to `./stdout/stdout_<yyyydddhhmmss>.log`)
- `--debug` (writes full-frame debug log to `<schedule_dir>/debug_yi-corr/debug_<epoch>.log`)
- `--mkxml`

Advanced hidden:

- `--process-index <N>`
- `--delay-reference ant1|ant2`
- `--compact-logs`

## Performance tuning notes

1. Use release build only.
2. Keep `--cpu` near physical/logical core availability.
3. Usually start with automatic `--chunk-frames` and `--pipeline-depth`; tune only if needed.
4. Storage matters:
   - NVMe > SATA SSD > HDD/RAID
   - first run may be slower due to cold page cache
5. Very large `.raw` (10s of GB to TB) is supported as streaming I/O; full-RAM loading is not required.
6. Process-window prefetch reads only the required scan window, not the whole file.

## Scan window and end-time handling

Processing length is clamped safely by both requested window and available samples:

- `requested = max(process_length - total_skip, 0)`
- `processed = min(requested, available_from_file)`
- frame count uses floor to complete FFT frames only

This prevents overrun past scan/file end. Short reads are zero-padded when needed.

## CPU affinity (optional)

If present, the file below is used to pin rayon workers:

- `$CARGO_HOME/tmp/yi-corr-affinity.txt`
- or `~/.cargo/tmp/yi-corr-affinity.txt`

Format: one affinity set per line.

Examples:

```text
0-15
16-31
```

When multiple lines exist, one slot is claimed via lock file so concurrent runs can use different core sets.

## Example full command (2-station corr)

```bash
time target/release/yi-corr \
  --sc r24131a_cor.xml \
  --raw raw \
  --cor test/ \
  --cpu 18
```

## Troubleshooting

- `Input files not found`
  - Check station names in XML and file naming pattern.
  - Or pass `--ant1/--ant2` explicitly.
- `No band overlap`
  - Check stream frequency/rotation/sideband consistency between stations.
- Large run-to-run wall-time variation
  - Commonly due to storage/cache state, not algorithmic change.
- Need ACF and XCF together
  - Use `yi-corr` (not `yi-xcf` only).

## External RAID (3 Gbps) operation tip

If raw data is on slow external RAID, total time can be dominated by I/O.
Use staged continuous processing:

1. stage next scan from external RAID to local fast cache
2. run `yi-corr` on local cache
3. overlap 1 and 2 (one-scan-ahead)

Helper:

```bash
python3 tools/continuous_corr_runner.py \
  --schedule /path/to/obs.xml \
  --raw-src-dir /mnt/external/raw \
  --cache-dir /tmp/yi-corr-cache \
  --cor-dir /path/to/cor \
  --yi-corr-bin target/release/yi-corr \
  --cpu 18 --chunk-frames 65536 --pipeline-depth 8
```

Notes:

- This keeps correlation throughput near `max(stage I/O, yi-corr compute)`.
- With two 2-bit streams at 1024 Msps, required input is about 4.096 Gbps; a 3 Gbps link cannot achieve true real-time without reducing data rate.

## Phase Check (No Fixed Offset)

When you want to test phase consistency without fitting a constant offset `C`,
use:

- `phi_pred = wrap180(360 * f_peak * delta_tau)`
- `residual = wrap180(phi_obs - phi_pred)`

Helper script:

```bash
python3 tools/phase_no_offset_check.py \
  --yi yi-corr-freq.txt \
  --ref gico3-freq.txt \
  --debug-dir /path/to/debug_yi-corr \
  --min-snr 20
```

Notes:

- By default, rows where yi/ref peak frequency differs are skipped.
- `delta_tau` is derived from each `debug_*.log` seek diagnostics.
- If debug logs are unavailable, you can use a fixed sample offset:

```bash
python3 tools/phase_no_offset_check.py \
  --yi yi-corr-freq.txt \
  --ref gico3-freq.txt \
  --delta-samples 5.3 \
  --allow-freq-mismatch
```

## License

MIT
