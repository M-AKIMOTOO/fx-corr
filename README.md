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
- `--coarse <s>`
- `--delay <samples>`
- `--rate <Hz>`
- `--resdelay <samples>`
- `--resrate <Hz>`
- `--rotation ...`
- `--skip <s>`
- `--length <s>`

Performance:

- `--cpu <N>`
- `--chunk-frames <N>`
- `--pipeline-depth <N>`

Diagnostics:

- `--debug` (writes per-frame debug log)
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

## License

MIT
