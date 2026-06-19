# fx-corr

FX correlator toolkit for VLBI-style dual-station processing.

## Development Purpose and Research Context

本プロジェクトの長期的な目的は、VLBI 観測および結合型干渉計観測に利用できる GPU ベースの高速 FX 相関器を開発することである。GPU 相関器は、多局・広帯域・長時間観測に対して高い処理スループットを実現できる一方で、実装初期から GPU 固有の並列化、メモリ転送、カーネル分割、数値再現性の問題を同時に扱う必要がある。そのため、本リポジトリでは GPU 実装に進む前段階として、まず CPU 上で動作する相関器を開発し、信号処理モデル、XML スケジュール解釈、遅延補償、周波数グリッド整合、出力形式を検証している。

CPU 相関器は単なる試作ではなく、GPU 相関器開発の基準実装として位置づけている。具体的には、観測 XML に記述されたアンテナ座標、天体座標、サンプリング周波数、観測周波数、sideband、rotation、clock delay/rate、process epoch/skip/length を用いて、外部から恣意的な固定値を与えずに相関処理を行う。CPU 実装で幾何学的遅延、地球自転による fringe rate、地球公転速度に由来する一次補正、clock delay/rate、integer/fractional delay correction、band alignment を明示的に扱うことで、GPU 実装時にどの処理を FFT 前、FFT 後、積分前、積分後に配置すべきかを切り分ける。

研究費獲得に向けた開発計画では、この CPU 相関器を以下の役割に使う。

- GPU 相関器の数値検証用 reference implementation として利用する。
- 既存相関器および fringe search ツールとの比較により、delay/rate/phase の符号規約と周波数規約を確定する。
- XML schedule だけから相関処理パラメータを再現できることを示し、観測運用への接続性を確認する。
- CPU 実装で処理時間、I/O 帯域、FFT 負荷、積分処理の比率を測定し、GPU 化で高速化すべき部分を定量化する。
- GPU 実装後も CPU 版を regression test として残し、GPU kernel の変更が科学的な delay/rate/phase 結果を壊していないか確認する。

現在の開発段階では、CPU 相関器により、2-bit packed raw data の復号、per-station shuffle、LSB/USB 正規化、FFT、frequency rotation、ACF/XCF 出力、clock delay/rate、幾何学的 delay/rate/accel、地球速度一次補正を含む相関処理を実装している。これにより、GPU 相関器開発に入る前に、観測モデルとデータ読み取りモデルの不整合を CPU 上で発見・修正できる状態を作っている。

今後の GPU 相関器開発では、CPU 版で確定した処理を段階的に GPU へ移植する。第一段階では raw decode、sideband normalization、integer delay shift、FFT 入力整形を GPU 化する。第二段階では batched FFT と fractional delay phase correction、frequency rotation、XCF/ACF accumulation を GPU 上で実装する。第三段階では multi-baseline/multi-station 処理、streaming I/O、観測 scan の連続処理を GPU pipeline として統合する。CPU 版はこの各段階で比較対象として使い、GPU 化による高速化と科学的再現性の両方を評価する。

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

Install the Cargo subcommand once:

```bash
cargo install --path .
```

After that, plain `cargo fx-corr` runs the project install command and then installs the versioned `yi-corr` binary:

```bash
cargo fx-corr
# runs: cargo install --path .
# also installs yi-corr-vX.Y.Z into `$CARGO_HOME/bin` or `~/.cargo/bin`
```

The versioned binary step is also available explicitly:

```bash
cargo fx-corr install-versioned
```

Cargo build-target names cannot contain dots, and `cargo install --path .` only
installs declared Cargo bin targets. The versioned binary is therefore produced
by copying `target/release/yi-corr` after the normal release build.

## Quick Start

### 1) Create a template XML

```bash
target/release/yi-corr --mkxml
```

This writes `example.xml` in the current directory.

### 2) Run correlation

```bash
target/release/yi-corr \
  --sc test.xml \
  --raw raw \
  --cor test/
```

### 3) Run per mode

```bash
# ACF only
target/release/yi-acf --sc test.xml --raw raw --cor test/

# XCF only
target/release/yi-xcf --sc test.xml --raw raw --cor test/

# ACF + XCF
target/release/yi-corr --sc test.xml --raw raw --cor test/

# Phased-array (writes outputs to current directory)
target/release/yi-phasedarray --sc test.xml --raw raw
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
  - optional `<inband>N</inband>` power-of-two sub-band output split
  - optional `<pulsar>` pulse-phase folding model
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

### In-band `.cor` splitting

XML `<stream><inband>N</inband></stream>` splits the already-correlated positive-frequency band into `N` `.cor` files. `N` must be a power of two. For a 512 MHz band and `<inband>4</inband>`, bins are written as `[0:128)`, `[128:256)`, `[256:384)`, and `[384:512)`.

Split files add `.chNbwBW` before `.cor`, for example `<ANT1>_<ANT2>_<TAG>_<LABEL>.ch1bw512.cor`. Each split `.cor` header uses the sub-band low frequency as observing frequency and the effective sampling speed corresponding to the sub-band bandwidth. No additional delay correction or bandwidth synthesis is applied.

### Quick-look fringe plots

`yi-xcf` and `yi-corr` support `--fringe <S>` to write quick-look products every `S` seconds. The implementation accumulates the already fringe-stopped XCF frequency spectrum, writes a spectrum plot/TSV, then performs a lightweight 1D inverse FFT across frequency to write a lag-fringe plot/TSV. This is intended for fast inspection, not a full 2D delay/rate fringe search.

Example:

```bash
yi-corr --sc scan.xml --raw ../raw --cor cor_out --fringe 60
```

For each interval, files are written under `<cor>/fringe_yicorr_ql/<schedule-basename>/<process-epoch>/` with names containing `_fringe_spectrum.tsv`, `_fringe_spectrum_amp.png`, `_fringe_lag.tsv`, and `_fringe_lag_amp.png`.

### Pulsar folding `.cor` output

XML `<stream><pulsar>...</pulsar></stream>` enables pulse-phase folding for `.cor` spectral products. The pulsar model belongs to the stream, while `<process>` continues to define the scan epoch, skip, length, source, and stations. The first implementation uses topocentric phase folding; barycentric correction is reserved for a future mode.

Example:

```xml
<stream>
  <label>yi-corr-pulsar</label>
  <frequency>+6600e+6</frequency>
  <fft>1024</fft>
  <output>1</output>
  <inband>4</inband>
  <pulsar>
    <epoch>2025/302 08:15:00</epoch>
    <period>0.0335028583</period>
    <pdot>0.0</pdot>
    <pddot>0.0</pddot>
    <bins>32</bins>
    <dm>0.0</dm>
    <dedisperse>false</dedisperse>
    <time-correction>topocentric</time-correction>
  </pulsar>
</stream>
```

`<epoch>` is the pulse phase-zero epoch. If omitted, the process epoch is used. `<period>` is seconds, `<pdot>` is s/s, `<pddot>` is s/s^2, and `<bins>` is the number of pulse phase bins. When `<dedisperse>true</dedisperse>` and `<dm>` is non-zero, each frequency bin is assigned to a pulse phase bin after applying the cold-plasma dispersion delay relative to the top of the processed band; ACF and XCF are then accumulated in that corrected pulse bin.

Pulse-bin output adds `.pbinNN` before `.cor`, for example `<ANT1>_<ANT2>_<TAG>_<LABEL>.pbin00.cor`. With `<inband>4</inband>`, names become `<ANT1>_<ANT2>_<TAG>_<LABEL>.ch1bw512.pbin00.cor`.

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
- `--resacel <Hz/s>`
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

#### 0.5-sample boundary problem fixed in `3.0.0`

`fx-corr 3.0.0` fixes a phase/delay-branch instability that appeared when the
post-read-align residual delay crossed a half-sample boundary.  This issue was
most visible in short integrations of YAMAGU32--YAMAGU34 validation data near
`t = 56 s` and again near `t = 166 s`: `fringe`/`frinZ` could jump from one
residual-delay branch to the adjacent branch by almost exactly one sample, and
`frinZ --search deep --add` could report a single-point phase excursion of
about 30 degrees.

The important point is that the half-sample itself was not uncorrected.  A
fractional delay of `+0.5 sample` or `-0.5 sample` is a valid phase slope and
can be corrected in the frequency domain.  The bug was caused by mixing three
conceptually different operations:

1. **read alignment**: moving the raw-file read start by a large integer number
   of samples, such as `2047` samples;
2. **frame-local integer delay shift**: splitting the remaining per-frame
   residual delay with `round(delay * fs)` and applying the integer part to the
   decoded FFT input window; and
3. **fractional phase correction**: applying the remaining sub-sample delay as
   a complex phase slope across FFT bins.

Before the fix, the sector reader recomputed the read-align integer delay for
each 1-second sector.  In addition, the sector-local residual delay was still
split into integer and fractional parts by
`split_delay_to_integer_and_fractional()`.  This looked safe because the usual
split is mathematically valid for a single frame:

```text
integer_samples = round(delay * fs)
fractional      = delay - integer_samples / fs
```

However, it creates a discontinuity at the `±0.5 sample` boundary when the
split is used inside a sector/integration path.  For example, with residual
correction on antenna 2:

```text
residual = +0.499999998 sample
  tau2 = -0.499999998 sample
  int2 =  0
  frac2 = -0.499999998 sample

residual = +0.500000008 sample
  tau2 = -0.500000008 sample
  int2 = -1
  frac2 = +0.499999992 sample
```

Thus two adjacent FFT frames only `1 us` apart could use different integer
sample branches.  When those frames were accumulated into the same 1-second,
2-second, or 4-second correlation sector, the sector contained a mixture of two
integer-delay conventions.  The result was a narrow phase disturbance at the
half-sample crossing.  In validation, this appeared as a branch hop such as:

```text
Res-Delay:  -1.28 sample branch  ->  -0.28 sample branch
```

A second problem occurred at the boundary between adjacent 1-second sectors.
The old sector read-align path recomputed the integer read offset independently
for each sector:

```text
sector_d_seek_samples = round(full_rel(t_sector_ref) * fs)
```

When `full_rel` crossed `N + 0.5` samples, the sector read offset changed from
`N` to `N + 1`.  This caused the sector-local residual branch to fold from
`+0.5 sample` to `-0.5 sample` even though the physical full relative delay was
continuous:

```text
old behavior:
  sector 6 end:   full_rel = 2047.501323 sample, residual = +0.501323 sample
  sector 7 start: full_rel = 2047.501323 sample, residual = -0.498677 sample
```

This was not a source or atmosphere effect; it was a bookkeeping artifact from
changing the raw read origin at a half-sample boundary.

The `3.0.0` correction makes the read-align and residual-delay paths stable:

1. **The process/window read-align integer delay is fixed.**
   `sector_d_seek_samples` is no longer recomputed from the middle of each
   sector.  The same `d_seek_samples` chosen at the process/window start is
   used for all sectors.  Consequently, the residual delay does not fold at
   sector boundaries:

   ```text
   new behavior:
     sector 6 end:   residual = +0.501323 sample
     sector 7 start: residual = +0.501323 sample
   ```

2. **The sector-local residual is not split again into an integer shift and a
   fractional part.**  In the sector read-aligned path, the residual correction
   is carried entirely by the XCF phase slope:

   ```text
   int1 = 0
   int2 = 0
   frac1 = tau1
   frac2 = tau2
   ```

   Therefore the same half-sample crossing becomes continuous:

   ```text
   residual = +0.499999998 sample -> int=(0,0), frac2=-0.499999998
   residual = +0.500000008 sample -> int=(0,0), frac2=-0.500000008
   residual = +0.501323111 sample -> int=(0,0), frac2=-0.501323111
   ```

3. **The XCF phase correction uses the residual fractional delay.**
   The large read-align delay is already represented by the decoded sample
   window.  Applying it again as a baseband phase slope decorrelates the XCF.
   After the integer sample shift, only the remaining sub-sample term is used:

   ```text
   phase_delay1 = frac1
   phase_delay2 = frac2
   ```

   The carrier/fringe phase may keep a continuous delay branch for phase origin
   bookkeeping, but the per-bin XCF phase slope must stay on the fractional
   residual left after integer delay tracking.

After these changes, validation over 200 seconds showed that both `fringe` and
`frinZ --search deep --add` stayed on the same residual-delay branch through the
previously problematic crossings.  Around the original `t = 56 s` crossing, for
example, the results changed from a one-sample branch hop to a stable branch:

```text
before:
  08:15:55  Res-Delay ~= -1.26 sample
  08:15:56  Res-Delay ~= -0.70 sample or -0.28 sample
  08:15:57  Res-Delay ~= -0.28 sample

after:
  08:15:54  Res-Delay ~= -1.28 sample
  08:15:56  Res-Delay ~= -1.31 sample
  08:15:58  Res-Delay ~= -1.31 sample
```

The practical rule for future GPU and CPU implementations is:

- Use raw read alignment only as a stable coarse alignment for the processing
  window.
- Do not update the raw read origin at every 1-second sector unless the output
  phase convention explicitly accounts for that change.
- Do not reintroduce a per-frame integer residual split in the sector
  read-aligned XCF path.
- Carry the time-varying residual delay with the frequency-domain phase slope.
- Treat `0.5 sample` as a diagnostic boundary: if a phase jump, delay-branch
  hop, or one-point amplitude anomaly occurs exactly there, first check integer
  delay branch bookkeeping before changing geometric or clock models.


### G9.62 Maser Validation and Station/Grid-Origin Correction in `3.0.1`

`fx-corr 3.0.1` adds and documents a station/backend frequency-grid-origin
correction that was identified during validation with the narrow-band G9.62
maser data observed on `2025/302 07:00:00`.

This validation was important because a continuum source and a narrow maser
line stress different parts of an FX correlator.  A continuum source has signal
over a broad range of FFT bins, so a one-bin frequency-grid error can be partly
hidden by fringe search over frequency and residual rate.  A maser line is
narrow in frequency; if two stations place the same spectral line in adjacent
FFT bins, the XCF product

```text
C12[k] = X1[k] * conj(X2[k])
```

does not multiply the two line peaks at the same `k`.  The correlation amplitude
then drops even if the geometric delay and fringe-rate model are otherwise
reasonable.  For this reason, the G9.62 maser was used as a sensitive diagnostic
for whether each station's real-FFT output was being mapped onto the `.cor` XML
frequency grid with the correct bin origin.

The schedule used for this test had:

```text
source        : G9.62
epoch         : 2025/302 07:00:00
sampling      : 1024 MHz
FFT           : 1048576
bandwidth     : 512 MHz
XML frequency : 6600 MHz
channel width : 1024 MHz / 1048576 = 0.0009765625 MHz/bin
```

Thus an output bin near `69865` corresponds to a baseband frequency near

```text
69865 * 0.0009765625 MHz = 68.2275390625 MHz
```

The essential test was not based on residual delay or residual fringe rate.  It
was based on ACF spectral-line placement.  For a real astronomical maser seen by
two stations with the same XML reference frequency, sideband, rotation, and
bandwidth, the line must appear in the same output frequency bin in the two ACF
spectra.  Therefore, the measured quantity was

```text
Delta k = k_station - k_reference
```

where `k_reference` was taken from YAMAGU32.

#### HITACH32 one-bin grid-origin offset

Before the `3.0.1` correction, G9.62 ACF peak measurements showed that
YAMAGU32 and HITACH32 did not place the maser in the same bin:

```text
YAMAGU32 yi-corr ACF:
  sector peaks     = 69865, 69865, 69865, 69865, 69866
  integrated peak  = 69865
  centroid         ~= 69864.999

HITACH32 yi-corr ACF before correction:
  sector peaks     = 69864, 69864, 69863, 69864, 69865
  integrated peak  = 69864
  centroid         ~= 69864.020
```

The difference was approximately

```text
Delta k = k_HITACH32 - k_YAMAGU32 ~= -1 bin
```

This was not caused by HITACH32 being antenna 2, nor by the antenna to which the
read-align residual delay was applied.  The baseline order was reversed so that
HITACH32 became antenna 1 and YAMAGU32 became antenna 2.  The same HITACH32
low-bin behavior remained.  This eliminated the hypothesis that the effect came
from `ant2`-specific read alignment, integer-delay correction, or XCF residual
delay handling.

The same dataset also contained a YAMAGU32--YAMAGU34 short-baseline pair.  These
two stations use the same ADS3000-style format and the same shuffle family
(`24,25,...`).  Without any station offset, their ACF peaks were aligned:

```text
YAMAGU32--YAMAGU34 yi-corr ACF:
  sector peaks = 69865/69865, 69865/69865, 69865/69865,
                 69865/69865, 69866/69866
```

This demonstrated that the general LSB handling, `.cor` writer, XML frequency
definition, and YAMAGU32/YAMAGU34 shuffle handling were not the origin of the
HITACH32 one-bin error.  The remaining interpretation was a station/backend
frequency-grid-origin convention difference for HITACH32.

`3.0.1` therefore introduces a station/backend grid-origin correction:

```text
frequency_grid_origin_offset_bins(HITACH32) = -1
```

In the implementation, the output XML-grid bin is mapped to the station raw FFT
bin as:

```text
raw_idx = out_idx - rotation_bins + station/grid offset convention
```

Operationally, for HITACH32 the `-1` setting shifts the HITACH32 spectral line
one output bin upward, aligning the G9.62 line with YAMAGU32.  This is not a
rotation-frequency correction.  In the G9.62 validation data, HITACH32 had
`rotation = 0 MHz`; therefore the effect cannot be explained by `round`, `floor`,
or `ceil` of a non-integer rotation frequency.  The correction is a
station/backend frequency-grid-origin convention correction.

After the correction, the ACF debug log showed the intended behavior:

```text
YAMAGU32:
  station_offset=0  env_offset=0  total_offset=0
  peak ~= 69865/69866

HITACH32:
  station_offset=-1 env_offset=0  total_offset=-1
  peak ~= 69865/69866
```

The correction also worked after reversing the baseline order.  HITACH32 was
aligned whether it appeared as antenna 1 or antenna 2.  This is the key reason
the correction is stored as a station/backend property rather than as an
`ant2` diagnostic offset.

#### Why this is not an arbitrary manual offset

The offset was introduced by the following reproducible decision procedure:

1. Use a narrow spectral-line source with a high-SNR ACF peak.
2. Process the same scan with identical XML frequency, sampling rate, FFT
   length, bandwidth, sideband, and rotation for the stations being compared.
3. Measure the ACF peak and/or centroid bin for a reference station.
4. Measure the ACF peak and/or centroid bin for the target station.
5. Reverse the baseline order and verify that the offset follows the station,
   not `ant1` or `ant2`.
6. Compare with a station pair that uses the same recorder/format convention
   and verify that no offset is needed there.
7. Apply the smallest integer grid-origin correction that removes the stable
   station-specific ACF bin difference.
8. Verify that XCF sector-by-sector peak amplitudes recover and that unrelated
   station pairs are not changed.

The correction is therefore not a fitted science parameter and is not derived
from fringe phase, residual delay, or residual rate.  It is a station/backend
frequency-index convention needed to map decoded real-FFT spectra onto the
common XML `.cor` frequency grid.

The correction should be treated like other station/backend metadata such as
bit level map, shuffle map, sideband, and clock delay.  For the validated G9.62
setup:

```text
YAMAGU32:
  frequency_grid_origin_offset_bins = 0

YAMAGU34:
  frequency_grid_origin_offset_bins = 0

HITACH32:
  frequency_grid_origin_offset_bins = -1
```

#### XCF recovery after the HITACH32 correction

Before the HITACH32 grid-origin correction, the G9.62 maser XCF amplitude was
strongly suppressed because the YAMAGU32 and HITACH32 line peaks were not
multiplied at the same FFT bin.  After applying the station correction, the
sector-by-sector XCF peak amplitudes became comparable to GICO3:

```text
sector 0: gico3 1.529e-05 / yi-corr 1.482e-05  -> 0.97
sector 1: gico3 1.470e-05 / yi-corr 1.327e-05  -> 0.90
sector 2: gico3 1.413e-05 / yi-corr 1.398e-05  -> 0.99
sector 3: gico3 1.518e-05 / yi-corr 1.426e-05  -> 0.94
sector 4: gico3 1.347e-05 / yi-corr 1.301e-05  -> 0.97
```

The sector-to-sector phase did not have to match GICO3 exactly because yi-corr
and GICO3 do not necessarily use identical delay-model conventions, phase
origins, or residual-rate definitions.  The important validation result was
that the line was placed in the correct frequency bin and the per-sector
coherent amplitude recovered to the GICO3 level.

#### YAMAGU32--YAMAGU34 amplitude issue was a clock-setting error

During the same validation, the YAMAGU32--YAMAGU34 short-baseline result initially
showed about a factor-of-two smaller yi-corr XCF amplitude than GICO3, even
though the ACF peak bins were aligned.  This was initially suspicious because it
looked like a possible decoder, normalization, or ACF/XCF accumulation problem.

The observed ratio was:

```text
YAMAGU32--YAMAGU34 frinZ result before clock correction:
  yi-corr Amp ~= 2393 %
  gico3 Amp after excluding 511--512 MHz edge feature ~= 4589 %
  ratio ~= 0.52
```

The ACF power of YAMAGU34 in yi-corr was also much smaller than expected.  A
word-aligned read-start trial did not change the YAMAGU34 ACF power, ruling out
the hypothesis that a non-word-aligned packed-bit read start was the cause.

The actual cause was that the YAMAGU34 clock delay had not been included in the
XML.  After adding the YAMAGU34 clock delay,

```text
YAMAGU34 clock delay = +1.720000e-6 s
YAMAGU34 clock rate  = 0
```

the YAMAGU34 ACF power increased from about `1.2e10--1.3e10` to about
`4.4e10--4.9e10`, and the XCF amplitude recovered:

```text
YAMAGU32--YAMAGU34 yi-corr after clock correction:
  Amp  = 4706.4 %
  SNR  = 2740.1
  freq = 68.2275391 MHz

YAMAGU32--YAMAGU34 gico3 with 511--512 MHz edge feature excluded:
  Amp  ~= 4589 %
  SNR  ~= 2675
  freq ~= 68.228516 MHz
```

Therefore, the YAMAGU32--YAMAGU34 factor-of-two amplitude problem was not a
yi-corr normalization bug.  It was a schedule/model configuration error: the
station clock term for YAMAGU34 was missing.  This is an important operational
lesson for future validation: if ACF bins are aligned but XCF amplitude is
suppressed, confirm station clock delay/rate before changing decoder or
normalization code.

#### Consequences for continuum and maser correlation

The G9.62 validation separated two classes of problems:

```text
HITACH32:
  problem  : maser line appeared about one FFT bin too low in yi-corr
  cause    : station/backend frequency-grid-origin convention
  solution : frequency_grid_origin_offset_bins = -1

YAMAGU34:
  problem  : short-baseline XCF amplitude was about half of GICO3
  cause    : missing station clock delay in XML
  solution : add YAMAGU34 clock delay = +1.72 us
```

With these corrections, yi-corr can process both continuum and narrow-band maser
data without ad hoc per-run bin offsets.  The remaining station-specific values
are not arbitrary science-tuning parameters; they are part of the instrumental
configuration needed to interpret the packed raw data and station clocks:

```text
YAMAGU32:
  clock delay/rate = 0 / 0
  grid-origin offset = 0

YAMAGU34:
  clock delay/rate = +1.72e-6 / 0
  grid-origin offset = 0

HITACH32:
  clock delay/rate from XML schedule
  grid-origin offset = -1
```

For continuum sources, a one-bin grid-origin error may be less obvious because
the signal is broadband.  For narrow maser lines, the same error is immediately
visible as ACF bin disagreement and XCF amplitude loss.  The maser validation
therefore provides a stricter test of the frequency-grid convention than a
continuum-only validation.

The practical validation rule after `3.0.1` is:

- ACF line peaks should agree between stations that observe the same RF band.
- Reversing baseline order must not move a station-specific ACF offset to the
  other antenna.
- Same-format station pairs, such as YAMAGU32--YAMAGU34, should not need the
  HITACH32 grid-origin correction.
- XCF sector peak amplitudes should be compared sector-by-sector before drawing
  conclusions from an integrated spectrum, because different correlators may
  use different phase origins or residual-rate conventions.
- Missing or wrong station clock terms can mimic a correlator amplitude bug,
  especially for high-SNR narrow-band maser data.

For future GPU implementation, this result means that the GPU pipeline must not
hard-code a single universal FFT-bin origin for all backends.  The station or
backend metadata must include the same frequency-grid-origin convention used by
the CPU reference implementation.  The CPU implementation in `3.0.1` is
therefore the reference for both the signal-processing order and the validated
station/backend convention.


### G9.62 Maser Validation, Station/Grid-Origin Correction, and Clock-Delay Lessons in `3.0.1`

`fx-corr 3.0.1` was produced after a focused validation campaign using the
G9.62 maser data observed on `2025/302 07:00:00`.  The purpose of this campaign
was not only to increase the SNR of one data set, but to identify which
remaining discrepancies between `yi-corr` and GICO3 were caused by real
instrumental metadata and which were caused by correlator implementation bugs.
This distinction is essential for using `yi-corr` as the CPU reference
implementation for the future GPU correlator: the GPU version must inherit
station/backend conventions and clock metadata, but it must not inherit
ad-hoc, data-set-specific tuning values.

The investigation separated two independent issues:

```text
1. HITACH32 placed the G9.62 maser line about one FFT bin lower than YAMAGU32
   in yi-corr, even though both stations used the same RF/grid definition.

2. YAMAGU32--YAMAGU34 initially showed about a factor-of-two lower yi-corr
   XCF amplitude than GICO3, even though their ACF line bins were aligned.
```

The first issue was solved by introducing a station/backend frequency-grid-origin
correction for HITACH32.  The second was solved by adding the missing YAMAGU34
clock delay to the XML schedule.  These are both instrumental/model parameters,
not arbitrary science-fitting knobs.

#### Why the G9.62 maser is a stricter test than a continuum source

Continuum and maser sources stress different parts of an FX correlator.  A
continuum source has signal over many FFT bins.  If one station's spectrum is
mis-registered by one bin, broadband fringe search can still find a strong
solution, especially when the search allows residual frequency/rate and when the
bandpass has broad structure.  A narrow maser line is much less forgiving.  If
station 1 puts the line at bin `k` and station 2 puts the same line at `k-1`,
then the XCF product at a fixed output bin

```text
C12[k] = X1[k] * conj(X2[k])
```

multiplies the peak of one station by the wing/noise of the other station.  The
ACF spectra may each show a strong maser line, but the XCF amplitude is reduced
because the two line peaks are not multiplied at the same frequency index.  This
made G9.62 a sensitive diagnostic for the correctness of:

- packed raw decoding;
- station-specific shuffle and level maps;
- LSB-to-USB normalization;
- real-FFT positive-frequency indexing;
- rotation/bin shifting onto the XML grid;
- `.cor` output-bin convention;
- station/backend-specific frequency-grid origin.

The validation scan used:

```text
source        : G9.62
process epoch : 2025/302 07:00:00
sampling      : 1024 MHz
FFT length    : 1048576
bandwidth     : 512 MHz
XML frequency : 6600 MHz
sideband      : LSB
rotation      : 0 MHz for the G9.62 YAMAGU32/HITACH32/YAMAGU34 tests
```

The spectral channel width was therefore:

```text
df = fs / FFT
   = 1024 MHz / 1048576
   = 0.0009765625 MHz/bin
   = 0.9765625 kHz/bin
```

The observed maser line near output bin `69865` corresponds to a baseband
frequency of approximately

```text
69865 * 0.0009765625 MHz = 68.2275390625 MHz.
```

This frequency is the one reported by `frinZ` near `+68.2275391 MHz` or
`+68.2285156 MHz`, depending on the exact peak-bin/parabolic-search result.

#### What was measured: ACF line placement, not residual fringe phase

The HITACH32 correction was not derived from the XCF phase, residual delay, or
residual fringe rate.  The primary measurement was the position of the maser
line in each station's ACF spectrum.  This was deliberate: ACF line placement is
insensitive to baseline phase conventions and to the particular fringe-search
reference used by `frinZ` or GICO3.  If two stations observe the same RF band and
are mapped onto the same XML `.cor` grid, the same astronomical line must appear
at the same ACF output bin.

For a station `S`, define a measured line location from the ACF spectrum:

```text
k_S = peak or centroid bin of the G9.62 maser line in ACF(S,S)
```

and compare it to the reference station YAMAGU32:

```text
Delta k_S = k_S - k_YAMAGU32.
```

For the G9.62 validation, `Delta k` was stable at approximately `-1 bin` for
HITACH32 before correction.  That is the calculation that motivated the
station/backend grid-origin correction.  It is not a delay calculation.  It is
not computed from baseline geometry.  It is a frequency-index convention check.

#### Initial symptom: HITACH32 was approximately one bin lower

Before the `3.0.1` correction, yi-corr placed the HITACH32 maser line about one
bin lower than YAMAGU32 in the same G9.62 scan.  A representative ACF comparison
was:

```text
YAMAGU32 yi-corr ACF:
  sector peaks     = 69865, 69865, 69865, 69865, 69866
  integrated peak  = 69865
  centroid         ~= 69864.999

HITACH32 yi-corr ACF before correction:
  sector peaks     = 69864, 69864, 69863, 69864, 69865
  integrated peak  = 69864
  centroid         ~= 69864.020
```

Thus:

```text
Delta k_HITACH32 ~= 69864.020 - 69864.999 ~= -0.98 bin.
```

GICO3 did not show the same yi-corr low-bin placement for HITACH32.  In the
GICO3 ACF, the HITACH32 maser centroid was around `69865.38`, i.e. on the same
side as the YAMAGU32 maser line.  The problem was therefore not that HITACH32
really observed a different sky frequency.  The problem was how yi-corr mapped
HITACH32's decoded real-FFT bins onto the output XML grid.

#### Eliminated explanations

The following possible explanations were tested and rejected during the
investigation.

1. **It was not caused by the antenna-2 path.**

   The baseline was reversed so that HITACH32 became antenna 1 and YAMAGU32
   became antenna 2.  Without the station correction, the low-bin behavior still
   followed HITACH32.  Therefore the effect was not caused by applying the
   residual delay to antenna 2, by `ant2` read alignment, or by an `ant2` XCF
   accumulation path.

2. **It was not a global LSB or `.cor` writer problem.**

   The same G9.62 data included a YAMAGU32--YAMAGU34 pair.  Both stations used
   the same ADS3000-style format and the same shuffle family (`24,25,...`).
   Without any station grid-origin offset, their ACF peaks were aligned:

   ```text
   YAMAGU32--YAMAGU34 yi-corr ACF:
     sector peaks = 69865/69865, 69865/69865, 69865/69865,
                    69865/69865, 69866/69866
   ```

   Therefore the general LSB-to-USB normalization, `.cor` output indexing,
   YAMAGU32/YAMAGU34 shuffle handling, and XML frequency definition were not
   globally wrong.

3. **It was not a non-integer rotation-rounding problem.**

   In the G9.62 validation used for the HITACH32 decision, the station rotation
   was `0 MHz`.  Therefore `round`, `floor`, or `ceil` of
   `rotation_hz / df` cannot explain the one-bin offset.  The correction is not
   a hidden rotation correction.

4. **It was not solved by bit-level or parity trials.**

   Several diagnostic trials changed parity/bit interpretation or applied
   one-off XCF bin offsets.  These were useful for identifying the direction of
   the mismatch, but they were not acceptable as final explanations.  The final
   correction had to follow the station/backend, not the baseline order or an
   environment variable.

5. **It was not a `[0, FFT/2)` versus `(0, FFT/2]` global convention error.**

   A global endpoint convention error would affect the same-format
   YAMAGU32--YAMAGU34 pair and other validation data.  Instead, the offset was
   tied to HITACH32's backend/format convention.

#### Meaning of `frequency_grid_origin_offset_bins = -1`

`frequency_grid_origin_offset_bins` is a station/backend metadata value used
when placing a station's real-FFT spectrum onto the common XML output grid.  It
is conceptually separate from:

- observing RF;
- XML `<frequency>`;
- sideband;
- rotation frequency;
- geometric delay;
- clock delay;
- residual fringe rate.

For a normal station with no extra grid-origin correction, the output-grid index
is mapped to the raw FFT index approximately as:

```text
raw_idx = out_idx - rotation_bins
```

With a station/backend grid-origin correction, the mapping is effectively:

```text
raw_idx = out_idx - rotation_bins + frequency_grid_origin_offset_bins
```

For HITACH32:

```text
frequency_grid_origin_offset_bins = -1
```

so the mapping becomes:

```text
raw_idx = out_idx - rotation_bins - 1.
```

Operationally, this means that a raw HITACH32 spectral feature that previously
appeared one output bin too low is placed one output bin higher on the XML grid.
For the G9.62 line, this moved HITACH32 from the `69864/69863` side to the
`69865/69866` side, matching YAMAGU32.

The sign is easy to misunderstand.  The stored value is `-1`, but the visible
line in the output spectrum moves toward a larger output-bin number because the
output bin now reads the raw spectrum from one bin lower:

```text
extra = -1  ->  raw_idx = out_idx - 1
             ->  raw peak at 69864 appears at output bin 69865
```

#### Why this is not an arbitrary manual offset

The correction was introduced by a reproducible decision procedure:

1. Use a high-SNR narrow spectral-line source.
2. Process the same scan with identical XML frequency, sampling rate, FFT
   length, bandwidth, sideband, and rotation for the stations being compared.
3. Measure ACF peak and centroid bins for a reference station.
4. Measure ACF peak and centroid bins for the target station.
5. Reverse the baseline order and verify that the offset follows the station,
   not `ant1` or `ant2`.
6. Compare with a same-format station pair and verify that no offset is needed
   there.
7. Apply the smallest integer grid-origin correction that removes the stable
   station-specific ACF bin difference.
8. Verify that XCF sector-by-sector peak amplitudes recover.
9. Verify that unrelated stations are unchanged.

The result for the validated G9.62 setup is:

```text
YAMAGU32:
  frequency_grid_origin_offset_bins = 0

YAMAGU34:
  frequency_grid_origin_offset_bins = 0

HITACH32:
  frequency_grid_origin_offset_bins = -1
```

This is the same type of information as a station's bit-level map, shuffle map,
sideband convention, or clock model.  It should be treated as instrumental
metadata, not as a free parameter adjusted to maximize a science result.

#### Debug logging in `3.0.1`

`3.0.1` prints the station/backend offset explicitly in ACF peak debug output.
For HITACH32, a typical log line has:

```text
[acfpeakdbg] ... ant=HITACH32 ... station_offset=-1 env_offset=0 total_offset=-1 peak=69865 ...
```

For YAMAGU32 and YAMAGU34, the corresponding fields are:

```text
station_offset=0 env_offset=0 total_offset=0
```

The three fields have different meanings:

```text
station_offset : built-in station/backend correction
                 e.g. HITACH32 = -1

env_offset     : optional diagnostic offset from environment variables
                 used only for experiments, not for normal processing

total_offset   : station_offset + env_offset
                 the actual offset applied in the FFT-to-XML-grid mapping
```

Normal validation and production processing should use `env_offset=0`.  If a
non-zero environment offset is used in a diagnostic run, it should be reported
explicitly and not mixed with production validation products.

#### XCF recovery after the HITACH32 correction

After applying the HITACH32 station correction, the per-sector G9.62 XCF peak
amplitudes became comparable to GICO3:

```text
sector 0: gico3 1.529e-05 / yi-corr 1.482e-05  -> 0.97
sector 1: gico3 1.470e-05 / yi-corr 1.327e-05  -> 0.90
sector 2: gico3 1.413e-05 / yi-corr 1.398e-05  -> 0.99
sector 3: gico3 1.518e-05 / yi-corr 1.426e-05  -> 0.94
sector 4: gico3 1.347e-05 / yi-corr 1.301e-05  -> 0.97
```

This is the key validation result.  The sector-to-sector phase difference
between yi-corr and GICO3 did not have to be identical because the two programs
may use different delay-model details, phase origins, and residual-rate
conventions.  The important result is that the same narrow spectral line is now
placed in the same bin and the per-sector coherent amplitude is restored to the
GICO3 level.

When checking this, compare per-sector peaks before interpreting a fully
integrated spectrum.  If two correlators have different residual phase/rate
conventions, the integrated spectrum can differ even when each individual sector
has the correct amplitude.  Sector-by-sector XCF comparison is therefore a more
direct diagnostic for frequency-grid alignment.

#### YAMAGU32--YAMAGU34 amplitude issue: missing YAMAGU34 clock delay

The YAMAGU32--YAMAGU34 short-baseline test initially appeared to reveal another
correlator amplitude problem.  The ACF bins were aligned, but yi-corr's XCF
amplitude was about half of GICO3:

```text
YAMAGU32--YAMAGU34, before clock correction:
  yi-corr Amp ~= 2392.7 %
  yi-corr SNR ~= 1394.9

GICO3, after excluding the 511--512 MHz edge feature:
  Amp ~= 4589.3 %
  SNR ~= 2675.5
```

This ratio was close to the square-root of the YAMAGU34 ACF power discrepancy,
which made the problem look like a possible decoder or normalization bug.  A
trial that forced read alignment to a 32-bit word boundary changed the input
start but did not change the YAMAGU34 ACF power.  This ruled out the hypothesis
that a non-word-aligned packed-bit read start was the cause.

The actual cause was simpler and more important operationally: the YAMAGU34
clock delay had not been included in the XML schedule.  After adding:

```text
YAMAGU34 clock delay = +1.720000e-6 s
YAMAGU34 clock rate  = 0
```

the read-align delay for the YAMAGU32--YAMAGU34 scan changed from about
`234.9 samples` to about `1996.2 samples`, and the XCF amplitude recovered:

```text
YAMAGU32--YAMAGU34, after YAMAGU34 clock correction:
  yi-corr Amp  = 4706.4 %
  yi-corr SNR  = 2740.1
  freq         = +68.2275391 MHz
  residual rate ~= +0.013650 Hz
```

The YAMAGU34 ACF debug power also increased strongly:

```text
before clock correction:
  YAMAGU34 ACF debug power ~= 1.2e10--1.3e10

after clock correction:
  YAMAGU34 ACF debug power ~= 4.4e10--4.9e10
```

The lesson is that a missing station clock term can mimic a correlator
normalization or decoder problem, especially for a high-SNR narrow-band maser.
If ACF line bins are aligned but XCF amplitude is unexpectedly low, station
clock delay/rate should be checked before changing decode, normalization, or
FFT-grid code.

#### Why the 511--512 MHz GICO3 feature was excluded in the YAMAGU32--YAMAGU34 comparison

In the YAMAGU32--YAMAGU34 GICO3 output, an extremely strong feature near the
Nyquist edge (`~511.999 MHz`) could be selected by an unconstrained frequency
search.  This feature is not the G9.62 maser line at `~68.228 MHz`.  Therefore,
for the maser comparison the GICO3 search was repeated with the edge region
masked:

```text
frinZ --in cor/YAMAGU32_YAMAGU34_2025302070000_gico3.cor \
      --plot --fre --search deep --rfi 511,512
```

After this mask, GICO3 selected the same maser-frequency region as yi-corr.  The
comparison of the YAMAGU32--YAMAGU34 amplitude should use this RFI/edge-masked
GICO3 result, not the unmasked `511.999 MHz` solution.

#### Final operational interpretation after `3.0.1`

The G9.62 investigation supports the following operational interpretation:

```text
YAMAGU32:
  clock delay/rate from XML, usually 0/0 for the tested schedule
  frequency_grid_origin_offset_bins = 0

YAMAGU34:
  clock delay/rate must include the measured station clock term
  validated value in this test: +1.720000e-6 s / 0
  frequency_grid_origin_offset_bins = 0

HITACH32:
  clock delay/rate from XML schedule
  frequency_grid_origin_offset_bins = -1
```

With these values, `yi-corr` can process both continuum and narrow-band maser
scans without per-run hand-tuned bin offsets.  The remaining non-zero values are
station/backend or station-clock metadata.  They should be documented in the
schedule/configuration layer and carried into the GPU implementation.

#### Recommended validation checklist for future scans

For each new station/backend/channel combination, especially before using the
GPU correlator for science processing, run the following checks.

1. **ACF line placement**

   For a narrow maser or injected tone, compare ACF peak/centroid bins between
   stations.  A stable integer-bin difference that follows a station after
   baseline reversal indicates a station/backend grid-origin convention.

2. **Baseline reversal**

   Reverse station order in the XML process.  A real station/backend offset
   follows the station.  An `ant1`/`ant2` implementation bug follows the antenna
   index.

3. **Same-format control pair**

   Compare two stations with the same recorder/format convention.  If the
   same-format pair is aligned without correction, do not introduce a global
   offset.

4. **XCF sector amplitude**

   Compare XCF peak amplitude sector-by-sector before comparing fully integrated
   spectra.  Phase/rate convention differences can affect integrated products
   even when per-sector amplitudes are correct.

5. **Clock delay/rate sanity check**

   Confirm XML `<clock>` values and clock epochs.  A missing clock delay can
   suppress XCF amplitude while leaving ACF line bins apparently reasonable.

6. **RFI and edge-bin masking**

   Ensure the frequency-search tool is selecting the intended astronomical line,
   not a Nyquist-edge or band-edge artifact.

7. **No environment offsets in production**

   Environment variables that alter grid offsets are useful diagnostics, but
   final validation products should use station/backend metadata with
   `env_offset=0`.

#### Implications for the GPU correlator

The GPU implementation should not assume that all stations share a universal
FFT-bin origin after real FFT.  The CPU implementation in `3.0.1` defines the
reference behavior:

- decode packed raw samples using the station's level and shuffle map;
- normalize sideband using absolute sample parity;
- FFT real-valued samples;
- shift spectra onto the XML grid using rotation plus station/backend
  grid-origin metadata;
- apply XCF phase correction on the same frequency grid used for accumulation;
- use station clock delay/rate from XML rather than substituting empirical
  post-correlation offsets.

For GPU development, `frequency_grid_origin_offset_bins` should be part of the
same metadata bundle as shuffle, level map, sideband, and clock terms.  It
should be applied before or during the GPU-side FFT-to-grid gather/scatter step,
not as a later science-product correction.  The G9.62 maser test should remain a
regression test because it can reveal one-bin mapping mistakes that continuum
validation may not expose.

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

#### 見かけ上 0.05 Hz の固定 residual rate が合った理由

YAMAGU32/USUDA64 の検証では、fringe search の residual fringe rate として約 `-0.05 Hz` が繰り返し現れた。この値は経験的な補正値として導入したものではなく、GICO3 に合わせるための固定オフセットでもない。相関器の幾何学的遅延モデルに物理項が不足していることを示す診断値だった。

観測された residual delay は時間に対してほぼ直線的に変化していた。例えば、その傾きが

```text
-0.0064 sample/s
```

程度であれば、`fs = 1024 MHz`、`f_obs = 8192 MHz` では

```text
residual_rate_hz = delay_slope_sample_per_s * (f_obs / fs)
                 = -0.0064 * (8192 / 1024)
                 = -0.0512 Hz
```

となる。つまり `0.05 Hz` は、相関器が出力した複素スペクトルに残っていた一次の delay-rate 成分を観測周波数で表したものに相当する。`frinZ` 固有の効果ではない。`frinZ` は書き出された複素スペクトルを fringe search しているだけで、アンテナ座標、天体座標、時刻モデルは参照しない。

不足していたモデル項は、天体方向に沿った地球の太陽系重心速度による一次補正である。ただし、この補正を絶対 read-align delay にそのまま掛けると、短い scan でも 100 sample 以上の定数 delay shift を生む。そこで CPU 基準実装では、process epoch を基準にした differential correction として扱う。

```text
tau_bary(t) = tau_geocentric(t) * (1 - V_earth(t) . s / c)
tau_model(t) = tau_geocentric(t_ref) + (tau_bary(t) - tau_bary(t_ref))
```

ここで `V_earth` は ERFA `eraEpv00` から得る地球速度、`s` は天体方向単位ベクトル、`c` は光速、`t_ref` は process epoch である。この形にすると、raw read-align を決める絶対 delay の基準を不要に動かさず、地球速度に由来する rate/accel 成分だけを fringe tracking に入れられる。この観測では、この物理的に計算される differential term により、モデル上の rate が約 `0.05 Hz` 分だけ動き、残っていた residual rate の主成分を説明できた。

短時間の scan ではこの補正量がほぼ一定に見えるため、`0.05 Hz` の固定値が合ったように見える。しかし実際には固定値ではなく、観測 epoch と天体方向に依存する地球速度ベクトルから決まる量である。GPU 相関器では `0.05 Hz` を hard-code してはいけない。CPU 基準実装と同じく、地球自転による fringe rate、clock delay/rate、地球速度の一次補正をモデルとして持つ必要がある。

### XML Frequency, Rotation, Sideband, Clock, and Baseline Semantics

Several GICO3-compatibility details are more fundamental than Hermitian
completion. These rules define the frequency and time reference before FFT
spectra are formed.

#### XML `<frequency>` is the correlation reference

The stream `<frequency>` is the output/correlation reference frequency. It is
not replaced by each station sampled-band edge. All spectra are shifted onto
this XML frequency grid before `.cor` output is written.

For each antenna:

```text
data_low_mhz = xml_frequency_mhz + rotation_hz / 1e6
xml_grid_low_mhz = data_low_mhz - rotation_hz / 1e6
                 = xml_frequency_mhz
```

Example:

```text
YAMAGU32: frequency=8192 MHz, rotation=0 MHz
  data band: 8192..8704 MHz
  XML grid:  8192..8704 MHz

USUDA64:  frequency=8192 MHz, rotation=-112 MHz
  data band: 8080..8592 MHz
  XML grid after rotation: 8192..8704 MHz
```

This means `<rotation>` describes how the sampled band is moved onto the XML
frequency grid. In normal `<rotation>` mode, XCF integration uses only XML-grid
bins where both antennas have real observed signal after the shift. The `.cor`
header observing frequency remains the XML `<frequency>` value.

`<rotation2>` is an explicit opt-in for the old Hermitian-completion behavior.
When `<rotation2>` is present in a station `<special>` block, its value is used
as that station rotation and missing image-band bins may be reconstructed from
real-FFT Hermitian symmetry. Use `<rotation2>` only for validation cases where
that image-band correlation is intentionally desired.

#### Rotation has two roles

Rotation affects both:

1. frequency placement of the sampled spectrum
2. fringe/phase tracking associated with that frequency move

The integer-bin shift is:

```text
rotation_bins = round(rotation_hz / df)
df = sampling_rate / fft
```

The exact shifted part is handled by moving FFT bins. Any non-integer residual
is treated as a fringe/phase term:

```text
rotation_residual_hz = rotation_delta_hz - rotation_bins * df
rotation_fringe_hz  = -rotation_residual_hz
```

For exact grid-aligned rotations such as `112 MHz` with `df=0.125 MHz`, the
residual is zero.

#### All internal correlation is USB-normalized

The correlator normalizes inputs into the USB domain before correlation.

- XML/CLI `sideband=USB`: no conversion
- XML/CLI `sideband=LSB`: apply LSB-to-USB conversion
- omitted sideband defaults to USB

For real sample streams, LSB-to-USB conversion is implemented as alternating
sample sign inversion:

```text
x_usb[n] = (-1)^n x_lsb[n]
```

The parity `n` is based on the absolute raw sample index, not on the local
index inside each FFT frame. This matters when a scan starts at a nonzero sample
or after read alignment.

#### Shuffle default and station-specific shuffle

The default shuffle map is reverse bit order:

```text
31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,
15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0
```

XML `<shuffle key="...">` overrides this per station. The display order follows
GICO3-style external bit numbering, while the decoder converts it to the
internal LSB-first word layout used by the packed raw reader.

#### Clock epoch is the delay reference epoch

For XML clocks:

```xml
<clock key="U"><epoch>2026/026 02:56:00</epoch><delay>1.0</delay><rate>0</rate></clock>
```

`delay` is interpreted at that clock `<epoch>`. If the process epoch differs,
clock delay is propagated to the process epoch using clock rate:

```text
clock_delay_at_process = delay_at_clock_epoch
                       + clock_rate * (process_epoch - clock_epoch)
```

The baseline-relative clock term is then:

```text
clock_relative_delay = clock2_delay - clock1_delay
clock_relative_rate  = clock2_rate  - clock1_rate
```


#### Earth orientation and delay-model reference

For schedule XML runs, Earth-orientation and delay-model reference parameters
come from the XML schedule. CLI values are only fallbacks for non-XML direct
runs. Add these optional top-level nodes under `<schedule>` when EOP values are
known:

```xml
<eop><dut1>0.0</dut1><tt-utc>69.184</tt-utc><xp>0.0</xp><yp>0.0</yp></eop>
<delay-model><time-offset>0.0</time-offset></delay-model>
```

If omitted, yi-corr uses DUT1=0, TT-UTC=69.184, xp=0, yp=0, and
time-offset=0.0 s; no empirical GICO3-matching delay offset is applied by default.

#### EOP file operation plan

For detection-oriented VLBI correlation, EOP is an a priori model: it only has
to keep the signal inside the delay/rate search window and avoid decorrelation
over the requested coherent integration time. It is not the final geodetic EOP
solution. A higher-quality EOP model improves long integrations and reduces the
load on fringe fitting, but yi-corr should still be able to run in an approximate
mode when the purpose is first detection.

The preferred precise/operational EOP source is IERS rapid/standard
`finals2000A.data`:

```text
https://data.iers.org/products/eop/rapid/standard/finals2000A.data
```

Downloaded EOP tables are time-dependent input data, not Rust source files. Do
not place `finals2000A.data` under `src/`; `src/eop.rs` should contain only the
parser, interpolation, cache-selection, and optional-update logic. A practical
local layout is:

```text
data/eop/finals2000A.data                  # current local default
data/eop/archive/finals2000A-YYYYMMDD.data # downloaded snapshots
```

Old snapshots are important because they preserve the exact Bulletin A prediction
state used by an earlier correlation run. Automatic network access, if
implemented, must be opt-in only, for example with `YI_EOP_AUTO_DOWNLOAD=1`, and
must not silently replace the EOP used for a run. The run log should record the
EOP file path, snapshot name or timestamp, selected MJD rows, interpolated
DUT1/xp/yp values, TT-UTC, and whether Bulletin A prediction or Bulletin B final
values were used.

Implemented EOP modes:

```text
YI_EOP_MODE=auto    use XML/file EOP when available; otherwise run approximate + warning
YI_EOP_MODE=strict  require XML/file EOP; fail if unavailable or outside table range
YI_EOP_MODE=none    force approximate EOP: DUT1=0, xp=0, yp=0, TT-UTC from fallback
```

Implemented precedence:

```text
1. XML <eop> explicit values, when present
2. XML <eop file="..."> or YI_EOP_FILE, interpolated from finals2000A.data
3. data/eop/finals2000A.data, when present
4. CLI fallback values --dut1, --tt-utc, --xp, --yp
5. Environment overrides YI_DUT1_S, YI_TT_UTC_S, YI_XP_ARCSEC, YI_YP_ARCSEC
6. Approximate EOP only when YI_EOP_MODE=auto/none allows it
```

For `finals2000A.data`, use Bulletin B values when available for past dates and
fall back to Bulletin A rapid/predicted values for recent or future dates. In
first-detection workflows, missing or stale EOP should be a logged warning rather
than a hard error unless strict mode is selected.

#### Source coordinates: J2000 header vs delay-model coordinates

Source RA/Dec in XML are treated as J2000 catalog coordinates. The `.cor`
header keeps these original J2000 values so downstream tools report the same
catalog position as GICO3/frinZ.

For geometric delay tracking, the internal model currently precesses the J2000
coordinates to the observation date before evaluating geometric delay/rate. A
fixed-J2000 delay-model trial in `2.0.9` caused severe decorrelation on
YAMAGU32/USUDA64 validation data and was reverted in `2.0.10`.

```text
XML/header coordinates: J2000 catalog RA/Dec
Delay-model coordinates: precessed-to-date RA/Dec
```

#### `.cor` clock fields are metadata

The `.cor` file header carries station clock delay/rate from XML `<clock>` for
GICO3-compatible header display. The written complex XCF spectrum is already the
result of yi-corr correlation processing; `fringe`/`frinZ` residual delay/rate
searches operate on that complex spectrum and should not be interpreted as
applying clock metadata from the `.cor` header.

#### Process baseline order and closure

For schedule XML, station keys define the participating antennas. The optional
`<process><baseline>...</baseline></process>` list defines the baseline order.
For multiple stations, requested baselines must form closure over the process
station set. This keeps multi-baseline correlation products internally
consistent.

#### Large read-align delay uses zero padding, not shorter integration

When clock/geometric delay shifts one antenna by a large amount, the correlator
keeps the requested process length. Missing data outside the available raw
window are zero-filled. Therefore a 5-second process remains a 5-second `.cor`
sector even if only 4 seconds contain overlapping raw samples after read
alignment. This matches the GICO3 behavior where nonexistent signal components
are represented by zeros rather than shortening the integration time.

### Frequency Shift, Hermitian Completion, and XCF Phase Correction

GICO3-compatible frequency handling follows this conceptual order:

1. decode packed samples
2. normalize sideband (`LSB -> USB` when requested)
3. apply fringe/rotation terms
4. FFT
5. shift spectra onto the XML frequency grid
6. apply fringe-rotation phase-slope correction
7. integrate ACF/XCF products

The XML `<frequency>` value is the correlation reference frequency. Per-antenna
`<rotation>` shifts the sampled band onto this XML grid. For example, if
YAMAGU32 covers `8192..8704 MHz` and USUDA64 covers `8080..8592 MHz`, then
`USUDA64 rotation=-112 MHz` shifts the USUDA64 spectrum onto the same
`8192..8704 MHz` XML grid.

Because the input samples are real-valued before FFT, the FFT spectrum has
Hermitian symmetry:

```text
X[N-k] = conj(X[k])
```

`yi-corr` stores and outputs the positive-frequency side. When a frequency
shift asks for a bin outside the stored positive-frequency range, the missing
bin is reconstructed using this Hermitian relation.

For an output/XML grid bin `k_grid` and integer shift `r_bins`:

```text
k_raw = k_grid - r_bins
```

If `k_raw` is outside the positive-frequency half:

```text
X_raw[k_raw] = conj(X_raw[N - k_raw])
```

This is sufficient for ACF power spectra because:

```text
|X[k]|^2 = |conj(X[k])|^2
```

For XCF, however, phase matters:

```text
C12[k] = X1[k] * conj(X2[k])
```

Therefore the fractional-delay phase correction must use the same frequency
axis as GICO3. Since GICO3 applies phase-slope correction after frequency
shift, the bin index used for this correction is the shifted XML/grid bin, not
the raw signed FFT bin.

The XCF fractional-delay correction is:

```text
phase_corr[k] = exp{-j 2*pi*df*(k1_grid*eps1 - k2_grid*eps2)}
```

where:

- `df = sampling_rate / fft`
- `eps1`, `eps2` are the fractional residual delays for antenna 1 and 2
- `k1_grid`, `k2_grid` are bin indices after frequency shift onto the XML grid

This detail is important for the non-overlap/image part of a shifted real FFT.
Using raw signed bins after Hermitian completion gives correct-looking ACF
power but breaks XCF phase coherence in the image region. `fx-corr 2.0.4`
uses shifted XML/grid bins for XCF phase-slope correction.

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

<!-- yi-corr-v302-grid-origin-notes:start -->
### `3.0.2` validation update: HITACH32 frequency-grid origin is a Hz-level backend convention

Version `3.0.2` clarifies and generalizes the HITACH32 frequency-grid correction.  The correction must not be interpreted as a fixed `-1 bin` shift.  It is a small station/backend frequency-origin offset expressed in physical frequency units:

```text
HITACH32 frequency_grid_origin_offset_hz = -976.5625 Hz
```

The number of FFT bins corresponding to this offset depends on the FFT resolution:

```text
offset_bins = round(offset_hz / df_hz)
df_hz       = fs_hz / fft_len
```

For the YI data discussed here, `fs = 1024 MHz`.

At high spectral resolution:

```text
FFT = 1048576
df  = 1024e6 / 1048576 = 976.5625 Hz/bin

HITACH32 offset = -976.5625 Hz = -1 bin
```

At low spectral resolution:

```text
FFT = 8192
df  = 1024e6 / 8192 = 125000 Hz/bin

HITACH32 offset = -976.5625 Hz = -0.0078125 bin
round(...)      = 0 bin
```

Therefore the correct rule is not:

```text
HITACH32 offset = always -1 bin
```

but instead:

```text
HITACH32 offset = -976.5625 Hz
```

This distinction is essential.  A fixed `-1 bin` correction is appropriate for `FFT=1048576`, where one bin is 976.5625 Hz, but it would become a 125 kHz shift at `FFT=8192`.  Such a large artificial shift can decorrelate low-resolution continuum data.  The correction must therefore be stored as a station/backend frequency offset and converted to bins only after the FFT length is known.

#### How the correction was identified

The correction was first identified using the G9.62 maser scan.  The maser line is narrow, so the auto-correlation spectra provide a direct check of the frequency-bin placement.  In this scan, YAMAGU32 and HITACH32 used the same RF setup, the same bandwidth, the same nominal sideband, zero rotation, and the same output frequency grid.  Therefore the maser line should appear at the same `.cor` frequency bin after decoding and FFT-to-XML-grid mapping.

Before the correction, the HITACH32 ACF line appeared approximately one high-resolution bin lower than the YAMAGU32 line.  The offset followed HITACH32 when the baseline order was reversed, so it was not an `ant2`-specific read-align, delay, or phase-slope effect.  The same G9.62 data also included a same-format YAMAGU32--YAMAGU34 control baseline; that pair did not require this offset.  This separated the problem from a general LSB convention, `.cor` writer convention, or FFT endpoint convention, and pointed to a HITACH32 station/backend frequency-grid origin convention.

With the `3.0.2` Hz-based correction, G9.62 at `FFT=1048576` gives:

```text
YAMAGU32:
  station_offset = 0
  env_offset     = 0
  total_offset   = 0
  ACF peak       = 69865/69866

HITACH32:
  station_offset = -1
  env_offset     = 0
  total_offset   = -1
  ACF peak       = 69864/69865/69866
```

The remaining sector-to-sector one-bin scatter is expected for a narrow line sampled near a channel maximum with finite integration.  The important point is that HITACH32 is no longer systematically displaced by one bin relative to YAMAGU32.

The corrected G9.62 cross-power is recovered with frequency search:

```bash
frinZ --in <corfile> --freq --search deep --len 5
```

Typical result:

```text
Amp       ≈ 1431 %
SNR       ≈ 917
Frequency ≈ +68.2275391 MHz
```

This confirms that the corrected grid placement recovers the narrow-band maser signal.

#### Continuum validation at `FFT=8192`

Ordinary continuum fringe searches often use a shorter FFT.  For `FFT=8192`, the HITACH32 offset of `-976.5625 Hz` rounds to zero bins.  Thus low-resolution continuum processing is not forced to apply a spurious 125 kHz shift.

For NRAO530 at `FFT=8192`, yi-corr recovers a normal continuum fringe comparable to GICO3:

```text
NRAO530, FFT=8192, yi-corr:
  Amp ≈ 0.755 %
  SNR ≈ 733

NRAO530, FFT=8192, GICO3:
  Amp ≈ 0.724 %
  SNR ≈ 709
```

This explains why earlier low-resolution continuum tests did not appear to need a `-1 bin` correction: at this FFT length, the physically correct HITACH32 offset is smaller than one coarse channel and rounds to zero bins.

#### Continuum validation at `FFT=1048576`

VLBI delay searches may also use `FFT=1048576`.  In that case, the same `-976.5625 Hz` correction becomes exactly one bin.  This was tested with the continuum source NRAO530.

With the default HITACH32 correction enabled:

```text
NRAO530, FFT=1048576, HITACH32 correction enabled:
  Amp ≈ 0.759 %
  SNR ≈ 980.5
```

When the correction was deliberately canceled with:

```bash
YI_ANT2_GRID_OFFSET=1
```

the same scan dropped to:

```text
NRAO530, FFT=1048576, HITACH32 correction canceled:
  Amp ≈ 0.015 %
  SNR ≈ 19.3
```

This is the decisive validation.  The HITACH32 correction is not a maser-specific adjustment.  It is also required for high-resolution continuum processing.  Therefore it should be treated as station/backend metadata, not as a hand-tuned operation for a particular source.

#### Operational interpretation

The correct behavior is:

```text
HITACH32 correction_hz = -976.5625 Hz
offset_bins = round(correction_hz / (fs_hz / fft_len))
```

Consequently:

```text
FFT=1048576:
  df = 976.5625 Hz/bin
  offset_bins = -1

FFT=8192:
  df = 125000 Hz/bin
  offset_bins = 0
```

This resolves the apparent contradiction between the old continuum tests and the high-resolution maser/VLBI tests.  Low-resolution continuum processing did not need an apparent `-1 bin` correction because the physical offset rounded to zero bins.  High-resolution spectral-line and delay-search processing does need the correction because the same physical offset becomes one FFT bin.

#### Maser and continuum should be validated with different fringe-search modes

A narrow-band maser should be evaluated with frequency search enabled:

```bash
frinZ --in <corfile> --freq --search deep --len 5
```

If a narrow-band maser `.cor` file is evaluated only with a broad-band delay/rate search, the signal may appear weak or may be assigned an unphysical residual delay.  This is an analysis-mode issue, not a correlator failure.

For continuum sources such as NRAO530, the standard delay/rate search remains appropriate:

```bash
frinZ --in <corfile> --search deep
```

This distinction is important when validating the same correlator using both continuum calibrators and narrow-band maser targets.

#### Diagnostic overrides

The environment variables such as:

```text
YI_ANT2_GRID_OFFSET
```

are diagnostic overrides.  They are useful for controlled tests, for example canceling the built-in HITACH32 correction to prove that the correction is required.  They should not be treated as normal production parameters.

The production correction should be derived from station/backend metadata:

```text
station_grid_origin_offset_hz("HITACH32") = -976.5625 Hz
```

and converted to FFT bins internally.

#### Implication for future GPU or streaming implementations

Any future GPU implementation must preserve this correction in physical frequency units.  The correction must not be hard-coded as a fixed bin shift unless the FFT length is fixed.  The recommended metadata path is:

```text
station/backend name
  -> frequency_grid_origin_offset_hz
  -> df_hz = fs_hz / fft_len
  -> offset_bins = round(offset_hz / df_hz)
  -> FFT-to-output-grid mapping
```

This keeps the behavior consistent across low-resolution continuum processing, high-resolution continuum delay searches, and narrow-band maser processing.
<!-- yi-corr-v302-grid-origin-notes:end -->

## Compatibility Change Log

This section summarizes the GICO3-compatibility fixes made during the 1.1.x
series. Versions not listed here were intermediate builds or unrelated local
state.

| Version | Summary |
|---|---|
| `3.1.12` | Reused normalized `.cor` output buffers for normal `yi-corr` sector writes, avoiding per-sector/per-inband allocation during ACF/XCF/phased product serialization. |
| `3.1.11` | Reused `PackedSampleReader` instances in the `yi-corr` sector reader and skipped seeks when sector reads are byte-contiguous, reducing raw I/O open/seek overhead. |
| `3.1.10` | Precomputed real-FFT XML-grid bin maps for normal `yi-corr` accumulation, removing per-bin modulo/Hermitian mapping from the ACF/XCF hot loop while keeping arbitrary FFT lengths. |
| `3.1.9` | Fixed the normal `yi-corr` accumulation kernel from XML/CLI state before the frame loop and cached FFT debug mode outside the hot path. |
| `3.1.8` | Fused normal `yi-corr` FFT-bin accumulation so ACF/XCF update directly from mapped FFT bins without materializing XML-grid spectra unless pulsar folding or FFT debug needs them. |
| `3.1.7` | Reduced hot-loop allocation in `yi-corr` by reusing XML-grid FFT scratch buffers per worker and shifting bit-offset input reads in-place. |
| `3.1.6` | Grouped `--fringe <seconds>` quick-look outputs under `<cor>/fringe_yicorr_ql/<schedule-basename>/<process-epoch>/` so repeated process windows from one schedule stay together. |
| `3.1.5` | Changed `--fringe` to `--fringe <seconds>` quick-look output for `yi-xcf`/`yi-corr`: each interval writes integrated XCF spectrum PNG/TSV plus a lightweight 1D frequency-IFFT lag-fringe PNG/TSV. |
| `3.1.4` | Added XML `<stream><pulsar>...</pulsar></stream>` pulse-phase folding for `.cor` spectral products. Pulsar handling lives in `pulsar.rs`; optional DM correction assigns each frequency bin to a corrected pulse phase bin before ACF/XCF accumulation, and output files use `.pbinNN.cor` suffixes. |
| `3.1.3` | Made `yi-phasedarray` apply residual delay/fringe phase on the same XML-grid bins used by `yi-corr`, so phased raw synthesis differs from correlation only at the final beamforming/output step. |
| `3.1.2` | Added `--resacel <Hz/s>` to apply residual fringe acceleration from fringe fitting to the delay model, useful for long `yi-phasedarray` integrations when residual phase curvature remains. |
| `3.1.1` | Updated `yi-phasedarray` beamforming to use the same post-FFT XML-grid alignment as `yi-corr` before phased summing, so phased raw and phased `.cor` products follow the validated delay/grid tracking path. |
| `3.1.0` | Added XML `<stream><inband>N</inband></stream>` output splitting for pulsar and sub-band correlation workflows. `N` must be a power of two; split `.cor` files use `.chNbwBW.cor` names and carry sub-band observing frequency / effective bandwidth metadata. |
| `3.0.2` | Reinterpreted the HITACH32 grid-origin correction as a physical frequency offset (`-976.5625 Hz`) rather than a fixed `-1` bin shift. This gives `-1` bin at `FFT=1048576` but `0` bins at `FFT=8192`, preserving low-resolution continuum processing while fixing high-resolution maser and VLBI delay-search modes. Validated with G9.62 maser and NRAO530 continuum; canceling the correction at `FFT=1048576` reduced the NRAO530 SNR from about 980 to about 19. |
| `3.0.1` | Added validated station/backend frequency-grid-origin correction for HITACH32 (`-1` bin), ACF peak debug logging of station/env/total grid offsets, and documented the G9.62 maser validation. The same investigation showed that the YAMAGU32--YAMAGU34 factor-of-two amplitude loss was caused by a missing YAMAGU34 clock delay, not by ACF/XCF normalization. |
| `3.0.0` | Major delay-path stabilization release. Fixed the 0.5-sample boundary problem by keeping process/window read alignment fixed, avoiding residual re-splitting into per-frame integer shifts in the sector read-aligned path, and carrying the time-varying residual with the XCF phase slope. This removes sector-boundary phase jumps and one-sample residual-delay branch hopping seen by `fringe`/`frinZ`. |
| `3.0.1` | Added validated station/backend frequency-grid-origin correction for HITACH32 (`-1` bin), ACF peak debug logging of station/env/total grid offsets, and detailed G9.62 maser validation notes. The same investigation showed that the YAMAGU32--YAMAGU34 factor-of-two amplitude loss was caused by a missing YAMAGU34 clock delay (`+1.72 us`), not by ACF/XCF normalization. |
| `1.1.1` | Core GICO3-compatible semantics were stabilized: station-key baseline order, multi-baseline closure, reverse default shuffle, XML clock epoch as delay reference, requested-length zero padding for large read-align delays, `.cor` output generation, and ACF/XCF file naming. |
| `1.1.8` | Added explicit logging for read-align residual delay and the residual-delay correction target antenna. |
| `1.1.9` | Trial change to delay phase reference. This reduced SNR on validation data and was not kept as the final direction. |
| `1.1.21` | Added GICO3-style frequency-grid handling: XML `<frequency>` is the output reference, `<rotation>` shifts each sampled band onto that XML grid, and post-FFT bin shifting with Hermitian completion was introduced. ACF became GICO3-like, but XCF image-region phase was still wrong. |
| `1.1.22` | Fixed rotation/fringe LO interpretation: carrier phase uses the raw sampled band frequency, while `.cor` observing frequency remains the XML grid frequency. XCF amplitude recovered. |
| `1.1.23` | Trial XCF fractional-delay correction using raw signed FFT bins. This later proved wrong for the Hermitian-completed image region. |
| `1.1.24` | All correlation is USB-normalized; LSB-to-USB sign alternation now uses absolute raw-sample parity instead of per-FFT-block local parity. |
| `1.1.25` | `.cor` file headers now carry station clock delay/rate from XML `<clock>` metadata for GICO3-compatible headers. This metadata does not change the written XCF spectrum. |
| `1.1.26` | Trial source-coordinate-frame change during validation. |
| `1.1.27` | Reverted the `1.1.26` source-frame trial while other delay-model checks were still in progress. |
| `1.1.28` | Trial residual-delay sign change. Validation showed lower SNR on YAMAGU32/USUDA64, so this was reverted. |
| `1.1.29` | Reverted the `1.1.28` residual-delay sign trial and restored the higher-SNR convention. |
| `1.1.30` | Fixed XCF fractional-delay phase correction to use shifted XML/grid bins after frequency shift, matching GICO3 order. |
| `2.0.0` | Baseline release after GICO3-compatible frequency/rotation/sideband/clock semantics, Hermitian completion, and XCF grid-bin phase correction were consolidated. Superseded by `2.0.1` because intermediate local rebuilds reused the same version string while performance fixes were still in progress. |
| `2.0.1` | Disabled the continuous-decode-offset experiment and restored scratch-buffer reuse, but still performed per-bin trigonometric XCF phase correction and could stall before the first sector completed. |
| `2.0.2` | Invalid performance recovery attempt; do not use for validation. |
| `2.0.3` | Reverted runtime source to the known `src20260602B/src` integration path after `.cor` sector output stalled in 2.0.0-2.0.2. |
| `2.0.4` | Split rotation semantics: `<rotation>` uses only real overlapping signal bins for XCF, while `<rotation2>` opts into Hermitian-completion/image-band correlation. |
| `2.0.5` | Residual read-align delay is corrected fully in the frequency domain; frame-local integer zero-fill shifts are no longer used for the small post-seek residual. |
| `2.0.6` | Versioned install now reads the package version from `Cargo.toml` at runtime, so `cargo fx-corr` creates the new `yi-corr-vX.Y.Z` even while replacing an older `cargo-fx-corr`. |
| `2.0.7` | Reverted the full frequency-domain residual-delay experiment after validation decorrelated YAMAGU32/USUDA64; keeps absolute-sample LSB parity handling. |
| `2.0.8` | XCF fractional-delay phase correction now uses each antenna raw/physical FFT bin (`xml_bin - rotation_bins`) while output indexing remains on the XML grid. |
| `2.0.9` | Trial fixed-J2000 delay-model change. Validation decorrelated YAMAGU32/USUDA64 and moved the residual delay to thousands of samples; do not use for validation. |
| `2.0.10` | Reverted the `2.0.9` fixed-J2000 delay-model trial and restored precessed-to-date geometric delay tracking while keeping J2000 coordinates in `.cor` headers. |
| `2.0.11` | Trial GMST/GAST coordinate-frame change. Validation worsened YAMAGU32/USUDA64 residual delay/rate, so this was rejected. |
| `2.0.12` | Reverted the `2.0.11` GMST/GAST trial and restored the previous geometric model. Focus returns to raw read/decode ordering and delay application inside the correlator. |
| `2.0.13` | XCF paths now decode the integer-delay-shifted FFT window directly from the packed raw chunk instead of applying per-frame zero-fill shifts after decode. This removes artificial frame-boundary discontinuities from residual read-align correction. |
| `2.0.14` | Trial explicit-shuffle parse change. Validation decorrelated YAMAGU32/USUDA64, so this was rejected. |
| `2.0.15` | Reverted the `2.0.14` shuffle trial and restored the previous external-to-internal shuffle conversion. |
| `2.0.16` | Added per-frame delay diagnostics under `--debug`: full relative delay, post-read-align residual, antenna tau values, integer shifts, and fractional-delay terms are now printed for sampled frames and written for every frame in the debug log. |
| `2.0.17` | XCF phase-bin diagnostics are now printed in normal logs whenever XCF products are written, so rotation raw-vs-XML bin-index mistakes can be checked without debug mode. |
| `2.0.18` | Trial XCF phase-start change to the XML grid origin. Validation did not improve YAMAGU32/USUDA64 residual delay/rate, so it was reverted in `2.0.19`. |
| `2.0.19` | Normal yi-corr logs now include sampled per-frame delay-correction values for the first and last XCF sectors, making residual, antenna tau, integer shift, and fractional terms visible without debug mode. |
| `2.0.20` | Per-frame geometric delay-rate residual is no longer converted into changing integer raw-window shifts. The integer raw shift is fixed at the initial post-read-align residual, while the time-varying residual is carried by the XCF frequency-domain phase slope. |
| `2.0.30` | Rejected delay/frequency-path experiment: direct RF-overlap bin selection, antenna-based geometric delay table, positive residual convention, and positive-baseband XCF phase. Validation collapsed YAMAGU32/USUDA64 XCF, so this path was not kept. |
| `2.0.31` | Rejected sector-reader trial applied on top of the bad `2.0.30` frequency/delay path; produced only ~896 overlap bins and decorrelated the XCF. |
| `2.0.32` | Restored the stable `2.0.19` rotation/delay/frequency path, then added sector-based bit-exact raw reads. Each 1-second sector now re-computes integer read-align from `full_rel` while the frame delay model uses the same sector `d_seek`, preventing geometric-rate integer shifts from accumulating across long scans. |
| `2.0.33` | Rejected GMST trial. It worsened YAMAGU32/USUDA64 residual delay and rate, so the apparent-sidereal geometric path was restored. |
| `2.0.34` | Rejected no-fringe-rate-stopping trial. Leaving geometric rate entirely in the XCF phase decorrelated 10-second integrations on YAMAGU32/USUDA64, so this mode is not suitable as the default correlator path. |
| `2.0.35` | Restored the high-SNR `2.0.32` sector read-align and per-frame geometric delay correction after the rejected `2.0.33`/`2.0.34` trials. The remaining sign investigation should target fractional phase/carrier phase, not the integer raw read-align direction. |
| `2.0.36` | Rejected trial: reversing only the fractional-delay XCF phase-slope sign kept residual delay/rate nearly unchanged and cut YAMAGU32/USUDA64 amplitude roughly in half. |
| `2.0.37` | Rejected trial: reversing only carrier/fringe `fr_lo` phase sign fully decorrelated YAMAGU32/USUDA64 XCF. |
| `2.0.38` | Trial: restore all phase signs from the high-SNR path and change sector read-align reference from the first frame of each 1-second sector to the sector midpoint. Validation showed essentially identical residual delay/rate, so the sector reference point is not the cause. |
| `2.0.39` | Confirmed time-reference trial: evaluating the delay model 1.0 second later moved residual delay by about 94.6 samples while preserving high SNR, proving the remaining offset is a delay-model time reference issue. |
| `2.0.40` | Confirmed tuned model-time offset: +0.713 s reduced YAMAGU32/USUDA64 residual delay to ~0.03 sample and residual rate to ~1 mHz while preserving high SNR. |
| `2.0.41` | Promoted the +0.713 s model-time offset to a named constant, extended the geometric-delay table by the offset duration to avoid endpoint clamping, and logs the applied model-time offset. |
| `2.0.42` | Rejected trial: removing sector integer `d_seek` from carrier/fringe phase preserved delay but made residual rate much worse, so the absolute carrier phase must use full model delay. |
| `2.0.43` | Restored the confirmed `2.0.41` full-model-delay carrier phase after rejecting `2.0.42`. |
| `2.0.44` | Fixed exact-duration bookkeeping: requested integer-second windows now round to the nearest FFT-frame count when appropriate, and debug end-time display no longer floors 1799.999999-style floating-point durations to the previous second. |
| `2.0.45` | Reduced progress-log spam by throttling Synth I/O reader reports to 10-second intervals and terminating the carriage-return progress line before the final I/O summary. |
| `2.0.46` | Disabled intermediate Synth I/O reader reports in normal logs; progress now stays on one carriage-return line and only the final I/O summary is printed. |
| `2.0.47` | Added first/last-sector XCF phase-reference diagnostics showing model time, XML reference frequency, raw/XML phase bins, and fringe-mix phase so phase-origin drift can be separated from delay/rate residuals. |
| `2.0.48` | Rejected trial: changing carrier/fringe phase reference from raw sampled-band low frequency to XML-grid low frequency severely reduced YAMAGU32/USUDA64 SNR and reintroduced a large residual rate in deep search. |
| `2.0.49` | Restored the confirmed raw sampled-band low-frequency carrier/fringe phase reference after rejecting `2.0.48`. |
| `2.0.50` | Reorganized schedule stdout for multi-baseline closure runs: common schedule/stream parameters, baseline delay/band summary, and compact antenna parameters are printed separately to reduce repeated wide tables. |
| `2.0.51` | Unified schedule stdout format for single- and multi-baseline runs; all runs now use Schedule/Stream, Baseline, and Antenna parameter sections with compact long fields. |
| `2.0.52` | Replaced ASCII bordered schedule/baseline/antenna parameter tables with plain sectioned key/value output, with antenna parameters grouped per antenna. |
| `2.0.53` | Added `--model-time-offset` so the empirical delay-model time reference can be swept at runtime, and logs its geometric-rate sample equivalent for diagnosing the 0.713 s offset. |
| `2.0.54` | Added EOP-aware geometric-delay controls (`--dut1`, `--tt-utc`, `--xp`, `--yp`) and uses separate UT1/TT epochs plus optional polar motion in the delay model for long-baseline GICO3 comparisons. |
| `2.0.55` | Moved EOP and delay-model time-offset control into schedule XML (`<eop>`, `<delay-model>`) for XML-driven correlation runs; CLI values are now only non-XML fallbacks. |
| `2.0.56` | Removed empirical delay-model time-offset defaults: schedule XML and non-XML fallback now default to time-offset=0.0 unless explicitly specified. |
| `2.0.57` | Rejected trial: full IAU 2006/2000A J2000-to-ITRS vector transform moved YAMAGU32/USUDA64 geometric delay by about -32 samples and worsened the residual delay. |
| `2.0.58` | Reverted the rejected `2.0.57` source-vector transform and restored the previous precessed-to-date apparent-sidereal delay path for continued timing investigation. |

Practical validation notes:

- ACF can look correct even when XCF image-region phase is wrong, because ACF uses power and is insensitive to complex conjugation phase.
- XCF must preserve phase after frequency shift; otherwise non-overlap/image-band bins lose coherence and SNR drops.
- `fringe`/`frinZ` residual delay/rate are measured from the complex spectrum itself. They should not be treated as corrections read from `.cor` clock metadata.
- For connected-element interferometer data, a non-zero residual rate can also come from applying VLBI geometric fringe tracking to data that should be correlated with zero geometric rate.

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
  --sc test.xml \
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
