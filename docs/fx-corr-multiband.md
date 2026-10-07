# Two-band correlation and calibrator transfer

`yi-corr` 3.9.0 can jointly use the simultaneous 6600–7112 MHz and
8192–8704 MHz observations. Each station/band has its own RAW file.
The occupied bandwidth is 512 + 512 = 1024 MHz; the unobserved gap is
1080 MHz and the RF extent is 2104 MHz. Each band retains its actual RF
coordinates. The reference frequency for these bands is 7652 MHz.

The intended 110 m baseline workflow is to estimate residual delay/rate and
instrumental phase on a strong calibrator, then apply those solutions across
the observation. Target noise is never fringe-fitted in this mode. Geometric
delay/fringe tracking still uses each target's XML coordinates and epoch.
The residual calibration assumes the same instrumental setup and sufficiently
stable propagation and clock errors over the transfer interval.

For independent noise, equal per-band sensitivity, the same coherent
integration time, and usable full bands, doubling the bandwidth lowers thermal
noise to `1/sqrt(2) = 0.707` of one band and increases S/N by `sqrt(2)`.
This follows the [NRAO baseline sensitivity equation](https://science.nrao.edu/facilities/vlba/docs/manuals/oss2022B/bsln-sens).
Longer coherent integration can further lower the noise if the transferred
phase model remains valid. RFI, unequal noise, lost common channels, gain
differences, source spectra, or calibration errors change the achieved gain.

## Run with a strong calibrator

Prepare two XML schedules with the same ordered baseline, source/process
entries, epochs, skips, lengths, sampling, FFT and final output rate. Set the
first XML `<frequency>` to `6600000000`, the second to `8192000000`.
For each 512 MHz real-sampled band, terminal speed is `1024000000` sample/s.
Keep each band's actual station clock, sideband, shuffle and rotation values.
Each XML uses its first stream, as in the existing single-band workflow.
Only `inband=1`, unfolded two-station correlation is supported here.

```bash
target/release/yi-corr \
  --sc low.xml --raw raw/low --cor output \
  --multiband-schedule high.xml \
  --multiband-raw-directory raw/high \
  --multiband-calibrator CALIBRATOR_NAME \
  --cpu 16
```

Replace `CALIBRATOR_NAME` with the XML `<object>` spelling. Multiple names
can be comma separated. The calibrator scans are processed first, even when
they occur later in the XML or `--process-index` selects only a target.
Directory/file naming follows the normal station/epoch RAW lookup.
For explicitly named files, use `--ant1/--ant2` for the low band and
`--multiband-ant1/--multiband-ant2` for the high band; `--raw` is still required
by the existing CLI. Explicit paths refer to the input of every selected scan,
so use directories for observations with different files per scan.

The first calibrator solution defaults to `--multiband-phase-mode auto`:
fit common delay/rate using the sum of the independent band powers, estimate
their relative constant IF phase, then keep that phase offset fixed for
subsequent phase-connected solutions. This first window obtains delay
information from the individual band widths; unknown independent IF phases
cannot provide the RF-gap delay precision until they are calibrated.
If phase-cal has already removed IF phase offsets, select
`--multiband-phase-mode connected`. Instrumental IF offsets must be removed
before combining IF phases, as described by the
[AIPS calibration guidance](https://www.aips.nrao.edu/CookHTML/CookBookse138.html).

Short calibrator correlations use 0.05 s integrations by default, with
10 s solution windows. Both bands contribute to the common fit. The final
pass rereads calibrator RAW and applies corrections to each FFT cross product
before the requested XML integration; target RAW is read only for the final
pass. Long final integrations therefore retain the phase correction that
would be lost by averaging uncorrected fringe rates first.
FFT plans and search workspaces are reused, rate trials run in parallel,
frequency-gap bins have zero weight, and only a solution window is buffered.

Targets use the nearest calibrator solution's linear model, extrapolated to
the target time. One calibrator scan can thus cover the entire observation.
With several scans/windows the nearest model is selected by reference time;
no interpolation or global acceleration fit is performed. A transition can
show a phase jump if adjacent calibrator models disagree. Add
`--multiband-calibration-max-gap SECONDS` to impose a maximum transfer distance.
The default has no time limit. `calibration.tsv` records all calibrator models
on Unix time; each target's `solutions.tsv` records the selected models in
seconds since its process epoch, including extrapolated reference times.

## Search controls and limits

| Option | Default | Meaning |
|---|---:|---|
| `--multiband-window` | 10 s | Minimum solution span; the last window includes a short remainder |
| `--multiband-solve-integration` | 0.05 s | Requested integration for calibrator/solution-pass data |
| `--multiband-delay-window-ns` | 100 ns | Symmetric residual delay search half-width |
| `--multiband-rate-window-hz` | 0.25 Hz | Symmetric residual fringe-rate half-width at the reference RF |
| `--multiband-min-coherence` | 0.05 | Minimum fitted coherent amplitude divided by summed visibility magnitudes |

Coherence is a diagnostic and rejection threshold, not a false-alarm
probability or a measured detection S/N. Choose strong calibrators and inspect
their solutions. Increase the search windows if a solution reaches a boundary.
Short integration must resolve the searched rate, with at least two rows in a
window. The delay range must stay below the channel-spacing ambiguity
`1/(2*df)`. RF separation must be an integer multiple of `df`; FFT 8192 at
1.024 GHz satisfies this for the two observation bands. Separated bands have
delay ambiguity sidelobes: a sufficient calibrator S/N and a realistic residual
delay window are necessary to choose the correct peak.

The correction is a nondispersive common delay/rate plus constant per-band
phase. It does not fit ionospheric dispersion, bandpass phase ripple,
acceleration, RFI flags or time-variable IF offsets. Residual calibration is
applied spectrally after the existing RAW sample alignment. Keep residual
delays (including drift when extrapolating) small compared with the FFT time
span; a post-FFT correction cannot recover sample overlap already lost in
misaligned FFT windows. Repair larger clock errors in the XML before fitting.
`--gain-phasecal`, phased validation, model sweep, fringe QL, `--band`, and
pulsar folding cannot be combined with this workflow. `yi-phasedarray` continues
to operate on individual bands.

Omitting `--multiband-calibrator` instead fits every scan's joint data, including
its own initial IF-phase bootstrap. That mode is for sources with sufficient
S/N; use calibrator transfer for faint targets.

## Products and coherent averaging

Each `output/multiband/scanNNNN` contains:

- `solve-band1/2`: short uncorrected XCF `.cor` files for scans that were fitted.
- `solutions.tsv`: applied residual delay/rate/common phase and fixed IF phases.
- `band1/2`: corrected native `.cor` products, compatible with existing readers.
- `joint.mbcor`: both corrected cross spectra, original frequency grids and time headers.

The native files and joint file deliberately preserve spectral information.
The existing `.cor` header cannot represent a gap between two IFs; ordinary
`.cor`/frinZ readers cannot directly read `.mbcor`. Export with the included
streaming Python tool (standard library only):

```bash
python3 tools/multiband_visibility.py \
  output/multiband/scan0001/joint.mbcor --average combined.tsv

# Further coherent time averaging within this scan:
python3 tools/multiband_visibility.py \
  output/multiband/scan0001/joint.mbcor --average scan-mean.tsv --time-average

# All spectral samples with their physical RF frequencies:
python3 tools/multiband_visibility.py \
  output/multiband/scan0001/joint.mbcor --spectra spectra.tsv
```

The complex mean defaults to equal per-channel noise weights. It combines
real and imaginary components before computing amplitude. Do not average
visibility magnitudes to detect faint sources: their noise has a positive bias.
`--band-scales SCALE1 SCALE2` supplies amplitude calibration before averaging;
`--band-weights W1 W2` supplies inverse per-channel noise variances **after**
those scales. For time averaging, weights also include integration duration.
The mean is a continuum statistic for a phase-centred source; resolved sources,
offset positions and spectral features should retain their per-channel data.
Exact zero bins, used by yi-corr for channels outside the common station band,
are excluded; `usable_bandwidth_hz` reports the included width. No automatic
RFI rejection, bandpass calibration or measured noise weighting is implied.
TSV times retain the native `.cor` timestamp convention.

### `.mbcor` version 1 layout

All numeric fields are little endian. Fixed header is 544 bytes:

| Offset | Type | Value |
|---:|---|---|
| 0 | 8 bytes | `YIMBCOR\0` |
| 8 | u32 | version 1 |
| 12 | u32 | two bands |
| 16 | f64 | occupied bandwidth in Hz |
| 24 | f64 | reference RF in Hz |
| 32 | 256 bytes | original low-band `.cor` file header |
| 288 | 256 bytes | original high-band `.cor` file header |

Each time record contains the low-band 128-byte `.cor` sector header and its
`Nfft/2` complex f32 channels, followed by the high-band header and channels.
The two time grids must match exactly. RF for channel `k` is the corresponding
embedded header frequency plus `k * sampling/Nfft`. Native sector counts and
FFT sizes determine total size. Joint output is written to `.part`, validated
and renamed only after successful completion.

## Validation

Unit tests use the actual 512 MHz observation bands and RF gap, injected
nonzero delay/rate, unknown IF phase and transfer between epochs. The RAW
end-to-end test introduces a clock error, fits a calibrator, transfers to a
different target epoch without a target solution pass, checks corrected complex
phase and joint RF metadata, and exercises both direct and mapped FFT kernels
plus the exporter. Real dual-band observation validation remains necessary to
measure the achieved sensitivity and calibration stability.
