use crate::utils::DynError;
use clap::Parser;
use std::f64::consts::PI;
use std::path::PathBuf;

const K_BOLTZMANN: f64 = 1.380_649e-23; // J/K
const SI_TO_JY: f64 = 1.0e26; // multiply W/m^2/Hz to express in Jy
pub const DEFAULT_FFT: usize = 16384;
pub const DEFAULT_MODEL_TIME_OFFSET_S: f64 = 0.0;
pub const DEFAULT_SHUFFLE_IN: [usize; 32] = [
    31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8,
    7, 6, 5, 4, 3, 2, 1, 0,
];

#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about = "yi toolkit for FX correlation and phased-array processing",
    long_about = "yi-corr: auto + cross correlation\nyi-phasedarray: phased-array synthesis and spectrum plots",
    infer_long_args = true,
    arg_required_else_help = true
)]
pub struct Args {
    // Schedule
    #[arg(
        long,
        visible_alias = "sc",
        value_name = "XML",
        help = "Observation schedule XML file (.xml)"
    )]
    pub schedule: Option<PathBuf>,

    #[arg(long, help = "Write example schedule XML to example.xml and exit")]
    pub mkxml: bool,

    // Input files
    #[arg(
        long = "raw-directory",
        visible_alias = "data",
        visible_alias = "raw",
        value_name = "DIR",
        help = "Directory containing input raw files"
    )]
    pub raw_directory: Option<PathBuf>,

    #[arg(long, value_name = "FILE", help = "Input raw file for antenna 1")]
    pub ant1: Option<PathBuf>,

    #[arg(long, value_name = "FILE", help = "Input raw file for antenna 2")]
    pub ant2: Option<PathBuf>,

    // Output files
    #[arg(
        long = "cor-directory",
        visible_alias = "output",
        visible_alias = "cor",
        value_name = "DIR",
        help = "Directory for output .cor/.raw products"
    )]
    pub cor_directory: Option<PathBuf>,

    #[arg(
        long,
        value_name = "NAME",
        default_value = "YAMAGU66",
        help = "Virtual-station name used by yi-phasedarray output"
    )]
    pub phased_name: String,

    #[arg(
        long,
        default_value_t = false,
        help = "Also write phased/input ACFs, input-pair XCF, and diagnostic plots"
    )]
    pub phased_diagnostics: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Run unattended three-station phased validation and write a complete NPZ"
    )]
    pub phased_validation: bool,

    #[arg(
        long,
        value_name = "NPZ",
        help = "Output NPZ path for --phased-validation"
    )]
    pub phased_validation_npz: Option<PathBuf>,

    #[arg(
        long,
        default_value_t = false,
        help = "Automatically solve gain-source Y32-Y34 group delay and phase with frinZ peak+acel, update L clock, and synthesize all scans"
    )]
    pub gain_phasecal: bool,

    #[arg(
        long,
        value_name = "SOURCE",
        help = "Schedule process/object name identifying the gain calibrator"
    )]
    pub gain_source: Option<String>,

    #[arg(
        long,
        value_name = "SECTORS",
        default_value_t = 1,
        help = "frinZ acel fitting window in correlator sectors; each gain scan must be divisible by it"
    )]
    pub gain_fringe_length: usize,

    #[arg(
        long,
        value_name = "XML",
        help = "Generated schedule with gain-derived L clock (default: OUTPUT/SCHEDULE-STEM_CORRECTED-KEY_phasecal.xml)"
    )]
    pub phasecal_schedule_output: Option<PathBuf>,

    #[arg(
        long,
        visible_alias = "gain-uncalibrated-sc",
        value_name = "XML",
        help = "Uncalibrated gain schedule used to write Y32-Y34 phasecal-off/on comparison .cor files"
    )]
    pub gain_uncalibrated_schedule: Option<PathBuf>,

    #[arg(
        long,
        value_name = "DIR",
        help = "Directory for gain phasecal-off/on .cor files (default: OUTPUT/gain_correlation)"
    )]
    pub gain_correlation_directory: Option<PathBuf>,

    #[arg(
        long,
        value_name = "KEY",
        default_value = "K",
        help = "Reference-station key which must remain unchanged in gain comparison XMLs"
    )]
    pub gain_reference_key: String,

    #[arg(
        long,
        value_name = "KEY",
        default_value = "L",
        help = "Station key receiving the gain delay calibration"
    )]
    pub gain_corrected_key: String,

    // Internal orchestration controls. They are never exposed on the command
    // line, but let yi-phasedarray reuse the yi-corr path without writing the
    // two ACF products or colliding with the normal stream label.
    #[arg(skip)]
    pub xcf_only: bool,

    #[arg(skip)]
    pub cor_label_override: Option<String>,

    // Correlation processing parameters
    #[arg(
        long,
        value_name = "HMS",
        help = "Source right ascension (HHhMMmSS.s or degrees)"
    )]
    pub ra: Option<String>,

    #[arg(
        long,
        value_name = "DMS",
        allow_hyphen_values = true,
        help = "Source declination (DDdMMmSS.s or degrees)"
    )]
    pub dec: Option<String>,

    #[arg(
        long,
        value_name = "YYYY/DDD HH:MM:SS",
        help = "Observation epoch (UTC)"
    )]
    pub epoch: Option<String>,

    #[arg(
        long,
        value_name = "MHZ",
        default_value_t = 0.0,
        help = "Reference observing frequency [MHz]"
    )]
    pub obsfreq: f64,

    #[arg(
        long,
        value_name = "MSPS",
        default_value_t = 1024.0,
        help = "Sampling rate [Msps]"
    )]
    pub sampling: f64,

    #[arg(long, value_name = "POINTS", default_value_t = DEFAULT_FFT, help = "FFT size [samples/frame]")]
    pub fft: usize,

    #[arg(
        long,
        value_name = "START:END",
        help = "Correlate/output only this band within the observed band [MHz], e.g. 0:8"
    )]
    pub band: Option<String>,

    #[arg(
        long = "bin",
        visible_alias = "bit",
        value_name = "BITS",
        num_args = 1..,
        help = "Quantization bits per sample [bits] (e.g. \"2\" or \"ant1:2 ant2:4\")"
    )]
    pub bin: Vec<String>,

    #[arg(
        long,
        value_name = "LEVEL",
        allow_negative_numbers = true,
        num_args = 1..,
        help = "Quantization levels [a.u.] (e.g. \"ant1:-1.5,0.5,-0.5,1.5\" or space-separated values)"
    )]
    pub level: Vec<String>,

    #[arg(
        long = "shuffle",
        value_name = "BIT",
        num_args = 1..=32,
        help = "Bit shuffle map [bit index] (e.g. \"ant1:31,30...0\" or 32 space-separated indices)"
    )]
    pub shuffle_in: Vec<String>,

    #[arg(long, value_name = "LSB|USB", num_args = 1.., help = "Input sampler sideband")]
    pub sideband: Vec<String>,

    #[arg(
        long,
        value_name = "S",
        allow_hyphen_values = true,
        help = "Fixed coarse delay [s]"
    )]
    pub coarse: Option<f64>,

    #[arg(
        long,
        value_name = "SAMPLES",
        allow_hyphen_values = true,
        default_value_t = 0.0,
        help = "Residual delay applied to ant2 [samples]"
    )]
    pub delay: f64,

    #[arg(
        long,
        value_name = "HZ",
        default_value_t = 0.0,
        allow_hyphen_values = true,
        help = "Residual fringe rate correction [Hz]"
    )]
    pub rate: f64,

    #[arg(
        long,
        value_name = "S",
        allow_hyphen_values = true,
        default_value_t = DEFAULT_MODEL_TIME_OFFSET_S,
        help = "Non-XML fallback: delay-model evaluation time offset from frame midpoint [s]"
    )]
    pub model_time_offset: f64,

    #[arg(
        long,
        value_name = "S",
        allow_hyphen_values = true,
        default_value_t = 0.0,
        help = "Non-XML fallback: Earth orientation DUT1 = UT1-UTC [s]"
    )]
    pub dut1: f64,

    #[arg(
        long = "tt-utc",
        value_name = "S",
        allow_hyphen_values = true,
        default_value_t = 69.184,
        help = "Non-XML fallback: TT-UTC used for IAU precession/nutation/sidereal time [s]"
    )]
    pub tt_utc: f64,

    #[arg(
        long,
        value_name = "ARCSEC",
        allow_hyphen_values = true,
        default_value_t = 0.0,
        help = "Non-XML fallback: polar motion x coordinate [arcsec]"
    )]
    pub xp: f64,

    #[arg(
        long,
        value_name = "ARCSEC",
        allow_hyphen_values = true,
        default_value_t = 0.0,
        help = "Non-XML fallback: polar motion y coordinate [arcsec]"
    )]
    pub yp: f64,

    #[arg(
        long,
        value_name = "SAMPLES",
        default_value_t = 0.0,
        allow_hyphen_values = true,
        help = "Additional residual delay [samples]"
    )]
    pub resdelay: f64,

    #[arg(
        long,
        value_name = "HZ",
        default_value_t = 0.0,
        allow_hyphen_values = true,
        help = "Additional residual rate [Hz]"
    )]
    pub resrate: f64,

    #[arg(
        long,
        value_name = "HZ/S",
        default_value_t = 0.0,
        allow_hyphen_values = true,
        help = "Additional residual fringe acceleration [Hz/s]"
    )]
    pub resacel: f64,

    #[arg(long, value_name = "HZ", num_args = 1.., allow_hyphen_values = true, help = "Per-antenna rotation frequency [Hz] (e.g. \"ant1:343000000 ant2:0\")")]
    pub rotation: Vec<String>,

    #[arg(
        long,
        value_name = "S",
        default_value_t = 0.0,
        help = "Start offset from input head [s]"
    )]
    pub skip: f64,

    #[arg(long, value_name = "S", help = "Processing duration [s]")]
    pub length: Option<f64>,

    // Phased-array analysis parameters
    #[arg(long, value_name = "K", num_args = 1.., help = "System temperature [K]")]
    pub tsys: Vec<String>,

    #[arg(long, value_name = "M", num_args = 1.., help = "Antenna diameter [m]")]
    pub diameter: Vec<String>,

    #[arg(long, value_name = "RATIO", num_args = 1.., help = "Antenna aperture efficiency [0..1]")]
    pub eta: Vec<String>,

    #[arg(long, value_name = "SCALE", num_args = 1.., help = "Per-antenna amplitude scale [a.u.]")]
    pub gain: Vec<String>,

    #[arg(long, value_name = "JY", num_args = 1.., help = "SEFD [Jy]")]
    pub sefd: Vec<String>,

    #[arg(
        long,
        value_name = "S",
        help = "Write quick-look XCF spectrum and 1D lag-fringe plots every S seconds"
    )]
    pub fringe: Option<f64>,

    // Other
    #[arg(
        long,
        value_name = "N",
        help = "Number of compute threads [default: logical CPU count - 2, min 1]"
    )]
    pub cpu: Option<usize>,

    #[arg(
        long,
        value_name = "N",
        help = "Reader chunk size [frames] (auto when omitted)"
    )]
    pub chunk_frames: Option<usize>,

    #[arg(
        long = "razoku5bay",
        help = "Force small I/O chunking at 0.5-second intervals (overrides --chunk-frames)"
    )]
    pub razoku5bay: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Read the two antenna inputs concurrently for USB-attached storage; unrelated to USB signal sideband"
    )]
    pub usb: bool,

    #[arg(
        long,
        value_name = "N",
        help = "Reader-worker queue depth [chunks] (auto when omitted)"
    )]
    pub pipeline_depth: Option<usize>,

    #[arg(long, help = "Enable diagnostic debug logging to file (all frames)")]
    pub debug: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Write one-pass full geometric model comparison diagnostics without changing correlation numerics"
    )]
    pub model_diagnostics: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Run an unattended five-case full-correlation model sweep into subdirectories"
    )]
    pub model_sweep: bool,

    #[arg(
        long = "stdout",
        default_value_t = false,
        help = "Write yi-corr runtime stdout log to ./stdout/stdout_<yyyydddhhmmss>.log"
    )]
    pub stdout: bool,

    #[arg(long, hide = true, value_name = "N")]
    pub process_index: Option<usize>,

    #[arg(long, hide = true)]
    pub ant1_name_override: Option<String>,

    #[arg(long, hide = true)]
    pub ant2_name_override: Option<String>,

    #[arg(long, hide = true, default_value = "ant2", value_parser = clap::builder::PossibleValuesParser::new(["ant1", "ant2"]))]
    pub delay_reference: String,

    #[arg(long, hide = true, default_value_t = false)]
    pub compact_logs: bool,
}

pub fn resolve_per_antenna_config<T: Clone>(
    args_list: &[String],
    default_val: T,
    parser: impl Fn(&str) -> Result<T, DynError>,
) -> Result<(T, T), DynError> {
    resolve_per_antenna_config_with_defaults(args_list, default_val.clone(), default_val, parser)
}

pub fn resolve_per_antenna_config_with_defaults<T: Clone>(
    args_list: &[String],
    default_ant1: T,
    default_ant2: T,
    parser: impl Fn(&str) -> Result<T, DynError>,
) -> Result<(T, T), DynError> {
    let mut ant1_val = default_ant1;
    let mut ant2_val = default_ant2;
    for entry in args_list {
        let entry = entry.trim();

        // Important for comma-separated values such as shuffle maps and level maps:
        //
        //   ant2:0,1,2,...,31
        //
        // must be passed to parser() as the full string "0,1,2,...,31".
        // Do not split this by comma before stripping the ant2: prefix.
        if let Some(val_str) = entry.strip_prefix("ant1:") {
            ant1_val = parser(val_str.trim())?;
            continue;
        }
        if let Some(val_str) = entry.strip_prefix("ant2:") {
            ant2_val = parser(val_str.trim())?;
            continue;
        }

        if entry.contains("ant1:") || entry.contains("ant2:") {
            let parts: Vec<&str> = entry
                .split(|c: char| c == ' ')
                .filter(|s| !s.is_empty())
                .collect();
            for part in parts {
                if let Some(val_str) = part.strip_prefix("ant1:") {
                    ant1_val = parser(val_str.trim())?;
                } else if let Some(val_str) = part.strip_prefix("ant2:") {
                    ant2_val = parser(val_str.trim())?;
                } else {
                    let v = parser(part.trim())?;
                    ant1_val = v.clone();
                    ant2_val = v;
                }
            }
        } else {
            let v = parser(entry)?;
            ant1_val = v.clone();
            ant2_val = v;
        }
    }
    Ok((ant1_val, ant2_val))
}

pub fn parse_levels(bit: usize, list: &str) -> Result<Vec<f64>, DynError> {
    if bit == 0 {
        return Err("Bit depth must be at least 1".into());
    }
    let expected = 1usize << bit;
    let levels = list
        .split(',')
        .map(|v| v.trim().parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    if levels.len() != expected {
        return Err(format!("Expected {expected} levels, received {}", levels.len()).into());
    }
    Ok(levels)
}

pub fn parse_shuffle(list: &str) -> Result<Vec<usize>, DynError> {
    let values_external: Vec<usize> = list
        .split(',')
        .map(|v| v.trim().parse::<usize>())
        .collect::<Result<Vec<_>, _>>()?;
    if values_external.len() != 32 {
        return Err("Shuffle map must contain exactly 32 entries".into());
    }
    let mut sorted = values_external.clone();
    sorted.sort_unstable();
    for (expected, &found) in (0usize..32).zip(sorted.iter()) {
        if expected != found {
            return Err("Shuffle map must be a permutation of 0..31".into());
        }
    }
    let mut values_internal = vec![0usize; 32];
    for (idx_msb_to_lsb, input_bit) in values_external.into_iter().enumerate() {
        let out_bit_lsb = 31 - idx_msb_to_lsb;
        values_internal[out_bit_lsb] = input_bit;
    }
    Ok(values_internal)
}

pub fn resolve_weight(
    tsys: f64,
    gain: f64,
    sefd: Option<f64>,
    diameter: Option<f64>,
    eta: f64,
    label: &str,
) -> Result<(f64, f64, Option<f64>, Option<f64>), DynError> {
    if !tsys.is_finite() || tsys <= 0.0 {
        return Err(format!("{label} Tsys must be positive and finite").into());
    }
    if let Some(sefd_v) = sefd {
        if !sefd_v.is_finite() || sefd_v <= 0.0 {
            return Err(format!("{label} SEFD must be positive").into());
        }
        let w = 1.0 / sefd_v;
        Ok((w, tsys * w, Some(sefd_v), None))
    } else if let Some(dia_v) = diameter {
        if !dia_v.is_finite() || !eta.is_finite() || dia_v <= 0.0 || eta <= 0.0 {
            return Err(format!("{label} diameter/eta must be positive").into());
        }
        let geom_area = PI * (dia_v / 2.0).powi(2);
        let eff_area = eta * geom_area;
        let sefd_si = (2.0 * K_BOLTZMANN * tsys) / eff_area;
        let sefd_jy = sefd_si * SI_TO_JY;
        let w = 1.0 / sefd_jy;
        Ok((w, tsys * w, Some(sefd_jy), Some(eff_area)))
    } else {
        if !gain.is_finite() || gain <= 0.0 {
            return Err(format!("{label} gain must be positive").into());
        }
        Ok((gain / tsys, gain, None, None))
    }
}

#[cfg(test)]
mod tests {
    use super::{resolve_weight, Args};
    use clap::Parser;
    use std::path::PathBuf;

    #[test]
    fn usb_storage_flag_parses_with_mkxml_and_defaults_off() {
        let default_args = Args::try_parse_from(["yi-corr", "--mkxml"]).unwrap();
        assert!(!default_args.usb);

        let usb_args = Args::try_parse_from(["yi-corr", "--mkxml", "--usb"]).unwrap();
        assert!(usb_args.usb);
        assert!(usb_args.mkxml);
    }

    #[test]
    fn usb_storage_flag_is_distinct_from_usb_signal_sideband() {
        let sideband_only =
            Args::try_parse_from(["yi-corr", "--mkxml", "--sideband", "USB"]).unwrap();
        assert!(!sideband_only.usb);
        assert_eq!(sideband_only.sideband, vec!["USB".to_string()]);

        let both =
            Args::try_parse_from(["yi-corr", "--mkxml", "--sideband", "USB", "--usb"]).unwrap();
        assert!(both.usb);
        assert_eq!(both.sideband, vec!["USB".to_string()]);
    }

    #[test]
    fn model_diagnostic_and_unattended_sweep_flags_parse() {
        let default_args = Args::try_parse_from(["yi-corr", "--mkxml"]).unwrap();
        assert!(!default_args.model_diagnostics);
        assert!(!default_args.model_sweep);

        let diagnostic =
            Args::try_parse_from(["yi-corr", "--mkxml", "--model-diagnostics"]).unwrap();
        assert!(diagnostic.model_diagnostics);
        assert!(!diagnostic.model_sweep);

        let sweep = Args::try_parse_from(["yi-corr", "--mkxml", "--model-sweep"]).unwrap();
        assert!(sweep.model_sweep);
    }

    #[test]
    fn phased_output_options_have_safe_defaults() {
        let args = Args::try_parse_from(["yi-phasedarray", "--mkxml"]).unwrap();
        assert_eq!(args.phased_name, "YAMAGU66");
        assert!(!args.phased_diagnostics);
        assert!(!args.phased_validation);
        assert!(args.phased_validation_npz.is_none());
        assert!(!args.gain_phasecal);
        assert!(args.gain_source.is_none());
        assert_eq!(args.gain_fringe_length, 1);
        assert!(args.phasecal_schedule_output.is_none());
        assert!(args.gain_uncalibrated_schedule.is_none());
        assert!(args.gain_correlation_directory.is_none());
        assert_eq!(args.gain_reference_key, "K");
        assert_eq!(args.gain_corrected_key, "L");
        assert!(!args.xcf_only);
        assert!(args.cor_label_override.is_none());

        let args = Args::try_parse_from([
            "yi-phasedarray",
            "--mkxml",
            "--phased-name",
            "ARRAY1",
            "--phased-diagnostics",
        ])
        .unwrap();
        assert_eq!(args.phased_name, "ARRAY1");
        assert!(args.phased_diagnostics);

        let args = Args::try_parse_from([
            "yi-phasedarray",
            "--mkxml",
            "--phased-validation",
            "--phased-validation-npz",
            "validation.npz",
        ])
        .unwrap();
        assert!(args.phased_validation);
        assert_eq!(
            args.phased_validation_npz,
            Some(PathBuf::from("validation.npz"))
        );

        let args = Args::try_parse_from([
            "yi-phasedarray",
            "--mkxml",
            "--gain-uncalibrated-sc",
            "gain-original.xml",
            "--gain-correlation-directory",
            "gain-cor",
            "--gain-reference-key",
            "K",
            "--gain-corrected-key",
            "L",
        ])
        .unwrap();
        assert_eq!(
            args.gain_uncalibrated_schedule,
            Some(PathBuf::from("gain-original.xml"))
        );
        assert_eq!(
            args.gain_correlation_directory,
            Some(PathBuf::from("gain-cor"))
        );
        assert_eq!(args.gain_reference_key, "K");
        assert_eq!(args.gain_corrected_key, "L");

        let args = Args::try_parse_from([
            "yi-phasedarray",
            "--mkxml",
            "--gain-phasecal",
            "--gain-source",
            "J1924-2914",
            "--gain-fringe-length",
            "2",
            "--phasecal-schedule-output",
            "phasecal.xml",
        ])
        .unwrap();
        assert!(args.gain_phasecal);
        assert_eq!(args.gain_source.as_deref(), Some("J1924-2914"));
        assert_eq!(args.gain_fringe_length, 2);
        assert_eq!(
            args.phasecal_schedule_output,
            Some(PathBuf::from("phasecal.xml"))
        );
    }

    #[test]
    fn phased_weights_reject_nonphysical_inputs() {
        assert!(resolve_weight(0.0, 1.0, None, None, 0.65, "A1").is_err());
        assert!(resolve_weight(-1.0, 1.0, None, None, 0.65, "A1").is_err());
        assert!(resolve_weight(100.0, f64::NAN, None, None, 0.65, "A1").is_err());
    }
}
