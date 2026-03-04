use crate::utils::DynError;
use clap::Parser;
use std::f64::consts::PI;
use std::path::PathBuf;

const K_BOLTZMANN: f64 = 1.380_649e-23; // J/K
const SI_TO_JY: f64 = 1.0e26; // multiply W/m^2/Hz to express in Jy
pub const DEFAULT_FFT: usize = 16384;
pub const DEFAULT_SHUFFLE_IN: [usize; 32] = [
    31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8,
    7, 6, 5, 4, 3, 2, 1, 0,
];

#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about = "yi toolkit for FX correlation and phased-array processing",
    long_about = "yi-acf: auto-correlation only\nyi-xcf: cross-correlation only\nyi-corr: auto + cross correlation\nyi-phasedarray: phased-array synthesis and spectrum plots",
    infer_long_args = true,
    arg_required_else_help = true
)]
pub struct Args {
    // Schedule
    #[arg(long, visible_alias = "sc", value_name = "XML", help = "Observation schedule XML file (.xml)")]
    pub schedule: Option<PathBuf>,

    #[arg(long, help = "Write example schedule XML to example.xml and exit")]
    pub mkxml: bool,

    // Input files
    #[arg(long = "raw-directory", visible_alias = "data", visible_alias = "raw", value_name = "DIR", help = "Directory containing input raw files")]
    pub raw_directory: Option<PathBuf>,

    #[arg(long, value_name = "FILE", help = "Input raw file for antenna 1")]
    pub ant1: Option<PathBuf>,

    #[arg(long, value_name = "FILE", help = "Input raw file for antenna 2")]
    pub ant2: Option<PathBuf>,

    // Output files
    #[arg(long = "cor-directory", visible_alias = "output", visible_alias = "cor", value_name = "DIR", help = "Directory for output .cor/.raw products")]
    pub cor_directory: Option<PathBuf>,

    // Correlation processing parameters
    #[arg(long, value_name = "HMS", help = "Source right ascension (HHhMMmSS.s or degrees)")]
    pub ra: Option<String>,

    #[arg(long, value_name = "DMS", allow_hyphen_values = true, help = "Source declination (DDdMMmSS.s or degrees)")]
    pub dec: Option<String>,

    #[arg(long, value_name = "YYYY/DDD HH:MM:SS", help = "Observation epoch (UTC)")]
    pub epoch: Option<String>,

    #[arg(long, value_name = "MHZ", default_value_t = 0.0, help = "Reference observing frequency [MHz]")]
    pub obsfreq: f64,

    #[arg(long, value_name = "MSPS", default_value_t = 1024.0, help = "Sampling rate [Msps]")]
    pub sampling: f64,

    #[arg(long, value_name = "POINTS", default_value_t = DEFAULT_FFT, help = "FFT size [samples/frame]")]
    pub fft: usize,

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

    #[arg(long, value_name = "LSB|USB", num_args = 1.., default_value = "LSB", help = "Input sampler sideband")]
    pub sideband: Vec<String>,

    #[arg(long, value_name = "S", allow_hyphen_values = true, help = "Fixed coarse delay [s]")]
    pub coarse: Option<f64>,

    #[arg(long, value_name = "SAMPLES", allow_hyphen_values = true, default_value_t = 0.0, help = "Residual delay applied to ant2 [samples]")]
    pub delay: f64,

    #[arg(long, value_name = "HZ", default_value_t = 0.0, allow_hyphen_values = true, help = "Residual fringe rate correction [Hz]")]
    pub rate: f64,

    #[arg(long, value_name = "SAMPLES", default_value_t = 0.0, allow_hyphen_values = true, help = "Additional residual delay [samples]")]
    pub resdelay: f64,

    #[arg(long, value_name = "HZ", default_value_t = 0.0, allow_hyphen_values = true, help = "Additional residual rate [Hz]")]
    pub resrate: f64,

    #[arg(long, value_name = "HZ", num_args = 1.., allow_hyphen_values = true, help = "Per-antenna rotation frequency [Hz] (e.g. \"ant1:343000000 ant2:0\")")]
    pub rotation: Vec<String>,

    #[arg(long, value_name = "S", default_value_t = 0.0, help = "Start offset from input head [s]")]
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

    #[arg(long, help = "Enable fringe search in yi-phasedarray mode only")]
    pub fringe: bool,

    // Other
    #[arg(long, value_name = "N", help = "Number of compute threads [default: logical CPU count - 2, min 1]")]
    pub cpu: Option<usize>,

    #[arg(long, value_name = "N", help = "Reader chunk size [frames] (auto when omitted)")]
    pub chunk_frames: Option<usize>,

    #[arg(
        long = "razoku5bay",
        help = "Force small I/O chunking at 0.5-second intervals (overrides --chunk-frames)"
    )]
    pub razoku5bay: bool,

    #[arg(long, value_name = "N", help = "Reader-worker queue depth [chunks] (auto when omitted)")]
    pub pipeline_depth: Option<usize>,

    #[arg(long, help = "Enable diagnostic debug logging to file (all frames)")]
    pub debug: bool,

    #[arg(
        long = "stdout",
        default_value_t = false,
        help = "Write yi-corr runtime stdout log to ./stdout/stdout_<yyyydddhhmmss>.log"
    )]
    pub stdout: bool,

    #[arg(long, hide = true, value_name = "N")]
    pub process_index: Option<usize>,

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
    let mut ant1_val = default_val.clone();
    let mut ant2_val = default_val;
    for entry in args_list {
        if entry.contains("ant1:") || entry.contains("ant2:") {
            let parts: Vec<&str> = entry.split(|c: char| c == ',' || c == ' ').filter(|s| !s.is_empty()).collect();
            for part in parts {
                if let Some(val_str) = part.strip_prefix("ant1:") { ant1_val = parser(val_str.trim())?; }
                else if let Some(val_str) = part.strip_prefix("ant2:") { ant2_val = parser(val_str.trim())?; }
                else { let v = parser(part.trim())?; ant1_val = v.clone(); ant2_val = v; }
            }
        } else {
            let v = parser(entry.trim())?; ant1_val = v.clone(); ant2_val = v;
        }
    }
    Ok((ant1_val, ant2_val))
}

pub fn parse_levels(bit: usize, list: &str) -> Result<Vec<f64>, DynError> {
    if bit == 0 { return Err("Bit depth must be at least 1".into()); }
    let expected = 1usize << bit;
    let levels = list.split(',').map(|v| v.trim().parse::<f64>()).collect::<Result<Vec<_>, _>>()?;
    if levels.len() != expected { return Err(format!("Expected {expected} levels, received {}", levels.len()).into()); }
    Ok(levels)
}

pub fn parse_shuffle(list: &str) -> Result<Vec<usize>, DynError> {
    let values_external: Vec<usize> = list.split(',').map(|v| v.trim().parse::<usize>()).collect::<Result<Vec<_>, _>>()?;
    if values_external.len() != 32 { return Err("Shuffle map must contain exactly 32 entries".into()); }
    let mut sorted = values_external.clone(); sorted.sort_unstable();
    for (expected, &found) in (0usize..32).zip(sorted.iter()) { if expected != found { return Err("Shuffle map must be a permutation of 0..31".into()); } }
    let mut values_internal = vec![0usize; 32];
    for (idx_msb_to_lsb, input_bit) in values_external.into_iter().enumerate() { let out_bit_lsb = 31 - idx_msb_to_lsb; values_internal[out_bit_lsb] = input_bit; }
    Ok(values_internal)
}

pub fn resolve_weight(tsys: f64, gain: f64, sefd: Option<f64>, diameter: Option<f64>, eta: f64, label: &str) -> Result<(f64, f64, Option<f64>, Option<f64>), DynError> {
    if tsys == 0.0 { return Err(format!("{label} Tsys must be non-zero").into()); }
    if let Some(sefd_v) = sefd {
        if sefd_v <= 0.0 { return Err(format!("{label} SEFD must be positive").into()); }
        let w = 1.0 / sefd_v; Ok((w, tsys * w, Some(sefd_v), None))
    } else if let Some(dia_v) = diameter {
        if dia_v <= 0.0 || eta <= 0.0 { return Err(format!("{label} diameter/eta must be positive").into()); }
        let geom_area = PI * (dia_v / 2.0).powi(2); let eff_area = eta * geom_area;
        let sefd_si = (2.0 * K_BOLTZMANN * tsys) / eff_area; let sefd_jy = sefd_si * SI_TO_JY;
        let w = 1.0 / sefd_jy; Ok((w, tsys * w, Some(sefd_jy), Some(eff_area)))
    } else {
        if gain <= 0.0 { return Err(format!("{label} gain must be positive").into()); }
        Ok((gain / tsys, gain, None, None))
    }
}
