use crate::utils::DynError;

pub fn finalize_auto_spectrum(
    accumulated_auto_spec: &mut [f64],
    power_norm: f64,
) -> Result<Vec<f64>, DynError> {
    let mut auto_mag_spectrum: Vec<f64> = accumulated_auto_spec
        .iter()
        .map(|&value| {
            if power_norm > 0.0 {
                value / power_norm
            } else {
                value
            }
        })
        .collect();

    if !auto_mag_spectrum.is_empty() {
        auto_mag_spectrum[0] = 0.0;
    }

    Ok(auto_mag_spectrum)
}
