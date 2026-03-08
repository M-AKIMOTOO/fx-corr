use crate::utils::DynError;
use plotters::prelude::PathElement;
use plotters::prelude::*;

pub use plotters::prelude::{RGBColor, BLUE, GREEN, RED};

const PLOT_FONT_SCALE: f64 = 1.2;

fn scaled_font_size(base: i32) -> i32 {
    ((base as f64) * PLOT_FONT_SCALE).round() as i32
}

fn scaled_area_size(base: i32) -> i32 {
    ((base as f64) * PLOT_FONT_SCALE).round() as i32
}

pub fn plot_series_f64_x(
    x_vals: &[f64],
    data: &[f64],
    _title: &str,
    filename: &str,
    x_label: &str,
    y_label: &str,
    y_range: Option<(f64, f64)>,
    label: &str,
) -> Result<(), DynError> {
    if x_vals.len() != data.len() {
        return Err("X-value vector length does not match data length".into());
    }
    if x_vals.is_empty() {
        return Err("No data points to plot".into());
    }

    let root = BitMapBackend::new(filename, (1280, 720)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_min = *x_vals
        .iter()
        .min_by(|a, b| a.total_cmp(b))
        .ok_or("Failed to determine minimum x value")?;
    let x_max = *x_vals
        .iter()
        .max_by(|a, b| a.total_cmp(b))
        .ok_or("Failed to determine maximum x value")?;

    let (min_val, max_val) = if let Some(range) = y_range {
        (range.0, range.1)
    } else {
        (
            data.iter().cloned().fold(f64::INFINITY, f64::min),
            data.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        )
    };

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(scaled_area_size(40))
        .y_label_area_size(scaled_area_size(60))
        .build_cartesian_2d(x_min..x_max, min_val..max_val)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .label_style(("sans-serif", scaled_font_size(20)).into_font())
        .axis_desc_style(("sans-serif", scaled_font_size(24)).into_font())
        .light_line_style(WHITE.mix(0.0))
        .draw()?;

    chart
        .draw_series(LineSeries::new(
            x_vals.iter().zip(data.iter()).map(|(x, y)| (*x, *y)),
            &BLUE,
        ))
        .map(|s| {
            s.label(label)
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], BLUE))
        })?;

    chart
        .configure_series_labels()
        .border_style(BLACK)
        .background_style(&WHITE.mix(0.8))
        .label_font(("sans-serif", scaled_font_size(20)).into_font())
        .draw()?;

    root.present()?;
    Ok(())
}

pub fn plot_multi_series_f64_x(
    x_vals: &[f64],
    series: &[(&[f64], &RGBColor, &str)],
    _title: &str,
    filename: &str,
    x_label: &str,
    y_label: &str,
    y_range: Option<(f64, f64)>,
) -> Result<(), DynError> {
    if series.is_empty() {
        return Err("No series provided to plot".into());
    }
    for (data_series, _, _) in series.iter() {
        if data_series.len() != x_vals.len() {
            return Err("X-value vector length does not match data length".into());
        }
    }
    if x_vals.is_empty() {
        return Err("No data points to plot".into());
    }

    let root = BitMapBackend::new(filename, (1280, 720)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_min = *x_vals
        .iter()
        .min_by(|a, b| a.total_cmp(b))
        .ok_or("Failed to determine minimum x value")?;
    let x_max = *x_vals
        .iter()
        .max_by(|a, b| a.total_cmp(b))
        .ok_or("Failed to determine maximum x value")?;

    let (min_val, max_val) = if let Some(range) = y_range {
        (range.0, range.1)
    } else {
        let mut min_val = f64::INFINITY;
        let mut max_val = f64::NEG_INFINITY;
        for (data_series, _, _) in series.iter() {
            min_val = min_val.min(data_series.iter().cloned().fold(f64::INFINITY, f64::min));
            max_val = max_val.max(
                data_series
                    .iter()
                    .cloned()
                    .fold(f64::NEG_INFINITY, f64::max),
            );
        }
        (min_val, max_val)
    };

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(scaled_area_size(40))
        .y_label_area_size(scaled_area_size(60))
        .build_cartesian_2d(x_min..x_max, min_val..max_val)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .label_style(("sans-serif", scaled_font_size(20)).into_font())
        .axis_desc_style(("sans-serif", scaled_font_size(24)).into_font())
        .light_line_style(WHITE.mix(0.0))
        .draw()?;

    for (data_series, color, label) in series.iter() {
        chart
            .draw_series(LineSeries::new(
                x_vals.iter().zip(data_series.iter()).map(|(x, y)| (*x, *y)),
                *color,
            ))?
            .label(*label)
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], *color));
    }

    chart
        .configure_series_labels()
        .border_style(BLACK)
        .background_style(&WHITE.mix(0.8))
        .label_font(("sans-serif", scaled_font_size(20)).into_font())
        .draw()?;

    root.present()?;
    Ok(())
}

pub fn plot_series_with_x(
    x_vals: &[i32],
    series: &[(&[f64], &RGBColor)],
    _title: &str,
    filename: &str,
    x_label: &str,
    y_label: &str,
    y_range: Option<(f64, f64)>,
    highlight_points: Option<&[(i32, f64)]>,
) -> Result<(), DynError> {
    if x_vals.is_empty() {
        return Err("No data points to plot".into());
    }
    if series.is_empty() {
        return Err("No series provided to plot".into());
    }
    for (data_series, _) in series.iter() {
        if x_vals.len() != data_series.len() {
            return Err(
                "X-value vector length does not match data length for one of the series".into(),
            );
        }
    }

    let root = BitMapBackend::new(filename, (1280, 720)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_min = *x_vals
        .iter()
        .min()
        .ok_or("Failed to determine minimum lag value")?;
    let x_max = *x_vals
        .iter()
        .max()
        .ok_or("Failed to determine maximum lag value")?;

    let (min_val, max_val) = if let Some(range) = y_range {
        range
    } else {
        let mut min_val = f64::INFINITY;
        let mut max_val = f64::NEG_INFINITY;
        for (data_series, _) in series.iter() {
            min_val = min_val.min(data_series.iter().cloned().fold(f64::INFINITY, f64::min));
            max_val = max_val.max(
                data_series
                    .iter()
                    .cloned()
                    .fold(f64::NEG_INFINITY, f64::max),
            );
        }
        (min_val, max_val)
    };

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(scaled_area_size(40))
        .y_label_area_size(scaled_area_size(60))
        .build_cartesian_2d(x_min..(x_max + 1), min_val..max_val)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .label_style(("sans-serif", scaled_font_size(20)).into_font())
        .axis_desc_style(("sans-serif", scaled_font_size(24)).into_font())
        .light_line_style(WHITE.mix(0.0))
        .draw()?;

    for (data_series, color) in series.iter() {
        chart.draw_series(LineSeries::new(
            x_vals.iter().zip(data_series.iter()).map(|(x, y)| (*x, *y)),
            color,
        ))?;
    }

    if let Some(points) = highlight_points {
        chart.draw_series(
            points
                .iter()
                .map(|&(x, y)| Circle::new((x, y), 5, RED.filled())),
        )?;
    }

    root.present()?;
    Ok(())
}

#[allow(dead_code)]
pub fn plot_scatter_complex(
    real_data: &[f64],
    imag_data: &[f64],
    _title: &str,
    filename: &str,
    x_label: &str,
    y_label: &str,
) -> Result<(), DynError> {
    if real_data.len() != imag_data.len() {
        return Err("Real and imaginary data lengths must match".into());
    }
    if real_data.is_empty() {
        return Err("No data points to plot".into());
    }

    let root = BitMapBackend::new(filename, (1280, 720)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut max_abs_val: f64 = 0.0;
    for &val in real_data.iter() {
        max_abs_val = max_abs_val.max(val.abs());
    }
    for &val in imag_data.iter() {
        max_abs_val = max_abs_val.max(val.abs());
    }

    let plot_limit = max_abs_val * 1.1;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(scaled_area_size(40))
        .y_label_area_size(scaled_area_size(60))
        .build_cartesian_2d(-plot_limit..plot_limit, -plot_limit..plot_limit)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .label_style(("sans-serif", scaled_font_size(20)).into_font())
        .axis_desc_style(("sans-serif", scaled_font_size(24)).into_font())
        .light_line_style(WHITE.mix(0.0))
        .draw()?;

    chart.draw_series(
        real_data
            .iter()
            .zip(imag_data.iter())
            .map(|(&r, &i)| Circle::new((r, i), 2, BLUE.filled())),
    )?;

    root.present()?;
    println!("[plot] Wrote scatter plot to {}", filename);
    Ok(())
}

#[allow(dead_code)]
pub fn plot_combined_real_imag_dual_axis(
    x_vals: &[i32],
    real_data: &[f64],
    imag_data: &[f64],
    title: &str,
    filename: &str,
    x_label: &str,
    y_label_real: &str,
    y_label_imag: &str,
) -> Result<(), DynError> {
    if x_vals.len() != real_data.len() || x_vals.len() != imag_data.len() {
        return Err("X-value, real, and imaginary data lengths must match".into());
    }
    if x_vals.is_empty() {
        return Err("No data points to plot".into());
    }

    let root = BitMapBackend::new(filename, (1280, 720)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_min = *x_vals
        .iter()
        .min()
        .ok_or("Failed to determine minimum lag value")?;
    let x_max = *x_vals
        .iter()
        .max()
        .ok_or("Failed to determine maximum lag value")?;

    let real_min = real_data.iter().cloned().fold(f64::INFINITY, f64::min);
    let real_max = real_data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let imag_min = imag_data.iter().cloned().fold(f64::INFINITY, f64::min);
    let imag_max = imag_data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", scaled_font_size(40)).into_font())
        .margin(10)
        .x_label_area_size(scaled_area_size(40))
        .y_label_area_size(scaled_area_size(60))
        .build_cartesian_2d(x_min..(x_max + 1), real_min..real_max)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label_real)
        .label_style(("sans-serif", scaled_font_size(20)).into_font())
        .axis_desc_style(("sans-serif", scaled_font_size(24)).into_font())
        .draw()?;

    chart.draw_series(LineSeries::new(
        x_vals.iter().zip(real_data.iter()).map(|(x, y)| (*x, *y)),
        &BLUE,
    ))?;

    let mut secondary_chart = chart.set_secondary_coord(x_min..(x_max + 1), imag_min..imag_max);

    secondary_chart
        .configure_secondary_axes()
        .y_desc(y_label_imag)
        .axis_desc_style(("sans-serif", scaled_font_size(24)).into_font())
        .draw()?;

    secondary_chart.draw_series(LineSeries::new(
        x_vals.iter().zip(imag_data.iter()).map(|(x, y)| (*x, *y)),
        &GREEN,
    ))?;

    root.present()?;
    println!("[plot] Wrote dual-axis plot to {}", filename);
    Ok(())
}
