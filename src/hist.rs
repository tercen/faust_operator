//! The threshold histograms: one PNG per selected marker, as FAUST's `.plotMarkerHistograms`
//! draws them for the R operator — the marker's expression over at most 1000 cells per sample,
//! pooled, trimmed to the pooled 1 %/99 % quantiles, binned at the Freedman–Diaconis width
//! rounded by R's `pretty`, with every sample's thresholds as red dashed lines and the marker
//! as title. The R operator wrote one file per gate number; here all gates of a marker are on
//! the one plot, so the table has one row per marker and joins the row relation on `.ri`.
//!
//! Memory: the pooled sample is ≤ 1000 × samples doubles and the canvas is one 800 × 600 RGB
//! buffer, reused; the PNGs are tens of KB each.
use anyhow::{Result, anyhow};
use plotters::prelude::*;
use plotters::style::register_font;
use std::sync::Once;

pub const WIDTH: u32 = 800;
pub const HEIGHT: u32 = 600;
pub const CELLS_PER_SAMPLE: usize = 1000;

static FONT: Once = Once::new();

fn ensure_font() {
    FONT.call_once(|| {
        register_font(
            "sans-serif",
            FontStyle::Normal,
            include_bytes!("../fonts/DejaVuSans.ttf"),
        )
        .ok();
    });
}

/// xorshift64*: enough for picking which cells to draw.
fn next(state: &mut u64) -> u64 {
    *state ^= *state >> 12;
    *state ^= *state << 25;
    *state ^= *state >> 27;
    state.wrapping_mul(0x2545_F491_4F6C_DD1D)
}

/// At most [`CELLS_PER_SAMPLE`] values of every sample, pooled (R: `sample(seq(n), min(1000, n))`
/// per experimental unit, then `rbind`). Non-finite values are dropped.
pub fn pooled_sample(per_sample: &[&[f64]], seed: u64) -> Vec<f64> {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let mut out = Vec::new();
    for v in per_sample {
        if v.len() <= CELLS_PER_SAMPLE {
            out.extend(v.iter().copied().filter(|x| x.is_finite()));
        } else {
            // partial Fisher–Yates on an index vector
            let mut idx: Vec<u32> = (0..v.len() as u32).collect();
            for k in 0..CELLS_PER_SAMPLE {
                let r = k + (next(&mut state) % (v.len() - k) as u64) as usize;
                idx.swap(k, r);
                let x = v[idx[k] as usize];
                if x.is_finite() {
                    out.push(x);
                }
            }
        }
    }
    out
}

/// R's `quantile(type = 7)` on sorted data.
fn quantile7(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    let h = (n - 1) as f64 * p;
    let lo = h.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    sorted[lo] + (h - lo as f64) * (sorted[hi] - sorted[lo])
}

/// Keep the values within the 1 % and 99 % quantiles (inclusive), sorted.
pub fn trim_1_99(values: &[f64]) -> Vec<f64> {
    let mut s: Vec<f64> = values.iter().copied().filter(|x| x.is_finite()).collect();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (lo, hi) = (quantile7(&s, 0.01), quantile7(&s, 0.99));
    s.into_iter().filter(|&x| x >= lo && x <= hi).collect()
}

/// `grDevices::nclass.FD` on sorted data.
fn nclass_fd(sorted: &[f64]) -> usize {
    let n = sorted.len();
    if n < 2 {
        return 1;
    }
    let mut h = 2.0 * (quantile7(sorted, 0.75) - quantile7(sorted, 0.25));
    if h == 0.0 {
        // R widens the quantile pair until the spread is non-zero
        let mut al = 0.25;
        while h == 0.0 && al >= 1.0 / 512.0 {
            al /= 2.0;
            h = (quantile7(sorted, 1.0 - al) - quantile7(sorted, al)) * (1.0 - 2.0 * al);
        }
    }
    if h > 0.0 {
        let range = sorted[n - 1] - sorted[0];
        ((range / (h * (n as f64).powf(-1.0 / 3.0))).ceil() as usize).max(1)
    } else {
        1
    }
}

/// The unit R's `pretty(c(lo, up), n)` chooses (the core of `R_pretty`, `min.n = 1`,
/// `shrink.sml = 0.75`, the default biases); the bin width `.getHistogram` uses.
pub fn pretty_unit(lo: f64, up: f64, ndiv: usize) -> f64 {
    let (h, h5) = (1.5f64, 0.5 + 1.5 * 1.5);
    let shrink_sml = 0.75;
    let dx = up - lo;
    let mut cell;
    let i_small;
    if dx == 0.0 && up == 0.0 {
        cell = 1.0;
        i_small = true;
    } else {
        cell = lo.abs().max(up.abs());
        let mut u = 1.0
            + if h5 >= 1.5 * h + 0.5 {
                1.0 / (1.0 + h)
            } else {
                1.5 / (1.0 + h5)
            };
        u *= (ndiv.max(1)) as f64 * f64::EPSILON;
        i_small = dx < cell * u * 3.0;
    }
    if i_small {
        if cell > 10.0 {
            cell = 9.0 + cell / 10.0;
        }
        cell *= shrink_sml;
    } else {
        cell = dx;
        if ndiv > 1 {
            cell /= ndiv as f64;
        }
    }
    if cell < 20.0 * f64::MIN_POSITIVE {
        cell = 20.0 * f64::MIN_POSITIVE;
    }
    let base = 10f64.powf(cell.log10().floor());
    let mut unit = base;
    if 2.0 * base - cell < h * (cell - unit) {
        unit = 2.0 * base;
        if 5.0 * base - cell < h5 * (cell - unit) {
            unit = 5.0 * base;
            if 10.0 * base - cell < h * (cell - unit) {
                unit = 10.0 * base;
            }
        }
    }
    unit
}

/// ggplot2 `geom_histogram(binwidth = w)` default bins: edges at `w/2 + k·w`, right-closed.
/// Returns `(edges, counts)`.
fn bins(sorted: &[f64], w: f64) -> (Vec<f64>, Vec<u32>) {
    let (lo, hi) = (sorted[0], sorted[sorted.len() - 1]);
    let boundary = w / 2.0;
    let origin = boundary + ((lo - boundary) / w).floor() * w;
    let n = (((hi - origin) / w).floor() as usize) + 1;
    let edges: Vec<f64> = (0..=n).map(|k| origin + k as f64 * w).collect();
    let mut counts = vec![0u32; n];
    for &x in sorted {
        // (a, b]: a value on an edge goes to the bin on its left
        let mut k = ((x - origin) / w).ceil() as isize - 1;
        if k < 0 {
            k = 0;
        }
        let k = (k as usize).min(n - 1);
        counts[k] += 1;
    }
    (edges, counts)
}

/// Render one marker's histogram with its gates as a PNG.
pub fn render(marker: &str, trimmed_sorted: &[f64], gates: &[f64]) -> Result<Vec<u8>> {
    ensure_font();
    if trimmed_sorted.is_empty() {
        return Err(anyhow!("no finite values for {marker}"));
    }
    let w = pretty_unit(
        trimmed_sorted[0],
        trimmed_sorted[trimmed_sorted.len() - 1],
        nclass_fd(trimmed_sorted),
    );
    let (edges, counts) = bins(trimmed_sorted, w);
    let x_min = edges[0].min(gates.iter().copied().fold(f64::INFINITY, f64::min));
    let x_max = edges[edges.len() - 1].max(gates.iter().copied().fold(f64::NEG_INFINITY, f64::max));
    let pad = (x_max - x_min).max(w) * 0.03;
    let (x_min, x_max) = (x_min - pad, x_max + pad);
    let y_max = (*counts.iter().max().unwrap_or(&1)).max(1) as f64 * 1.05;

    let mut buf = vec![255u8; (WIDTH * HEIGHT * 3) as usize];
    {
        let root = BitMapBackend::with_buffer(&mut buf, (WIDTH, HEIGHT)).into_drawing_area();
        root.fill(&WHITE).map_err(|e| anyhow!("{e}"))?;
        let mut chart = ChartBuilder::on(&root)
            .caption(marker, ("sans-serif", 30))
            .margin(18)
            .x_label_area_size(48)
            .y_label_area_size(70)
            .build_cartesian_2d(x_min..x_max, 0f64..y_max)
            .map_err(|e| anyhow!("{e}"))?;
        // theme_bw: white panel, grey grid, black border
        chart
            .configure_mesh()
            .x_desc(marker)
            .y_desc("count")
            .axis_desc_style(("sans-serif", 20))
            .label_style(("sans-serif", 16))
            .y_label_formatter(&|v| format!("{}", v.round() as i64))
            .light_line_style(TRANSPARENT)
            .bold_line_style(RGBColor(225, 225, 225))
            .axis_style(BLACK)
            .draw()
            .map_err(|e| anyhow!("{e}"))?;
        let fill = RGBColor(89, 89, 89).filled();
        chart
            .draw_series(
                counts
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| **c > 0)
                    .map(|(k, &c)| {
                        Rectangle::new([(edges[k], 0.0), (edges[k + 1], c as f64)], fill)
                    }),
            )
            .map_err(|e| anyhow!("{e}"))?;
        for &g in gates {
            if g.is_finite() {
                chart
                    .draw_series(DashedLineSeries::new(
                        [(g, 0.0), (g, y_max)],
                        10,
                        6,
                        RED.stroke_width(2),
                    ))
                    .map_err(|e| anyhow!("{e}"))?;
            }
        }
        root.present().map_err(|e| anyhow!("{e}"))?;
    }

    let mut out = Vec::with_capacity(64 << 10);
    {
        let mut enc = png::Encoder::new(&mut out, WIDTH, HEIGHT);
        enc.set_color(png::ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        let mut wr = enc.write_header()?;
        wr.write_image_data(&buf)?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pretty_matches_r_on_known_cases() {
        // R 4.6.1: diff(pretty(c(lo, up), n))[1]
        assert_eq!(pretty_unit(0.0, 1.0, 7), 0.2);
        assert_eq!(pretty_unit(0.0, 100.0, 5), 20.0);
        assert_eq!(pretty_unit(1.0, 1.0, 5), 1.0);
        assert!((pretty_unit(-1.3, 4.7, 23) - 0.2).abs() < 1e-12);
        assert!((pretty_unit(0.5, 3.2, 12) - 0.2).abs() < 1e-12);
    }

    #[test]
    fn a_png_comes_out() {
        let v: Vec<f64> = (0..5000)
            .map(|i| ((i * 7919) % 1000) as f64 / 100.0)
            .collect();
        let t = trim_1_99(&v);
        let png = render("CD3", &t, &[2.5, 7.0]).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        assert!(png.len() > 5_000);
    }

    #[test]
    fn pooling_caps_each_sample() {
        let a: Vec<f64> = (0..5000).map(|i| i as f64).collect();
        let b: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let p = pooled_sample(&[&a, &b], 7);
        assert_eq!(p.len(), 1010);
    }
}
