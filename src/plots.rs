//! The diagnostic histograms (`R/plotMarkerHistograms.R`): for every selected marker and every
//! threshold index, a histogram of up to 1,000 cells per sample (between the 1st and 99th
//! percentile) with one dashed red line per sample at that sample's threshold. One PNG per
//! (marker, threshold), 6 x 6 inches at 100 dpi, named `hist_<marker>_ab_<k>.png` like R's.
//!
//! Returned as the rows of the R operator's third table (`tim::plot_file_to_df`): file name,
//! MIME type, MD5 of the bytes, base64 content.
use anyhow::{Result, anyhow};
use base64::Engine;
use md5::{Digest, Md5};
use plotters::prelude::*;

use crate::output::Png;
use faust::pipeline::FaustRun;

/// The same deterministic draw as the embedding jitter: FAUST's `sample()` is unseeded here.
struct Pick(u64);
impl Pick {
    fn next(&mut self, n: usize) -> usize {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x % n as u64) as usize
    }
}

/// Freedman–Diaconis class count (`grDevices::nclass.FD`), then `pretty` breaks, as R.
fn fd_binwidth(x: &[f64]) -> f64 {
    let mut s = x.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = s.len() as f64;
    let q = |p: f64| faust::scamp::r_quantile7(&s, p);
    let mut h = 2.0 * (q(0.75) - q(0.25));
    if h == 0.0 {
        // nclass.FD falls back to the MAD-based width
        let med = q(0.5);
        let mut d: Vec<f64> = s.iter().map(|v| (v - med).abs()).collect();
        d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        h = 2.0 * 1.4826 * faust::scamp::r_quantile7(&d, 0.5);
    }
    let range = s[s.len() - 1] - s[0];
    let nclass = if h > 0.0 {
        (range / (h * n.powf(-1.0 / 3.0))).ceil().max(1.0)
    } else {
        1.0
    };
    // pretty(): a "nice" step near range / nclass (1, 2, 5 × 10^k)
    let raw = range / nclass;
    if raw <= 0.0 || !raw.is_finite() {
        return 1.0;
    }
    let k = 10f64.powf(raw.log10().floor());
    let f = raw / k;
    let nice = if f < 1.5 {
        1.0
    } else if f < 3.0 {
        2.0
    } else if f < 7.0 {
        5.0
    } else {
        10.0
    };
    nice * k
}

pub fn marker_histograms(run: &FaustRun, names: &[String], seed: u64) -> Result<Vec<Png>> {
    let n_samples = run.selected_exprs.len();
    let mut pick = Pick(seed.max(1) ^ 0x2545_F491_4F6C_DD1D);
    let mut out = Vec::new();
    for (j, name) in names.iter().enumerate() {
        // 1,000 random cells per sample
        let mut values: Vec<f64> = Vec::new();
        for g in 0..n_samples {
            let col = &run.selected_exprs[g][j];
            let n = col.len();
            if n <= 1000 {
                values.extend_from_slice(col);
            } else {
                let mut chosen = vec![false; n];
                let mut c = 0;
                while c < 1000 {
                    let i = pick.next(n);
                    if !chosen[i] {
                        chosen[i] = true;
                        values.push(col[i]);
                        c += 1;
                    }
                }
            }
        }
        if values.is_empty() {
            continue;
        }
        let mut s = values.clone();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let lo = faust::scamp::r_quantile7(&s, 0.01);
        let hi = faust::scamp::r_quantile7(&s, 0.99);
        let hist: Vec<f64> = values.iter().copied().filter(|v| *v >= lo && *v <= hi).collect();
        if hist.len() < 2 {
            continue;
        }
        let bw = fd_binwidth(&hist);
        let xmin = hist.iter().copied().fold(f64::INFINITY, f64::min);
        let xmax = hist.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let n_bins = ((xmax - xmin) / bw).ceil().max(1.0) as usize + 1;
        let mut counts = vec![0u32; n_bins];
        for v in &hist {
            let b = (((v - xmin) / bw).floor() as usize).min(n_bins - 1);
            counts[b] += 1;
        }
        let n_thresholds = run.thresholds[0].get(j).map_or(0, Vec::len);
        for k in 0..n_thresholds {
            let gates: Vec<f64> = (0..n_samples)
                .filter_map(|g| run.thresholds[g].get(j).and_then(|t| t.get(k)).copied())
                .collect();
            let png = draw(name, &counts, xmin, bw, &gates)?;
            let mut hasher = Md5::new();
            hasher.update(&png);
            out.push(Png {
                filename: format!("hist_{name}_ab_{}.png", k + 1),
                checksum: format!("{:x}", hasher.finalize()),
                content_base64: base64::engine::general_purpose::STANDARD.encode(&png),
            });
        }
    }
    Ok(out)
}

fn draw(title: &str, counts: &[u32], xmin: f64, bw: f64, gates: &[f64]) -> Result<Vec<u8>> {
    const W: u32 = 600;
    const H: u32 = 600;
    let mut buf = vec![0u8; (W * H * 3) as usize];
    {
        let root = BitMapBackend::with_buffer(&mut buf, (W, H)).into_drawing_area();
        root.fill(&WHITE).map_err(|e| anyhow!("{e}"))?;
        let xmax = xmin + bw * counts.len() as f64;
        let ymax = counts.iter().copied().max().unwrap_or(1).max(1) as f64 * 1.05;
        let x_lo = gates.iter().copied().fold(xmin, f64::min);
        let x_hi = gates.iter().copied().fold(xmax, f64::max);
        let mut chart = ChartBuilder::on(&root)
            .caption(title, ("sans-serif", 28))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(x_lo..x_hi, 0.0..ymax)
            .map_err(|e| anyhow!("{e}"))?;
        chart
            .configure_mesh()
            .x_desc("x")
            .y_desc("count")
            .draw()
            .map_err(|e| anyhow!("{e}"))?;
        chart
            .draw_series(counts.iter().enumerate().map(|(i, &c)| {
                let x0 = xmin + bw * i as f64;
                Rectangle::new([(x0, 0.0), (x0 + bw, c as f64)], RGBColor(89, 89, 89).filled())
            }))
            .map_err(|e| anyhow!("{e}"))?;
        for &g in gates {
            chart
                .draw_series(DashedLineSeries::new(
                    [(g, 0.0), (g, ymax)],
                    6,
                    4,
                    RED.stroke_width(2),
                ))
                .map_err(|e| anyhow!("{e}"))?;
        }
        root.present().map_err(|e| anyhow!("{e}"))?;
    }
    let mut png = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut png, W, H);
        enc.set_color(png::ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        let mut writer = enc.write_header().map_err(|e| anyhow!("{e}"))?;
        writer.write_image_data(&buf).map_err(|e| anyhow!("{e}"))?;
    }
    Ok(png)
}
