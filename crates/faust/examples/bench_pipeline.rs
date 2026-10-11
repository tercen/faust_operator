//! Time the whole FAUST run on synthetic data shaped like `dev/profile_faust.R` (the R
//! baseline): `cargo run --release -p faust --example bench_pipeline -- <cells per sample> <markers> <samples> [threads]`.
use std::time::Instant;

use faust::pipeline::{FaustSettings, Sample, run_faust};

fn peak_mb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1)?.parse::<u64>().ok())
        })
        .unwrap_or(0)
        / 1024
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(50_000);
    let d: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(10);
    let g: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2);
    let threads: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut unif = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let markers: Vec<String> = (1..=d).map(|j| format!("M{j}")).collect();
    let samples: Vec<Sample> = (0..g)
        .map(|s| {
            let cols = (0..d)
                .map(|j| {
                    let modes = 1 + (j % 3);
                    (0..n)
                        .map(|_| {
                            let u = unif();
                            let comp = ((u * modes as f64) as usize).min(modes - 1);
                            let u1 = unif().max(f64::MIN_POSITIVE);
                            let u2 = unif();
                            let z = (-2.0 * u1.ln()).sqrt()
                                * (2.0 * std::f64::consts::PI * u2).cos();
                            comp as f64 * 2.5 + s as f64 * 0.1 + 0.45 * z
                        })
                        .collect()
                })
                .collect();
            Sample {
                name: format!("s{s}"),
                cols,
            }
        })
        .collect();
    let settings = FaustSettings {
        selection_quantile: 0.5,
        depth_score_threshold: 0.01,
        seed: 123,
        threads,
    };
    let t = Instant::now();
    let mut last = Instant::now();
    let run = run_faust(&samples, &markers, &settings, &mut |m| {
        eprintln!("[{:7.1}s +{:6.1}s] {m}", t.elapsed().as_secs_f64(), last.elapsed().as_secs_f64());
        last = Instant::now();
    })
    .unwrap();
    println!(
        "pipeline: {g} x {n} cells x {d} markers, threads {threads}: {:.1} s, peak {} MB, {} markers selected, {} phenotypes",
        t.elapsed().as_secs_f64(),
        peak_mb(),
        run.selected.len(),
        run.phenotypes.len()
    );
}
