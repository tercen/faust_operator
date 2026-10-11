//! Time phase 1 (noise + exhaustive annotation forest) on synthetic data shaped like
//! `dev/profile_faust.R`: `cargo run --release -p faust --example bench_forest -- 50000 10`.
use std::time::Instant;

use faust::forest::{SearchParams, grow_annotation_forest};
use faust::rrng::forest_seed;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(50_000);
    let d: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(10);
    // xorshift + Box-Muller mixture: marker j has 1 + (j % 3) modes
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut unif = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut cols = Vec::new();
    for j in 0..d {
        let modes = 1 + (j % 3);
        let mut c = Vec::with_capacity(n);
        for _ in 0..n {
            let u = unif();
            let comp = ((u * modes as f64) as usize).min(modes - 1);
            let u1 = unif().max(f64::MIN_POSITIVE);
            let u2 = unif();
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            c.push(comp as f64 * 2.5 + 0.45 * z);
        }
        cols.push(c);
    }
    let params = SearchParams::annotation_forest(3);
    let t = Instant::now();
    let gates = grow_annotation_forest(&cols, &params, forest_seed(123), 4.0, 2);
    let peak = std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1)?.parse::<u64>().ok())
        })
        .unwrap_or(0);
    println!(
        "forest: {n} cells x {d} markers -> {} gates in {:.2} s, peak {} MB",
        gates.len(),
        t.elapsed().as_secs_f64(),
        peak / 1024
    );
}
