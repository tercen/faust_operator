//! The noise FAUST adds to the data before growing a forest or running scamp.
//!
//! Port of `src/addNoiseToDataMatrix.cpp`: `addDistributionNoiseToDataVector` (per column:
//! Gaussian jitter scaled by `scaleGaussian`, and a uniform component between repeated
//! values — read the source for the exact rule, including the `unordered_map` count of ties)
//! and `addNoiseToDataMatrix`, which seeds column `k` with `seed + k + 1` when a seed is set
//! (the thread count does not change the result; both branches derive the same per-column
//! seeds). The R-exposed `addNoiseToDataVector(x, scaleGaussian, rSeed)` is the single-column
//! entry point and is what the fixtures were made with.
//!
//! Parity: `tests/parity/noise.rs` against `*_noise_s123.txt` / `*_noise_s124.txt`.
use crate::rng::Mt19937;

/// `addNoiseToDataVector(x, scale_gaussian, seed)`: the noised copy of `x` (not sorted).
pub fn add_noise_to_vector(x: &[f64], scale_gaussian: f64, seed: u32) -> Vec<f64> {
    let _ = (x, scale_gaussian, seed);
    todo!("port addNoiseToDataMatrix.cpp")
}

/// `addNoiseToDataMatrix`: every column noised with its own seed derived from `seed`
/// (`seed + k + 1` for column `k` when `seed > 0`). Returns the noised columns and the seed
/// the C++ leaves in `randomSeed` afterwards (it is incremented per column and reused by the
/// caller for the search).
pub fn add_noise_to_matrix(cols: &[Vec<f64>], scale_gaussian: f64, seed: u64) -> (Vec<Vec<f64>>, u64) {
    let _ = (cols, scale_gaussian, seed);
    todo!()
}

/// A helper the search uses: a fresh `Mt19937` from a 64-bit seed the way the C++ constructs
/// it (`std::mt19937 mt(rSeed)` with an `unsigned long long` narrowed to the engine's 32-bit
/// `result_type`).
pub fn engine_from_seed(seed: u64) -> Mt19937 {
    Mt19937::new(seed as u32)
}
