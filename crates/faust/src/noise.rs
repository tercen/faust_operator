//! The noise FAUST adds to the data before growing a forest or running scamp.
//!
//! Port of `src/addNoiseToDataMatrix.cpp`. `addNoiseToDataVector(x, scaleGaussian, rSeed)`
//! (the R-exposed single-column entry point the fixtures were made with) seeds one
//! `std::mt19937 mt(rSeed)` and runs `addDistributionNoiseToDataVector` twice on it:
//!
//! 1. **uniform pass** (`unifNoise = true`): if every value is unique the vector is returned
//!    as is (no draw); otherwise every value that occurs more than once (an `unordered_map`
//!    count — the map is only looked up, never iterated, so its order is irrelevant) is
//!    replaced by `cppRunif(low, high)` where, with `u` the sorted unique values and `i` the
//!    position of the value in `u`: `low = (u[i-1] + u[i]) / 2` (or `u[0]` at the bottom),
//!    `high = (u[i] + u[i+1]) / 2` (or `u[last]` at the top). Single occurrences keep their
//!    value and draw nothing. Row order is kept.
//! 2. **Gaussian pass** (`unifNoise = false`) on the result of 1: every value `x` at position
//!    `i` in the (new) unique vector draws `cppRnorm(x, sd)` with `sd = dHigh/s` at the bottom,
//!    `dLow/s` at the top, `(dHigh + dLow) / (2 s)` inside, `s = scaleGaussian` (FAUST uses 4).
//!    `cppRnorm` constructs a fresh `std::normal_distribution` per call, so the polar method's
//!    saved second deviate is never reused: every draw costs two canonicals (four 32-bit
//!    outputs) per polar iteration.
//!
//! `addNoiseToDataMatrix(dm, numThreads, scaleGaussian, randomSeed&)`: with `randomSeed > 0`
//! column `k` (0-based) is noised with seed `randomSeed + k + 1`, and `randomSeed` is left at
//! `randomSeed + ncol` for the caller. The threaded branch (`numThreads > 2`) hands out the
//! same `currentSeed += 1` sequence in column order, so its result is identical to the serial
//! one; only the `randomSeed == 0` case (a `std::random_device` per column) is irreproducible
//! and is not supported here.
//!
//! Parity: `tests/parity_noise.rs` against `*_noise_s123.txt` / `*_noise_s124.txt`.
use crate::rng::{Mt19937, Normal, uniform_real};
use std::collections::HashMap;

/// Hash key for a double the way `std::unordered_map<double, long>` groups them: `==` on the
/// values, so `-0.0` and `0.0` are one key.
fn key(v: f64) -> u64 {
    if v == 0.0 {
        0.0f64.to_bits()
    } else {
        v.to_bits()
    }
}

/// `std::lower_bound(uniq.begin(), uniq.end(), obs)` on the sorted unique values.
fn lower_index(uniq: &[f64], obs: f64) -> usize {
    uniq.partition_point(|&u| u < obs)
}

/// `addDistributionNoiseToDataVector(dataVector, unifNoise, scaleGaussian, twistRef)`.
fn add_distribution_noise(
    data: &[f64],
    unif_noise: bool,
    scale_gaussian: f64,
    rng: &mut Mt19937,
) -> Vec<f64> {
    let mut uniq = data.to_vec();
    // std::sort on doubles: a total order on non-NaN inputs.
    uniq.sort_by(|a, b| a.partial_cmp(b).expect("noise: NaN in data"));
    uniq.dedup();
    let un = uniq.len();
    let dvn = data.len();
    if un == 0 {
        return Vec::new();
    }
    let d_un = un - 1;
    let mut result = vec![0.0; dvn];

    if unif_noise {
        if un == dvn {
            // no ties: the original vector, no draw.
            return data.to_vec();
        }
        let mut count: HashMap<u64, i64> = HashMap::with_capacity(un);
        for &v in data {
            *count.entry(key(v)).or_insert(0) += 1;
        }
        for (k, &obs) in data.iter().enumerate() {
            let new_val = if count[&key(obs)] > 1 {
                let li = lower_index(&uniq, obs);
                let (low, high) = if li == 0 {
                    (uniq[li], (uniq[li] + uniq[li + 1]) / 2.0)
                } else if li == d_un {
                    ((uniq[li] + uniq[li - 1]) / 2.0, uniq[li])
                } else {
                    (
                        (uniq[li] + uniq[li - 1]) / 2.0,
                        (uniq[li] + uniq[li + 1]) / 2.0,
                    )
                };
                uniform_real(low, high, rng)
            } else {
                obs
            };
            result[k] = new_val;
        }
    } else {
        // whp order-preserving Gaussian noise centred at each observation, sds from the gaps
        // between unique observations. (The C++ reads uniqVec[1] when there is a single
        // unique value — out of bounds; we take that gap as 0.)
        let mut d_low = vec![0.0; un];
        let mut d_high = vec![0.0; un];
        for j in 0..un {
            if j == 0 {
                d_low[j] = 0.0;
                d_high[j] = if un > 1 { uniq[j + 1] - uniq[j] } else { 0.0 };
            } else if j == d_un {
                d_low[j] = uniq[j] - uniq[j - 1];
                d_high[j] = 0.0;
            } else {
                d_low[j] = uniq[j] - uniq[j - 1];
                d_high[j] = uniq[j + 1] - uniq[j];
            }
        }
        for (k, &obs) in data.iter().enumerate() {
            let li = lower_index(&uniq, obs);
            let sd = if li == 0 {
                d_high[li] / scale_gaussian
            } else if li == d_un {
                d_low[li] / scale_gaussian
            } else {
                (d_high[li] + d_low[li]) / (scale_gaussian * 2.0)
            };
            // cppRnorm: a fresh std::normal_distribution<> per draw.
            result[k] = Normal::new().sample(obs, sd, rng);
        }
    }
    result
}

/// `addNoiseToDataVector(x, scale_gaussian, seed)`: the noised copy of `x` (not sorted).
pub fn add_noise_to_vector(x: &[f64], scale_gaussian: f64, seed: u32) -> Vec<f64> {
    let mut mt = Mt19937::new(seed);
    let with_unif = add_distribution_noise(x, true, scale_gaussian, &mut mt);
    add_distribution_noise(&with_unif, false, scale_gaussian, &mut mt)
}

/// `addNoiseToDataMatrix`: every column noised with its own seed derived from `seed`
/// (`seed + k + 1` for column `k` when `seed > 0`). Returns the noised columns and the seed
/// the C++ leaves in `randomSeed` afterwards (it is incremented per column and reused by the
/// caller for the search). The threaded branch of the C++ yields the same columns.
///
/// `seed == 0` means "a `std::random_device` per column" in the C++ and cannot be reproduced;
/// it is treated like `seed = 0` with the same `+ k + 1` derivation (a deterministic stand-in)
/// and the returned seed is 0, as the C++ leaves it.
pub fn add_noise_to_matrix(
    cols: &[Vec<f64>],
    scale_gaussian: f64,
    seed: u64,
) -> (Vec<Vec<f64>>, u64) {
    let mut current = seed;
    let noised = cols
        .iter()
        .map(|col| {
            current += 1;
            // std::mt19937 mt(unsigned long long) keeps the low 32 bits.
            add_noise_to_vector(col, scale_gaussian, current as u32)
        })
        .collect();
    let seed_after = if seed > 0 { current } else { 0 };
    (noised, seed_after)
}

/// A helper the search uses: a fresh `Mt19937` from a 64-bit seed the way the C++ constructs
/// it (`std::mt19937 mt(rSeed)` with an `unsigned long long` narrowed to the engine's 32-bit
/// `result_type`).
pub fn engine_from_seed(seed: u64) -> Mt19937 {
    Mt19937::new(seed as u32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_ties_uniform_pass_draws_nothing() {
        // With unique values the uniform pass returns the input and the Gaussian pass starts
        // on a fresh engine: the first draw equals a direct cppRnorm on seed 123.
        let x = [1.0, 2.0, 4.0];
        let got = add_noise_to_vector(&x, 4.0, 123);
        let mut mt = Mt19937::new(123);
        let first = Normal::new().sample(1.0, (2.0 - 1.0) / 4.0, &mut mt);
        assert_eq!(got[0], first);
    }

    #[test]
    fn ties_are_replaced_then_jittered() {
        let x = [1.0, 1.0, 3.0];
        let mut mt = Mt19937::new(5);
        let u0 = uniform_real(1.0, 2.0, &mut mt);
        let u1 = uniform_real(1.0, 2.0, &mut mt);
        let after_unif = [u0, u1, 3.0];
        let expect = add_distribution_noise(&after_unif, false, 4.0, &mut mt);
        assert_eq!(add_noise_to_vector(&x, 4.0, 5), expect);
    }

    #[test]
    fn matrix_seeds_columns_in_order() {
        let cols = vec![vec![1.0, 1.0, 2.0], vec![3.0, 5.0, 5.0]];
        let (noised, after) = add_noise_to_matrix(&cols, 4.0, 10);
        assert_eq!(after, 12);
        assert_eq!(noised[0], add_noise_to_vector(&cols[0], 4.0, 11));
        assert_eq!(noised[1], add_noise_to_vector(&cols[1], 4.0, 12));
        let (_, after0) = add_noise_to_matrix(&cols, 4.0, 0);
        assert_eq!(after0, 0);
    }

    #[test]
    fn negative_zero_counts_as_zero() {
        let x = [0.0, -0.0, 1.0];
        let mut mt = Mt19937::new(1);
        let uniq_count = {
            let mut c: HashMap<u64, i64> = HashMap::new();
            for &v in &x {
                *c.entry(key(v)).or_insert(0) += 1;
            }
            c[&key(0.0)]
        };
        assert_eq!(uniq_count, 2);
        let out = add_distribution_noise(&x, true, 4.0, &mut mt);
        assert_ne!(out[0], 0.0);
        assert_eq!(out[2], 1.0);
    }
}
