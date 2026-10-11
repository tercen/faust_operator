//! libstdc++ (GCC 9.4, the compiler of the reference image) random number machinery, so the
//! noise FAUST adds and the draws scamp makes are the same numbers R's C++ produced.
//!
//! What has to match, bit for bit:
//! - `std::mt19937` seeded from an integer (`seed_seq` is NOT used: FAUST constructs
//!   `std::mt19937 mt(seed)` directly, so this is the plain Knuth-style init with 1812433253).
//! - `std::normal_distribution<double>` — libstdc++'s Marsaglia polar method **with the saved
//!   second value** (`_M_saved_available`), and `generate_canonical<double, 53>` which takes
//!   **two** 32-bit draws per double.
//! - `std::uniform_real_distribution<double>` — `a + (b - a) * generate_canonical`.
//! - `std::bernoulli_distribution` — `generate_canonical < p`.
//! - `std::uniform_int_distribution<int>` — libstdc++'s algorithm for a 32-bit URBG: the
//!   downscaling branch with rejection (`__urange < __urngrange`), the exact one in
//!   `bits/uniform_int_dist.h` of GCC 9.
//! - `std::shuffle` — libstdc++'s implementation draws **two** indices per uniform when the
//!   range allows (`__gen_two_uniform_ints`); port that loop, not a textbook Fisher–Yates.
//!
//! Parity: `tests/parity/noise.rs` exercises all of these through the noise fixtures.

/// `std::mt19937` as libstdc++ implements it.
#[derive(Debug, Clone)]
pub struct Mt19937 {
    state: [u32; 624],
    index: usize,
}

impl Mt19937 {
    pub fn new(seed: u32) -> Self {
        let _ = seed;
        todo!()
    }
    /// The next 32-bit output.
    pub fn next_u32(&mut self) -> u32 {
        todo!()
    }
}

/// `std::normal_distribution<double>` with libstdc++'s cached second value.
#[derive(Debug, Clone, Default)]
pub struct Normal {
    saved: Option<f64>,
}

impl Normal {
    pub fn sample(&mut self, mean: f64, sd: f64, rng: &mut Mt19937) -> f64 {
        let _ = (mean, sd, rng);
        todo!()
    }
}

/// `std::generate_canonical<double, 53, mt19937>`: two draws per value.
pub fn generate_canonical(rng: &mut Mt19937) -> f64 {
    let _ = rng;
    todo!()
}

/// `std::uniform_real_distribution<double>(a, b)(rng)`.
pub fn uniform_real(a: f64, b: f64, rng: &mut Mt19937) -> f64 {
    let _ = (a, b, rng);
    todo!()
}

/// `std::bernoulli_distribution(p)(rng)`.
pub fn bernoulli(p: f64, rng: &mut Mt19937) -> bool {
    let _ = (p, rng);
    todo!()
}

/// `std::uniform_int_distribution<int>(a, b)(rng)`, inclusive.
pub fn uniform_int(a: i64, b: i64, rng: &mut Mt19937) -> i64 {
    let _ = (a, b, rng);
    todo!()
}

/// `std::shuffle(v.begin(), v.end(), rng)` as libstdc++ does it.
pub fn shuffle<T>(v: &mut [T], rng: &mut Mt19937) {
    let _ = (v, rng);
    todo!()
}
