//! The one R random draw FAUST makes: `set.seed(seedValue); ceiling(runif(1, 10, 1e8))` in
//! `growAnnotationForest`, which is the seed handed to the C++ forest search.
//!
//! R's default generator is Mersenne-Twister with "Inversion" normals; `set.seed(seed)`
//! scrambles the seed through 50 rounds of the LCG `69069 * s + 1`, fills the 625-word table
//! with the same LCG, and marks the twister as "generate on next call" (`dummy[0] = 624`).
//! `unif_rand()` returns `fixup(genrand * 2.3283064365386963e-10)`. R code (GPL-2-or-later),
//! `src/main/RNG.c` and `src/nmath/...`; the same port exists in flowsom-rs.
const N: usize = 624;
const M: usize = 397;

pub struct RMersenneTwister {
    mt: [u32; N],
    mti: usize,
}

impl RMersenneTwister {
    /// `set.seed(seed)` with the default kind.
    pub fn set_seed(seed: u32) -> Self {
        let mut s = seed;
        for _ in 0..50 {
            s = s.wrapping_mul(69069).wrapping_add(1);
        }
        let mut table = [0u32; N + 1];
        for slot in table.iter_mut() {
            s = s.wrapping_mul(69069).wrapping_add(1);
            *slot = s;
        }
        // FixupSeeds: dummy[0] = 624 → mti = N
        let mut mt = [0u32; N];
        mt.copy_from_slice(&table[1..]);
        Self { mt, mti: N }
    }

    fn genrand(&mut self) -> f64 {
        const MATRIX_A: u32 = 0x9908_b0df;
        const UPPER: u32 = 0x8000_0000;
        const LOWER: u32 = 0x7fff_ffff;
        let mag01 = [0u32, MATRIX_A];
        if self.mti >= N {
            let mut kk = 0;
            while kk < N - M {
                let y = (self.mt[kk] & UPPER) | (self.mt[kk + 1] & LOWER);
                self.mt[kk] = self.mt[kk + M] ^ (y >> 1) ^ mag01[(y & 1) as usize];
                kk += 1;
            }
            while kk < N - 1 {
                let y = (self.mt[kk] & UPPER) | (self.mt[kk + 1] & LOWER);
                self.mt[kk] = self.mt[kk + M - N] ^ (y >> 1) ^ mag01[(y & 1) as usize];
                kk += 1;
            }
            let y = (self.mt[N - 1] & UPPER) | (self.mt[0] & LOWER);
            self.mt[N - 1] = self.mt[M - 1] ^ (y >> 1) ^ mag01[(y & 1) as usize];
            self.mti = 0;
        }
        let mut y = self.mt[self.mti];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        f64::from(y) * 2.328_306_436_538_696_3e-10
    }

    /// `unif_rand()`.
    pub fn unif_rand(&mut self) -> f64 {
        let v = self.genrand();
        // fixup: ensure in (0, 1)
        const I2_32M1: f64 = 2.328_306_437_080_797e-10;
        if v <= 0.0 {
            0.5 * I2_32M1
        } else if 1.0 - v <= 0.0 {
            1.0 - 0.5 * I2_32M1
        } else {
            v
        }
    }

    /// `runif(1, min, max)`: `min + (max - min) * unif_rand()`.
    pub fn runif(&mut self, min: f64, max: f64) -> f64 {
        let u = self.unif_rand();
        min + (max - min) * u
    }
}

/// `set.seed(seed); ceiling(runif(1, min = 10, max = 100000000))` — the forest seed FAUST
/// derives from `seedValue` for the (single) forest iteration of every experimental unit.
pub fn forest_seed(seed_value: u32) -> u64 {
    let mut r = RMersenneTwister::set_seed(seed_value);
    r.runif(10.0, 100_000_000.0).ceil() as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forest_seed_matches_r() {
        // R 4.0.4: set.seed(123); ceiling(runif(1, min = 10, max = 1e8)) → 28757760
        assert_eq!(forest_seed(123), 28_757_760);
    }

    #[test]
    fn set_seed_42_first_uniforms() {
        // R: set.seed(42); runif(3) → 0.9148060435 0.9370754133 0.2861395348
        let mut r = RMersenneTwister::set_seed(42);
        let u: Vec<f64> = (0..3).map(|_| r.unif_rand()).collect();
        assert!((u[0] - 0.914806043496355).abs() < 1e-12, "{}", u[0]);
        assert!((u[1] - 0.937075413297862).abs() < 1e-12, "{}", u[1]);
        assert!((u[2] - 0.286139534786344).abs() < 1e-12, "{}", u[2]);
    }
}
