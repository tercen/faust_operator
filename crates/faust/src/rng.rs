//! libstdc++ (GCC 9.4, the compiler of the reference image) random number machinery, so the
//! noise FAUST adds and the draws scamp makes are the same numbers R's C++ produced.
//!
//! What has to match, bit for bit:
//! - `std::mt19937` seeded from an integer (`seed_seq` is NOT used: FAUST constructs
//!   `std::mt19937 mt(seed)` directly, so this is the plain Knuth-style init with 1812433253).
//!   The engine's `result_type` is `uint_fast32_t` but every step is reduced mod 2^32, so a
//!   wider seed is simply truncated to its low 32 bits.
//! - `std::generate_canonical<double, 53>`: `k = ceil(53 / 32) = 2` draws, summed as
//!   `r1 + r2 * 2^32` **in double** (that addition rounds) and divided by `2^64`; a result of
//!   `1.0` is replaced by `nextafter(1, 0)`.
//! - `std::normal_distribution<double>` — libstdc++'s Marsaglia polar method **with the saved
//!   second value** (`_M_saved_available`). NB FAUST's `cppRnorm` constructs a fresh
//!   distribution per call, so in the noise the saved value is never reused ([`Normal::new`]
//!   per draw); scamp code that keeps one distribution alive reuses it.
//! - `std::uniform_real_distribution<double>` — `canonical * (b - a) + a`.
//! - `std::bernoulli_distribution` — `canonical < p` (the adaptor's range is `[0, 1)`).
//! - `std::uniform_int_distribution` — GCC 9's algorithm for a 32-bit URBG
//!   (`bits/uniform_int_dist.h`): the *downscaling* branch with rejection when the requested
//!   range is narrower than the engine's (`scaling = urngrange / uerange`, reject `>= past`),
//!   the pass-through when equal, the *upscaling* recursion when wider. GCC 11+ replaced the
//!   downscaling branch with Lemire's method, which gives **different** numbers: verified by
//!   compiling the same program with GCC 9.4 (the reference image) and GCC 13.3.
//! - `std::shuffle` — GCC 9's `bits/stl_algo.h`: while `urngrange / n >= n` it draws **two**
//!   swap positions per engine call (`__gen_two_uniform_ints`, a `uniform_int_distribution`
//!   over `unsigned long` = 64-bit arithmetic), otherwise one `uniform_int_distribution<
//!   ptrdiff_t>` draw per element. Not a textbook Fisher–Yates.
//!
//! Oracle: every constant in the unit tests below was printed by a C++ program compiled with
//! `g++ (Ubuntu 9.4.0-1ubuntu1~20.04) 9.4.0` inside the `tercen/faust_operator:0.1.2` image.
//! Parity against R itself: `tests/parity_noise.rs` through the noise fixtures.

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;

/// `std::mt19937` as libstdc++ implements it.
#[derive(Debug, Clone)]
pub struct Mt19937 {
    state: [u32; N],
    index: usize,
}

impl Mt19937 {
    /// `std::mt19937 mt(seed)`: `x[0] = seed`, `x[i] = 1812433253 * (x[i-1] ^ (x[i-1] >> 30)) + i`
    /// (mod 2^32).
    pub fn new(seed: u32) -> Self {
        let mut state = [0u32; N];
        state[0] = seed;
        for i in 1..N {
            let prev = state[i - 1];
            state[i] = 1_812_433_253u32
                .wrapping_mul(prev ^ (prev >> 30))
                .wrapping_add(i as u32);
        }
        Self { state, index: N }
    }

    /// The default-constructed engine (`std::mt19937{}`), seed 5489.
    pub fn default_seeded() -> Self {
        Self::new(5489)
    }

    fn gen_rand_all(&mut self) {
        for k in 0..N {
            let y = (self.state[k] & UPPER_MASK) | (self.state[(k + 1) % N] & LOWER_MASK);
            let mut yy = self.state[(k + M) % N] ^ (y >> 1);
            if y & 1 != 0 {
                yy ^= MATRIX_A;
            }
            self.state[k] = yy;
        }
        self.index = 0;
    }

    /// The next 32-bit output.
    pub fn next_u32(&mut self) -> u32 {
        if self.index >= N {
            self.gen_rand_all();
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// `mt.discard(z)`.
    pub fn discard(&mut self, z: u64) {
        for _ in 0..z {
            self.next_u32();
        }
    }

    /// `std::mt19937::min()`.
    pub const MIN: u32 = 0;
    /// `std::mt19937::max()`.
    pub const MAX: u32 = u32::MAX;
}

/// `std::generate_canonical<double, 53, mt19937>`: two draws per value.
pub fn generate_canonical(rng: &mut Mt19937) -> f64 {
    // __r = max - min + 1 = 2^32; __k = 2; the sum is accumulated in double, so the second
    // term's addition rounds exactly as libstdc++'s does.
    let r = 4_294_967_296.0_f64;
    let mut sum = 0.0_f64;
    let mut tmp = 1.0_f64;
    for _ in 0..2 {
        sum += f64::from(rng.next_u32()) * tmp;
        tmp *= r;
    }
    let ret = sum / tmp;
    if ret >= 1.0 {
        // std::nextafter(1.0, 0.0)
        f64::from_bits(1.0_f64.to_bits() - 1)
    } else {
        ret
    }
}

/// `std::normal_distribution<double>` with libstdc++'s cached second value.
#[derive(Debug, Clone, Default)]
pub struct Normal {
    saved: Option<f64>,
}

impl Normal {
    /// A fresh distribution: no saved deviate (what FAUST's `cppRnorm` builds on every call).
    pub fn new() -> Self {
        Self::default()
    }

    /// `std::normal_distribution<>(mean, sd)(rng)` — Marsaglia polar, the unused deviate kept.
    pub fn sample(&mut self, mean: f64, sd: f64, rng: &mut Mt19937) -> f64 {
        let ret = match self.saved.take() {
            Some(v) => v,
            None => {
                let (mut x, mut y, mut r2);
                loop {
                    x = 2.0 * generate_canonical(rng) - 1.0;
                    y = 2.0 * generate_canonical(rng) - 1.0;
                    r2 = x * x + y * y;
                    if !(r2 > 1.0 || r2 == 0.0) {
                        break;
                    }
                }
                let mult = (-2.0 * r2.ln() / r2).sqrt();
                self.saved = Some(x * mult);
                y * mult
            }
        };
        ret * sd + mean
    }

    /// `reset()`: forget the saved deviate.
    pub fn reset(&mut self) {
        self.saved = None;
    }
}

/// `std::uniform_real_distribution<double>(a, b)(rng)`.
pub fn uniform_real(a: f64, b: f64, rng: &mut Mt19937) -> f64 {
    generate_canonical(rng) * (b - a) + a
}

/// `std::bernoulli_distribution(p)(rng)`.
pub fn bernoulli(p: f64, rng: &mut Mt19937) -> bool {
    // (__aurng() - 0) < __p * (1 - 0)
    generate_canonical(rng) < p * 1.0
}

/// `std::uniform_int_distribution<I>(0, urange)` on a 32-bit engine, GCC 9, in the unsigned
/// common type (u64 here: identical to the u32 arithmetic of `<int>` whenever the values fit,
/// and what `<unsigned long>` / `<ptrdiff_t>` actually use).
fn uniform_uint_gcc9(urange: u64, rng: &mut Mt19937) -> u64 {
    const URNGRANGE: u64 = u32::MAX as u64; // max - min
    if URNGRANGE > urange {
        // downscaling
        let uerange = urange + 1; // urange can be zero
        let scaling = URNGRANGE / uerange;
        let past = uerange * scaling;
        loop {
            let ret = u64::from(rng.next_u32());
            if ret < past {
                return ret / scaling;
            }
        }
    } else if URNGRANGE < urange {
        // upscaling: __ret = tmp + urng(); reject on overflow past urange or wrap-around
        let uerngrange = URNGRANGE + 1;
        loop {
            let tmp = uerngrange.wrapping_mul(uniform_uint_gcc9(urange / uerngrange, rng));
            let ret = tmp.wrapping_add(u64::from(rng.next_u32()));
            if !(ret > urange || ret < tmp) {
                return ret;
            }
        }
    } else {
        u64::from(rng.next_u32())
    }
}

/// `std::uniform_int_distribution<int>(a, b)(rng)`, inclusive. `a <= b` required.
pub fn uniform_int(a: i64, b: i64, rng: &mut Mt19937) -> i64 {
    assert!(a <= b, "uniform_int: a > b");
    let urange = (b as u64).wrapping_sub(a as u64);
    let ret = uniform_uint_gcc9(urange, rng);
    (ret.wrapping_add(a as u64)) as i64
}

/// `__gen_two_uniform_ints(b0, b1, g)`: one `uniform_int_distribution<unsigned long>{0, b0*b1-1}`
/// draw split into `(x / b1, x % b1)`.
fn gen_two_uniform_ints(b0: u64, b1: u64, rng: &mut Mt19937) -> (u64, u64) {
    let x = uniform_uint_gcc9(b0 * b1 - 1, rng);
    (x / b1, x % b1)
}

/// `std::shuffle(v.begin(), v.end(), rng)` as libstdc++ (GCC 9) does it.
pub fn shuffle<T>(v: &mut [T], rng: &mut Mt19937) {
    let n = v.len();
    if n == 0 {
        return;
    }
    let urngrange = u32::MAX as u64;
    let urange = n as u64;
    if urngrange / urange >= urange {
        let mut i = 1usize;
        // An even element count means an odd number of elements to swap: do one up front.
        if urange.is_multiple_of(2) {
            let d = uniform_uint_gcc9(1, rng) as usize;
            v.swap(i, d);
            i += 1;
        }
        while i != n {
            let swap_range = i as u64 + 1;
            let (p0, p1) = gen_two_uniform_ints(swap_range, swap_range + 1, rng);
            v.swap(i, p0 as usize);
            i += 1;
            v.swap(i, p1 as usize);
            i += 1;
        }
        return;
    }
    for i in 1..n {
        let d = uniform_uint_gcc9(i as u64, rng) as usize;
        v.swap(i, d);
    }
}

#[cfg(test)]
mod tests {
    // Literals are the oracle's %.17g output, kept verbatim on purpose.
    #![allow(clippy::excessive_precision)]
    // All expected values: g++ (Ubuntu 9.4.0-1ubuntu1~20.04) 9.4.0, libstdc++ of the
    // tercen/faust_operator:0.1.2 image. The mt19937 / canonical / normal / uniform_real /
    // bernoulli values are identical under g++ 13.3.0; uniform_int for ranges 65535, 1000000,
    // 2147483647 and the n = 1001 shuffle differ there (Lemire's method).
    use super::*;

    fn first5(seed: u32) -> Vec<u32> {
        let mut mt = Mt19937::new(seed);
        (0..5).map(|_| mt.next_u32()).collect()
    }

    #[test]
    fn mt19937_first_outputs() {
        assert_eq!(
            first5(5489),
            [3499211612, 581869302, 3890346734, 3586334585, 545404204]
        );
        assert_eq!(
            first5(1),
            [1791095845, 4282876139, 3093770124, 4005303368, 491263]
        );
        assert_eq!(
            first5(123),
            [2991312382, 3062119789, 1228959102, 1840268610, 974319580]
        );
        assert_eq!(
            first5(124),
            [455545294, 2996933276, 3201775633, 4252590215, 2458067906]
        );
        assert_eq!(
            first5(0),
            [2357136044, 2546248239, 3071714933, 3626093760, 2588848963]
        );
        assert_eq!(
            first5(u32::MAX),
            [419326371, 479346978, 3918654476, 2416749639, 3388880820]
        );
    }

    #[test]
    fn mt19937_10000th_output() {
        // The C++ standard's own check: the 10000th output of the default engine.
        let mut mt = Mt19937::default_seeded();
        mt.discard(9999);
        assert_eq!(mt.next_u32(), 4123659995);
    }

    #[test]
    fn mt19937_wide_seed_is_truncated() {
        // std::mt19937 mt(4294967296ULL + 123ULL) == std::mt19937 mt(123)
        let seed: u64 = 4_294_967_296 + 123;
        assert_eq!(Mt19937::new(seed as u32).next_u32(), 2991312382);
    }

    #[test]
    fn canonical_two_draws() {
        let mut mt = Mt19937::new(123);
        let got: Vec<f64> = (0..4).map(|_| generate_canonical(&mut mt)).collect();
        assert_eq!(
            got,
            [
                0.71295532158028085,
                0.42847092502893397,
                0.69088485143772593,
                0.71915030888079545
            ]
        );
    }

    #[test]
    fn normal_fresh_distribution_each_call() {
        let mut mt = Mt19937::new(123);
        let got: Vec<f64> = (0..4)
            .map(|_| Normal::new().sample(0.5, 2.0, &mut mt))
            .collect();
        assert_eq!(
            got,
            [
                -0.63921939894421231,
                2.7217429826816089,
                3.542721613813582,
                3.6906021215853562
            ]
        );
    }

    #[test]
    fn normal_reuses_saved_deviate() {
        let mut mt = Mt19937::new(123);
        let mut d = Normal::new();
        let got: Vec<f64> = (0..4).map(|_| d.sample(0.0, 1.0, &mut mt)).collect();
        assert_eq!(
            got,
            [
                -0.56960969947210616,
                1.6958337120310407,
                1.1108714913408044,
                0.96759407127432384
            ]
        );
    }

    #[test]
    fn uniform_real_values() {
        let mut mt = Mt19937::new(124);
        let got: Vec<f64> = (0..4).map(|_| uniform_real(-1.5, 2.25, &mut mt)).collect();
        assert_eq!(
            got,
            [
                1.1166671387380322,
                2.2129999392306243,
                0.71371710204895944,
                0.60156003868170149
            ]
        );
    }

    #[test]
    fn bernoulli_values() {
        let mut mt = Mt19937::new(7);
        let got: Vec<i32> = (0..16)
            .map(|_| i32::from(bernoulli(0.3, &mut mt)))
            .collect();
        assert_eq!(got, [1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]);
    }

    fn uint8(a: i64, b: i64) -> Vec<i64> {
        let mut mt = Mt19937::new(42);
        (0..8).map(|_| uniform_int(a, b, &mut mt)).collect()
    }

    #[test]
    fn uniform_int_gcc9_downscaling() {
        assert_eq!(uint8(0, 0), [0; 8]);
        assert_eq!(uint8(0, 1), [0, 1, 1, 0, 1, 1, 1, 1]);
        assert_eq!(uint8(0, 6), [2, 5, 6, 1, 5, 5, 4, 4]);
        assert_eq!(uint8(0, 99), [37, 79, 95, 18, 73, 77, 59, 59]);
        assert_eq!(uint8(0, 999), [374, 796, 950, 183, 731, 779, 598, 596]);
        // These three differ under GCC 11+ (Lemire): GCC 9.4 values.
        assert_eq!(
            uint8(0, 65535),
            [24546, 52203, 62306, 12021, 47972, 51098, 39234, 39115]
        );
        assert_eq!(
            uint8(0, 1_000_000),
            [
                374624, 796722, 950928, 183476, 732158, 779866, 598793, 596984
            ]
        );
        assert_eq!(
            uint8(0, 2_147_483_647),
            [
                1608637542, 787846414, 670094950, 1914837113, 669991378, 429389014, 249467210,
                1972458954
            ]
        );
        assert_eq!(uint8(-5, 5), [-1, 3, 5, -3, 3, 3, 1, 1]);
    }

    #[test]
    fn uniform_int_passthrough_and_upscaling_consume_draws() {
        // range == 2^32 - 1: the raw draw.
        let mut mt = Mt19937::new(5489);
        assert_eq!(uniform_int(0, u32::MAX as i64, &mut mt), 3499211612);
        // wider than the engine: at least two draws, in range.
        let mut mt = Mt19937::new(5489);
        let v = uniform_int(0, 1 << 40, &mut mt);
        assert!((0..=(1i64 << 40)).contains(&v));
    }

    fn shuf(n: usize) -> (Vec<i64>, u32) {
        let mut mt = Mt19937::new(99);
        let mut v: Vec<i64> = (0..n as i64).collect();
        shuffle(&mut v, &mut mt);
        v.truncate(12);
        (v, mt.next_u32())
    }

    #[test]
    fn shuffle_gcc9() {
        assert_eq!(shuf(1), (vec![0], 2887414401));
        assert_eq!(shuf(2), (vec![0, 1], 3103284259));
        assert_eq!(shuf(3), (vec![0, 2, 1], 3103284259));
        assert_eq!(shuf(7), (vec![0, 2, 5, 1, 4, 3, 6], 3573189800));
        assert_eq!(shuf(8), (vec![3, 1, 5, 0, 2, 6, 7, 4], 3545474761));
        assert_eq!(
            shuf(100),
            (
                vec![75, 89, 5, 18, 9, 47, 72, 61, 15, 51, 85, 42],
                951667555
            )
        );
        // Differs under GCC 11+ (Lemire inside __gen_two_uniform_ints): GCC 9.4 values.
        assert_eq!(
            shuf(1001),
            (
                vec![496, 29, 228, 520, 436, 723, 763, 764, 11, 575, 749, 131],
                1436834425
            )
        );
        // n*n > 2^32 - 1: the one-draw-per-element branch.
        assert_eq!(
            shuf(70000),
            (
                vec![
                    18435, 36529, 62519, 13358, 7052, 28407, 67709, 27484, 65803, 50148, 41098,
                    57347
                ],
                2625892722
            )
        );
    }

    #[test]
    fn shuffle_empty_is_noop() {
        let mut mt = Mt19937::new(99);
        let mut v: Vec<i64> = vec![];
        shuffle(&mut v, &mut mt);
        assert_eq!(mt.next_u32(), Mt19937::new(99).next_u32());
    }
}
