//! Hartigan's dip test and FAUST's tabulated p-values.
//!
//! Port of `src/cppdip.cpp` (the dip statistic, from the R `diptest` package's C code) and of
//! `src/singleDip.cpp`, `src/doubleDip.cpp`, `src/tripleDip.cpp`, which look the statistic up in
//! null tables FAUST ships as C++ array literals (`_dip_<k>`, `_dip2v_<k>_`, `_dip3v_<k>_`). The
//! tables are embedded here as binary data (`include_bytes!`), generated from the C++ sources by
//! `dev/extract_dip_tables.py` (run it with `python3 -I dev/extract_dip_tables.py`), so the
//! p-values are the same numbers R returns.
//!
//! Table files (`src/dip_tables/*.f64le`, little-endian IEEE f64, row-major, no header):
//!
//! | file               | shape    | C++ source                                      |
//! |--------------------|----------|-------------------------------------------------|
//! | `dip1_crit.f64le`  | 21 × 26  | `_dip1_crit_vals_` (diptest `qDiptab`)          |
//! | `dip1_rows.f64le`  | 21       | `_dip1_row_names_` (sample sizes 4 … 72000)     |
//! | `dip1_cols.f64le`  | 26       | `_dip1_col_names_` (probabilities 0 … 1)        |
//! | `dip2_crit.f64le`  | 500 × 100| `_dip2_critical_values_` (row N−1 ↔ sample size N) |
//! | `dip3_crit.f64le`  | 400 × 100| `_dip3_critical_values_` (row N−1 ↔ sample size N) |
//!
//! They are decoded once, on first use (`OnceLock`).
//!
//! All inputs are **sorted ascending**; the R-exposed functions take sorted vectors too
//! (`faust:::singleDip(sort(x))`).
//!
//! Deviation from the C++: `cppdip` accumulates the GCM/LCM distance `d`/`dx` in `long double`
//! (x87 80-bit on x86-64 Linux). Those values are only ever *compared* (`dx >= d`,
//! `d < *dip`), never written to the result, and Rust has no 80-bit float, so they are `f64`
//! here. A difference would need two candidate distances within one f64 ulp of each other;
//! none of the fixtures (n up to 20000) shows one.
//!
//! Parity: `tests/parity_dip.rs` — `singleDip`, `doubleDip`, `tripleDip` on the fixtures in
//! `tests/parity/fixtures/*_singleDip.txt` etc. must match to the last bit.

use std::sync::OnceLock;

// ---------------------------------------------------------------------------------------------
// Embedded null tables
// ---------------------------------------------------------------------------------------------

const DIP1_ROWS_N: usize = 21;
const DIP1_COLS_N: usize = 26;
const DIP2_ROWS_N: usize = 500;
const DIP3_ROWS_N: usize = 400;
/// Columns of the double/triple tables: critical values at p = 1.00, 0.99, …, 0.01.
const DIPK_COLS_N: usize = 100;

static DIP1_CRIT_BYTES: &[u8] = include_bytes!("dip_tables/dip1_crit.f64le");
static DIP1_ROWS_BYTES: &[u8] = include_bytes!("dip_tables/dip1_rows.f64le");
static DIP1_COLS_BYTES: &[u8] = include_bytes!("dip_tables/dip1_cols.f64le");
static DIP2_CRIT_BYTES: &[u8] = include_bytes!("dip_tables/dip2_crit.f64le");
static DIP3_CRIT_BYTES: &[u8] = include_bytes!("dip_tables/dip3_crit.f64le");

const _: () = assert!(DIP1_CRIT_BYTES.len() == DIP1_ROWS_N * DIP1_COLS_N * 8);
const _: () = assert!(DIP1_ROWS_BYTES.len() == DIP1_ROWS_N * 8);
const _: () = assert!(DIP1_COLS_BYTES.len() == DIP1_COLS_N * 8);
const _: () = assert!(DIP2_CRIT_BYTES.len() == DIP2_ROWS_N * DIPK_COLS_N * 8);
const _: () = assert!(DIP3_CRIT_BYTES.len() == DIP3_ROWS_N * DIPK_COLS_N * 8);

fn decode_f64le(bytes: &[u8]) -> Vec<f64> {
    bytes
        .chunks_exact(8)
        .map(|c| {
            let mut b = [0u8; 8];
            b.copy_from_slice(c);
            f64::from_le_bytes(b)
        })
        .collect()
}

/// `qDiptab` as FAUST ships it: critical values (rows = sample sizes, cols = probabilities).
struct Dip1Table {
    /// 21 × 26, row-major.
    crit: Vec<f64>,
    /// The 21 tabulated sample sizes (`_dip1_row_names_`).
    rows: Vec<i64>,
    /// The 26 tabulated probabilities (`_dip1_col_names_`).
    cols: Vec<f64>,
}

fn dip1_table() -> &'static Dip1Table {
    static T: OnceLock<Dip1Table> = OnceLock::new();
    T.get_or_init(|| Dip1Table {
        crit: decode_f64le(DIP1_CRIT_BYTES),
        rows: decode_f64le(DIP1_ROWS_BYTES)
            .into_iter()
            .map(|v| v as i64)
            .collect(),
        cols: decode_f64le(DIP1_COLS_BYTES),
    })
}

/// `_dip2_critical_values_`: 500 × 100, row-major.
fn dip2_table() -> &'static [f64] {
    static T: OnceLock<Vec<f64>> = OnceLock::new();
    T.get_or_init(|| decode_f64le(DIP2_CRIT_BYTES))
}

/// `_dip3_critical_values_`: 400 × 100, row-major.
fn dip3_table() -> &'static [f64] {
    static T: OnceLock<Vec<f64>> = OnceLock::new();
    T.get_or_init(|| decode_f64le(DIP3_CRIT_BYTES))
}

// ---------------------------------------------------------------------------------------------
// The dip statistic (cppdip.cpp == diptest's diptst, C translation)
// ---------------------------------------------------------------------------------------------

/// Scratch arrays for [`dip_statistic_with`]. 1-based like the C (`--gcm;` etc.): index 0 is
/// unused, so every array has `n + 1` slots.
#[derive(Default)]
pub struct DipWorkspace {
    gcm: Vec<usize>,
    lcm: Vec<usize>,
    mn: Vec<usize>,
    mj: Vec<usize>,
}

impl DipWorkspace {
    pub fn new() -> Self {
        Self::default()
    }

    fn reserve(&mut self, n: usize) {
        for v in [&mut self.gcm, &mut self.lcm, &mut self.mn, &mut self.mj] {
            if v.len() < n + 1 {
                v.resize(n + 1, 0);
            }
        }
    }
}

/// Hartigan's dip statistic of sorted `x`.
///
/// `cppdip` with `min_is_0 = 0` (FAUST always passes 0). Returns `0.0` where the C++ returns
/// early with `ifault != 0` and leaves the caller's zero-initialised `dip` untouched (empty or
/// unsorted input).
pub fn dip_statistic(x: &[f64]) -> f64 {
    dip_statistic_with(x, &mut DipWorkspace::new())
}

/// [`dip_statistic`] with caller-provided scratch space (the double/triple dips call it
/// O(n²) times).
pub fn dip_statistic_with(x0: &[f64], ws: &mut DipWorkspace) -> f64 {
    let n = x0.len();
    if n == 0 {
        return 0.0; // ifault = 1
    }
    // x[1..=n], as in the C after `--x`.
    let x = |i: usize| x0[i - 1];

    // Check that X is sorted --- if not, return with ifault = 2.
    for k in 2..=n {
        if x(k) < x(k - 1) {
            return 0.0;
        }
    }

    let mut low: usize = 1;
    let mut high: usize = n;
    let two_n = (2 * n) as f64;

    // M.Maechler speedup: work with (2n * dip) everywhere but the very end.
    let mut dip: f64 = 1.0; // min_is_0 == 0
    if n < 2 || x(n) == x(1) {
        return dip / two_n;
    }

    ws.reserve(n);
    let DipWorkspace { gcm, lcm, mn, mj } = ws;

    // Indices mn[1..n] over which combination is necessary for the convex MINORANT (GCM) fit.
    mn[1] = 1;
    for j in 2..=n {
        mn[j] = j - 1;
        loop {
            let mnj = mn[j];
            let mnmnj = mn[mnj];
            if mnj == 1
                || (x(j) - x(mnj)) * ((mnj - mnmnj) as f64)
                    < (x(mnj) - x(mnmnj)) * ((j - mnj) as f64)
            {
                break;
            }
            mn[j] = mnmnj;
        }
    }

    // Indices mj[1..n] over which combination is necessary for the concave MAJORANT (LCM) fit.
    mj[n] = n;
    for k in (1..n).rev() {
        mj[k] = k + 1;
        loop {
            let mjk = mj[k];
            let mjmjk = mj[mjk];
            // (mjk - mjmjk) and (k - mjk) are negative ints in the C; keep them signed.
            if mjk == n
                || (x(k) - x(mjk)) * ((mjk as i64 - mjmjk as i64) as f64)
                    < (x(mjk) - x(mjmjk)) * ((k as i64 - mjk as i64) as f64)
            {
                break;
            }
            mj[k] = mjmjk;
        }
    }

    // ----------------------- Start the cycling. -------------------------------
    loop {
        // Collect the change points for the GCM from HIGH to LOW.
        gcm[1] = high;
        let mut i = 1;
        while gcm[i] > low {
            gcm[i + 1] = mn[gcm[i]];
            i += 1;
        }
        let l_gcm = i; // relevant_length(GCM)
        let mut ig = l_gcm;
        let mut ix = ig - 1; // ix, ig are counters for the convex minorant.

        // Collect the change points for the LCM from LOW to HIGH.
        lcm[1] = low;
        let mut i = 1;
        while lcm[i] < high {
            lcm[i + 1] = mj[lcm[i]];
            i += 1;
        }
        let l_lcm = i; // relevant_length(LCM)
        let mut ih = l_lcm;
        let mut iv = 2; // iv, ih are counters for the concave majorant.

        // Find the largest distance greater than 'DIP' between the GCM and the LCM from LOW
        // to HIGH.  (`long double` in the C++; f64 here — see the module doc.)
        let mut d: f64 = 0.0;
        if l_gcm != 2 || l_lcm != 2 {
            loop {
                let gcmix = gcm[ix];
                let lcmiv = lcm[iv];
                let dx: f64;
                if gcmix > lcmiv {
                    // If the next point of either the GCM or LCM is from the LCM, calculate
                    // the distance here.
                    let gcmi1 = gcm[ix + 1];
                    dx = ((lcmiv as f64 - gcmi1 as f64) + 1.0)
                        - (x(lcmiv) - x(gcmi1)) * ((gcmix - gcmi1) as f64) / (x(gcmix) - x(gcmi1));
                    iv += 1;
                    if dx >= d {
                        d = dx;
                        ig = ix + 1;
                        ih = iv - 1;
                    }
                } else {
                    // If the next point of either the GCM or LCM is from the GCM, calculate
                    // the distance here.
                    let lcmiv1 = lcm[iv - 1];
                    dx = (x(gcmix) - x(lcmiv1)) * ((lcmiv - lcmiv1) as f64)
                        / (x(lcmiv) - x(lcmiv1))
                        - (gcmix as f64 - lcmiv1 as f64 - 1.0);
                    ix -= 1;
                    if dx >= d {
                        d = dx;
                        ig = ix + 1;
                        ih = iv;
                    }
                }
                if ix < 1 {
                    ix = 1;
                }
                if iv > l_lcm {
                    iv = l_lcm;
                }
                if gcm[ix] == lcm[iv] {
                    break;
                }
            }
        } else {
            // l_gcm or l_lcm == 2
            d = 1.0; // min_is_0 == 0
        }

        if d < dip {
            break;
        }

        // Calculate the DIPs for the current LOW and HIGH.

        // The DIP for the convex minorant.
        let mut dip_l = 0.0f64;
        for j in ig..l_gcm {
            let mut max_t = 1.0f64;
            let jb = gcm[j + 1];
            let je = gcm[j];
            if je - jb > 1 && x(je) != x(jb) {
                let c = ((je - jb) as f64) / (x(je) - x(jb));
                for jj in jb..=je {
                    let t = ((jj - jb + 1) as f64) - (x(jj) - x(jb)) * c;
                    if max_t < t {
                        max_t = t;
                    }
                }
            }
            if dip_l < max_t {
                dip_l = max_t;
            }
        }

        // The DIP for the concave majorant.
        let mut dip_u = 0.0f64;
        for j in ih..l_lcm {
            let mut max_t = 1.0f64;
            let jb = lcm[j];
            let je = lcm[j + 1];
            if je - jb > 1 && x(je) != x(jb) {
                let c = ((je - jb) as f64) / (x(je) - x(jb));
                for jj in jb..=je {
                    // (jj - jb - 1) is -1 at jj == jb in the C; keep it signed.
                    let t = (x(jj) - x(jb)) * c - (jj as f64 - jb as f64 - 1.0);
                    if max_t < t {
                        max_t = t;
                    }
                }
            }
            if dip_u < max_t {
                dip_u = max_t;
            }
        }

        // Determine the current maximum.
        let dipnew = if dip_u > dip_l { dip_u } else { dip_l };
        if dip < dipnew {
            dip = dipnew;
        }

        // The following if-clause is NECESSARY (may loop infinitely otherwise)!
        if low == gcm[ig] && high == lcm[ih] {
            break;
        }
        low = gcm[ig];
        high = lcm[ih];
    }

    dip / two_n
}

// ---------------------------------------------------------------------------------------------
// cppApprox.cpp — restricted port of R's approx (linear interpolation, rule = 2)
// ---------------------------------------------------------------------------------------------

/// Linear interpolation of `(x, y)` at `v` with constant extrapolation outside `x`'s range
/// (`cppApprox`, one value). `x` must be non-empty and non-decreasing.
fn approx1(x: &[f64], y: &[f64], v: f64) -> f64 {
    let nxy = x.len();
    let mut i = 0usize;
    let mut j = nxy - 1;
    let cv = v;
    if cv < x[i] {
        return y[0];
    }
    if cv > x[j] {
        return y[nxy - 1];
    }
    // find the correct interval by bisection
    while i + 1 < j {
        // x[i] <= v <= x[j]
        let ij = (i + j) / 2;
        if cv < x[ij] {
            j = ij;
        } else {
            i = ij;
        }
    }
    // provably have i == j-1
    if cv == x[j] {
        y[j]
    } else if cv == x[i] {
        y[i]
    } else {
        y[i] + (y[j] - y[i]) * ((cv - x[i]) / (x[j] - x[i]))
    }
}

// ---------------------------------------------------------------------------------------------
// singleDip.cpp
// ---------------------------------------------------------------------------------------------

/// p-value of unimodality for sorted `x` (FAUST `singleDip`).
///
/// `sqrt(n) * dip` interpolated in diptest's `qDiptab` between the two bracketing tabulated
/// sample sizes (linear in n, after scaling each row by `sqrt` of its n), then linear in the
/// probability columns. Returns `1.0` for `n <= 3`; for `n >= 72000` the 72000 row is used as
/// the asymptote.
pub fn single_dip(x: &[f64]) -> f64 {
    let big_n = x.len();
    let dip = dip_statistic(x);
    let s_dip = dip * (big_n as f64).sqrt();

    if big_n <= 3 {
        return 1.0;
    }
    let t = dip1_table();
    let n_i = big_n as i64;
    let (n0, i0, n1, i1, f_n);
    if n_i >= 72000 {
        // sample size greater than values in qDiptab: use 72000 as an asymptotic value
        n0 = 72000i64;
        i0 = 20usize;
        n1 = 72000i64;
        i1 = 20usize;
        f_n = 0.0f64;
    } else {
        // The C++ scans all 20 brackets and takes the first that holds; n in [4, 72000) always
        // falls in exactly one.
        let mut found = None;
        for i in 0..20 {
            if found.is_none() && n_i < t.rows[i + 1] && n_i >= t.rows[i] {
                found = Some(i);
            }
        }
        let i = match found {
            Some(i) => i,
            None => return 1.0, // unreachable for 4 <= n < 72000
        };
        n0 = t.rows[i];
        i0 = i;
        n1 = t.rows[i + 1];
        i1 = i + 1;
        f_n = ((n_i - n0) as f64) / ((n1 - n0) as f64);
    }
    let y_0 = &t.crit[i0 * DIP1_COLS_N..(i0 + 1) * DIP1_COLS_N];
    let y_1 = &t.crit[i1 * DIP1_COLS_N..(i1 + 1) * DIP1_COLS_N];
    let sq0 = (n0 as f64).sqrt();
    let sq1 = (n1 as f64).sqrt();
    let mut d_v = [0.0f64; DIP1_COLS_N];
    for ((d, &y0), &y1) in d_v.iter_mut().zip(y_0).zip(y_1) {
        let a = y0 * sq0;
        let b = y1 * sq1;
        *d = a + (f_n * (b - a));
    }
    1.0 - approx1(&d_v, &t.cols, s_dip)
}

// ---------------------------------------------------------------------------------------------
// doubleDip.cpp / tripleDip.cpp
// ---------------------------------------------------------------------------------------------

/// The `cvloc` lookup shared by `doubleDip` and `tripleDip`: `crit_values` is the 100-column
/// row for the sample size; returns `cvloc / 100`.
///
/// Quirk kept from the C++: when no bracket matches (`cvloc == -1`, which needs a NaN) the
/// result is `-0.01`. Rows for tiny n are all 1s, so `crit_values[99] <= 1.0` and the p-value is
/// `0.0` there.
fn tabulated_p(crit_values: &[f64], current_dip: f64) -> f64 {
    let mut cvloc: i64 = -1;
    if crit_values[99] <= current_dip {
        cvloc = 0;
    } else if crit_values[0] >= current_dip {
        cvloc = 100;
    } else {
        for i in 0..99 {
            if cvloc == -1 && crit_values[i + 1] > current_dip && crit_values[i] <= current_dip {
                cvloc = 100 - i as i64 - 1;
                break;
            }
        }
    }
    (cvloc as f64) / 100.0
}

/// p-value of "at most two modes" for sorted `x` (FAUST `doubleDip`); used for n ≤ 500.
///
/// Minimises over split points `i in 2..n-2` the larger of the two scaled part dips
/// (`dip(x[..i]) * i / sqrt(n)`, `dip(x[i..]) * (n-i) / sqrt(n)`), then looks the minimum up in
/// the 500 × 100 table. Returns `-1.0` for `n > 500` (unsupported in the C++); `1.0` for an
/// empty input (the C++ reads out of bounds there).
pub fn double_dip(x: &[f64]) -> f64 {
    let big_n = x.len();
    if big_n > DIP2_ROWS_N {
        return -1.0;
    }
    if big_n == 0 {
        return 1.0;
    }
    let sqrt_n = (big_n as f64).sqrt();
    let mut ws = DipWorkspace::new();
    let mut current_dip = 1.0f64;
    // for (int i = 2; i < (N-2); i++)
    for i in 2..big_n.saturating_sub(2) {
        let dip_1 = dip_statistic_with(&x[..i], &mut ws);
        let ratio_1 = (i as f64) / sqrt_n;
        let temp_dip_1 = dip_1 * ratio_1;

        let dip_2 = dip_statistic_with(&x[i..], &mut ws);
        let ratio_2 = ((big_n - i) as f64) / sqrt_n;
        let temp_dip_2 = dip_2 * ratio_2;

        current_dip = current_dip.min(temp_dip_1.max(temp_dip_2));
    }
    let row = &dip2_table()[(big_n - 1) * DIPK_COLS_N..big_n * DIPK_COLS_N];
    tabulated_p(row, current_dip)
}

/// p-value of "at most three modes" for sorted `x` (FAUST `tripleDip`); used for n ≤ 400.
///
/// Minimises over split pairs `3 <= i < n-9`, `i+3 <= j < n-5` the largest of the three scaled
/// part dips, then looks the minimum up in the 400 × 100 table. Returns `-1.0` for `n > 400`
/// (unsupported in the C++); `1.0` for an empty input (the C++ reads out of bounds there).
pub fn triple_dip(x: &[f64]) -> f64 {
    let big_n = x.len();
    if big_n > DIP3_ROWS_N {
        return -1.0;
    }
    if big_n == 0 {
        return 1.0;
    }
    let sqrt_n = (big_n as f64).sqrt();
    let mut ws = DipWorkspace::new();
    let mut current_dip = 1.0f64;
    // for (int i = 3; i < (N-9); i++) for (int j = (i+3); j < (N-5); j++)
    for i in 3..big_n.saturating_sub(9) {
        let dip_1 = dip_statistic_with(&x[..i], &mut ws);
        let ratio_1 = (i as f64) / sqrt_n;
        let temp_dip_1 = dip_1 * ratio_1;
        for j in (i + 3)..big_n.saturating_sub(5) {
            let dip_2 = dip_statistic_with(&x[i..j], &mut ws);
            let ratio_2 = ((j - i) as f64) / sqrt_n;
            let temp_dip_2 = dip_2 * ratio_2;

            let dip_3 = dip_statistic_with(&x[j..], &mut ws);
            let ratio_3 = ((big_n - j) as f64) / sqrt_n;
            let temp_dip_3 = dip_3 * ratio_3;

            current_dip = current_dip.min(temp_dip_1.max(temp_dip_2).max(temp_dip_3));
        }
    }
    let row = &dip3_table()[(big_n - 1) * DIPK_COLS_N..big_n * DIPK_COLS_N];
    tabulated_p(row, current_dip)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tables_decode() {
        let t = dip1_table();
        assert_eq!(t.rows[0], 4);
        assert_eq!(t.rows[20], 72000);
        assert_eq!(t.cols[0], 0.0);
        assert_eq!(t.cols[25], 1.0);
        assert_eq!(t.crit[0], 0.125);
        assert_eq!(dip2_table().len(), DIP2_ROWS_N * DIPK_COLS_N);
        assert_eq!(dip3_table().len(), DIP3_ROWS_N * DIPK_COLS_N);
        // Rows for tiny n are all ones.
        assert!(dip2_table()[..DIPK_COLS_N].iter().all(|&v| v == 1.0));
    }

    #[test]
    fn degenerate_inputs() {
        assert_eq!(dip_statistic(&[]), 0.0);
        assert_eq!(dip_statistic(&[1.0]), 0.5);
        assert_eq!(dip_statistic(&[2.0, 2.0, 2.0]), 1.0 / 6.0);
        assert_eq!(dip_statistic(&[2.0, 1.0]), 0.0); // unsorted: ifault = 2
        assert_eq!(single_dip(&[1.0, 2.0, 3.0]), 1.0);
        assert_eq!(double_dip(&vec![0.0; 501]), -1.0);
        assert_eq!(triple_dip(&vec![0.0; 401]), -1.0);
        // all-ones table rows: p = 0 (C++ quirk)
        assert_eq!(double_dip(&[1.0, 2.0, 3.0]), 0.0);
        assert_eq!(triple_dip(&[1.0, 2.0, 3.0]), 0.0);
    }

    #[test]
    fn approx_matches_r_rule_2() {
        let x = [0.0, 1.0, 2.0];
        let y = [10.0, 20.0, 30.0];
        assert_eq!(approx1(&x, &y, -1.0), 10.0);
        assert_eq!(approx1(&x, &y, 5.0), 30.0);
        assert_eq!(approx1(&x, &y, 1.0), 20.0);
        assert_eq!(approx1(&x, &y, 0.5), 15.0);
    }
}
