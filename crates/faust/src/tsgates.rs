//! Gate locations from the taut-string density and 1-D k-medians.
//!
//! Port of `src/cpPmden.cpp` + `src/tautstring.cpp` + `src/kkuiper.cpp` (Davies–Kovac taut
//! string density estimate with the Kuiper-based bandwidth search), `src/kMedDP.cpp` +
//! `src/fill_SMAWK.cpp` + `src/dynamic_prog.cpp` (exact 1-D k-medians by dynamic programming),
//! `src/tsGates.cpp` (`tsGates`: antimodes of the taut string, else k-medians boundaries) and
//! `src/tsModeEstimate.cpp`, `src/getTautStringApprox.cpp`.
//!
//! All inputs are **sorted ascending**. Outputs are the cut points **including the minimum and
//! the maximum** of the data, exactly as the C++ returns them.
//!
//! Parity: `tests/parity/tsgates.rs` against `*_tsGates0.txt`, `*_tsGates2.txt`,
//! `*_tsGates3.txt` and `*_tautString.txt`, bit for bit.

/// The taut string of sorted `x` (`getTautStringApprox`): one value per observation.
pub fn taut_string_approx(x: &[f64]) -> Vec<f64> {
    let _ = x;
    todo!("port cpPmden.cpp / tautstring.cpp / kkuiper.cpp")
}

/// Cut points for sorted `x`. `mode_prior` 0 lets the taut string decide (2 centres when it
/// sees one mode); `k > 0` forces `k` k-medians centres (FAUST passes 2/3/4 from the dip tests
/// when n ≤ 500).
pub fn ts_gates(x: &[f64], mode_prior: usize) -> Vec<f64> {
    let _ = (x, mode_prior);
    todo!("port tsGates.cpp")
}

/// 1-D k-medians by dynamic programming (`kMedDP`): the cluster index (0-based) of each
/// observation of sorted `x`.
pub fn k_med_dp(x: &[f64], k: usize) -> Vec<usize> {
    let _ = (x, k);
    todo!("port kMedDP.cpp + fill_SMAWK.cpp + dynamic_prog.cpp")
}

/// Mode locations the taut string estimates (`tsModeEstimate`).
pub fn ts_mode_estimate(x: &[f64]) -> Vec<f64> {
    let _ = x;
    todo!("port tsModeEstimate.cpp")
}
