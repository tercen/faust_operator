//! Hartigan's dip test and FAUST's tabulated p-values.
//!
//! Port of `src/cppdip.cpp` (the dip statistic, from the R `diptest` package's C code) and of
//! `src/singleDip.cpp`, `src/doubleDip.cpp`, `src/tripleDip.cpp`, which look the statistic up in
//! null tables FAUST ships as C++ array literals (`_dip2v_<k>_` …). The tables are embedded
//! here as binary data (`include_bytes!`), generated from the C++ sources by
//! `dev/extract_dip_tables.py`, so the p-values are the same numbers R returns.
//!
//! All inputs are **sorted ascending**; the R-exposed functions take sorted vectors too.
//!
//! Parity: `tests/parity/dip.rs` — `singleDip`, `doubleDip`, `tripleDip` on the fixtures in
//! `tests/parity/fixtures/*_singleDip.txt` etc. must match to the last bit.

/// The dip statistic of a sorted sample (`diptest`'s `diptst`), plus what FAUST keeps of it.
pub struct DipResult {
    pub dip: f64,
}

/// Hartigan's dip statistic of sorted `x`.
pub fn dip_statistic(x: &[f64]) -> f64 {
    let _ = x;
    todo!("port cppdip.cpp")
}

/// p-value of unimodality for sorted `x` (FAUST `singleDip`).
pub fn single_dip(x: &[f64]) -> f64 {
    let _ = x;
    todo!("port singleDip.cpp")
}

/// p-value of "at most two modes" for sorted `x` (FAUST `doubleDip`); used for n ≤ 500.
pub fn double_dip(x: &[f64]) -> f64 {
    let _ = x;
    todo!("port doubleDip.cpp")
}

/// p-value of "at most three modes" for sorted `x` (FAUST `tripleDip`); used for n ≤ 400.
pub fn triple_dip(x: &[f64]) -> f64 {
    let _ = x;
    todo!("port tripleDip.cpp")
}
