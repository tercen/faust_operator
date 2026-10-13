//! Bit-exact parity of `faust::tsgates` against values dumped from R FAUST 0.5.5
//! (`dev/make_fixtures.R`, `%.17g`).
//!
//! - `<tag>_<col>_tsGates{0,2,3}.txt` = `faust:::tsGates(sort(col), 0/2/3)`
//! - `<tag>_sub<n>_tsGates0.txt`      = `faust:::tsGates(<tag>_sub<n>_x, 0)`
//! - `<tag>_<col>_tautString.txt`     = `faust:::getTautStringApprox(sort(col))`, which R wrote
//!   with `as.character(<list>)`, i.e. `deparse` at 15 significant digits: that one is compared
//!   after rounding both sides to 15 significant digits.

use std::fs;
use std::path::PathBuf;

use faust::tsgates::{taut_string_pairs, ts_gates};

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/parity/fixtures")
}

fn read_vec(name: &str) -> Vec<f64> {
    let p = fixture_dir().join(name);
    let s = fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    s.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(|l| {
            l.parse::<f64>()
                .unwrap_or_else(|e| panic!("{name}: bad number {l:?}: {e}"))
        })
        .collect()
}

/// Column `col` of a header CSV, sorted ascending.
fn read_sorted_column(name: &str, col: &str) -> Vec<f64> {
    let p = fixture_dir().join(name);
    let s = fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    let mut lines = s.lines();
    let header = lines.next().expect("csv header");
    let j = header
        .split(',')
        .position(|h| h.trim() == col)
        .unwrap_or_else(|| panic!("{name}: no column {col}"));
    let mut v: Vec<f64> = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.split(',')
                .nth(j)
                .and_then(|f| f.trim().parse::<f64>().ok())
                .unwrap_or_else(|| panic!("{name}: bad row {l:?}"))
        })
        .collect();
    v.sort_by(f64::total_cmp);
    v
}

/// `as.character(list(stringX, stringY))`: two `c(...)` deparsed vectors, possibly wrapped.
fn read_deparsed_pair(name: &str) -> (Vec<f64>, Vec<f64>) {
    let p = fixture_dir().join(name);
    let s = fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    let parts: Vec<Vec<f64>> = s
        .split("c(")
        .skip(1)
        .map(|seg| {
            seg.split(')')
                .next()
                .unwrap_or("")
                .split(',')
                .map(str::trim)
                .filter(|t| !t.is_empty())
                .map(|t| {
                    t.parse::<f64>()
                        .unwrap_or_else(|e| panic!("{name}: bad number {t:?}: {e}"))
                })
                .collect()
        })
        .collect();
    assert_eq!(parts.len(), 2, "{name}: expected stringX and stringY");
    (parts[0].clone(), parts[1].clone())
}

fn round15(v: f64) -> f64 {
    format!("{v:.14e}").parse::<f64>().expect("round trip")
}

fn assert_bits_eq(label: &str, got: &[f64], want: &[f64]) {
    assert_eq!(
        got.len(),
        want.len(),
        "{label}: length {} vs expected {}\n got={got:.17?}\nwant={want:.17?}",
        got.len(),
        want.len()
    );
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        assert!(
            g.to_bits() == w.to_bits(),
            "{label}: first mismatch at index {i}: got {g:.17e} want {w:.17e} (diff {:.3e})",
            g - w
        );
    }
}

fn assert_15_eq(label: &str, got: &[f64], want: &[f64]) {
    assert_eq!(
        got.len(),
        want.len(),
        "{label}: length {} vs expected {}",
        got.len(),
        want.len()
    );
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        let g15 = round15(*g);
        assert!(
            g15.to_bits() == w.to_bits(),
            "{label}: first mismatch at index {i}: got {g:.17e} (15 digits {g15:.15e}) want {w:.17e}"
        );
    }
}

fn check_dataset(csv: &str, tag: &str, cols: &[&str]) {
    for col in cols {
        let x = read_sorted_column(csv, col);
        for prior in [0usize, 2, 3] {
            let name = format!("{tag}_{col}_tsGates{prior}.txt");
            let want = read_vec(&name);
            let got = ts_gates(&x, prior);
            assert_bits_eq(&name, &got, &want);
        }
        let name = format!("{tag}_{col}_tautString.txt");
        let (want_x, want_y) = read_deparsed_pair(&name);
        let (got_x, got_y) = taut_string_pairs(&x);
        assert_15_eq(&format!("{name} stringX"), &got_x, &want_x);
        assert_15_eq(&format!("{name} stringY"), &got_y, &want_y);
    }
}

fn check_subsets(tag: &str) {
    for n in [30usize, 100, 300, 400, 500] {
        let x = read_vec(&format!("{tag}_sub{n}_x.txt"));
        assert_eq!(x.len(), n);
        let name = format!("{tag}_sub{n}_tsGates0.txt");
        let want = read_vec(&name);
        let got = ts_gates(&x, 0);
        assert_bits_eq(&name, &got, &want);
    }
}

#[test]
fn test1_columns_bit_exact() {
    check_dataset("test1_matrix.csv", "test1", &["V1", "V2", "V3", "V4", "V5"]);
}

#[test]
fn syn1_columns_bit_exact() {
    check_dataset(
        "syn_s1.csv",
        "syn1",
        &["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"],
    );
}

#[test]
fn test1_subsets_bit_exact() {
    check_subsets("test1");
}

#[test]
fn syn1_subsets_bit_exact() {
    check_subsets("syn1");
}
