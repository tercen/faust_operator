//! Bit-exact parity of `faust::dip` against values dumped from R FAUST 0.5.5 (`%.17g`).
//!
//! Fixtures live in `tests/parity/fixtures/`:
//! - `test1_matrix.csv` (V1..V5, 500 rows) and `syn_s1.csv` (M1..M8, 20000 rows): each column
//!   sorted ascending is an input; `<tag>_<col>_singleDip.txt` / `_doubleDip.txt` hold the
//!   p-values.
//! - `<tag>_sub<n>_x.txt` (already sorted) with `_singleDip`, `_doubleDip`, `_tripleDip`
//!   (n ≤ 400 only).

use std::path::{Path, PathBuf};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/parity/fixtures")
}

fn read_vec(name: &str) -> Vec<f64> {
    let p = fixtures().join(name);
    let s = std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    s.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(|l| {
            l.parse::<f64>()
                .unwrap_or_else(|e| panic!("{name}: {l:?}: {e}"))
        })
        .collect()
}

fn read_scalar(name: &str) -> f64 {
    let v = read_vec(name);
    assert_eq!(v.len(), 1, "{name}: expected one value");
    v[0]
}

/// Columns of a header-ed CSV, each sorted ascending.
fn read_sorted_columns(name: &str) -> (Vec<String>, Vec<Vec<f64>>) {
    let p = fixtures().join(name);
    let s = std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    let mut lines = s.lines();
    let header: Vec<String> = lines
        .next()
        .expect("header")
        .split(',')
        .map(|h| h.trim().to_string())
        .collect();
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); header.len()];
    for line in lines.filter(|l| !l.trim().is_empty()) {
        for (j, tok) in line.split(',').enumerate() {
            cols[j].push(tok.trim().parse::<f64>().expect("number"));
        }
    }
    for c in &mut cols {
        c.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    }
    (header, cols)
}

#[track_caller]
fn check(label: &str, got: f64, want: f64) {
    assert!(
        got.to_bits() == want.to_bits(),
        "{label}: mismatch\n  rust = {got:.17e} ({:#018x})\n  R    = {want:.17e} ({:#018x})",
        got.to_bits(),
        want.to_bits()
    );
}

#[test]
fn single_dip_test1_columns() {
    let (names, cols) = read_sorted_columns("test1_matrix.csv");
    assert_eq!(cols[0].len(), 500);
    for (name, col) in names.iter().zip(&cols) {
        let want = read_scalar(&format!("test1_{name}_singleDip.txt"));
        check(
            &format!("test1 {name} singleDip"),
            faust::dip::single_dip(col),
            want,
        );
    }
}

#[test]
fn double_dip_test1_columns() {
    let (names, cols) = read_sorted_columns("test1_matrix.csv");
    for (name, col) in names.iter().zip(&cols) {
        let want = read_scalar(&format!("test1_{name}_doubleDip.txt"));
        check(
            &format!("test1 {name} doubleDip"),
            faust::dip::double_dip(col),
            want,
        );
    }
}

#[test]
fn single_dip_syn1_columns() {
    let (names, cols) = read_sorted_columns("syn_s1.csv");
    assert_eq!(cols[0].len(), 20000);
    for (name, col) in names.iter().zip(&cols) {
        let want = read_scalar(&format!("syn1_{name}_singleDip.txt"));
        check(
            &format!("syn1 {name} singleDip"),
            faust::dip::single_dip(col),
            want,
        );
    }
}

fn subsets(tag: &str) {
    for n in [30usize, 100, 300, 400, 500] {
        let x = read_vec(&format!("{tag}_sub{n}_x.txt"));
        assert_eq!(x.len(), n);
        assert!(
            x.windows(2).all(|w| w[0] <= w[1]),
            "{tag}_sub{n}_x not sorted"
        );

        let want = read_scalar(&format!("{tag}_sub{n}_singleDip.txt"));
        check(
            &format!("{tag} sub{n} singleDip"),
            faust::dip::single_dip(&x),
            want,
        );

        let want = read_scalar(&format!("{tag}_sub{n}_doubleDip.txt"));
        check(
            &format!("{tag} sub{n} doubleDip"),
            faust::dip::double_dip(&x),
            want,
        );

        if n <= 400 {
            let want = read_scalar(&format!("{tag}_sub{n}_tripleDip.txt"));
            check(
                &format!("{tag} sub{n} tripleDip"),
                faust::dip::triple_dip(&x),
                want,
            );
        }
    }
}

#[test]
fn dips_test1_subsets() {
    subsets("test1");
}

#[test]
fn dips_syn1_subsets() {
    subsets("syn1");
}
