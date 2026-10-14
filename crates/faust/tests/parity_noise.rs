//! Bit-exact parity of `faust::noise::add_noise_to_vector` against R FAUST 0.5.5
//! (`faust:::addNoiseToDataVector(col, 4, seed)`, compiled with GCC 9.4 libstdc++), dumped
//! with `%.17g` into `tests/parity/fixtures/<set>_<col>_noise_s<seed>.txt`.
//!
//! The R-exposed function takes the raw (unsorted) column and returns the noised values in
//! the same row order, so the fixture lines align with the CSV rows.
use faust::noise::add_noise_to_vector;
use std::fs;

mod common;
use common::fixtures;

/// Columns of a header-ed CSV, by header name, in file order.
fn read_csv(name: &str) -> (Vec<String>, Vec<Vec<f64>>) {
    let text = fs::read_to_string(fixtures().join(name)).expect("read csv fixture");
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<String> = lines
        .next()
        .expect("csv header")
        .split(',')
        .map(|s| s.trim().trim_matches('"').to_string())
        .collect();
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); header.len()];
    for line in lines {
        for (j, field) in line.split(',').enumerate() {
            cols[j].push(field.trim().parse::<f64>().expect("csv number"));
        }
    }
    (header, cols)
}

fn read_vector(name: &str) -> Vec<f64> {
    fs::read_to_string(fixtures().join(name))
        .unwrap_or_else(|e| panic!("read fixture {name}: {e}"))
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim().parse::<f64>().expect("fixture number"))
        .collect()
}

/// Asserts bit equality; on failure reports the first mismatching index with both values at
/// 17 significant digits.
fn assert_bit_equal(label: &str, got: &[f64], want: &[f64]) {
    assert_eq!(got.len(), want.len(), "{label}: length");
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        if g.to_bits() != w.to_bits() {
            panic!(
                "{label}: first mismatch at index {i}: got {g:.17e} (0x{:016x}) want {w:.17e} (0x{:016x})",
                g.to_bits(),
                w.to_bits()
            );
        }
    }
}

fn check_set(csv: &str, prefix: &str, seeds: &[u32]) {
    let (header, cols) = read_csv(csv);
    for (name, col) in header.iter().zip(&cols) {
        for &seed in seeds {
            let fixture = format!("{prefix}_{name}_noise_s{seed}.txt");
            let want = read_vector(&fixture);
            let got = add_noise_to_vector(col, 4.0, seed);
            assert_bit_equal(&fixture, &got, &want);
        }
    }
}

#[test]
fn test1_matrix_noise_matches_r_seeds_123_124() {
    check_set("test1_matrix.csv", "test1", &[123, 124]);
}

#[test]
fn syn_s1_noise_matches_r_seeds_123_124() {
    check_set("syn_s1.csv", "syn1", &[123, 124]);
}

#[test]
fn matrix_entry_point_reproduces_vector_fixtures() {
    // add_noise_to_matrix seeds column k with seed + k + 1: with seed 122 the first column
    // is the s123 fixture, and with seed 123 it is the s124 fixture.
    let (header, cols) = read_csv("test1_matrix.csv");
    let (noised, after) = faust::noise::add_noise_to_matrix(&cols[..1], 4.0, 122);
    assert_eq!(after, 123);
    let want = read_vector(&format!("test1_{}_noise_s123.txt", header[0]));
    assert_bit_equal("matrix seed 122 col 0", &noised[0], &want);
    let (noised, after) = faust::noise::add_noise_to_matrix(&cols[..1], 4.0, 123);
    assert_eq!(after, 124);
    let want = read_vector(&format!("test1_{}_noise_s124.txt", header[0]));
    assert_bit_equal("matrix seed 123 col 0", &noised[0], &want);
}

#[test]
fn ties_fixture_exercises_the_uniform_pass() {
    // The CSV fixtures have no repeated values, so they never reach the uniform-between-ties
    // rule. `ties_x.txt` (ties at the bottom, inside and at the top of the unique range, plus
    // a -0.0) was noised by faust 0.5.5 itself inside the tercen/faust_operator:0.1.2 image
    // (GCC 9.4): see dev/libstdcxx_oracle/make_ties_fixture.R.
    let x = read_vector("ties_x.txt");
    for seed in [123u32, 7] {
        let fixture = format!("ties_noise_s{seed}.txt");
        let want = read_vector(&fixture);
        let got = add_noise_to_vector(&x, 4.0, seed);
        assert_bit_equal(&fixture, &got, &want);
    }
}
