//! The annotation forest of one sample, gate for gate, against R's `growAnnotationForest`
//! (`gateData`, dumped with %.17g by `dev/make_fixtures.R`; FAUST 0.5.5, 1 thread, seed 123,
//! the operator's parameters).
//!
//! MUTATION CONTRACT — each of these turns a test RED: a different traversal order (the
//! `numCuts` sequence changes), a wrong seed derivation (every cut point moves), a wrong path
//! score or node size, a gate kept or dropped at the wrong cut-point count.
use std::path::PathBuf;

use faust::forest::{SearchParams, gate_data, grow_annotation_forest};
use faust::rrng::forest_seed;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/parity/fixtures")
}

fn read_matrix(name: &str) -> (Vec<String>, Vec<Vec<f64>>) {
    let text = std::fs::read_to_string(fixtures().join(name)).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines.next().unwrap().split(',').map(|s| s.to_string()).collect();
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); header.len()];
    for l in lines {
        if l.trim().is_empty() {
            continue;
        }
        for (j, v) in l.split(',').enumerate() {
            cols[j].push(v.parse().unwrap());
        }
    }
    (header, cols)
}

fn read_vec(path: PathBuf) -> Vec<f64> {
    match std::fs::read_to_string(&path) {
        Ok(t) => t
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.trim().parse().unwrap_or_else(|_| panic!("{}: {l}", path.display())))
            .collect(),
        Err(_) => Vec::new(),
    }
}

fn check(tag: &str, matrix: &str) {
    let (markers, cols) = read_matrix(matrix);
    let params = SearchParams::annotation_forest(3);
    let seed = forest_seed(123);
    let gates = grow_annotation_forest(&cols, &params, seed, 4.0, 2);
    let data = gate_data(&gates, cols.len());
    let dir = fixtures().join(format!("{tag}_forest_gateData"));
    let mut n_checked = 0;
    for (j, m) in markers.iter().enumerate() {
        let want = [
            ("cutPoints", &data[j].cut_points),
            ("numCuts", &data[j].num_cuts),
            ("cutDepth", &data[j].cut_depth),
            ("nodePathScore", &data[j].node_path_score),
            ("nodePopSize", &data[j].node_pop_size),
        ];
        for (what, got) in want {
            let exp = read_vec(dir.join(format!("{m}_{what}.txt")));
            assert_eq!(
                got.len(),
                exp.len(),
                "{tag} {m} {what}: {} gates values vs {} in R",
                got.len(),
                exp.len()
            );
            for (i, (g, e)) in got.iter().zip(&exp).enumerate() {
                assert!(
                    g.to_bits() == e.to_bits(),
                    "{tag} {m} {what}[{i}]: got {g:.17} want {e:.17}"
                );
            }
            n_checked += exp.len();
        }
    }
    assert!(n_checked > 0, "{tag}: no fixture values compared");
    eprintln!("{tag}: {} gates, {n_checked} values identical", gates.len());
}

#[test]
fn test1_forest_is_identical_to_r() {
    check("test1", "test1_matrix.csv");
}

#[test]
fn syn1_forest_is_identical_to_r() {
    check("syn1", "syn_s1.csv");
}

#[test]
fn syn2_forest_is_identical_to_r() {
    check("syn2", "syn_s2.csv");
}
