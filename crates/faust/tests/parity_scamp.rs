//! scamp on one sample against R FAUST's `scampAnnotation.csv` (the per-cell label strings
//! `cppNoisyScamp` returns), given R's own selected markers and standardised thresholds
//! (`gateData/root_selectedChannels`, `gateData/root_resList`). FAUST 0.5.5, 1 thread, seed
//! 123, float32 input, empirical channel bounds (restricted values on).
//!
//! MUTATION CONTRACT: a wrong `scale()`, a different random column choice, a different
//! candidate scoring or selection order, a wrong annotation quantile, or dropping the
//! restriction matrix changes labels.
use std::path::PathBuf;

use faust::scamp::{ScampParams, run_scamp};
use faust::select::{any_restricted, empirical_channel_bounds, restriction_matrix};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/parity/fixtures")
}

fn read_matrix_f32(name: &str) -> (Vec<String>, Vec<Vec<f64>>) {
    let text = std::fs::read_to_string(fixtures().join(name)).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split(',')
        .map(|s| s.trim_matches('"').to_string())
        .collect();
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); header.len()];
    for l in lines {
        if l.trim().is_empty() {
            continue;
        }
        for (j, v) in l.split(',').enumerate() {
            cols[j].push(v.parse::<f64>().unwrap() as f32 as f64);
        }
    }
    (header, cols)
}

fn read_lines(path: PathBuf) -> Vec<String> {
    std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("{}", path.display()))
        .lines()
        .map(|l| l.trim().trim_matches('"').to_string())
        .filter(|l| !l.is_empty())
        .collect()
}

/// `run_tag`: the full-run dump; `sample`: the unit to check; `matrices`: every sample of the
/// run (the channel bounds are experiment-wide), with the checked one named by `which`.
fn check(run_tag: &str, sample: &str, matrices: &[&str], which: usize) {
    let all: Vec<(Vec<String>, Vec<Vec<f64>>)> =
        matrices.iter().map(|m| read_matrix_f32(m)).collect();
    let markers = all[0].0.clone();
    let views: Vec<Vec<&[f64]>> = all
        .iter()
        .map(|(_, cols)| cols.iter().map(Vec::as_slice).collect())
        .collect();
    let bounds = empirical_channel_bounds(&views).unwrap();
    let cols = &all[which].1;
    let res_all = restriction_matrix(&views[which], &bounds);

    let gd = fixtures().join(format!("{run_tag}_faust/gateData"));
    let selected = read_lines(gd.join("root_selectedChannels.txt"));
    let idx: Vec<usize> = selected
        .iter()
        .map(|m| markers.iter().position(|x| x == m).unwrap())
        .collect();
    let raw: Vec<Vec<f64>> = idx.iter().map(|&j| cols[j].clone()).collect();
    let res: Vec<Vec<u8>> = idx.iter().map(|&j| res_all[j].clone()).collect();
    let thresholds: Vec<Vec<f64>> = selected
        .iter()
        .map(|m| {
            read_lines(gd.join(format!("root_resList/{m}/{sample}.txt")))
                .iter()
                .filter_map(|l| l.parse::<f64>().ok())
                .collect()
        })
        .collect();
    let want = read_lines(fixtures().join(format!(
        "{run_tag}_faust/sampleData/{sample}/scampAnnotation.csv"
    )));
    let restriction = if any_restricted(&res) {
        Some(res.as_slice())
    } else {
        None
    };
    let r = run_scamp(
        &raw,
        &selected,
        &thresholds,
        restriction,
        &ScampParams::default(),
    );
    assert_eq!(r.labels.len(), want.len(), "{run_tag}/{sample}: cell count");
    let mismatches: Vec<usize> = (0..want.len())
        .filter(|&i| r.labels[i] != want[i])
        .collect();
    if !mismatches.is_empty() {
        let scaled = faust::scamp::r_scale(&raw);
        let th = faust::scamp::thresholds_to_scaled(&raw, &scaled, &thresholds);
        eprintln!("scaled thresholds: {th:?}");
        for &i in &mismatches {
            let cl = r.clustering[i];
            let same: usize = r.clustering.iter().filter(|&&c| c == cl).count();
            let vals: Vec<String> = (0..raw.len())
                .map(|j| format!("{:.6}", scaled[j][i]))
                .collect();
            eprintln!(
                "cell {i}: cluster {cl} (size {same}) got {} want {} scaled {}",
                r.labels[i],
                want[i],
                vals.join(",")
            );
        }
    }
    if let Some(&i) = mismatches.first() {
        panic!(
            "{run_tag}/{sample}: {} of {} labels differ; first at cell {i}: got {} want {} \
             (candidates {}, residual rounds {})",
            mismatches.len(),
            want.len(),
            r.labels[i],
            want[i],
            r.n_candidates,
            r.residual_rounds
        );
    }
    eprintln!(
        "{run_tag}/{sample}: {} labels identical ({} candidates, {} residual rounds)",
        want.len(),
        r.n_candidates,
        r.residual_rounds
    );
}

#[test]
fn test1_scamp_labels_are_identical_to_r() {
    check("test1", "sample", &["test1_matrix.csv"], 0);
}

#[test]
fn syn_s1_scamp_labels_are_identical_to_r() {
    check("syn", "s1", &["syn_s1.csv", "syn_s2.csv"], 0);
}

#[test]
fn syn_s2_scamp_labels_are_identical_to_r() {
    check("syn", "s2", &["syn_s1.csv", "syn_s2.csv"], 1);
}
