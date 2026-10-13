//! The whole run (`pipeline::run_faust`) against R FAUST's final per-cell labels
//! (`sampleData/<unit>/faustAnnotation.csv`, printed as the embedding does) and its selected
//! markers and thresholds. FAUST 0.5.5, `faust()` with the operator's parameters, seed 123.
use std::path::PathBuf;

use faust::pipeline::{FaustSettings, Sample, run_faust};
use faust::select::embedding_label;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/parity/fixtures")
}

fn read_matrix(name: &str) -> (Vec<String>, Vec<Vec<f64>>) {
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
            cols[j].push(v.parse().unwrap());
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

fn check(run_tag: &str, units: &[&str], matrices: &[&str]) {
    let mut markers = Vec::new();
    let samples: Vec<Sample> = units
        .iter()
        .zip(matrices)
        .map(|(u, m)| {
            let (h, cols) = read_matrix(m);
            markers = h;
            Sample {
                name: u.to_string(),
                cols,
            }
        })
        .collect();
    let s = FaustSettings {
        selection_quantile: 0.5,
        depth_score_threshold: 0.01,
        seed: 123,
        threads: 2,
    };
    let mut log = Vec::new();
    let run = run_faust(&samples, &markers, &s, &mut |m| log.push(m.to_string())).unwrap();
    for l in &log {
        eprintln!("  {l}");
    }
    let gd = fixtures().join(format!("{run_tag}_faust/gateData"));
    let want_sel = read_lines(gd.join("root_selectedChannels.txt"));
    let got_sel: Vec<String> = run.selected.iter().map(|&j| markers[j].clone()).collect();
    assert_eq!(got_sel, want_sel, "{run_tag}: selected markers");
    for (u, unit) in units.iter().enumerate() {
        for (k, m) in want_sel.iter().enumerate() {
            let want: Vec<f64> = read_lines(gd.join(format!("root_resList/{m}/{unit}.txt")))
                .iter()
                .filter_map(|l| l.parse::<f64>().ok())
                .collect();
            let got = &run.thresholds[u][k];
            assert_eq!(
                got.len(),
                want.len(),
                "{run_tag}/{unit} {m}: threshold count"
            );
            for (g, w) in got.iter().zip(&want) {
                assert!(
                    g.to_bits() == w.to_bits(),
                    "{run_tag}/{unit} {m}: threshold {g:.17} vs {w:.17}"
                );
            }
        }
        let want: Vec<String> = read_lines(fixtures().join(format!(
            "{run_tag}_faust/sampleData/{unit}/faustAnnotation.csv"
        )))
        .iter()
        .map(|l| embedding_label(l))
        .collect();
        let got = &run.labels[u];
        assert_eq!(got.len(), want.len(), "{run_tag}/{unit}: cell count");
        let n_diff = got.iter().zip(&want).filter(|(g, w)| g != w).count();
        if n_diff > 0 {
            let i = got.iter().zip(&want).position(|(g, w)| g != w).unwrap();
            panic!(
                "{run_tag}/{unit}: {n_diff} of {} labels differ; first at {i}: got {} want {}",
                want.len(),
                got[i],
                want[i]
            );
        }
        eprintln!("{run_tag}/{unit}: {} labels identical", want.len());
    }
}

#[test]
fn test1_pipeline_is_identical_to_r() {
    check("test1", &["sample"], &["test1_matrix.csv"]);
}

#[test]
fn syn_pipeline_is_identical_to_r() {
    check("syn", &["s1", "s2"], &["syn_s1.csv", "syn_s2.csv"]);
}
