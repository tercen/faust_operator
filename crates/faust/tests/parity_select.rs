//! Parity of `faust::select` with R FAUST 0.5.5, bit-exact against the dumps in
//! `tests/parity/fixtures/` (R `%.17g`):
//!
//! - `*_forest_gateData` / `*_forest_parsed`: `growAnnotationForest` → `.parseAnnotationForest`
//!   on test1 (500×5), syn1 and syn2 (20000×8);
//! - `test1_faust/`, `syn_faust/`: every intermediate of a full `faust()` run (1 and 2 samples);
//! - `recon_synth/`: `.reconcileAnnotationBoundaries`, `.upConvert`/`.downConvert`,
//!   `.computeElbow`/`fdrtool::gcmlcm`/`pvt.isoMean` and R's numeric primitives driven directly
//!   (`dev/make_select_fixtures.R`), covering the branches the full runs never reach.
use std::collections::HashSet;
use std::path::{Path, PathBuf};

use faust::forest::Gate;
use faust::select::*;

mod common;
use common::fixtures;

fn lines(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .map(str::to_string)
        .collect()
}

fn parse_r(s: &str) -> f64 {
    match s {
        "NA" => f64::NAN,
        _ => s.parse().unwrap_or_else(|_| panic!("not a number: {s:?}")),
    }
}

fn is_r_number(s: &str) -> bool {
    s == "NA" || s.parse::<f64>().is_ok()
}

/// A dumped numeric vector; a leading `names` line (anything not numeric) is skipped.
fn num_vec(path: &Path) -> Vec<f64> {
    let ls = lines(path);
    let skip = usize::from(!ls.is_empty() && !is_r_number(&ls[0]));
    ls.iter().skip(skip).map(|s| parse_r(s)).collect()
}

fn usize_vec(path: &Path) -> Vec<usize> {
    num_vec(path).into_iter().map(|v| v as usize).collect()
}

/// `dumpObj` of a matrix: header `rowname,c1,...`, rows `name,v1,...`.
fn csv_matrix(path: &Path) -> (Vec<String>, Vec<String>, Vec<Vec<f64>>) {
    let ls = lines(path);
    let header: Vec<String> = ls[0].split(',').skip(1).map(str::to_string).collect();
    let mut rows = Vec::new();
    let mut values = Vec::new();
    for l in &ls[1..] {
        let mut it = l.split(',');
        rows.push(it.next().unwrap().to_string());
        values.push(it.map(parse_r).collect());
    }
    (header, rows, values)
}

/// The operator input matrices (`test1_matrix.csv`, `syn_s1.csv`): header of marker names, no
/// row names. FAUST's GatingSet stores float32, so the values are rounded through `f32`.
fn input_matrix_f32(path: &Path) -> (Vec<String>, Vec<Vec<f64>>) {
    let ls = lines(path);
    let header: Vec<String> = ls[0].split(',').map(str::to_string).collect();
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); header.len()];
    for l in &ls[1..] {
        for (j, v) in l.split(',').enumerate() {
            cols[j].push(v.parse::<f64>().unwrap() as f32 as f64);
        }
    }
    (header, cols)
}

/// `annotationMatrix.csv`: no header, integer levels, column-major result.
fn annotation_csv(path: &Path) -> Vec<Vec<u8>> {
    let ls = lines(path);
    let ncol = ls[0].split(',').count();
    let mut cols: Vec<Vec<u8>> = vec![Vec::new(); ncol];
    for l in &ls {
        for (j, v) in l.split(',').enumerate() {
            cols[j].push(v.parse().unwrap());
        }
    }
    cols
}

fn gate_data_dir(dir: &Path, markers: &[&str]) -> Vec<MarkerGateData> {
    markers
        .iter()
        .map(|m| MarkerGateData {
            cut_points: num_vec(&dir.join(format!("{m}_cutPoints.txt"))),
            num_cuts: usize_vec(&dir.join(format!("{m}_numCuts.txt"))),
            cut_depth: usize_vec(&dir.join(format!("{m}_cutDepth.txt"))),
            node_path_score: num_vec(&dir.join(format!("{m}_nodePathScore.txt"))),
            node_pop_size: usize_vec(&dir.join(format!("{m}_nodePopSize.txt"))),
        })
        .collect()
}

fn parsed_dir(dir: &Path, markers: &[&str]) -> Vec<ParsedMarker> {
    markers
        .iter()
        .map(|m| {
            let score = num_vec(&dir.join(m).join("channelScore.txt"));
            let gates = lines(&dir.join(m).join("gates.txt"));
            ParsedMarker {
                channel_score: score[0],
                gates: if gates == ["NA"] {
                    None
                } else {
                    Some(gates.iter().map(|s| parse_r(s)).collect())
                },
            }
        })
        .collect()
}

fn assert_exact(a: f64, b: f64, what: &str) {
    assert!(
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
        "{what}: rust {a:.17e} ({:#x}) vs R {b:.17e} ({:#x})",
        a.to_bits(),
        b.to_bits()
    );
}

fn assert_exact_vec(a: &[f64], b: &[f64], what: &str) {
    assert_eq!(a.len(), b.len(), "{what}: length");
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        assert_exact(*x, *y, &format!("{what}[{i}]"));
    }
}

fn assert_parsed_eq(got: &[ParsedMarker], want: &[ParsedMarker], markers: &[&str], what: &str) {
    assert_eq!(got.len(), want.len());
    for ((g, w), m) in got.iter().zip(want).zip(markers) {
        assert_exact(
            g.channel_score,
            w.channel_score,
            &format!("{what} {m} channelScore"),
        );
        match (&g.gates, &w.gates) {
            (None, None) => {}
            (Some(a), Some(b)) => assert_exact_vec(a, b, &format!("{what} {m} gates")),
            _ => panic!(
                "{what} {m}: gates NA mismatch: {:?} vs {:?}",
                g.gates, w.gates
            ),
        }
    }
}

const V: [&str; 5] = ["V1", "V2", "V3", "V4", "V5"];
const M: [&str; 8] = ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"];

fn strings(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

// ---------------------------------------------------------------------------------------------
// 1. parseAnnotationForest
// ---------------------------------------------------------------------------------------------

/// The forest dumps: `growAnnotationForest` with the operator parameters, parsed with
/// `effectiveChannelSizes = nrow` (no channel bounds). Both input layouts are checked.
#[test]
fn parse_annotation_forest_matches_r_forest_dumps() {
    for (tag, markers, n) in [
        ("test1", &V[..], 500usize),
        ("syn1", &M[..], 20000),
        ("syn2", &M[..], 20000),
    ] {
        let gd = gate_data_dir(&fixtures().join(format!("{tag}_forest_gateData")), markers);
        let want = parsed_dir(&fixtures().join(format!("{tag}_forest_parsed")), markers);
        let sizes = vec![n; markers.len()];
        if tag != "test1" {
            // the trimodal synthetic markers exercise the two-cut path
            assert!(
                gd.iter().any(|m| m.num_cuts.contains(&2)),
                "{tag}: expected some two-cut gates in the forest"
            );
        }
        let got = parse_gate_data(&gd, &sizes);
        assert_parsed_eq(&got, &want, markers, tag);

        // the same forest as `Gate`s, through the C++ gateData assembly
        let gates: Vec<Gate> = gates_from_gate_data(&gd);
        let round_trip = gate_data_from_gates(markers.len(), &gates, CUT_POINT_UPPER_BOUND);
        assert_eq!(round_trip, gd, "{tag}: gateData round trip");
        let got2 = parse_annotation_forest(&strings(markers), &gates, &sizes);
        assert_parsed_eq(&got2, &want, markers, &format!("{tag} via Gate"));
        // a gate with more cuts than the bound is dropped by the C++
        let mut extra = gates.clone();
        extra.push(Gate {
            column: 0,
            cut_points: vec![f64::NEG_INFINITY, 0.1, 0.2, 0.3, f64::INFINITY],
            depth: 1,
            path_score: 1.0,
            pop_size: 10,
        });
        assert_eq!(
            gate_data_from_gates(markers.len(), &extra, CUT_POINT_UPPER_BOUND),
            gd
        );
    }
}

/// `ePop` of every unit of a full run, from R's restriction matrices (`sum(resMat == 0)` per
/// column of the dumped `sampleData/*/resMat.csv`; `dev/make_channel_bounds_fixtures.R`
/// reproduces them from the bounds rule).
fn r_effective_sizes(run: &str) -> Vec<Vec<usize>> {
    match run {
        "test1_faust" => vec![vec![490; 5]],
        "syn_faust" => vec![
            vec![19665, 19666, 19675, 19699, 19675, 19666, 19687, 19691],
            vec![19674, 19672, 19675, 19673, 19639, 19685, 19673, 19664],
        ],
        _ => unreachable!(),
    }
}

/// `channelBoundsUsedByFAUST` (`channelBounds = ""`): per-sample 1%/99% quantiles of the
/// float32 expression, then 5%/95% across samples; the restriction matrices and the effective
/// sizes they imply. `syn_s2.csv` is not in the fixture set, so its per-sample quantiles come
/// from the R reference and its sizes from R's `resMat`.
#[test]
fn channel_bounds_and_effective_sizes_match_r() {
    let cb = fixtures().join("recon_synth/channel_bounds");
    let q_csv = |name: &str| -> (Vec<f64>, Vec<f64>) {
        let ls = lines(&cb.join(name));
        let row = |i: usize| ls[i].split(',').map(parse_r).collect::<Vec<f64>>();
        (row(1), row(2))
    };
    // test1: one sample
    let (_, cols) = input_matrix_f32(&fixtures().join("test1_matrix.csv"));
    let exprs: Vec<&[f64]> = cols.iter().map(Vec::as_slice).collect();
    let (q01, q99) = sample_quantile_bounds(&exprs);
    let want = q_csv("test1_sample_q01_q99.csv");
    assert_exact_vec(&q01, &want.0, "test1 q01");
    assert_exact_vec(&q99, &want.1, "test1 q99");
    let bounds = empirical_channel_bounds(std::slice::from_ref(&exprs)).unwrap();
    let (_, rows, vals) = csv_matrix(
        &fixtures().join("test1_faust/metaData/channelBoundsUsedByFAUST/allSamples.csv"),
    );
    assert_eq!(rows, ["Low", "High"]);
    assert_exact_vec(&bounds.low, &vals[0], "test1 Low");
    assert_exact_vec(&bounds.high, &vals[1], "test1 High");
    let res = restriction_matrix(&exprs, &bounds);
    assert!(any_restricted(&res));
    assert_eq!(effective_sizes(&res), r_effective_sizes("test1_faust")[0]);
    // syn: two samples
    let (_, cols) = input_matrix_f32(&fixtures().join("syn_s1.csv"));
    let exprs: Vec<&[f64]> = cols.iter().map(Vec::as_slice).collect();
    let s1 = sample_quantile_bounds(&exprs);
    let want = q_csv("syn_s1_q01_q99.csv");
    assert_exact_vec(&s1.0, &want.0, "syn s1 q01");
    assert_exact_vec(&s1.1, &want.1, "syn s1 q99");
    let s2 = q_csv("syn_s2_q01_q99.csv");
    let bounds = empirical_channel_bounds_from_quantiles(&[s1, s2]).unwrap();
    let (_, _, vals) =
        csv_matrix(&fixtures().join("syn_faust/metaData/channelBoundsUsedByFAUST/allSamples.csv"));
    assert_exact_vec(&bounds.low, &vals[0], "syn Low");
    assert_exact_vec(&bounds.high, &vals[1], "syn High");
    assert_eq!(
        effective_sizes(&restriction_matrix(&exprs, &bounds)),
        r_effective_sizes("syn_faust")[0]
    );
    // a tiny sample is skipped; nothing left is an error
    let tiny: Vec<&[f64]> = vec![&[1.0, 2.0]];
    assert!(empirical_channel_bounds(&[tiny]).is_err());
}

/// The forests of the full runs (`root_annF/gateData` → `root_pAnnF`), bit-exact with the
/// true effective sizes (490 of 500 cells for test1: the empirical bounds restrict 10 cells
/// per marker).
#[test]
fn parse_annotation_forest_matches_full_run_units() {
    for (run, units, markers) in [
        ("test1_faust", &["sample"][..], &V[..]),
        ("syn_faust", &["s1", "s2"][..], &M[..]),
    ] {
        for (u, unit) in units.iter().enumerate() {
            let unit_dir = fixtures().join(run).join("expUnitData").join(unit);
            let gd = gate_data_dir(&unit_dir.join("root_annF/gateData"), markers);
            let want = parsed_dir(&unit_dir.join("root_pAnnF"), markers);
            let got = parse_gate_data(&gd, &r_effective_sizes(run)[u]);
            assert_parsed_eq(&got, &want, markers, &format!("{run}/{unit}"));
        }
    }
}

#[test]
fn parse_annotation_forest_null_parse_is_na() {
    let gd = vec![MarkerGateData::default(); 3];
    let got = parse_gate_data(&gd, &[10, 10, 10]);
    assert!(
        got.iter()
            .all(|p| p.channel_score.is_nan() && p.gates.is_none())
    );
    assert_eq!(first_depths(&gd), vec![f64::INFINITY; 3]);
    // a marker that never gated in a non-empty forest scores 0
    let gates = vec![Gate {
        column: 1,
        cut_points: vec![0.0, 1.5, 3.0],
        depth: 1,
        path_score: 0.9,
        pop_size: 10,
    }];
    let got = parse_annotation_forest(&strings(&["a", "b", "c"]), &gates, &[10, 10, 10]);
    assert_eq!(got[0].channel_score, 0.0);
    assert_eq!(got[1].channel_score, 1.0);
    assert_eq!(got[1].gates, Some(vec![1.5]));
    assert_eq!(got[2].gates, None);
}

// ---------------------------------------------------------------------------------------------
// 2. selectChannels
// ---------------------------------------------------------------------------------------------

struct FullRun {
    dir: PathBuf,
    units: Vec<String>,
    markers: Vec<String>,
    parsed: Vec<Vec<ParsedMarker>>,
    score_mat: ScoreMat,
}

/// The full-run chain: this port's parse of the dumped forests (exact, see
/// `parse_annotation_forest_matches_full_run_units`) feeds everything downstream.
fn full_run(run: &str, units: &[&str], markers: &[&str]) -> FullRun {
    let dir = fixtures().join(run);
    let mut parsed = Vec::new();
    let mut depths = Vec::new();
    for (u, unit) in units.iter().enumerate() {
        let unit_dir = dir.join("expUnitData").join(unit);
        let gd = gate_data_dir(&unit_dir.join("root_annF/gateData"), markers);
        parsed.push(parse_gate_data(&gd, &r_effective_sizes(run)[u]));
        depths.push(first_depths(&gd));
    }
    let score_mat = ScoreMat::new(strings(units), strings(markers), &parsed, depths);
    FullRun {
        dir,
        units: strings(units),
        markers: strings(markers),
        parsed,
        score_mat,
    }
}

fn runs() -> Vec<FullRun> {
    vec![
        full_run("test1_faust", &["sample"], &V),
        full_run("syn_faust", &["s1", "s2"], &M),
    ]
}

#[test]
fn score_mat_and_select_channels_match_r() {
    for run in runs() {
        let (hdr, rows, vals) = csv_matrix(&run.dir.join("metaData/scoreMat.csv"));
        assert_eq!(hdr, run.markers);
        assert_eq!(rows, run.units);
        for (u, row) in vals.iter().enumerate() {
            assert_exact_vec(
                &run.score_mat.scores[u],
                row,
                &format!("scoreMat {}", run.units[u]),
            );
        }
        let (_, _, dvals) = csv_matrix(&run.dir.join("metaData/depthMat.csv"));
        for (u, row) in dvals.iter().enumerate() {
            assert_exact_vec(
                &run.score_mat.first_depth[u],
                row,
                &format!("depthMat {}", run.units[u]),
            );
        }
        let selected = select_channels(&run.score_mat, 0.01, 0.5).unwrap();
        let got: Vec<String> = selected.iter().map(|&j| run.markers[j].clone()).collect();
        assert_eq!(
            got,
            lines(&run.dir.join("metaData/initSelC.txt")),
            "initSelC"
        );
    }
}

#[test]
fn select_channels_rules() {
    // quantile over units, threshold, decreasing total score, stable ties, NA -> 0
    let parsed = vec![
        vec![
            ParsedMarker {
                channel_score: 0.5,
                gates: None,
            },
            ParsedMarker {
                channel_score: 1.0,
                gates: None,
            },
            ParsedMarker {
                channel_score: f64::NAN,
                gates: None,
            },
            ParsedMarker {
                channel_score: 0.5,
                gates: None,
            },
        ],
        vec![
            ParsedMarker {
                channel_score: 0.0,
                gates: None,
            },
            ParsedMarker {
                channel_score: 0.2,
                gates: None,
            },
            ParsedMarker {
                channel_score: 1.0,
                gates: None,
            },
            ParsedMarker {
                channel_score: 0.5,
                gates: None,
            },
        ],
    ];
    let sm = ScoreMat::new(
        strings(&["u1", "u2"]),
        strings(&["a", "b", "c", "d"]),
        &parsed,
        vec![vec![1.0; 4]; 2],
    );
    assert_eq!(sm.scores[0][2], 0.0);
    // medians: a 0.25, b 0.6, c 0.5, d 0.5; sums: a 0.5, b 1.2, c 1.0, d 1.0
    assert_eq!(select_channels(&sm, 0.3, 0.5).unwrap(), vec![1, 2, 3]);
    assert_eq!(select_channels(&sm, 0.01, 0.5).unwrap(), vec![1, 2, 3, 0]);
    assert!(select_channels(&sm, 0.7, 0.5).is_err());
}

// ---------------------------------------------------------------------------------------------
// 3. reconcileAnnotationBoundaries
// ---------------------------------------------------------------------------------------------

fn res_list_dir(dir: &Path, markers: &[String], units: &[String]) -> Vec<Vec<Vec<f64>>> {
    markers
        .iter()
        .map(|m| {
            units
                .iter()
                .map(|u| num_vec(&dir.join(m).join(format!("{u}.txt"))))
                .collect()
        })
        .collect()
}

fn assert_res_list_eq(got: &ResList, dir: &Path, what: &str) {
    let want = res_list_dir(dir, &got.marker_names(), &got.exp_units);
    for (m, w) in got.markers.iter().zip(&want) {
        for (u, unit) in got.exp_units.iter().enumerate() {
            assert_exact_vec(
                &m.per_unit[u],
                &w[u],
                &format!("{what} {} {unit}", m.marker),
            );
        }
    }
}

fn possibility_list(dir: &Path, markers: &[String]) -> Vec<Vec<(usize, usize)>> {
    markers
        .iter()
        .map(|m| {
            let ls = lines(&dir.join(format!("{m}.txt")));
            let names: Vec<usize> = ls[0].split(',').map(|s| s.parse().unwrap()).collect();
            names
                .into_iter()
                .zip(ls[1..].iter().map(|s| s.parse().unwrap()))
                .collect()
        })
        .collect()
}

#[test]
fn reconcile_matches_full_runs() {
    for run in runs() {
        let selected_idx = select_channels(&run.score_mat, 0.01, 0.5).unwrap();
        let selected: Vec<String> = selected_idx
            .iter()
            .map(|&j| run.markers[j].clone())
            .collect();
        let raw = make_gate_list(&run.parsed, &selected_idx);
        // rawGateList
        for (u, unit) in run.units.iter().enumerate() {
            for (k, m) in selected.iter().enumerate() {
                let want = num_vec(
                    &run.dir
                        .join("gateData/root_rawGateList")
                        .join(unit)
                        .join(format!("{m}.txt")),
                );
                assert_exact_vec(
                    raw[u][k].as_ref().unwrap(),
                    &want,
                    &format!("rawGateList {unit} {m}"),
                );
            }
        }
        let rec = reconcile_boundaries(&run.units, &selected, &raw, &ReconcileOptions::default())
            .unwrap();
        assert_eq!(
            rec.selected_channels,
            lines(&run.dir.join("gateData/root_selectedChannels.txt"))
        );
        assert_eq!(
            rec.possibilities,
            possibility_list(&run.dir.join("metaData/possibilityList"), &selected)
        );
        assert_res_list_eq(
            &rec.res_list,
            &run.dir.join("gateData/root_resListPrep"),
            "resListPrep",
        );
        // no selection list: resList is a copy of resListPrep
        assert_res_list_eq(
            &rec.res_list,
            &run.dir.join("gateData/root_resList"),
            "resList",
        );
    }
}

/// The synthetic scenario of `dev/make_select_fixtures.R`: 7 units in 3 imputation levels,
/// markers with differing gate counts, an outlier, NA units, a tie resolved by preference, a
/// forced marker.
fn synthetic_raw(dir: &Path, units: &[String], markers: &[String]) -> Vec<Vec<Option<Vec<f64>>>> {
    units
        .iter()
        .map(|u| {
            markers
                .iter()
                .map(|m| {
                    let v = num_vec(&dir.join("rawGateList").join(u).join(format!("{m}.txt")));
                    if v.len() == 1 && v[0].is_nan() {
                        None
                    } else {
                        Some(v)
                    }
                })
                .collect()
        })
        .collect()
}

#[test]
fn reconcile_matches_synthetic_scenarios() {
    let units = strings(&["u2", "u1", "u3", "u4", "u6", "u5", "u7"]);
    let imp_h = strings(&["A", "A", "A", "A", "B", "B", "C"]);
    let markers = strings(&["X", "Y", "Z", "W"]);
    for (tag, opts) in [
        (
            "recon_default",
            ReconcileOptions {
                imputation_hierarchy: Some(imp_h.clone()),
                ..Default::default()
            },
        ),
        (
            "recon_supervised",
            ReconcileOptions {
                imputation_hierarchy: Some(imp_h.clone()),
                preference: vec![("W".into(), 2)],
                force: vec![("X".into(), vec![7.5]), ("F".into(), vec![5.0, 6.0])],
            },
        ),
    ] {
        let dir = fixtures().join("recon_synth").join(tag);
        let raw = synthetic_raw(&dir, &units, &markers);
        let rec = reconcile_boundaries(&units, &markers, &raw, &opts).unwrap();
        assert_eq!(
            rec.selected_channels,
            lines(&dir.join("selectedChannels.txt")),
            "{tag}"
        );
        assert_eq!(
            rec.possibilities,
            possibility_list(&dir.join("possibilityList"), &markers),
            "{tag}"
        );
        assert_res_list_eq(&rec.res_list, &dir.join("resListPrep"), tag);
    }
    // the standard counts behind the two runs: ties -> min, preference -> 2
    let dir = fixtures().join("recon_synth/recon_default");
    let raw = synthetic_raw(&dir, &units, &markers);
    let rec = reconcile_boundaries(
        &units,
        &markers,
        &raw,
        &ReconcileOptions {
            imputation_hierarchy: Some(imp_h.clone()),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(rec.num_gates, vec![1, 2, 1, 1]);
    let rec = reconcile_boundaries(
        &units,
        &markers,
        &raw,
        &ReconcileOptions {
            imputation_hierarchy: Some(imp_h),
            preference: vec![("W".into(), 2), ("X".into(), 5)],
            force: vec![],
        },
    )
    .unwrap();
    assert_eq!(
        rec.num_gates,
        vec![1, 2, 1, 2],
        "unobserved preference falls back to the rule"
    );
}

#[test]
fn up_and_down_convert_match_r() {
    let dir = fixtures().join("recon_synth/convert");
    for case in [
        "up1", "up2", "up3", "up4", "down1", "down2", "down3", "down4", "down5",
    ] {
        let from = num_vec(&dir.join(format!("{case}_from.txt")));
        let to = num_vec(&dir.join(format!("{case}_to.txt")));
        let want = num_vec(&dir.join(format!("{case}_out.txt")));
        let got = if case.starts_with("up") {
            up_convert(&from, &to)
        } else {
            down_convert(&from, &to)
        }
        .unwrap();
        assert_exact_vec(&got, &want, case);
    }
    // R stops when the matched set does not have the right size (duplicate standard gates)
    assert!(down_convert(&[1.0, 1.0], &[1.0, 1.0]).is_err());
    assert_eq!(up_convert(&[1.0, 3.0], &[1.0]).unwrap(), vec![1.0]);
}

// ---------------------------------------------------------------------------------------------
// 4. gating, phenotype selection, count matrices, labels
// ---------------------------------------------------------------------------------------------

struct GatedRun {
    run: FullRun,
    rec: Reconciled,
    /// `annotationMatrix` per unit (computed where the expression input is available, else the
    /// R dump), column-major.
    ann: Vec<Vec<Vec<u8>>>,
}

fn gated_run(run: FullRun, inputs: &[Option<&str>]) -> GatedRun {
    let selected_idx = select_channels(&run.score_mat, 0.01, 0.5).unwrap();
    let selected: Vec<String> = selected_idx
        .iter()
        .map(|&j| run.markers[j].clone())
        .collect();
    let raw = make_gate_list(&run.parsed, &selected_idx);
    let rec =
        reconcile_boundaries(&run.units, &selected, &raw, &ReconcileOptions::default()).unwrap();
    let mut ann = Vec::new();
    for (u, unit) in run.units.iter().enumerate() {
        let want = annotation_csv(
            &run.dir
                .join("sampleData")
                .join(unit)
                .join("annotationMatrix.csv"),
        );
        match inputs[u] {
            Some(csv) => {
                let (hdr, cols) = input_matrix_f32(&fixtures().join(csv));
                let exprs: Vec<&[f64]> = selected
                    .iter()
                    .map(|m| cols[hdr.iter().position(|h| h == m).unwrap()].as_slice())
                    .collect();
                let got = annotation_matrix(&exprs, &rec.res_list.thresholds(u));
                assert_eq!(got, want, "annotationMatrix {unit}");
                ann.push(got);
            }
            None => ann.push(want),
        }
    }
    GatedRun { run, rec, ann }
}

fn gated_runs() -> Vec<GatedRun> {
    let mut v = runs().into_iter();
    vec![
        gated_run(v.next().unwrap(), &[Some("test1_matrix.csv")]),
        // syn_s2.csv is not in the fixture set: the R annotation matrix is the gating input
        gated_run(v.next().unwrap(), &[Some("syn_s1.csv"), None]),
    ]
}

fn name_summary(path: &Path) -> Vec<(String, usize)> {
    let ls = lines(path);
    ls[0]
        .split(',')
        .map(str::to_string)
        .zip(ls[1..].iter().map(|s| s.parse().unwrap()))
        .collect()
}

#[test]
fn phenotype_selection_matches_r() {
    for g in gated_runs() {
        let dir = &g.run.dir;
        let per_unit: Vec<Vec<String>> = g
            .run
            .units
            .iter()
            .map(|u| {
                let unit_dir = dir.join("expUnitData").join(u);
                let want = lines(&unit_dir.join("scampClusterLabels.txt"));
                let scamp_csv = dir.join("sampleData").join(u).join("scampAnnotation.csv");
                if scamp_csv.exists() {
                    // the scamp output of the sample: consensus of two identical clusterings
                    let labels = lines(&scamp_csv);
                    let cons = consensus_labels(&labels, &labels).unwrap();
                    assert_eq!(cons, labels);
                    assert_eq!(cluster_names(&cons), want, "scampClusterLabels {u}");
                }
                want
            })
            .collect();
        let sel = phenotype_selection(&per_unit, 0).unwrap();
        assert_eq!(
            sel.name_summary,
            name_summary(&dir.join("metaData/scampNameSummary.txt"))
        );
        assert_eq!(
            sel.elbow,
            num_vec(&dir.join("metaData/phenotypeElbowValue.txt"))[0]
        );
        assert_eq!(
            sel.selected,
            lines(&dir.join("metaData/scampClusterNames.txt"))
        );
        assert_eq!(
            sel.all,
            lines(&dir.join("metaData/allScampClusterNames.txt"))
        );
    }
}

#[test]
fn consensus_and_cluster_names_rules() {
    let run = strings(&["a~1~2~", "b~1~2~", "a~1~2~", "Uncertain"]);
    let max = strings(&["a~1~2~", "a~1~2~", "a~1~2~", "Uncertain"]);
    let cons = consensus_labels(&run, &max).unwrap();
    assert_eq!(
        cons,
        strings(&["a~1~2~", "Uncertain", "a~1~2~", "Uncertain"])
    );
    assert_eq!(cluster_names(&cons), strings(&["a~1~2~"]));
    assert!(consensus_labels(&strings(&["x", "y"]), &strings(&["y", "x"])).is_err());
    // nameOccuranceNum > 0 overrides the elbow
    let per_unit = vec![strings(&["p", "q"]), strings(&["p"]), strings(&["p", "r"])];
    let sel = phenotype_selection(&per_unit, 2).unwrap();
    assert_eq!(sel.selected, strings(&["p"]));
    assert_eq!(sel.all, strings(&["p", "q", "r"]));
    assert_eq!(
        sel.name_summary,
        vec![
            ("p".to_string(), 3),
            ("q".to_string(), 1),
            ("r".to_string(), 1)
        ]
    );
    assert_eq!(sel.elbow, 1.0, "fewer than 4 rows keep everything");
}

#[test]
fn elbow_gcm_and_isomean_match_r() {
    let dir = fixtures().join("recon_synth/elbow");
    for case in ["e1", "e2", "e3", "e4", "e5", "e6"] {
        let y = num_vec(&dir.join(format!("{case}_y.txt")));
        let x: Vec<f64> = (1..=y.len()).map(|v| v as f64).collect();
        let want = num_vec(&dir.join(format!("{case}_elbow.txt")))[0];
        assert_exact(compute_elbow(&x, &y), want, &format!("elbow {case}"));
        if y.len() >= 4 {
            let n = y.len();
            let lo = 2usize.max((0.05 * n as f64).ceil() as usize);
            let hi = (n - 1).min((0.95 * n as f64).floor() as usize);
            let (xk, yk) = gcm(&x[lo - 1..hi], &y[lo - 1..hi]);
            assert_exact_vec(
                &xk,
                &num_vec(&dir.join(format!("{case}_xknots.txt"))),
                &format!("{case} x.knots"),
            );
            assert_exact_vec(
                &yk,
                &num_vec(&dir.join(format!("{case}_yknots.txt"))),
                &format!("{case} y.knots"),
            );
        }
    }
    assert_exact_vec(
        &iso_mean(
            &[3.0, 1.0, 2.0, 2.0, 5.0, 4.0],
            &[1.0, 2.0, 1.0, 1.0, 3.0, 1.0],
        ),
        &num_vec(&dir.join("isomean1.txt")),
        "isomean1",
    );
    assert_exact_vec(
        &iso_mean(
            &[1.5, -2.0, -2.0, 0.25, 7.0, 7.0, 6.5],
            &[2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0],
        ),
        &num_vec(&dir.join("isomean2.txt")),
        "isomean2",
    );
}

#[test]
fn r_numeric_primitives_match_r() {
    let dir = fixtures().join("recon_synth/num");
    let input = |n: &str| num_vec(&dir.join("in").join(format!("{n}.txt")));
    let out = |n: &str| num_vec(&dir.join("out").join(format!("{n}.txt")));
    // long double accumulation: 1e16 cancels and leaves the 64-bit-mantissa remainder
    assert_exact(rstats::sum(&input("sum1")), out("sum1")[0], "sum1");
    assert_exact_vec(&rstats::cumsum(&input("sum1")), &out("cumsum1"), "cumsum1");
    assert_exact(rstats::sum(&input("sum2")), out("sum2")[0], "sum2");
    assert_exact(
        rstats::median(&input("med_even")),
        out("med_even")[0],
        "median even",
    );
    assert_exact(rstats::mad(&input("mad1")), out("mad1")[0], "mad");
    assert_exact(
        rstats::quantile7(&input("q1"), 0.5),
        out("q1_50")[0],
        "quantile 0.5",
    );
    assert_exact(
        rstats::quantile7(&input("q1"), 0.3),
        out("q1_30")[0],
        "quantile 0.3",
    );
    assert_exact(
        rstats::quantile7(&input("q1"), 0.95),
        out("q1_95")[0],
        "quantile 0.95",
    );
    assert_exact(
        rstats::weighted_mean(
            &[3.320695596908827, 3.7970430343892203, 3.8790180350111636],
            &[0.2, 0.3, 0.5],
        ),
        out("wm")[0],
        "weighted.mean",
    );
    assert_exact(QT_975_DF1, out("qt")[0], "qt(0.975, 1)");
    // naive double sums would differ on the first case
    assert_ne!(input("sum1").iter().sum::<f64>(), out("sum1")[0]);
}

#[test]
fn gate_sample_matches_faust_annotation() {
    for g in gated_runs() {
        let dir = &g.run.dir;
        let selected: HashSet<String> = lines(&dir.join("metaData/scampClusterNames.txt"))
            .into_iter()
            .collect();
        let all: HashSet<String> = lines(&dir.join("metaData/allScampClusterNames.txt"))
            .into_iter()
            .collect();
        let markers = g.rec.res_list.marker_names();
        for (u, unit) in g.run.units.iter().enumerate() {
            let gate_nums = g.rec.res_list.gate_nums(u);
            let got = gate_sample(&markers, &g.ann[u], &gate_nums, &selected);
            let want = lines(
                &dir.join("sampleData")
                    .join(unit)
                    .join("faustAnnotation.csv"),
            );
            assert_eq!(got.len(), want.len());
            assert!(got == want, "faustAnnotation {unit}");
            let exhaustive = dir
                .join("sampleData")
                .join(unit)
                .join("exhaustiveFaustAnnotation.csv");
            if exhaustive.exists() {
                let got = gate_sample(&markers, &g.ann[u], &gate_nums, &all);
                assert!(
                    got == lines(&exhaustive),
                    "exhaustiveFaustAnnotation {unit}"
                );
            }
        }
    }
}

fn assert_count_matrix_eq(got: &CountMatrix, path: &Path, what: &str) {
    let (hdr, rows, vals) = csv_matrix(path);
    assert_eq!(got.pretty, hdr, "{what} column names");
    assert_eq!(got.samples, rows, "{what} rows");
    for (s, row) in vals.iter().enumerate() {
        let counts: Vec<f64> = got.counts[s].iter().map(|&c| c as f64).collect();
        assert_eq!(counts, *row, "{what} {}", rows[s]);
    }
}

#[test]
fn count_matrices_and_label_maps_match_r() {
    for g in gated_runs() {
        let dir = &g.run.dir;
        let markers = g.rec.res_list.marker_names();
        let selected_names = lines(&dir.join("metaData/scampClusterNames.txt"));
        let all_names = lines(&dir.join("metaData/allScampClusterNames.txt"));
        let selected: HashSet<String> = selected_names.iter().cloned().collect();
        let all: HashSet<String> = all_names.iter().cloned().collect();
        let mut labels = Vec::new();
        let mut exhaustive = Vec::new();
        for u in 0..g.run.units.len() {
            let gate_nums = g.rec.res_list.gate_nums(u);
            labels.push(gate_sample(&markers, &g.ann[u], &gate_nums, &selected));
            exhaustive.push(gate_sample(&markers, &g.ann[u], &gate_nums, &all));
        }
        let cm = count_matrix(&g.run.units, &labels, &selected_names);
        assert_count_matrix_eq(&cm, &dir.join("faustCountMatrix.csv"), "faustCountMatrix");
        let (_, _, map) = csv_matrix_str(&dir.join("metaData/colNameMap.csv"));
        assert_eq!(
            cm.phenotypes,
            map.iter().map(|r| r[0].clone()).collect::<Vec<_>>()
        );
        assert_eq!(
            cm.pretty,
            map.iter().map(|r| r[1].clone()).collect::<Vec<_>>()
        );
        let ecm = count_matrix(&g.run.units, &exhaustive, &all_names);
        assert_count_matrix_eq(
            &ecm,
            &dir.join("exhaustiveFaustCountMatrix.csv"),
            "exhaustiveFaustCountMatrix",
        );
        let (_, _, emap) = csv_matrix_str(&dir.join("metaData/exhaustiveColNameMap.csv"));
        assert_eq!(
            ecm.pretty,
            emap.iter().map(|r| r[1].clone()).collect::<Vec<_>>()
        );
    }
}

fn csv_matrix_str(path: &Path) -> (Vec<String>, Vec<String>, Vec<Vec<String>>) {
    let ls = lines(path);
    let header: Vec<String> = ls[0].split(',').skip(1).map(str::to_string).collect();
    let mut rows = Vec::new();
    let mut values = Vec::new();
    for l in &ls[1..] {
        let mut it = l.split(',');
        rows.push(it.next().unwrap().to_string());
        values.push(it.map(str::to_string).collect());
    }
    (header, rows, values)
}

/// The `faustLabels` column of `makeAnnotationEmbedding` (the operator's `label`) is the
/// embedding pretty-print of `faustAnnotation`.
#[test]
fn embedding_labels_match_r() {
    let dir = fixtures().join("test1_faust");
    let fa = lines(&dir.join("sampleData/sample/faustAnnotation.csv"));
    let ls = lines(&fixtures().join("test1_annotationEmbedding.csv"));
    let header: Vec<&str> = ls[0].split(',').map(|s| s.trim_matches('"')).collect();
    let col = header.iter().position(|h| *h == "faustLabels").unwrap();
    assert_eq!(ls.len() - 1, fa.len());
    for (l, label) in ls[1..].iter().zip(&fa) {
        let want = l.split(',').nth(col).unwrap().trim_matches('"');
        assert_eq!(embedding_label(label), want);
    }
    assert_eq!(
        embedding_label("A~2~3~B~3~4~C~1~4~D~4~4~"),
        "ADimBMedHighC-DVeryBright"
    );
    assert_eq!(
        pretty_count_name("A~2~3~B~3~4~C~1~4~D~4~4~"),
        "ADimBMed+C-D++"
    );
    assert_eq!(pretty_count_name("A~3~5~"), "A~3~5~");
    assert_eq!(pretty_count_name(NO_PHENOTYPE), NO_PHENOTYPE);
    assert_eq!(level_label(2, 3), "Dim");
}

#[test]
fn annotation_matrix_and_label_rules() {
    let exprs: Vec<&[f64]> = vec![&[0.0, 1.0, 2.0, 3.0], &[5.0, 0.5, 1.5, 1.0]];
    let th: Vec<&[f64]> = vec![&[1.0, 2.5], &[1.0]];
    let ann = annotation_matrix(&exprs, &th);
    assert_eq!(ann, vec![vec![1, 2, 2, 3], vec![2, 1, 2, 2]]);
    let markers = strings(&["A", "B"]);
    assert_eq!(phenotype_label(&markers, &[3, 1], &[3, 2]), "A~3~3~B~1~2~");
    let allowed: HashSet<String> = ["A~2~3~B~1~2~".to_string()].into_iter().collect();
    assert_eq!(
        gate_sample(&markers, &ann, &[3, 2], &allowed),
        strings(&[NO_PHENOTYPE, "A~2~3~B~1~2~", NO_PHENOTYPE, NO_PHENOTYPE])
    );
}
