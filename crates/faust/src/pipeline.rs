//! The whole FAUST run on a set of samples: what the operator calls.
//!
//! The order is `faust()`'s with `annotationsApproved = TRUE`, `channelBounds = ""`,
//! `experimentalUnit = ""` (one experimental unit per sample), `nameOccuranceNum = 0`:
//!
//! 1. round the expression to float32 (the GatingSet stores single precision);
//! 2. empirical channel bounds across samples → per-sample restriction matrices;
//! 3. per sample: noise, exhaustive annotation forest, depth scores (`.growAnnForest`);
//! 4. marker selection (`.selectChannels`) and threshold standardisation
//!    (`.reconcileAnnotationBoundaries`);
//! 5. per sample: scamp on the selected markers with the standardised thresholds, the
//!    consensus labels and the cluster names (`.clusterExpUnitWithScamp`);
//! 6. phenotype selection across samples (elbow), then every sample's cells gated against
//!    the thresholds and labelled when their phenotype was selected (`gateSample`).
//!
//! Samples run in parallel in steps 3 and 5 (rayon); each sample's work is the single-thread
//! C++ path, so the result does not depend on the thread count.
use std::collections::HashSet;

use rayon::prelude::*;

use crate::forest::{SearchParams, grow_annotation_forest};
use crate::rrng::forest_seed;
use crate::scamp::{ScampParams, run_scamp};
use crate::select::{
    CUT_POINT_UPPER_BOUND, ParsedMarker, ReconcileOptions, ScoreMat, annotation_matrix,
    any_restricted, cluster_names, consensus_labels, effective_sizes, embedding_label,
    empirical_channel_bounds, first_depths, gate_data_from_gates, gate_sample, make_gate_list,
    parse_gate_data, phenotype_selection, reconcile_boundaries, restriction_matrix,
    select_channels,
};
use crate::{FaustError, Result};

/// One sample: all markers, column-major (`cols[j][cell]`).
pub struct Sample {
    pub name: String,
    pub cols: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct FaustSettings {
    pub selection_quantile: f64,
    pub depth_score_threshold: f64,
    /// FAUST `seedValue`.
    pub seed: u32,
    /// Worker threads for the per-sample loops; 0 = all cores.
    pub threads: usize,
}

/// What comes out, per sample (outer index) where applicable.
pub struct FaustRun {
    /// Indices of the selected markers, in FAUST's order (`selectedChannels`).
    pub selected: Vec<usize>,
    /// Expression (float32-rounded) of the selected markers per sample, column-major
    /// (`[g][j][cell]`).
    pub selected_exprs: Vec<Vec<Vec<f64>>>,
    /// Standardised thresholds per sample per selected marker (`resList`).
    pub thresholds: Vec<Vec<Vec<f64>>>,
    /// Number of annotation levels per selected marker (thresholds + 1).
    pub levels: Vec<usize>,
    /// Annotation level (1-based) of every cell on every selected marker (`[g][cell][j]`).
    pub ann: Vec<Vec<Vec<u8>>>,
    /// Final label per cell (`[g][cell]`), as the embedding prints it (`CD3+CD4-…`), or
    /// `0_0_0_0_0` when the cell's phenotype was not selected.
    pub labels: Vec<Vec<String>>,
    /// Depth score of every selected marker in every sample (`[g][j]`).
    pub score_mat: Vec<Vec<f64>>,
    /// Depth score of every marker in every sample (`[g][marker]`).
    pub full_score_mat: Vec<Vec<f64>>,
    /// Phenotypes kept, in the internal `marker~level~levels~` form.
    pub phenotypes: Vec<String>,
    /// The occurrence threshold the elbow rule chose.
    pub elbow: f64,
}

fn pool(threads: usize) -> Result<rayon::ThreadPool> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| FaustError::Invalid(format!("thread pool: {e}")))
}

/// Run FAUST. `progress` receives one line per stage.
pub fn run_faust(
    samples: &[Sample],
    markers: &[String],
    s: &FaustSettings,
    progress: &mut dyn FnMut(&str),
) -> Result<FaustRun> {
    if samples.is_empty() {
        return Err(FaustError::Invalid("no samples".into()));
    }
    let n_markers = markers.len();
    let n_samples = samples.len();
    let pool = pool(s.threads)?;

    // 1. float32, as the GatingSet stores it
    let data: Vec<Vec<Vec<f64>>> = samples
        .iter()
        .map(|smp| {
            smp.cols
                .iter()
                .map(|c| c.iter().map(|&x| x as f32 as f64).collect())
                .collect()
        })
        .collect();

    // 2. channel bounds and restriction matrices
    let views: Vec<Vec<&[f64]>> = data
        .iter()
        .map(|d| d.iter().map(Vec::as_slice).collect())
        .collect();
    let bounds = empirical_channel_bounds(&views)?;
    let res: Vec<Vec<Vec<u8>>> = views
        .iter()
        .map(|v| restriction_matrix(v, &bounds))
        .collect();

    // 3. forests + depth scores
    let seed = forest_seed(s.seed);
    progress(&format!(
        "annotation forests: {n_samples} samples x {n_markers} markers (seed {} -> {seed})",
        s.seed
    ));
    let params = SearchParams::annotation_forest(3);
    let forests: Vec<(Vec<ParsedMarker>, Vec<f64>)> = pool.install(|| {
        data.par_iter()
            .zip(res.par_iter())
            .map(|(d, r)| {
                let gates = grow_annotation_forest(
                    d,
                    if any_restricted(r) {
                        Some(r.as_slice())
                    } else {
                        None
                    },
                    &params,
                    seed,
                    4.0,
                    CUT_POINT_UPPER_BOUND,
                );
                let gd = gate_data_from_gates(n_markers, &gates, CUT_POINT_UPPER_BOUND);
                let parsed = parse_gate_data(&gd, &effective_sizes(r));
                let depths = first_depths(&gd);
                (parsed, depths)
            })
            .collect()
    });
    let parsed: Vec<_> = forests.iter().map(|f| f.0.clone()).collect();
    let depths: Vec<_> = forests.iter().map(|f| f.1.clone()).collect();
    let unit_names: Vec<String> = samples.iter().map(|smp| smp.name.clone()).collect();
    let score_mat = ScoreMat::new(unit_names.clone(), markers.to_vec(), &parsed, depths);

    // 4. selection + reconciliation
    let selected = select_channels(&score_mat, s.depth_score_threshold, s.selection_quantile)?;
    let selected_names: Vec<String> = selected.iter().map(|&j| markers[j].clone()).collect();
    progress(&format!(
        "selected {} of {n_markers} markers: {}",
        selected.len(),
        selected_names.join(", ")
    ));
    let raw = make_gate_list(&parsed, &selected);
    let rec = reconcile_boundaries(
        &unit_names,
        &selected_names,
        &raw,
        &ReconcileOptions::default(),
    )?;
    if rec.selected_channels != selected_names {
        return Err(FaustError::Invalid(
            "reconciliation changed the selected markers without a force list".into(),
        ));
    }
    let thresholds: Vec<Vec<Vec<f64>>> = (0..n_samples)
        .map(|u| {
            rec.res_list
                .thresholds(u)
                .iter()
                .map(|t| t.to_vec())
                .collect()
        })
        .collect();
    let gate_nums: Vec<usize> = rec.res_list.gate_nums(0);
    progress(&format!(
        "thresholds per marker: {}",
        selected_names
            .iter()
            .zip(&gate_nums)
            .map(|(m, n)| format!("{m}:{}", n - 1))
            .collect::<Vec<_>>()
            .join(" ")
    ));

    // 5. scamp per sample on the selected markers
    let selected_exprs: Vec<Vec<Vec<f64>>> = data
        .iter()
        .map(|d| selected.iter().map(|&j| d[j].clone()).collect())
        .collect();
    let selected_res: Vec<Vec<Vec<u8>>> = res
        .iter()
        .map(|r| selected.iter().map(|&j| r[j].clone()).collect())
        .collect();
    progress("scamp clustering per sample");
    let scamp_params = ScampParams {
        seed: s.seed as u64,
        ..ScampParams::default()
    };
    let per_unit_labels: Vec<Result<Vec<String>>> = pool.install(|| {
        (0..n_samples)
            .into_par_iter()
            .map(|u| {
                let r = run_scamp(
                    &selected_exprs[u],
                    &selected_names,
                    &thresholds[u],
                    if any_restricted(&selected_res[u]) {
                        Some(selected_res[u].as_slice())
                    } else {
                        None
                    },
                    &scamp_params,
                );
                // one iteration: run == max, the consensus keeps every label
                consensus_labels(&r.labels, &r.labels)
            })
            .collect()
    });
    let mut per_unit_names = Vec::with_capacity(n_samples);
    for l in per_unit_labels {
        per_unit_names.push(cluster_names(&l?));
    }

    // 6. phenotype selection and gating
    let sel = phenotype_selection(&per_unit_names, 0)?;
    progress(&format!(
        "{} phenotypes found, {} kept (occurrence threshold {:.1})",
        sel.all.len(),
        sel.selected.len(),
        sel.elbow
    ));
    let allowed: HashSet<String> = sel.selected.iter().cloned().collect();
    let mut ann = Vec::with_capacity(n_samples);
    let mut labels = Vec::with_capacity(n_samples);
    for u in 0..n_samples {
        let exprs: Vec<&[f64]> = selected_exprs[u].iter().map(Vec::as_slice).collect();
        let th: Vec<&[f64]> = thresholds[u].iter().map(Vec::as_slice).collect();
        let cols = annotation_matrix(&exprs, &th);
        let gated = gate_sample(&selected_names, &cols, &gate_nums, &allowed);
        labels.push(gated.iter().map(|l| embedding_label(l)).collect::<Vec<_>>());
        // cell-major levels for the embedding
        let n = cols.first().map_or(0, Vec::len);
        ann.push(
            (0..n)
                .map(|i| cols.iter().map(|c| c[i]).collect::<Vec<u8>>())
                .collect::<Vec<_>>(),
        );
    }
    let score_sel: Vec<Vec<f64>> = score_mat
        .scores
        .iter()
        .map(|row| selected.iter().map(|&j| row[j]).collect())
        .collect();
    Ok(FaustRun {
        selected,
        selected_exprs,
        thresholds,
        levels: gate_nums,
        ann,
        labels,
        score_mat: score_sel,
        full_score_mat: score_mat.scores.clone(),
        phenotypes: sel.selected,
        elbow: sel.elbow,
    })
}
