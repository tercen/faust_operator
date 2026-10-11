//! The whole FAUST run on a set of samples: what the operator calls.
//!
//! Phase 1 per sample (parallel): noise → exhaustive annotation forest → depth scores.
//! Across samples: marker selection, threshold standardisation. Phase 2 per sample
//! (parallel): scamp on the selected markers → per-cell labels. Across samples: phenotype
//! selection → final labels. Everything here is a composition of the module ports; the
//! numbers it produces are R's.
use crate::forest::{SearchParams, grow_annotation_forest};
use crate::rrng::forest_seed;
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
    /// Indices of the selected markers, in FAUST's order.
    pub selected: Vec<usize>,
    /// Expression of the selected markers per sample, column-major (`[g][j][cell]`).
    pub selected_exprs: Vec<Vec<Vec<f64>>>,
    /// Standardised thresholds per sample per selected marker (`resList`).
    pub thresholds: Vec<Vec<Vec<f64>>>,
    /// Number of annotation levels per selected marker (thresholds + 1).
    pub levels: Vec<usize>,
    /// Annotation level (1-based) of every cell on every selected marker (`[g][cell][j]`).
    pub ann: Vec<Vec<Vec<u8>>>,
    /// Final FAUST label per cell (`[g][cell]`): e.g. `CD3+CD4-CD8+`, or the "no phenotype"
    /// label `0_0_0` (one `0` per selected marker, joined by `_`).
    pub labels: Vec<Vec<String>>,
    /// Depth score of every selected marker in every sample (`[g][j]`), for the embedding.
    pub score_mat: Vec<Vec<f64>>,
    /// Depth score of every marker in every sample (`[g][marker]`).
    pub full_score_mat: Vec<Vec<f64>>,
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
    let params = SearchParams::annotation_forest(3);
    let seed = forest_seed(s.seed);
    progress(&format!(
        "phase 1: annotation forests on {} samples (seed {} → forest seed {seed})",
        samples.len(),
        s.seed
    ));
    let _gates_per_sample: Vec<_> = samples
        .iter()
        .map(|smp| grow_annotation_forest(&smp.cols, &params, seed, 4.0, 2))
        .collect();
    let _ = n_markers;
    todo!("wire select.rs (depth scores, selection, reconciliation) and scamp.rs")
}
