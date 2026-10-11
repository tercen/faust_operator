//! The R half of FAUST: depth scores, marker selection, threshold standardisation across
//! samples, gating cells against the standardised thresholds, and phenotype selection.
//!
//! Port of `R/parseAnnotationForest.R`, `R/selectChannels.R`,
//! `R/reconcileAnnotationBoundaries.R`, `src/gateSample.cpp`, `R/phenotypeSelection.R` and the
//! count/occurrence helpers in `R/discoverPhenotypes.R` / `R/clusterExpUnitsWithScamp.R`.
//!
//! Inputs are the per-sample annotation forests ([`crate::forest::Gate`] lists) and, later,
//! the per-sample scamp labels. Parity: `tests/parity/select.rs` against the dumps of a full
//! `faust()` run (`*_faust/metaData/scoreMat.csv`, `resList`, `selectedChannels`,
//! `expUnitData/*/root_pAnnF`, `sampleData/*/`).
use crate::forest::Gate;

/// `.parseAnnotationForest`: the depth score of every marker of one sample, and the gates it
/// settled on (the cut-point count with the best score and the median cut locations).
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedMarker {
    pub channel_score: f64,
    /// `NA` in R when the marker never gated: `None`.
    pub gates: Option<Vec<f64>>,
}

pub fn parse_annotation_forest(
    markers: &[String],
    gates: &[Gate],
    effective_sizes: &[usize],
) -> Vec<ParsedMarker> {
    let _ = (markers, gates, effective_sizes);
    todo!("port R/parseAnnotationForest.R")
}
