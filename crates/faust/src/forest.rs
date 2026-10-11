//! The annotation forest and the candidate-cluster search (shared by both FAUST phases).
//!
//! Port of `src/findCandidateClusters.cpp` + `src/candidateClusterSearch.cpp` (FAUST) — the
//! same code exists in scamp. Implemented by the operator author; the primitives it calls are
//! [`crate::dip`], [`crate::tsgates`], [`crate::noise`].
//!
//! Design (same results as the C++ on the single-thread path, less work):
//! - every column is sorted **once**; a node is a row bitmask and its sorted active values are
//!   a filtered walk of the global order (the C++ copies and re-sorts every active column at
//!   every node);
//! - gates are memoised by (column, row mask): the exhaustive forest reaches the same subset
//!   through every ordering of its markers, and `tsGates` is a pure function of the sorted
//!   subset when no subsampling is in play.
//!
//! Parity: `tests/parity/forest.rs` against the dumped `gateData` of `growAnnotationForest`.

/// One gate the search recorded (`gateInfo`).
#[derive(Debug, Clone, PartialEq)]
pub struct Gate {
    pub column: usize,
    /// Cut points including min and max.
    pub cut_points: Vec<f64>,
    pub depth: usize,
    pub path_score: f64,
    pub pop_size: usize,
}

#[derive(Debug, Clone)]
pub struct ForestParams {
    pub dip_threshold: f64,
    pub min_cluster_size: usize,
    pub max_search_depth: usize,
    pub max_gate_num: usize,
    pub repeats_allowed: bool,
    pub random_search: bool,
    pub max_clusters: usize,
}

/// Result of a search: candidate clusters (as row index lists) and/or gates.
#[derive(Debug, Default)]
pub struct SearchResults {
    pub candidates: Vec<Vec<usize>>,
    pub gates: Vec<Gate>,
    pub vacuous: bool,
}
