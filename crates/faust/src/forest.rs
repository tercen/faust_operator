//! The candidate-cluster search: the annotation forest (phase 1) and scamp's candidate
//! clusters (phase 2) are the same depth-first search with different parameters.
//!
//! Port of `candidateClusterSearch.cpp` + `findCandidateClusters.cpp` (FAUST 0.5.5; scamp
//! 0.5.0 carries the same code), **single-thread path** (`numThreadsTotal <= 2`), without the
//! restricted-value branches (the operator never sets `channelBounds`, so every restriction is
//! 0 and `useRestrictedValue` is false).
//!
//! What the C++ does at every node, and what this port does instead:
//!
//! | C++ | here |
//! |---|---|
//! | copy the active rows of **every** column into a fresh `vector<vector<double>>`, copy it again and sort each column | one global sort order per column; a node's sorted values are read off it through the node's row mask (or, for a small node, by sorting the subset) |
//! | `vector<bool>` row mask copied per child | the same, as a `Vec<bool>` (the traversal order must match the C++ exactly, so the stack discipline is kept) |
//! | `singleDip` + `tsGates` recomputed for the same cell subset reached through another marker ordering | memoised by `(column, hash of the active row set)` — exact, because both are pure functions of the sorted subset when no subsampling is in play |
//!
//! The gate list comes out in the **same order** as the C++ pushes it: R's
//! `.parseAnnotationForest` breaks a tie between cut-point counts by first appearance
//! (`unique()`), so the order is part of the result.
//!
//! Parity: `tests/parity_forest.rs` against the dumped `gateData` of `growAnnotationForest`.
use std::collections::HashMap;

use crate::dip::{double_dip, single_dip, triple_dip};
use crate::noise::{add_noise_to_matrix, engine_from_seed};
use crate::rng::{Mt19937, uniform_int};
use crate::tsgates::ts_gates;

/// One gate the search recorded (`gateInfo`).
#[derive(Debug, Clone, PartialEq)]
pub struct Gate {
    pub column: usize,
    /// Cut points **including the minimum and the maximum** of the node's values.
    pub cut_points: Vec<f64>,
    /// Depth of the gate (root gates are depth 1).
    pub depth: usize,
    /// The parent node's path score (product of `1 - dip p-value` along the path).
    pub path_score: f64,
    /// Cells in the node that was gated.
    pub pop_size: usize,
}

#[derive(Debug, Clone)]
pub struct SearchParams {
    /// `pValueThreshold`: a column is viable when its dip p-value is below this.
    pub dip_threshold: f64,
    /// `minimumClusterSize`: nodes smaller than this are not searched.
    pub min_cluster_size: usize,
    /// `maximumSearchDepth`: nodes deeper than this are not searched (root is depth 0).
    pub max_search_depth: usize,
    /// `maxGateNum` (`maximumAntimodeNumber`): a gating with more cut points is dropped.
    pub max_gate_num: usize,
    /// `allowRepeatedSplitting`.
    pub repeats_allowed: bool,
    /// `randomCandidateSearch`: gate one viable column chosen at random (scamp) instead of
    /// every viable column (annotation forest).
    pub random_search: bool,
    /// `maximumClusterNum` / `maximumGatingNum`: stop once this many candidate clusters
    /// (scamp) or gates (forest) have been recorded.
    pub max_clusters: usize,
    /// `annotationForestSearch`: record gates (true) or candidate clusters (false).
    pub annotation_forest: bool,
    /// `subSamplingThreshold`: above this many cells a node's gates come from subsamples.
    /// FAUST's default is 1e6; the subsample path is not ported yet (see [`search`]).
    pub sub_sample_threshold: usize,
}

impl SearchParams {
    /// The parameters FAUST's `.growForestForExpUnit` uses (`growAnnotationForest`).
    pub fn annotation_forest(annotation_forest_depth: usize) -> Self {
        Self {
            dip_threshold: 0.25,
            min_cluster_size: 25,
            max_search_depth: annotation_forest_depth.saturating_sub(1),
            max_gate_num: 100,
            repeats_allowed: false,
            random_search: false,
            max_clusters: usize::MAX,
            annotation_forest: true,
            sub_sample_threshold: 1_000_000,
        }
    }

    /// The parameters FAUST's `.clusterExpUnitWithScamp` passes to `scamp()` (its defaults:
    /// `maximumSearchDepth = 1e6`, `maximumClusterNum = 50 * ncol`, random candidate search).
    pub fn scamp(n_markers: usize) -> Self {
        Self {
            dip_threshold: 0.25,
            min_cluster_size: 25,
            max_search_depth: 1_000_000,
            max_gate_num: 100,
            repeats_allowed: false,
            random_search: true,
            max_clusters: 50 * n_markers,
            annotation_forest: false,
            sub_sample_threshold: 1_000_000,
        }
    }
}

/// Result of a search: candidate clusters (as sorted row index lists) and/or gates.
#[derive(Debug, Default)]
pub struct SearchResults {
    pub candidates: Vec<Vec<usize>>,
    pub gates: Vec<Gate>,
    /// The C++ `vacuousSearch`: no gate was found at all.
    pub vacuous: bool,
}

/// A node of the depth-first search (`subsetInfo`).
struct Node {
    rows: Vec<bool>,
    cols: Vec<bool>,
    depth: usize,
    path_score: f64,
}

/// Per-column state shared by every node.
struct Column<'a> {
    values: &'a [f64],
    /// Row indices sorted by value (stable, so ties keep row order; values with ties are
    /// equal, so the sorted values are the same whichever order the C++ sort left them in).
    order: Vec<u32>,
}

#[derive(Default)]
struct Memo {
    /// (column, subset hash) → dip p-value of the subset.
    dip: HashMap<(usize, u64), f64>,
    /// (column, subset hash) → cut points of the subset.
    gates: HashMap<(usize, u64), Vec<f64>>,
    hits: usize,
    misses: usize,
}

fn subset_hash(base_index: &[u32]) -> u64 {
    use std::hash::{BuildHasher, Hasher};
    let mut h = ahash::RandomState::with_seeds(1, 2, 3, 4).build_hasher();
    h.write_usize(base_index.len());
    for &i in base_index {
        h.write_u32(i);
    }
    h.finish()
}

/// Sorted active values of one column: walk the global order when the node is large,
/// sort the subset when it is small (the walk is O(n), the sort O(m log m)).
fn sorted_subset(col: &Column, rows: &[bool], base_index: &[u32], out: &mut Vec<f64>) {
    out.clear();
    let m = base_index.len();
    let n = col.values.len();
    if m == 0 {
        return;
    }
    let walk_cost = n as f64;
    let sort_cost = m as f64 * (m as f64).log2().max(1.0) * 2.0;
    if walk_cost <= sort_cost {
        out.extend(
            col.order
                .iter()
                .filter(|&&i| rows[i as usize])
                .map(|&i| col.values[i as usize]),
        );
    } else {
        out.extend(base_index.iter().map(|&i| col.values[i as usize]));
        out.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    }
}

/// `getCutPointsForViableColumn`: the C++ decides the mode prior from the dip tests only for
/// small nodes (≤ 500 cells), then asks the taut string for the gates.
fn cut_points_for_column(
    sorted: &[f64],
    dip_threshold: f64,
    sub_sample_threshold: usize,
) -> Vec<f64> {
    let mut mode_est = 0usize;
    if sorted.len() <= 500 {
        let dd = double_dip(sorted);
        if dd >= dip_threshold / 2.0 {
            mode_est = 2;
        } else if sorted.len() <= 400 {
            let td = triple_dip(sorted);
            mode_est = if td >= dip_threshold / 3.0 { 3 } else { 4 };
        } else {
            mode_est = 3;
        }
    }
    if sorted.len() < sub_sample_threshold {
        ts_gates(sorted, mode_est)
    } else {
        // subSampleTsGates: not ported yet. FAUST's default threshold is 1e6 cells per node;
        // the operator refuses larger samples until this lands (see lib.rs).
        unimplemented!("subsampled taut-string gates (nodes of >= {sub_sample_threshold} cells)")
    }
}

/// `candidateClusterSearch` on noised, column-major data, seeded like the C++ (`std::mt19937
/// mtRCCS(rndSeed)`).
#[allow(clippy::needless_range_loop)]
pub fn search(cols: &[Vec<f64>], params: &SearchParams, rnd_seed: u64) -> SearchResults {
    let n_cols = cols.len();
    let n_rows = cols.first().map_or(0, Vec::len);
    let columns: Vec<Column> = cols
        .iter()
        .map(|v| {
            let mut order: Vec<u32> = (0..v.len() as u32).collect();
            order.sort_by(|&a, &b| {
                v[a as usize]
                    .partial_cmp(&v[b as usize])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            Column { values: v, order }
        })
        .collect();
    let mut rng = engine_from_seed(rnd_seed);
    let mut memo = Memo::default();
    let mut results = SearchResults::default();
    let mut stack: Vec<Node> = Vec::new();

    let mut first = true;
    let mut node = Node {
        rows: vec![true; n_rows],
        cols: vec![true; n_cols],
        depth: 0,
        path_score: 1.0,
    };
    let mut base_index: Vec<u32> = Vec::with_capacity(n_rows);
    let mut sorted = vec![Vec::<f64>::new(); n_cols];
    let mut sorted_hash = vec![0u64; n_cols];

    loop {
        if first {
            first = false;
        } else if let Some(next) = stack.pop() {
            node = next;
        } else {
            break;
        }
        if node.depth > params.max_search_depth {
            continue;
        }
        base_index.clear();
        base_index.extend(
            node.rows
                .iter()
                .enumerate()
                .filter(|(_, b)| **b)
                .map(|(i, _)| i as u32),
        );
        if base_index.len() < params.min_cluster_size {
            continue;
        }
        let h = subset_hash(&base_index);
        // Sorted active values of every column (the C++ builds them for all columns).
        for (j, col) in columns.iter().enumerate() {
            sorted_subset(col, &node.rows, &base_index, &mut sorted[j]);
            sorted_hash[j] = h;
        }
        // getViableColumns
        let mut viable: Vec<usize> = Vec::new();
        let mut viable_scores: Vec<f64> = Vec::new();
        for j in 0..n_cols {
            if node.cols[j] {
                let p = match memo.dip.get(&(j, h)) {
                    Some(&p) => {
                        memo.hits += 1;
                        p
                    }
                    None => {
                        memo.misses += 1;
                        let p = single_dip(&sorted[j]);
                        memo.dip.insert((j, h), p);
                        p
                    }
                };
                if p < params.dip_threshold {
                    viable.push(j);
                    viable_scores.push(p);
                }
            }
        }
        if params.random_search && !viable.is_empty() {
            let k = uniform_int(0, viable.len() as i64 - 1, &mut rng) as usize;
            viable = vec![viable[k]];
            viable_scores = vec![viable_scores[k]];
        }
        if viable.is_empty() {
            if !params.annotation_forest {
                results
                    .candidates
                    .push(base_index.iter().map(|&i| i as usize).collect());
                if results.candidates.len() > params.max_clusters {
                    break;
                }
            }
        } else {
            // updateStackWithViableCols
            let node_depth = node.depth + 1;
            for (i, &col) in viable.iter().enumerate() {
                let cut_points = match memo.gates.get(&(col, h)) {
                    Some(c) => c.clone(),
                    None => {
                        let c = cut_points_for_column(
                            &sorted[col],
                            params.dip_threshold,
                            params.sub_sample_threshold,
                        );
                        memo.gates.insert((col, h), c.clone());
                        c
                    }
                };
                let mut r_cols = node.cols.clone();
                if params.repeats_allowed {
                    r_cols.iter_mut().for_each(|c| *c = true);
                }
                r_cols[col] = false;
                if cut_points.len() > params.max_gate_num {
                    continue;
                }
                results.gates.push(Gate {
                    column: col,
                    cut_points: cut_points.clone(),
                    depth: node_depth,
                    path_score: node.path_score,
                    pop_size: base_index.len(),
                });
                let child_path_score = node.path_score * (1.0 - viable_scores[i]);
                let values = columns[col].values;
                for j in 1..cut_points.len() {
                    let gate_low = cut_points[j - 1];
                    let gate_high = cut_points[j];
                    let mut r_rows = node.rows.clone();
                    for &k in &base_index {
                        let v = values[k as usize];
                        if v <= gate_low || v > gate_high {
                            r_rows[k as usize] = false;
                        }
                        if j == 1 && v == gate_low {
                            r_rows[k as usize] = true;
                        }
                    }
                    stack.push(Node {
                        rows: r_rows,
                        cols: r_cols.clone(),
                        depth: node_depth,
                        path_score: child_path_score,
                    });
                }
            }
        }
        if params.annotation_forest && results.gates.len() > params.max_clusters {
            break;
        }
    }
    if results.gates.is_empty() {
        // The C++ records one "null" gate {0, 0, {0,0}, 0} and one cluster of every row; the
        // null gate is filtered out of gateData by its (unsigned) `numGates - 2`.
        results.candidates.push((0..n_rows).collect());
        results.vacuous = params.random_search;
    }
    let _ = &mut memo;
    results
}

/// `cppGrowAnnotationForest` for one sample, single thread: noise the data with `seed`,
/// search exhaustively, and keep the gates with at most `cut_point_upper_bound` interior cut
/// points (R: `cutPointUpperBound = 2`, so gates with up to 4 cut points in the C++ count,
/// `numGateUB = cutPointUpperBound + 2`). Returns the gates and the seed the C++ leaves in
/// `currentRandomSeed` (the search is seeded with it + 1).
pub fn grow_annotation_forest(
    cols: &[Vec<f64>],
    params: &SearchParams,
    seed: u64,
    gaussian_scale: f64,
    cut_point_upper_bound: usize,
) -> Vec<Gate> {
    let (noised, seed_after) = add_noise_to_matrix(cols, gaussian_scale, seed);
    let res = search(&noised, params, seed_after.wrapping_add(1));
    let num_gate_ub = cut_point_upper_bound + 2;
    res.gates
        .into_iter()
        .filter(|g| g.cut_points.len() >= 2 && g.cut_points.len() - 2 <= num_gate_ub)
        .collect()
}

/// R's `gateData` layout for one sample: per marker, five flat vectors (`cutPoints`,
/// `numCuts`, `cutDepth`, `nodePathScore`, `nodePopSize`), the last four repeated once per
/// interior cut point. This is what `.parseAnnotationForest` reads.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GateData {
    pub cut_points: Vec<f64>,
    pub num_cuts: Vec<f64>,
    pub cut_depth: Vec<f64>,
    pub node_path_score: Vec<f64>,
    pub node_pop_size: Vec<f64>,
}

pub fn gate_data(gates: &[Gate], n_markers: usize) -> Vec<GateData> {
    let mut out = vec![GateData::default(); n_markers];
    for g in gates {
        let interior = g.cut_points.len().saturating_sub(2);
        let d = &mut out[g.column];
        for j in 1..g.cut_points.len().saturating_sub(1) {
            d.cut_points.push(g.cut_points[j]);
            d.num_cuts.push(interior as f64);
            d.cut_depth.push(g.depth as f64);
            d.node_path_score.push(g.path_score);
            d.node_pop_size.push(g.pop_size as f64);
        }
    }
    out
}

/// Expose the engine type for the scamp module.
pub type Engine = Mt19937;
