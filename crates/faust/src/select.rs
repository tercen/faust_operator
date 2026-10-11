//! The R half of FAUST: depth scores, marker selection, threshold standardisation across
//! samples, gating cells against the standardised thresholds, and phenotype selection.
//!
//! Port of `R/processChannelBounds.R` + `R/makeRestrictionMatrices.R` (the empirical channel
//! bounds and the restriction matrix the forest, scamp and the depth score all consume),
//! `R/parseAnnotationForest.R`, `R/selectChannels.R`, `R/reconcileAnnotationBoundaries.R`,
//! `R/mkAnnMats.R`, `src/gateSample.cpp`, `R/clusterExpUnitsWithScamp.R` (consensus labels,
//! the elbow, the name summary), `R/getFaustCountMatrix.R` /
//! `R/getExhaustiveFaustCountMatrix.R` and the label pretty-printing of
//! `R/makeAnnotationEmbedding.R`.
//!
//! Pipeline order, per experimental unit: [`empirical_channel_bounds`] (all samples) →
//! [`restriction_matrix`] → [`effective_sizes`] → forest → [`parse_annotation_forest`] →
//! [`ScoreMat`] / [`select_channels`] (all units) → [`make_gate_list`] →
//! [`reconcile_boundaries`] (all units) → [`annotation_matrix`] → scamp →
//! [`consensus_labels`] / [`cluster_names`] → [`phenotype_selection`] (all units) →
//! [`gate_sample`] → [`count_matrix`] / [`embedding_label`].
//!
//! Inputs are the per-sample annotation forests ([`crate::forest::Gate`] lists, or R's
//! flattened `gateData` layout — [`MarkerGateData`]) and, later, the per-sample scamp labels.
//! Parity: `tests/parity_select.rs` against the dumps of full `faust()` runs
//! (`*_faust/metaData/scoreMat.csv`, `resList`, `selectedChannels`, `expUnitData/*/root_pAnnF`,
//! `sampleData/*/`), the forest dumps, and synthetic reconciliation scenarios driven through
//! the R functions.
//!
//! FAUST reads the expression matrix back from a `GatingSet`, which stores **float32**: every
//! value the R pipeline sees is the input rounded to `f32` (`x as f32 as f64`); the thresholds,
//! quantiles and annotation levels are computed on those values.
//!
//! # R numeric semantics
//!
//! Every value is compared bit-exactly with the R dumps, so the arithmetic follows R:
//! `sum`, `mean` and `cumsum` accumulate in 80-bit `long double` (x86_64 R), emulated in
//! software by [`rstats::Ext80`]; `median`, `mad` and `quantile` (type 7) are ported from
//! `stats`. The remaining operations are plain IEEE doubles in the same order as the R code.
//!
//! # String ordering
//!
//! R sorts phenotype labels with the locale collation (`en_US.UTF-8` in the operator image).
//! All labels of one run have the same shape (`marker~level~levels~` per selected marker, in
//! the same marker order and with the same number of levels), so the collation reduces to a
//! digit-by-digit comparison at fixed positions, which is the byte order used here. Sample
//! names (`sort(unique(...))` of the imputation hierarchy) are also byte-ordered: the order
//! only decides tie-breaks among samples with identical gate counts.
use std::collections::{BTreeMap, HashMap, HashSet};

use crate::forest::Gate;
use crate::{FaustError, Result};

/// `cutPointUpperBound` of `growAnnotationForest` in FAUST: gates with more cut points are
/// dropped when the C++ assembles `gateData`.
pub const CUT_POINT_UPPER_BOUND: usize = 2;
/// The label of a cell whose gated annotation is not a selected phenotype.
pub const NO_PHENOTYPE: &str = "0_0_0_0_0";
/// The scamp label of a cell whose two clusterings disagree.
pub const UNCERTAIN: &str = "Uncertain";
/// `stats::qt(0.975, df = 1)` as R 4.0.4 computes it (`1 / tanpi(0.025)`), %.17g.
pub const QT_975_DF1: f64 = 12.706204736174694;

// ---------------------------------------------------------------------------------------------
// R numeric semantics
// ---------------------------------------------------------------------------------------------

/// R's numeric primitives, bit-exact with x86_64 R 4.0.
pub mod rstats {
    /// An x87 80-bit extended value (64-bit mantissa) — what R's `LDOUBLE` accumulators hold.
    ///
    /// Only what the sums need: exact conversion from `f64`, round-to-nearest-even addition,
    /// halving, and the final rounding to `f64` (the hardware's double rounding, reproduced).
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct Ext80 {
        neg: bool,
        /// value = mant * 2^exp; `mant` has its top bit set unless the value is zero.
        exp: i32,
        mant: u64,
    }

    impl Ext80 {
        pub const ZERO: Ext80 = Ext80 {
            neg: false,
            exp: 0,
            mant: 0,
        };

        /// Exact (every finite double is representable).
        pub fn from_f64(x: f64) -> Ext80 {
            let bits = x.to_bits();
            let neg = bits >> 63 != 0;
            let e = ((bits >> 52) & 0x7ff) as i32;
            let frac = bits & ((1u64 << 52) - 1);
            let (m, exp) = if e == 0 {
                if frac == 0 {
                    return Ext80 {
                        neg,
                        exp: 0,
                        mant: 0,
                    };
                }
                (frac, -1074)
            } else {
                (frac | (1u64 << 52), e - 1075)
            };
            let sh = m.leading_zeros() as i32;
            Ext80 {
                neg,
                exp: exp - sh,
                mant: m << sh,
            }
        }

        /// `(double) s`: round the 64-bit mantissa to 53 bits, nearest-even.
        pub fn to_f64(self) -> f64 {
            if self.mant == 0 {
                return if self.neg { -0.0 } else { 0.0 };
            }
            let mut m = self.mant >> 11;
            let low = self.mant & 0x7ff;
            let mut exp = self.exp + 11;
            if low > 0x400 || (low == 0x400 && (m & 1) == 1) {
                m += 1;
                if m == (1u64 << 53) {
                    m = 1u64 << 52;
                    exp += 1;
                }
            }
            // m in [2^52, 2^53): biased exponent of m * 2^exp is exp + 52 + 1023.
            let biased = exp + 52 + 1023;
            let v = if (1..=2046).contains(&biased) {
                f64::from_bits(((biased as u64) << 52) | (m & ((1u64 << 52) - 1)))
            } else {
                (m as f64) * 2f64.powi(exp)
            };
            if self.neg { -v } else { v }
        }

        /// Exact division by two.
        pub fn half(self) -> Ext80 {
            if self.mant == 0 {
                self
            } else {
                Ext80 {
                    exp: self.exp - 1,
                    ..self
                }
            }
        }
    }

    impl std::ops::Neg for Ext80 {
        type Output = Ext80;
        fn neg(self) -> Ext80 {
            Ext80 {
                neg: !self.neg,
                ..self
            }
        }
    }

    impl std::ops::Sub for Ext80 {
        type Output = Ext80;
        fn sub(self, o: Ext80) -> Ext80 {
            self + (-o)
        }
    }

    impl std::ops::Add for Ext80 {
        type Output = Ext80;
        /// Round-to-nearest-even addition with a 64-bit mantissa.
        fn add(self, o: Ext80) -> Ext80 {
            if self.mant == 0 {
                return o;
            }
            if o.mant == 0 {
                return self;
            }
            let (big, small) = if self.exp > o.exp || (self.exp == o.exp && self.mant >= o.mant) {
                (self, o)
            } else {
                (o, self)
            };
            let shift = (big.exp - small.exp) as u32;
            let big_m: i128 = (big.mant as i128) << 62;
            // The small operand is aligned under the big one; bits that fall off the i128 are
            // summarised by `sticky` (the true value is `small_m + frac`, 0 < frac < 1).
            let (mut small_m, sticky): (i128, bool) = if shift <= 62 {
                ((small.mant as i128) << (62 - shift), false)
            } else {
                let sh = shift - 62;
                if sh >= 64 {
                    (0, true)
                } else {
                    (
                        (small.mant >> sh) as i128,
                        small.mant & ((1u64 << sh) - 1) != 0,
                    )
                }
            };
            let same_sign = big.neg == small.neg;
            if !same_sign && sticky {
                // big - (small_m + frac) = (big - small_m - 1) + (1 - frac): keep a positive
                // remainder so the rounding rule below stays the same.
                small_m += 1;
            }
            let sum: i128 = if same_sign {
                big_m + small_m
            } else {
                big_m - small_m
            };
            if sum == 0 {
                return Ext80::ZERO;
            }
            let mag = sum as u128;
            let bits = 128 - mag.leading_zeros();
            if bits > 64 {
                let drop = bits - 64;
                let mut m = (mag >> drop) as u64;
                let low = mag & ((1u128 << drop) - 1);
                let half = 1u128 << (drop - 1);
                let mut exp = big.exp - 62 + drop as i32;
                if low > half || (low == half && (sticky || (m & 1) == 1)) {
                    m = m.wrapping_add(1);
                    if m == 0 {
                        m = 1u64 << 63;
                        exp += 1;
                    }
                }
                Ext80 {
                    neg: big.neg,
                    exp,
                    mant: m,
                }
            } else {
                let up = 64 - bits;
                Ext80 {
                    neg: big.neg,
                    exp: big.exp - 62 - up as i32,
                    mant: (mag << up) as u64,
                }
            }
        }
    }

    /// R's `sum` of doubles: sequential `long double` accumulation, rounded once at the end.
    /// Non-finite inputs fall back to plain double arithmetic (the result is ±Inf/NaN anyway).
    pub fn sum(xs: &[f64]) -> f64 {
        if xs.iter().any(|x| !x.is_finite()) {
            return xs.iter().sum();
        }
        xs.iter()
            .fold(Ext80::ZERO, |acc, &x| acc + Ext80::from_f64(x))
            .to_f64()
    }

    /// R's `cumsum` of doubles: the running `long double` sum, rounded at every step.
    pub fn cumsum(xs: &[f64]) -> Vec<f64> {
        if xs.iter().any(|x| !x.is_finite()) {
            let mut acc = 0.0;
            return xs
                .iter()
                .map(|&x| {
                    acc += x;
                    acc
                })
                .collect();
        }
        let mut acc = Ext80::ZERO;
        xs.iter()
            .map(|&x| {
                acc = acc + Ext80::from_f64(x);
                acc.to_f64()
            })
            .collect()
    }

    /// R's `mean(c(a, b))`: `long double` sum, divide, then one correction pass.
    pub fn mean2(a: f64, b: f64) -> f64 {
        if !a.is_finite() || !b.is_finite() {
            return (a + b) / 2.0;
        }
        let (ea, eb) = (Ext80::from_f64(a), Ext80::from_f64(b));
        let s = (ea + eb).half();
        let t = (ea - s) + (eb - s);
        (s + t.half()).to_f64()
    }

    /// `sort` for the finite doubles R feeds these helpers (NaN sorts last; R would drop it).
    pub fn sorted(xs: &[f64]) -> Vec<f64> {
        let mut v = xs.to_vec();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        v
    }

    /// `stats::median`; `NA` (NaN) for an empty vector.
    pub fn median(xs: &[f64]) -> f64 {
        let n = xs.len();
        if n == 0 {
            return f64::NAN;
        }
        let v = sorted(xs);
        let half = n.div_ceil(2);
        if n % 2 == 1 {
            v[half - 1]
        } else {
            mean2(v[half - 1], v[half])
        }
    }

    /// `stats::mad` with the default `constant = 1.4826` and `center = median(x)`.
    pub fn mad(xs: &[f64]) -> f64 {
        let center = median(xs);
        let dev: Vec<f64> = xs.iter().map(|x| (x - center).abs()).collect();
        1.4826 * median(&dev)
    }

    /// `stats::quantile(x, p)` type 7, one probability, no NAs.
    pub fn quantile7(xs: &[f64], p: f64) -> f64 {
        let n = xs.len();
        if n == 0 {
            return f64::NAN;
        }
        let v = sorted(xs);
        let index = 1.0 + (n as f64 - 1.0).max(0.0) * p;
        let lo = index.floor();
        let hi = index.ceil();
        let qs = v[lo as usize - 1];
        let xhi = v[hi as usize - 1];
        if index > lo && xhi != qs {
            let h = index - lo;
            (1.0 - h) * qs + h * xhi
        } else {
            qs
        }
    }

    /// `stats::weighted.mean(x, w)`: `sum((x * w)[w != 0]) / sum(w)`.
    pub fn weighted_mean(xs: &[f64], ws: &[f64]) -> f64 {
        let prod: Vec<f64> = xs
            .iter()
            .zip(ws)
            .filter(|&(_, &w)| w != 0.0)
            .map(|(x, w)| x * w)
            .collect();
        sum(&prod) / sum(ws)
    }
}

// ---------------------------------------------------------------------------------------------
// 1. parseAnnotationForest
// ---------------------------------------------------------------------------------------------

/// R's flattened `gateData` for one marker (`X_cutPoints`, `X_numCuts`, `X_cutDepth`,
/// `X_nodePathScore`, `X_nodePopSize`), as `cppGrowAnnotationForest.cpp` assembles it: every
/// gate with `k` cut points contributes `k` consecutive entries, `cut_points` holding the cuts
/// and the other four vectors repeating the gate's value `k` times.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MarkerGateData {
    pub cut_points: Vec<f64>,
    pub num_cuts: Vec<usize>,
    /// 1 = a gate of the root node (R's `cutDepth`).
    pub cut_depth: Vec<usize>,
    pub node_path_score: Vec<f64>,
    pub node_pop_size: Vec<usize>,
}

impl MarkerGateData {
    pub fn len(&self) -> usize {
        self.num_cuts.len()
    }
    pub fn is_empty(&self) -> bool {
        self.num_cuts.is_empty()
    }
    /// The gates as `(first entry index, number of cuts)` groups.
    fn gate_groups(&self) -> Vec<(usize, usize)> {
        let mut groups = Vec::new();
        let mut i = 0;
        while i < self.num_cuts.len() {
            let k = self.num_cuts[i].max(1);
            groups.push((i, k));
            i += k;
        }
        groups
    }
}

/// The C++ `gateData` assembly: gates in search order, each appended to its marker's five
/// vectors, gates with more than `cut_point_upper_bound` cut points dropped. `Gate.cut_points`
/// includes the column min and max (the C++ `gates[0]` and `gates[n-1]`), which are not cuts.
pub fn gate_data_from_gates(
    n_markers: usize,
    gates: &[Gate],
    cut_point_upper_bound: usize,
) -> Vec<MarkerGateData> {
    let mut out = vec![MarkerGateData::default(); n_markers];
    for g in gates {
        let num_cuts = g.cut_points.len().saturating_sub(2);
        if num_cuts > cut_point_upper_bound {
            continue;
        }
        let m = &mut out[g.column];
        for &cp in &g.cut_points[1..g.cut_points.len() - 1] {
            m.cut_points.push(cp);
            m.num_cuts.push(num_cuts);
            m.cut_depth.push(g.depth);
            m.node_path_score.push(g.path_score);
            m.node_pop_size.push(g.pop_size);
        }
    }
    out
}

/// The inverse of [`gate_data_from_gates`] up to the lost min/max (set to ±∞): one [`Gate`]
/// per group of entries, in `gateData` order marker by marker.
pub fn gates_from_gate_data(gate_data: &[MarkerGateData]) -> Vec<Gate> {
    let mut gates = Vec::new();
    for (column, m) in gate_data.iter().enumerate() {
        for (start, k) in m.gate_groups() {
            let mut cut_points = vec![f64::NEG_INFINITY];
            cut_points.extend_from_slice(&m.cut_points[start..start + k]);
            cut_points.push(f64::INFINITY);
            gates.push(Gate {
                column,
                cut_points,
                depth: m.cut_depth[start],
                path_score: m.node_path_score[start],
                pop_size: m.node_pop_size[start],
            });
        }
    }
    gates
}

/// `.parseAnnotationForest`: the depth score of every marker of one sample, and the gates it
/// settled on (the cut-point count with the best score and the node-weighted mean of the cut
/// locations).
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedMarker {
    /// Normalised to the best marker of the sample (1 for that marker). `NaN` reproduces R's
    /// `NA` of the null parse (a forest without a single gate); a marker that never gated in a
    /// non-empty forest scores 0.
    pub channel_score: f64,
    /// `NA` in R when the marker never gated: `None`.
    pub gates: Option<Vec<f64>>,
}

/// `.parseAnnotationForest(gateData, effectiveChannelSizes)` on [`Gate`]s. `markers` only fixes
/// the number of columns; `effective_sizes[j]` is the number of unrestricted cells of marker
/// `j` (`eChannelSize`, the full cell count without channel bounds).
pub fn parse_annotation_forest(
    markers: &[String],
    gates: &[Gate],
    effective_sizes: &[usize],
) -> Vec<ParsedMarker> {
    let gate_data = gate_data_from_gates(markers.len(), gates, CUT_POINT_UPPER_BOUND);
    parse_gate_data(&gate_data, effective_sizes)
}

/// `.parseAnnotationForest` on R's flattened layout.
///
/// R oddities kept: the best cut-point count is the first in order of appearance among the
/// tied counts (the comment says "pick min", the code does not sort); a depth that no node
/// uses gets weight `1/0 = Inf` but is never read; the final normalisation divides by the max
/// score even when it is 0 (all scores become NaN).
pub fn parse_gate_data(
    gate_data: &[MarkerGateData],
    effective_sizes: &[usize],
) -> Vec<ParsedMarker> {
    assert_eq!(gate_data.len(), effective_sizes.len());
    // maxDepth over all channels; -Inf (no gates anywhere) -> null parse.
    let Some(max_depth) = gate_data
        .iter()
        .filter_map(|m| m.cut_depth.iter().copied().max())
        .max()
    else {
        return gate_data
            .iter()
            .map(|_| ParsedMarker {
                channel_score: f64::NAN,
                gates: None,
            })
            .collect();
    };
    // .getEmpiricalDepthCounts: one count per gate (the stride trick of the R code).
    let mut depth_counts = vec![0.0f64; max_depth + 1];
    for m in gate_data {
        let mut uniq: Vec<usize> = m.num_cuts.clone();
        uniq.sort_unstable();
        uniq.dedup();
        for cut_num in uniq {
            for i in stride_lookups(&m.num_cuts, cut_num, 0) {
                let d = m.cut_depth[i];
                if (1..=max_depth).contains(&d) {
                    depth_counts[d] += 1.0;
                }
            }
        }
    }
    let depth_norm: Vec<f64> = depth_counts.iter().map(|c| 1.0 / c).collect();

    let mut parsed: Vec<ParsedMarker> = Vec::with_capacity(gate_data.len());
    for (m, &e_size) in gate_data.iter().zip(effective_sizes) {
        let e_size = e_size as f64;
        let mut uniq_num_cuts: Vec<usize> = Vec::new();
        for &k in &m.num_cuts {
            if !uniq_num_cuts.contains(&k) {
                uniq_num_cuts.push(k);
            }
        }
        if uniq_num_cuts.is_empty() {
            parsed.push(ParsedMarker {
                channel_score: 0.0,
                gates: None,
            });
            continue;
        }
        // .getScoreListForChannel
        let node_score_list: Vec<Vec<f64>> = uniq_num_cuts
            .iter()
            .map(|&cp| {
                stride_lookups(&m.num_cuts, cp, 0)
                    .into_iter()
                    .map(|i| {
                        let pop_weight = m.node_pop_size[i] as f64 / e_size;
                        let depth_weight = depth_norm[m.cut_depth[i]];
                        m.node_path_score[i] * pop_weight * depth_weight
                    })
                    .collect()
            })
            .collect();
        let cut_scores: Vec<f64> = node_score_list.iter().map(|s| rstats::sum(s)).collect();
        let channel_score = rstats::sum(&cut_scores);
        if channel_score > 0.0 {
            let best = cut_scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let sel = cut_scores.iter().position(|&s| s == best).unwrap();
            let sel_num = uniq_num_cuts[sel];
            let weight_sum = rstats::sum(&node_score_list[sel]);
            let node_weights: Vec<f64> = node_score_list[sel]
                .iter()
                .map(|s| s / weight_sum)
                .collect();
            let gates: Vec<f64> = (0..sel_num)
                .map(|cut| {
                    let locs: Vec<f64> = stride_lookups(&m.num_cuts, sel_num, cut)
                        .into_iter()
                        .map(|i| m.cut_points[i])
                        .collect();
                    rstats::weighted_mean(&locs, &node_weights)
                })
                .collect();
            parsed.push(ParsedMarker {
                channel_score,
                gates: Some(gates),
            });
        } else {
            parsed.push(ParsedMarker {
                channel_score: 0.0,
                gates: None,
            });
        }
    }
    let max_score = parsed
        .iter()
        .map(|p| p.channel_score)
        .fold(f64::NEG_INFINITY, f64::max);
    for p in &mut parsed {
        p.channel_score /= max_score;
    }
    parsed
}

/// `which(numCuts == cp)[seq(1 + offset, ., by = cp)]`: every `cp`-th index among the entries
/// with `cp` cuts, starting at the `offset`-th — one index per gate (`offset = 0`) or the
/// `offset`-th cut of every gate.
fn stride_lookups(num_cuts: &[usize], cp: usize, offset: usize) -> Vec<usize> {
    let all: Vec<usize> = (0..num_cuts.len()).filter(|&i| num_cuts[i] == cp).collect();
    all.into_iter().skip(offset).step_by(cp.max(1)).collect()
}

/// `firstDepth` of `.selectChannels`: the shallowest depth at which a marker gated, `Inf` when
/// it never did (`depthMat`).
pub fn first_depths(gate_data: &[MarkerGateData]) -> Vec<f64> {
    gate_data
        .iter()
        .map(|m| {
            m.cut_depth
                .iter()
                .copied()
                .min()
                .map_or(f64::INFINITY, |d| d as f64)
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// 2. selectChannels
// ---------------------------------------------------------------------------------------------

/// `scoreMat` and `depthMat`: experimental units × markers.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoreMat {
    pub exp_units: Vec<String>,
    pub markers: Vec<String>,
    /// `scores[unit][marker]`; R's `NA` (null parse) is already 0.
    pub scores: Vec<Vec<f64>>,
    /// `depthMat`: `first_depth[unit][marker]`, `Inf` when the marker never gated.
    pub first_depth: Vec<Vec<f64>>,
}

impl ScoreMat {
    /// From the parsed forests (`parsed[unit]`) and the per-unit `firstDepth` vectors
    /// ([`first_depths`]).
    pub fn new(
        exp_units: Vec<String>,
        markers: Vec<String>,
        parsed: &[Vec<ParsedMarker>],
        first_depth: Vec<Vec<f64>>,
    ) -> ScoreMat {
        assert_eq!(parsed.len(), exp_units.len());
        assert_eq!(first_depth.len(), exp_units.len());
        let scores = parsed
            .iter()
            .map(|unit| {
                assert_eq!(unit.len(), markers.len());
                unit.iter()
                    .map(|p| {
                        if p.channel_score.is_nan() {
                            0.0
                        } else {
                            p.channel_score
                        }
                    })
                    .collect()
            })
            .collect();
        ScoreMat {
            exp_units,
            markers,
            scores,
            first_depth,
        }
    }

    pub fn column(&self, j: usize) -> Vec<f64> {
        self.scores.iter().map(|row| row[j]).collect()
    }
}

/// `.selectChannels`: the markers whose `quantile(score, selectionQuantile)` across units is
/// `>= depthScoreThreshold`, ordered by their total score, decreasing (ties keep the marker
/// order: R's radix sort is stable). Returns marker indices. `depthMat` is computed by R but
/// does not take part in the selection.
pub fn select_channels(
    score_mat: &ScoreMat,
    depth_score_threshold: f64,
    selection_quantile: f64,
) -> Result<Vec<usize>> {
    let mut selected: Vec<(usize, f64)> = Vec::new();
    for j in 0..score_mat.markers.len() {
        let col = score_mat.column(j);
        let q = rstats::quantile7(&col, selection_quantile);
        if q >= depth_score_threshold {
            selected.push((j, rstats::sum(&col)));
        }
    }
    if selected.is_empty() {
        return Err(FaustError::Invalid(
            "No channels selected at current settings: increase selectionQuantile, decrease depthScoreThreshold".into(),
        ));
    }
    selected.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    Ok(selected.into_iter().map(|(j, _)| j).collect())
}

// ---------------------------------------------------------------------------------------------
// 3. reconcileAnnotationBoundaries
// ---------------------------------------------------------------------------------------------

/// `.makeGateList`: `raw[unit][k]` = the parsed gates of selected marker `selected[k]`
/// (`None` = R's `NA`).
pub fn make_gate_list(
    parsed: &[Vec<ParsedMarker>],
    selected: &[usize],
) -> Vec<Vec<Option<Vec<f64>>>> {
    parsed
        .iter()
        .map(|unit| selected.iter().map(|&j| unit[j].gates.clone()).collect())
        .collect()
}

/// The standardised thresholds of one marker, per experimental unit (sorted ascending).
#[derive(Debug, Clone, PartialEq)]
pub struct MarkerThresholds {
    pub marker: String,
    pub per_unit: Vec<Vec<f64>>,
}

/// `resList`: marker → experimental unit → thresholds. Markers are in `selectedChannels`
/// order (forced markers appended); units in the analysis-map order.
#[derive(Debug, Clone, PartialEq)]
pub struct ResList {
    pub exp_units: Vec<String>,
    pub markers: Vec<MarkerThresholds>,
}

impl ResList {
    pub fn marker_names(&self) -> Vec<String> {
        self.markers.iter().map(|m| m.marker.clone()).collect()
    }
    /// The thresholds of every marker for one unit (`scampAF`).
    pub fn thresholds(&self, unit: usize) -> Vec<&[f64]> {
        self.markers
            .iter()
            .map(|m| m.per_unit[unit].as_slice())
            .collect()
    }
    /// `gateNums`: thresholds + 1 per marker, for one unit.
    pub fn gate_nums(&self, unit: usize) -> Vec<usize> {
        self.markers
            .iter()
            .map(|m| m.per_unit[unit].len() + 1)
            .collect()
    }
}

/// The supervision FAUST accepts (`forceList`, `preferenceList`, the imputation hierarchy).
/// The operator runs the default path: everything empty.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ReconcileOptions {
    /// `impH` per experimental unit; `None` = one level for all units (`allSamples`).
    pub imputation_hierarchy: Option<Vec<String>>,
    /// `preferenceList`: marker → requested number of gates (used when observed).
    pub preference: Vec<(String, usize)>,
    /// `forceList`: marker → gates imposed on every unit (overwrites, or appends the marker).
    pub force: Vec<(String, Vec<f64>)>,
}

/// Output of [`reconcile_boundaries`].
#[derive(Debug, Clone, PartialEq)]
pub struct Reconciled {
    /// `resListPrep` (= `resList` without a selection list).
    pub res_list: ResList,
    /// `selectedChannels` (the input selection plus forced markers).
    pub selected_channels: Vec<String>,
    /// `numSel`: the standard number of gates per selected marker.
    pub num_gates: Vec<usize>,
    /// `possibilityList`: per selected marker, `(number of gates, number of units)` observed.
    pub possibilities: GateCountTable,
}

/// `possibilityList`: per marker, `(number of gates, number of units with that many)`.
pub type GateCountTable = Vec<Vec<(usize, usize)>>;

/// `.getConsistentNumberOfGatesForMarkers`: the most frequent non-zero gate count, the smallest
/// on ties; a preference wins when it was observed at all.
fn consistent_number_of_gates(
    num_gate: &[Vec<usize>],
    selected: &[String],
    preference: &[(String, usize)],
) -> Result<(Vec<usize>, GateCountTable)> {
    let mut num_sel = Vec::with_capacity(selected.len());
    let mut possibilities = Vec::with_capacity(selected.len());
    for (c, name) in selected.iter().enumerate() {
        let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
        for row in num_gate {
            if row[c] > 0 {
                *counts.entry(row[c]).or_insert(0) += 1;
            }
        }
        if counts.is_empty() {
            return Err(FaustError::Invalid(format!(
                "Invalid reconciliation: marker {name} has no gates in any experimental unit"
            )));
        }
        let column_max = *counts.values().max().unwrap();
        let standard = *counts.iter().find(|&(_, &n)| n == column_max).unwrap().0;
        let sel = match preference.iter().find(|(p, _)| p == name) {
            Some((_, want)) if counts.contains_key(want) => *want,
            _ => standard,
        };
        num_sel.push(sel);
        possibilities.push(counts.into_iter().collect());
    }
    Ok((num_sel, possibilities))
}

fn column_medians(rows: &[Vec<f64>], ncol: usize) -> Vec<f64> {
    (0..ncol)
        .map(|i| rstats::median(&rows.iter().map(|r| r[i]).collect::<Vec<_>>()))
        .collect()
}

fn column_mads(rows: &[Vec<f64>], ncol: usize) -> Vec<f64> {
    (0..ncol)
        .map(|i| rstats::mad(&rows.iter().map(|r| r[i]).collect::<Vec<_>>()))
        .collect()
}

/// `lowVal`/`highVal`: medians ∓ `qt(0.975, 1) * mad`; any zero MAD makes every MAD infinite
/// (R: `mgMAD <- Inf`, a scalar).
fn outlier_bounds(med: &[f64], mad: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let any_zero = mad.contains(&0.0);
    let mad = |i: usize| if any_zero { f64::INFINITY } else { mad[i] };
    let low = (0..med.len())
        .map(|i| med[i] - (QT_975_DF1 * mad(i)))
        .collect();
    let high = (0..med.len())
        .map(|i| med[i] + (QT_975_DF1 * mad(i)))
        .collect();
    (low, high)
}

/// `.reconcileAnnotationBoundaries` without the file system: `raw[unit][k]` are the parsed
/// gates of `selected[k]` ([`make_gate_list`]), units in analysis-map order.
///
/// Per marker: the standard gate count; units with that count supply the standard locations
/// (column medians after replacing MAD outliers by the median of the others); units with a
/// different count are converted to the standard by nearest-gate matching
/// ([`up_convert`] / [`down_convert`]); units without gates get the medians — first within
/// each level of the imputation hierarchy, then across the experiment.
///
/// R oddities kept: in the experiment-wide fallback the conversion direction test is inverted
/// (`gateNumber < length(modVals)` calls `.upConvert`, i.e. a unit with *more* gates than the
/// standard goes through the up-conversion, which then keeps the first of its gates nearest to
/// each standard gate), and only that branch clamps converted gates to the outlier bounds.
pub fn reconcile_boundaries(
    exp_units: &[String],
    selected: &[String],
    raw: &[Vec<Option<Vec<f64>>>],
    opts: &ReconcileOptions,
) -> Result<Reconciled> {
    let n_units = exp_units.len();
    assert_eq!(raw.len(), n_units);
    let num_gate: Vec<Vec<usize>> = raw
        .iter()
        .map(|unit| {
            assert_eq!(unit.len(), selected.len());
            unit.iter()
                .map(|g| g.as_ref().map_or(0, Vec::len))
                .collect()
        })
        .collect();
    let (num_sel, possibilities) =
        consistent_number_of_gates(&num_gate, selected, &opts.preference)?;

    let imp_h: Vec<String> = match &opts.imputation_hierarchy {
        Some(h) => {
            assert_eq!(h.len(), n_units);
            h.clone()
        }
        None => vec!["allSamples".to_string(); n_units],
    };
    let mut unique_ih: Vec<&String> = Vec::new();
    for h in &imp_h {
        if !unique_ih.contains(&h) {
            unique_ih.push(h);
        }
    }

    let mut markers: Vec<MarkerThresholds> = Vec::with_capacity(selected.len());
    for (c, name) in selected.iter().enumerate() {
        let gate_number = num_sel[c];
        let mut update: Vec<Option<Vec<f64>>> = vec![None; n_units];
        let matching: Vec<usize> = (0..n_units)
            .filter(|&u| num_gate[u][c] == gate_number)
            .collect();
        for ih in &unique_ih {
            // sort(unique(expUnits in this level)): byte order of the names.
            let mut all_in_ih: Vec<usize> = (0..n_units).filter(|&u| &imp_h[u] == *ih).collect();
            all_in_ih.sort_by(|&a, &b| exp_units[a].cmp(&exp_units[b]));
            let matching_in_ih: Vec<usize> = matching
                .iter()
                .copied()
                .filter(|u| all_in_ih.contains(u))
                .collect();
            if matching_in_ih.is_empty() {
                continue;
            }
            let mut rows: Vec<usize> = matching_in_ih.clone();
            let mut gate_matrix: Vec<Vec<f64>> = rows
                .iter()
                .map(|&u| rstats::sorted(raw[u][c].as_ref().unwrap()))
                .collect();
            let mut possible_gates: Vec<usize> =
                all_in_ih.iter().map(|&u| num_gate[u][c]).collect();
            possible_gates.sort_unstable();
            possible_gates.dedup();
            possible_gates.retain(|&g| g != 0 && g != gate_number);
            for (r, &u) in rows.iter().enumerate() {
                update[u] = Some(rstats::sorted(&gate_matrix[r]));
            }
            // outliers -> medians of the others
            let mg_med = column_medians(&gate_matrix, gate_number);
            let mg_mad = column_mads(&gate_matrix, gate_number);
            let (low, high) = outlier_bounds(&mg_med, &mg_mad);
            let mut amended: Vec<usize> = Vec::new();
            for i in 0..gate_number {
                let bad: Vec<usize> = (0..gate_matrix.len())
                    .filter(|&r| !(gate_matrix[r][i] >= low[i] && gate_matrix[r][i] <= high[i]))
                    .collect();
                if !bad.is_empty() {
                    amended.extend(bad.iter().map(|&r| rows[r]));
                    let good: Vec<f64> = (0..gate_matrix.len())
                        .filter(|r| !bad.contains(r))
                        .map(|r| gate_matrix[r][i])
                        .collect();
                    let new_val = rstats::median(&good);
                    for &r in &bad {
                        gate_matrix[r][i] = new_val;
                    }
                }
            }
            amended.sort_by(|&a, &b| exp_units[a].cmp(&exp_units[b]));
            amended.dedup();
            for &u in &amended {
                let r = rows.iter().position(|&x| x == u).unwrap();
                update[u] = Some(rstats::sorted(&gate_matrix[r]));
            }
            let mut final_vals = column_medians(&gate_matrix, gate_number);
            // units with another number of gates: convert to the standard
            for &gate_num in &possible_gates {
                let mod_units: Vec<usize> = (0..n_units)
                    .filter(|&u| num_gate[u][c] == gate_num && all_in_ih.contains(&u))
                    .collect();
                for u in mod_units {
                    let mod_vals = raw[u][c].as_ref().unwrap();
                    let new_vals = if gate_num < final_vals.len() {
                        up_convert(mod_vals, &final_vals)?
                    } else {
                        down_convert(mod_vals, &final_vals)?
                    };
                    update[u] = Some(rstats::sorted(&new_vals));
                }
            }
            // units without gates: the medians of the (now standard) units of the level
            let na_units: Vec<usize> = (0..n_units)
                .filter(|&u| update[u].is_none() && all_in_ih.contains(&u))
                .collect();
            if !na_units.is_empty() {
                rows = all_in_ih
                    .iter()
                    .copied()
                    .filter(|u| !na_units.contains(u))
                    .collect();
                gate_matrix = rows
                    .iter()
                    .map(|&u| rstats::sorted(update[u].as_ref().unwrap()))
                    .collect();
                final_vals = column_medians(&gate_matrix, gate_number);
                for u in na_units {
                    update[u] = Some(rstats::sorted(&final_vals));
                }
            }
        }
        if update.iter().any(Option::is_none) {
            // standardise across the whole experiment
            let gate_matrix: Vec<Vec<f64>> = (0..n_units)
                .filter_map(|u| update[u].as_ref().map(|g| rstats::sorted(g)))
                .collect();
            let final_vals = column_medians(&gate_matrix, gate_number);
            let mg_mad = column_mads(&gate_matrix, gate_number);
            let (low, high) = outlier_bounds(&final_vals, &mg_mad);
            let na_units: Vec<usize> = (0..n_units).filter(|&u| update[u].is_none()).collect();
            for u in na_units {
                match &raw[u][c] {
                    Some(mod_vals) => {
                        // R: `if (gateNumber < length(modVals)) .upConvert else .downConvert`
                        let new_vals = if gate_number < mod_vals.len() {
                            up_convert(mod_vals, &final_vals)?
                        } else {
                            down_convert(mod_vals, &final_vals)?
                        };
                        let mut new_vals = rstats::sorted(&new_vals);
                        for k in 0..new_vals.len() {
                            if new_vals[k] <= low[k] || new_vals[k] >= high[k] {
                                new_vals[k] = final_vals[k];
                            }
                        }
                        update[u] = Some(rstats::sorted(&new_vals));
                    }
                    None => update[u] = Some(rstats::sorted(&final_vals)),
                }
            }
        }
        markers.push(MarkerThresholds {
            marker: name.clone(),
            per_unit: update.into_iter().map(Option::unwrap).collect(),
        });
    }

    let mut selected_out: Vec<String> = selected.to_vec();
    for (forced, gates) in &opts.force {
        let per_unit = vec![gates.clone(); n_units];
        if let Some(m) = markers.iter_mut().find(|m| &m.marker == forced) {
            m.per_unit = per_unit;
        } else {
            markers.push(MarkerThresholds {
                marker: forced.clone(),
                per_unit,
            });
            selected_out.push(forced.clone());
        }
    }
    Ok(Reconciled {
        res_list: ResList {
            exp_units: exp_units.to_vec(),
            markers,
        },
        selected_channels: selected_out,
        num_gates: num_sel,
        possibilities,
    })
}

fn nearest_index(to: &[f64], gate: f64) -> (usize, f64) {
    let diffs: Vec<f64> = to.iter().map(|t| (t - gate).abs()).collect();
    let min = diffs.iter().copied().fold(f64::INFINITY, f64::min);
    (diffs.iter().position(|&d| d == min).unwrap(), min)
}

fn sorted_unique(xs: &[f64]) -> Vec<f64> {
    let mut uniq: Vec<f64> = Vec::new();
    for &x in xs {
        if !uniq.contains(&x) {
            uniq.push(x);
        }
    }
    rstats::sorted(&uniq)
}

/// `.upConvert(fromGates, toGates)`: each from-gate claims its nearest to-gate (the lower one
/// on ties); when several claim the same, the first from-gate is kept; unclaimed to-gates are
/// copied over. Errors when the result does not have `to.len()` distinct values.
pub fn up_convert(from: &[f64], to: &[f64]) -> Result<Vec<f64>> {
    let drop_inds: Vec<usize> = from.iter().map(|&g| nearest_index(to, g).0).collect();
    let mut uniq_drops = drop_inds.clone();
    uniq_drops.sort_unstable();
    uniq_drops.dedup();
    let mut final_og: Vec<f64> = if uniq_drops.len() != drop_inds.len() {
        uniq_drops
            .iter()
            .map(|&idx| from[drop_inds.iter().position(|&d| d == idx).unwrap()])
            .collect()
    } else {
        from.to_vec()
    };
    final_og.extend(
        (0..to.len())
            .filter(|i| !uniq_drops.contains(i))
            .map(|i| to[i]),
    );
    let final_og = sorted_unique(&final_og);
    if final_og.len() != to.len() {
        return Err(FaustError::Invalid("Upconversion error.".into()));
    }
    Ok(final_og)
}

/// `.downConvert(fromGates, toGates)`: each from-gate is matched to its nearest to-gate (the
/// lower one on ties); per to-gate the closest from-gate is kept (the last on ties); to-gates
/// nobody matched are copied over.
pub fn down_convert(from: &[f64], to: &[f64]) -> Result<Vec<f64>> {
    let matched: Vec<(usize, f64)> = from.iter().map(|&g| nearest_index(to, g)).collect();
    let mut uniq: Vec<usize> = matched.iter().map(|m| m.0).collect();
    uniq.sort_unstable();
    uniq.dedup();
    let mut final_og: Vec<f64> = Vec::new();
    for &idx in &uniq {
        let sub: Vec<usize> = (0..from.len()).filter(|&k| matched[k].0 == idx).collect();
        let min_diff = sub
            .iter()
            .map(|&k| matched[k].1)
            .fold(f64::INFINITY, f64::min);
        let keep = *sub.iter().rfind(|&&k| matched[k].1 == min_diff).unwrap();
        final_og.push(from[keep]);
    }
    if uniq.len() != to.len() {
        final_og.extend((0..to.len()).filter(|i| !uniq.contains(i)).map(|i| to[i]));
    }
    let final_og = sorted_unique(&final_og);
    if final_og.len() != to.len() {
        return Err(FaustError::Invalid("Downconversion error.".into()));
    }
    Ok(final_og)
}

// ---------------------------------------------------------------------------------------------
// 4. gating: annotation matrix, gateSample, consensus, phenotype selection, count matrix
// ---------------------------------------------------------------------------------------------

/// `.mkAnnMats`: the annotation level of every cell for every selected marker, column-major:
/// `1 + #{thresholds t : x >= t}`. `exprs[k]` is the expression column of selected marker
/// `k` **as FAUST sees it** (the GatingSet stores float32: round the input with `as f32`).
pub fn annotation_matrix(exprs: &[&[f64]], thresholds: &[&[f64]]) -> Vec<Vec<u8>> {
    assert_eq!(exprs.len(), thresholds.len());
    exprs
        .iter()
        .zip(thresholds)
        .map(|(col, gates)| {
            col.iter()
                .map(|&x| 1 + gates.iter().filter(|&&g| x >= g).count() as u8)
                .collect()
        })
        .collect()
}

/// The internal phenotype label of one cell: `marker~level~gateNum~` per marker, concatenated
/// (`gateSample.cpp` and scamp's `assignNumLabelToScoreVector`).
pub fn phenotype_label(markers: &[String], levels: &[u8], gate_nums: &[usize]) -> String {
    let mut s = String::new();
    for ((m, &l), &n) in markers.iter().zip(levels).zip(gate_nums) {
        s.push_str(m);
        s.push('~');
        s.push_str(&l.to_string());
        s.push('~');
        s.push_str(&n.to_string());
        s.push('~');
    }
    s
}

/// `gateSample`: the label of every cell, or [`NO_PHENOTYPE`] when it is not in `allowed`
/// (the selected phenotypes for `faustAnnotation`, every observed one for
/// `exhaustiveFaustAnnotation`). `ann_cols` is the output of [`annotation_matrix`].
pub fn gate_sample(
    markers: &[String],
    ann_cols: &[Vec<u8>],
    gate_nums: &[usize],
    allowed: &HashSet<String>,
) -> Vec<String> {
    let n_cells = ann_cols.first().map_or(0, Vec::len);
    let mut levels = vec![0u8; markers.len()];
    (0..n_cells)
        .map(|r| {
            for (k, col) in ann_cols.iter().enumerate() {
                levels[k] = col[r];
            }
            let label = phenotype_label(markers, &levels, gate_nums);
            if allowed.contains(&label) {
                label
            } else {
                NO_PHENOTYPE.to_string()
            }
        })
        .collect()
}

/// The scamp consensus of `.clusterExpUnitWithScamp`: a cell keeps its label when the
/// `runClustering` and `maxClustering` labels agree, else [`UNCERTAIN`]; an error when no cell
/// agrees (R kills the job).
pub fn consensus_labels(run: &[String], max: &[String]) -> Result<Vec<String>> {
    assert_eq!(run.len(), max.len());
    let out: Vec<String> = run
        .iter()
        .zip(max)
        .map(|(r, m)| {
            if r == m {
                m.clone()
            } else {
                UNCERTAIN.to_string()
            }
        })
        .collect();
    if !out.iter().any(|l| l != UNCERTAIN) {
        return Err(FaustError::Invalid(
            "SCAMP reached no consensus: increase numScampIter until consensus is reached".into(),
        ));
    }
    Ok(out)
}

/// `scampClusterLabels` of one unit: the sorted distinct labels without [`UNCERTAIN`].
pub fn cluster_names(labels: &[String]) -> Vec<String> {
    let mut names: Vec<String> = labels
        .iter()
        .filter(|l| l.as_str() != UNCERTAIN)
        .cloned()
        .collect();
    names.sort();
    names.dedup();
    names
}

/// `fdrtool:::pvt.isoMean` (`C_isomean`): weighted isotonic regression by pool-adjacent
/// violators, pooling on `>=`.
pub fn iso_mean(y: &[f64], w: &[f64]) -> Vec<f64> {
    let nn = y.len();
    if nn <= 1 {
        return y.to_vec();
    }
    let mut k = vec![0usize; nn];
    let mut gew = vec![0.0f64; nn];
    let mut ghat = vec![0.0f64; nn];
    let mut c = 0usize;
    gew[0] = w[0];
    ghat[0] = y[0];
    for j in 1..nn {
        c += 1;
        k[c] = j;
        gew[c] = w[j];
        ghat[c] = y[j];
        while ghat[c - 1] >= ghat[c] {
            let neu = gew[c] + gew[c - 1];
            ghat[c - 1] += (gew[c] / neu) * (ghat[c] - ghat[c - 1]);
            gew[c - 1] = neu;
            c -= 1;
            if c == 0 {
                break;
            }
        }
    }
    // expand the pooled blocks (block c covers positions k[c]..nn)
    let mut n = nn;
    loop {
        let v = ghat[c];
        for item in ghat.iter_mut().take(n).skip(k[c]) {
            *item = v;
        }
        n = k[c];
        if c == 0 {
            break;
        }
        c -= 1;
    }
    ghat
}

/// `fdrtool::gcmlcm(x, y, type = "gcm")`: the knots of the greatest convex minorant.
/// Returns `(x.knots, y.knots)`.
pub fn gcm(x: &[f64], y: &[f64]) -> (Vec<f64>, Vec<f64>) {
    assert_eq!(x.len(), y.len());
    let n = x.len();
    if n < 2 {
        return (x.to_vec(), y.to_vec());
    }
    let dx: Vec<f64> = (1..n).map(|i| x[i] - x[i - 1]).collect();
    let dy: Vec<f64> = (1..n).map(|i| y[i] - y[i - 1]).collect();
    let raw: Vec<f64> = dy
        .iter()
        .zip(&dx)
        .map(|(dy, dx)| {
            let s = dy / dx;
            if s == f64::INFINITY {
                f64::MAX
            } else if s == f64::NEG_INFINITY {
                -f64::MAX
            } else {
                s
            }
        })
        .collect();
    let slope = iso_mean(&raw, &dx);
    // keep = !duplicated(slope): the first occurrence of every distinct value, anywhere.
    let mut seen: Vec<f64> = Vec::new();
    let keep: Vec<bool> = slope
        .iter()
        .map(|&s| {
            if seen.contains(&s) {
                false
            } else {
                seen.push(s);
                true
            }
        })
        .collect();
    let mut x_knots: Vec<f64> = (0..n - 1).filter(|&i| keep[i]).map(|i| x[i]).collect();
    x_knots.push(x[n - 1]);
    let dx_knots: Vec<f64> = (1..x_knots.len())
        .map(|i| x_knots[i] - x_knots[i - 1])
        .collect();
    let slope_knots: Vec<f64> = (0..n - 1).filter(|&i| keep[i]).map(|i| slope[i]).collect();
    let prods: Vec<f64> = dx_knots
        .iter()
        .zip(&slope_knots)
        .map(|(d, s)| d * s)
        .collect();
    let mut y_knots = vec![y[0]];
    y_knots.extend(rstats::cumsum(&prods).into_iter().map(|c| y[0] + c));
    (x_knots, y_knots)
}

/// `.computeElbow(phenoDF)`: `x` the occurrence counts `1..max`, `y` the number of phenotypes
/// seen at least that often. Fewer than 4 rows keep every phenotype (elbow 1); otherwise the
/// rows between `max(2, ceil(0.05 n))` and `min(n - 1, floor(0.95 n))` are kept, and the elbow
/// is the GCM knot farthest from the chord through the first and last kept rows (first on
/// ties).
pub fn compute_elbow(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len());
    let n = x.len();
    if n < 4 {
        return 1.0;
    }
    let lower = 2usize.max((0.05 * n as f64).ceil() as usize);
    let upper = (n - 1).min((0.95 * n as f64).floor() as usize);
    let ex = &x[lower - 1..upper];
    let ey = &y[lower - 1..upper];
    let (x1, y1) = (ex[0], ey[0]);
    let (x2, y2) = (ex[ex.len() - 1], ey[ey.len() - 1]);
    let m = (y2 - y1) / (x2 - x1);
    let c = y1 - (m * x1);
    let line_min = |px: f64, py: f64| (py - m * px - c).abs() / (1.0 + m * m).sqrt();
    let (xk, yk) = gcm(ex, ey);
    let cand: Vec<f64> = xk
        .iter()
        .zip(&yk)
        .map(|(&px, &py)| line_min(px, py))
        .collect();
    let best = cand.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    xk[cand.iter().position(|&d| d == best).unwrap()]
}

/// Output of [`phenotype_selection`].
#[derive(Debug, Clone, PartialEq)]
pub struct PhenotypeSelection {
    /// `scampNameSummary`: every phenotype (sorted) and the number of units it was found in.
    pub name_summary: Vec<(String, usize)>,
    /// `phenotypeElbowValue`: the automatic occurrence threshold (used when
    /// `nameOccuranceNum == 0`).
    pub elbow: f64,
    /// `scampClusterNames`: the phenotypes kept for gating.
    pub selected: Vec<String>,
    /// `allScampClusterNames`: every phenotype (for the exhaustive gating).
    pub all: Vec<String>,
}

/// The phenotype filter of `.clusterExpUnitsWithScamp`: `per_unit_names[u]` are the
/// [`cluster_names`] of unit `u`; a phenotype is kept when it occurs in at least
/// `name_occurrence_num` units (`> 0`), or at least the elbow of the occurrence curve.
pub fn phenotype_selection(
    per_unit_names: &[Vec<String>],
    name_occurrence_num: usize,
) -> Result<PhenotypeSelection> {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for unit in per_unit_names {
        for name in unit {
            *counts.entry(name.clone()).or_insert(0) += 1;
        }
    }
    let name_summary: Vec<(String, usize)> = counts.into_iter().collect();
    let max_count = name_summary.iter().map(|(_, n)| *n).max().ok_or_else(|| {
        FaustError::Invalid("no scamp phenotypes found in any experimental unit".into())
    })?;
    let x: Vec<f64> = (1..=max_count).map(|v| v as f64).collect();
    let y: Vec<f64> = (1..=max_count)
        .map(|v| name_summary.iter().filter(|(_, n)| *n >= v).count() as f64)
        .collect();
    let elbow = compute_elbow(&x, &y);
    let threshold = if name_occurrence_num > 0 {
        name_occurrence_num as f64
    } else {
        elbow
    };
    let selected = name_summary
        .iter()
        .filter(|(_, n)| *n as f64 >= threshold)
        .map(|(s, _)| s.clone())
        .collect();
    let all = name_summary.iter().map(|(s, _)| s.clone()).collect();
    Ok(PhenotypeSelection {
        name_summary,
        elbow,
        selected,
        all,
    })
}

/// `faustCountMatrix` / `exhaustiveFaustCountMatrix`: samples × phenotypes (+ [`NO_PHENOTYPE`]).
#[derive(Debug, Clone, PartialEq)]
pub struct CountMatrix {
    pub samples: Vec<String>,
    /// Internal labels (`faustColNames`), [`NO_PHENOTYPE`] last.
    pub phenotypes: Vec<String>,
    /// `newColNames`: the readable labels ([`pretty_count_name`]).
    pub pretty: Vec<String>,
    /// `counts[sample][phenotype]`.
    pub counts: Vec<Vec<u64>>,
}

/// `.getFaustCountMatrix`: count the gated labels ([`gate_sample`] output) of every sample over
/// `phenotypes` plus [`NO_PHENOTYPE`].
pub fn count_matrix(
    samples: &[String],
    labels: &[Vec<String>],
    phenotypes: &[String],
) -> CountMatrix {
    assert_eq!(samples.len(), labels.len());
    let mut cols: Vec<String> = phenotypes.to_vec();
    cols.push(NO_PHENOTYPE.to_string());
    let index: HashMap<&str, usize> = cols
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let counts = labels
        .iter()
        .map(|sample| {
            let mut row = vec![0u64; cols.len()];
            for l in sample {
                if let Some(&i) = index.get(l.as_str()) {
                    row[i] += 1;
                }
            }
            row
        })
        .collect();
    CountMatrix {
        samples: samples.to_vec(),
        pretty: cols.iter().map(|s| pretty_count_name(s)).collect(),
        phenotypes: cols,
        counts,
    }
}

/// Split `marker~level~levels~…` into `(marker, level, levels)` units; `None` when the string
/// is not of that shape (e.g. [`NO_PHENOTYPE`]).
fn label_units(label: &str) -> Option<Vec<(&str, u32, u32)>> {
    let parts: Vec<&str> = label.split('~').collect();
    if parts.len() < 4 || parts.len() % 3 != 1 || !parts[parts.len() - 1].is_empty() {
        return None;
    }
    let mut units = Vec::with_capacity(parts.len() / 3);
    for chunk in parts[..parts.len() - 1].chunks(3) {
        let level = chunk[1].parse().ok()?;
        let levels = chunk[2].parse().ok()?;
        units.push((chunk[0], level, levels));
    }
    Some(units)
}

fn map_label(label: &str, map: impl Fn(u32, u32) -> Option<&'static str>) -> String {
    let Some(units) = label_units(label) else {
        return label.to_string();
    };
    let mut out = String::with_capacity(label.len());
    for (marker, level, levels) in units {
        out.push_str(marker);
        match map(level, levels) {
            Some(s) => out.push_str(s),
            None => out.push_str(&format!("~{level}~{levels}~")),
        }
    }
    out
}

/// The `gsub` chain of `.getFaustCountMatrix` (`newColNames`): `~1~d~` → `-`, `~2~2~` → `+`,
/// `~3~3~` → `Bright`, `~4~4~` → `++`, `~2~3~` → `Dim`, `~2~4~` → `Med-`, `~3~4~` → `Med+`;
/// anything else is left as is.
pub fn pretty_count_name(label: &str) -> String {
    map_label(label, |level, levels| match (level, levels) {
        (1, _) => Some("-"),
        (2, 2) => Some("+"),
        (3, 3) => Some("Bright"),
        (4, 4) => Some("++"),
        (2, 3) => Some("Dim"),
        (2, 4) => Some("Med-"),
        (3, 4) => Some("Med+"),
        _ => None,
    })
}

/// The `gsub` chain of `.prepareSamplesForEmbed` (`faustLabels` of the annotation embedding,
/// the operator's `label` column): `-`/`+` for two levels, `-`/`Dim`/`Bright` for three,
/// `-`/`MedLow`/`MedHigh`/`VeryBright` for four; anything else is left as is.
pub fn embedding_label(label: &str) -> String {
    map_label(label, |level, levels| match (level, levels) {
        (1, 2) | (1, 3) | (1, 4) => Some("-"),
        (2, 2) => Some("+"),
        (2, 3) => Some("Dim"),
        (3, 3) => Some("Bright"),
        (2, 4) => Some("MedLow"),
        (3, 4) => Some("MedHigh"),
        (4, 4) => Some("VeryBright"),
        _ => None,
    })
}

/// The per-marker annotation strings of the embedding (`*_faust_annotation`, the operator's
/// `*_gate` columns): the level of a cell for a marker with `levels` annotation levels.
pub fn level_label(level: u8, levels: usize) -> String {
    match (levels, level) {
        (2, 1) | (3, 1) | (4, 1) => "-".into(),
        (2, 2) => "+".into(),
        (3, 2) => "Dim".into(),
        (3, 3) => "Bright".into(),
        (4, 2) => "MedLow".into(),
        (4, 3) => "MedHigh".into(),
        (4, 4) => "VeryBright".into(),
        _ => level.to_string(),
    }
}

// ---------------------------------------------------------------------------------------------
// 0. channel bounds and restriction matrices (the input of the forest, scamp and the scores)
// ---------------------------------------------------------------------------------------------

/// `channelBoundsUsedByFAUST`: the `Low`/`High` restriction bound of every marker.
///
/// The operator leaves `channelBounds = ""`, so FAUST derives them empirically
/// ([`empirical_channel_bounds`]). Cells outside the bounds are *restricted*: the annotation
/// forest and scamp drop them from the one-dimensional density estimates of that marker
/// (`anyValueRestricted = TRUE`, `resValMatrix`), and the depth score divides node sizes by the
/// number of unrestricted cells (`effectiveChannelSizes`, [`effective_sizes`]).
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelBounds {
    pub low: Vec<f64>,
    pub high: Vec<f64>,
}

/// First stage of `.processChannelBounds`: the 1% and 99% quantiles (type 7) of every marker
/// of one sample, on the expression values as FAUST holds them (float32-rounded).
pub fn sample_quantile_bounds(exprs: &[&[f64]]) -> (Vec<f64>, Vec<f64>) {
    let q01 = exprs.iter().map(|c| rstats::quantile7(c, 0.01)).collect();
    let q99 = exprs.iter().map(|c| rstats::quantile7(c, 0.99)).collect();
    (q01, q99)
}

/// Second stage: across the samples, `Low` is the 5% quantile of the per-sample 1% quantiles
/// and `High` the 95% quantile of the per-sample 99% quantiles.
pub fn empirical_channel_bounds_from_quantiles(
    per_sample: &[(Vec<f64>, Vec<f64>)],
) -> Result<ChannelBounds> {
    let Some(first) = per_sample.first() else {
        return Err(FaustError::Invalid(
            "channel bounds: no sample with more than 25 cells".into(),
        ));
    };
    let n_markers = first.0.len();
    let low = (0..n_markers)
        .map(|j| rstats::quantile7(&per_sample.iter().map(|s| s.0[j]).collect::<Vec<_>>(), 0.05))
        .collect();
    let high = (0..n_markers)
        .map(|j| rstats::quantile7(&per_sample.iter().map(|s| s.1[j]).collect::<Vec<_>>(), 0.95))
        .collect();
    Ok(ChannelBounds { low, high })
}

/// `.processChannelBounds` with `channelBounds = ""`: `samples[s][j]` is the expression column
/// of marker `j` in sample `s`; samples with 25 cells or fewer are skipped.
pub fn empirical_channel_bounds(samples: &[Vec<&[f64]>]) -> Result<ChannelBounds> {
    let per_sample: Vec<(Vec<f64>, Vec<f64>)> = samples
        .iter()
        .filter(|s| s.first().is_some_and(|c| c.len() > 25))
        .map(|s| sample_quantile_bounds(s))
        .collect();
    empirical_channel_bounds_from_quantiles(&per_sample)
}

/// `.makeRestrictionMatrices` (`resMat`, column-major): 0 for a free cell, 1 where
/// `x <= Low`, 2 where `x >= High` (a cell at both bounds is 2).
pub fn restriction_matrix(exprs: &[&[f64]], bounds: &ChannelBounds) -> Vec<Vec<u8>> {
    exprs
        .iter()
        .enumerate()
        .map(|(j, col)| {
            col.iter()
                .map(|&x| {
                    if x >= bounds.high[j] {
                        2
                    } else if x <= bounds.low[j] {
                        1
                    } else {
                        0
                    }
                })
                .collect()
        })
        .collect()
}

/// `ePop` of `.growForestForExpUnit`: the unrestricted cells of every marker, the
/// `effective_sizes` of [`parse_annotation_forest`].
pub fn effective_sizes(res: &[Vec<u8>]) -> Vec<usize> {
    res.iter()
        .map(|col| col.iter().filter(|&&r| r == 0).count())
        .collect()
}

/// `resFlag` of the forest / scamp calls: whether any cell of any marker is restricted.
pub fn any_restricted(res: &[Vec<u8>]) -> bool {
    res.iter().any(|col| col.iter().any(|&r| r > 0))
}
