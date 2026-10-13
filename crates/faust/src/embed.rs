//! The annotation embedding's input matrix (`R/makeAnnotationEmbedding.R`,
//! `.prepareSamplesForEmbed` / `.prepareSampleForEmbed`, `vizType = "uniform"`).
//!
//! Every cell's expression on the selected markers is winsorised at the 1st/99th percentile
//! (ties at the bounds jittered), standardised **inside its phenotype** (the cells sharing the
//! same annotation vector), then scaled so each phenotype's extent fits between the midpoints
//! of its neighbouring annotation levels, and shifted to its level's landmark (`level * 10` on
//! the "uniform" layout). UMAP of that matrix (n_neighbors 15, min_dist 0.2) is the embedding.
//!
//! The R code draws unseeded `rnorm` jitters, so its coordinates are not reproducible run to
//! run; this port draws them from a seeded generator, so they are. Coordinates therefore differ
//! from the R operator's by design; the layout logic is the same.
use crate::scamp::r_quantile7;

/// A tiny deterministic normal generator for the jitters (xorshift + Box–Muller). Not the
/// reference's RNG: see the module doc.
pub struct Jitter {
    state: u64,
}

impl Jitter {
    pub fn new(seed: u64) -> Self {
        Self {
            state: seed.max(1) ^ 0x9E37_79B9_7F4A_7C15,
        }
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn normal(&mut self, sd: f64) -> f64 {
        let u1 = self.uniform().max(f64::MIN_POSITIVE);
        let u2 = self.uniform();
        sd * (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

fn quantile7_unsorted(v: &[f64], p: f64) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    r_quantile7(&s, p)
}

/// `.windsorizeExpressionForEmbed`.
fn winsorize_for_embed(v: &[f64], rng: &mut Jitter) -> Vec<f64> {
    let lo = quantile7_unsorted(v, 0.01);
    let hi = quantile7_unsorted(v, 0.99);
    v.iter()
        .map(|&x| {
            if x <= lo {
                lo + rng.normal(0.01)
            } else if x >= hi {
                hi + rng.normal(0.01)
            } else {
                x
            }
        })
        .collect()
}

/// `.windsorizeAndScaleExpression`: winsorised and mapped to [0, 1] (the `_Windsorized`
/// columns of R's output, kept for completeness).
pub fn winsorize_and_scale(v: &[f64]) -> Vec<f64> {
    let lo = quantile7_unsorted(v, 0.01);
    let hi = quantile7_unsorted(v, 0.99);
    let w: Vec<f64> = v.iter().map(|&x| x.clamp(lo, hi)).collect();
    let mn = w.iter().copied().fold(f64::INFINITY, f64::min);
    let mx = w.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    w.iter().map(|x| (x - mn) / (mx - mn)).collect()
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

fn sd(v: &[f64]) -> f64 {
    let m = mean(v);
    (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (v.len() as f64 - 1.0)).sqrt()
}

/// One sample's input to the embedding.
pub struct EmbedSample<'a> {
    /// Expression of the **selected** markers, column-major (`cols[j][cell]`).
    pub exprs: &'a [Vec<f64>],
    /// Annotation level (1-based) of every cell on every selected marker: `ann[cell][j]`.
    pub ann: &'a [Vec<u8>],
}

/// The prepared matrix for every cell of every sample (samples concatenated), cell-major:
/// `rows[cell][j]` over the selected markers, plus which markers UMAP should use.
pub struct Prepared {
    pub rows: Vec<Vec<f64>>,
    /// Indices (into the selected markers) kept for the embedding: those whose
    /// 0.75-quantile depth score across samples is ≥ 0.01.
    pub embedding_markers: Vec<usize>,
}

/// `.prepareSamplesForEmbed(vizType = "uniform")` + the marker filter of
/// `makeAnnotationEmbedding`.
///
/// `levels[j]` is the number of annotation levels of selected marker `j` (thresholds + 1);
/// `score_mat[s][j]` the depth score of marker `j` in sample `s`.
pub fn prepare(
    samples: &[EmbedSample],
    levels: &[usize],
    score_mat: &[Vec<f64>],
    seed: u64,
) -> Prepared {
    let p = levels.len();
    let n_samples = score_mat.len().max(1) as f64;
    let mut rng = Jitter::new(seed);
    // annotLevelList: level k sits at k * 10 (uniform layout).
    let centers: Vec<Vec<f64>> = levels
        .iter()
        .map(|&k| (1..=k).map(|l| l as f64 * 10.0).collect())
        .collect();
    // targetThresholds: midpoints between consecutive centres.
    let targets: Vec<Vec<f64>> = centers
        .iter()
        .map(|c| c.windows(2).map(|w| (w[0] + w[1]) / 2.0).collect())
        .collect();
    // quantile weights from the raw score sums
    let raw_sums: Vec<f64> = (0..p)
        .map(|j| score_mat.iter().map(|r| r[j]).sum())
        .collect();
    let total: f64 = raw_sums.iter().sum();
    let mut q_low = vec![0.0; p];
    let mut acc = 0.0;
    for j in 0..p {
        // quantileScoreWeights = c(0, w[-last]) then cumsum * 0.025
        let w = if j == 0 { 0.0 } else { raw_sums[j - 1] / total };
        acc += w;
        q_low[j] = acc * 0.025;
    }
    let q_high: Vec<f64> = q_low.iter().map(|q| 1.0 - q).collect();
    let _ = n_samples;

    let mut rows: Vec<Vec<f64>> = Vec::new();
    for s in samples {
        let n = s.ann.len();
        // winsorised expression per marker (column-major)
        let wins: Vec<Vec<f64>> = s
            .exprs
            .iter()
            .map(|c| winsorize_for_embed(c, &mut rng))
            .collect();
        // phenotype key per cell
        let keys: Vec<Vec<u8>> = s.ann.to_vec();
        let mut uniq: Vec<Vec<u8>> = Vec::new();
        for k in &keys {
            if !uniq.contains(k) {
                uniq.push(k.clone());
            }
        }
        // standardise inside each phenotype; record per (marker, level) the cluster mins/maxs
        let mut scaled = vec![vec![0.0f64; p]; n];
        let mut rec_low: Vec<Vec<Vec<f64>>> = levels.iter().map(|&k| vec![Vec::new(); k]).collect();
        let mut rec_high: Vec<Vec<Vec<f64>>> =
            levels.iter().map(|&k| vec![Vec::new(); k]).collect();
        for u in &uniq {
            let idx: Vec<usize> = (0..n).filter(|&i| keys[i] == *u).collect();
            for j in 0..p {
                let col: Vec<f64> = idx.iter().map(|&i| wins[j][i]).collect();
                let z: Vec<f64> = if col.len() == 1 {
                    vec![rng.normal(1e-7)]
                } else {
                    let m = mean(&col);
                    let s = sd(&col);
                    col.iter()
                        .map(|x| if s == 0.0 { x - m } else { (x - m) / s })
                        .collect()
                };
                let mn = z.iter().copied().fold(f64::INFINITY, f64::min);
                let mx = z.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let lvl = (u[j] as usize).clamp(1, levels[j]) - 1;
                rec_low[j][lvl].push(mn);
                rec_high[j][lvl].push(mx);
                for (k, &i) in idx.iter().enumerate() {
                    scaled[i][j] = z[k];
                }
            }
        }
        // scaling factors per (marker, level)
        let mut factors: Vec<Vec<f64>> = levels.iter().map(|&k| vec![1.0; k]).collect();
        for j in 0..p {
            let k = levels[j];
            for lvl in 0..k {
                let min_v = if rec_low[j][lvl].is_empty() {
                    f64::NAN
                } else {
                    quantile7_unsorted(&rec_low[j][lvl], q_low[j])
                };
                let max_v = if rec_high[j][lvl].is_empty() {
                    f64::NAN
                } else {
                    quantile7_unsorted(&rec_high[j][lvl], q_high[j])
                };
                let translate = centers[j][lvl];
                let f = if k == 1 {
                    1.0
                } else if lvl == 0 {
                    (targets[j][0] - translate) / max_v
                } else if lvl == k - 1 {
                    (targets[j][k - 2] - translate) / min_v
                } else {
                    let hi = (targets[j][lvl] - translate) / max_v;
                    let lo = (targets[j][lvl - 1] - translate) / min_v;
                    hi.min(lo)
                };
                factors[j][lvl] = if f.is_finite() { f } else { 1.0 };
            }
        }
        // apply: scale, jitter, translate
        for i in 0..n {
            let mut row = vec![0.0; p];
            for j in 0..p {
                let lvl = (keys[i][j] as usize).clamp(1, levels[j]) - 1;
                row[j] = factors[j][lvl] * scaled[i][j] + rng.normal(1e-4) + centers[j][lvl];
            }
            rows.push(row);
        }
    }
    // embedding markers: quantile 0.75 of the depth score across samples ≥ 0.01
    let embedding_markers: Vec<usize> = (0..p)
        .filter(|&j| {
            let col: Vec<f64> = score_mat.iter().map(|r| r[j]).collect();
            if col.is_empty() {
                return true;
            }
            quantile7_unsorted(&col, 0.75) >= 0.01
        })
        .collect();
    Prepared {
        rows,
        embedding_markers,
    }
}
