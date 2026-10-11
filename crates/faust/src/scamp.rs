//! scamp: selective clustering annotated using modes of projections (RGLab/scamp 0.5.0), as
//! FAUST's `.clusterExpUnitWithScamp` calls it: `numberIterations = 1`, `useAnnForest = TRUE`
//! with the reconciled thresholds as annotation boundaries, random candidate search,
//! `minimumClusterSize = 25`, `maximumClusterNum = 50 * ncol`, no restricted values, no time
//! budget, single thread.
//!
//! Port of the R `scamp()` prologue (`scale()` and the threshold transfer into the scaled
//! space), `cppNoisyScamp.cpp` (one iteration, so the vote is trivial), `scamp.cpp`
//! (noise → candidate search → selection → residual search loop → annotation),
//! `selectCandidateClusters.cpp` (+ `absSampleLmoments.cpp`), `annotateCluster.cpp`,
//! `rQuantile.cpp`.
//!
//! The per-cell result is the label string `cppNoisyScamp` returns with forest values:
//! `"<marker>~<score+1>~<boundaries+1>~"` per marker, `"Uncertain"` when empty. FAUST's R
//! turns that into `CD3+CD4-…` with the thresholds ([`crate::select`]).
use crate::dip::single_dip;
use crate::forest::{SearchParams, search};
use crate::noise::add_noise_to_matrix;
use crate::select::rstats;

/// R's `scale(x)`: centre by `colMeans`, divide by `sqrt(sum(v^2) / (n - 1))` of the centred
/// column, both sums accumulated in long double as R does.
pub fn r_scale(cols: &[Vec<f64>]) -> Vec<Vec<f64>> {
    cols.iter()
        .map(|c| {
            let n = c.len();
            if n == 0 {
                return Vec::new();
            }
            // R's colMeans and sum() accumulate in 80-bit long double: rstats::sum emulates it.
            let mean = rstats::sum(c) / n as f64;
            let centred: Vec<f64> = c.iter().map(|v| v - mean).collect();
            let sq: Vec<f64> = centred.iter().map(|v| v * v).collect();
            let ss = rstats::sum(&sq);
            let sd = (ss / (n.max(2) - 1) as f64).sqrt();
            centred.iter().map(|v| v / sd).collect()
        })
        .collect()
}

/// R `quantile(x, probs, type = 7)` for one probability (as `scamp()` calls it on a scaled
/// column): `index = 1 + (n-1) p`, interpolate between the order statistics when `index` is
/// not integral and the two differ.
pub fn r_quantile7(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    let index = 1.0 + ((n - 1) as f64).max(0.0) * p;
    let lo = index.floor();
    let hi = index.ceil();
    let qs = sorted[lo as usize - 1];
    if index > lo && sorted[hi as usize - 1] != qs {
        let h = index - lo;
        (1.0 - h) * qs + h * sorted[hi as usize - 1]
    } else {
        qs
    }
}

/// The annotation thresholds moved into the scaled space, as `scamp()` does: for each
/// threshold `val` of marker `i`, `probVal = #(x <= val) / n` on the raw column and the
/// threshold becomes `quantile(scaled_i, probVal, type = 7)`.
pub fn thresholds_to_scaled(
    raw: &[Vec<f64>],
    scaled: &[Vec<f64>],
    thresholds: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    raw.iter()
        .zip(scaled)
        .zip(thresholds)
        .map(|((r, s), th)| {
            let n = r.len() as f64;
            let mut sorted = s.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            th.iter()
                .map(|&val| {
                    let prob = r.iter().filter(|&&x| x <= val).count() as f64 / n;
                    rstats::quantile7(&sorted, prob)
                })
                .collect()
        })
        .collect()
}

/// `rQuantile(dataVector, probs)` from scamp (sorts, then the type-7 rule with a 0-based
/// index `(n - 1) p`).
pub fn r_quantile_cpp(data: &[f64], probs: &[f64]) -> Vec<f64> {
    let mut d = data.to_vec();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = d.len() as f64;
    probs
        .iter()
        .map(|&p| {
            let index = (n - 1.0) * p;
            let lo = index.floor();
            let hi = index.ceil();
            let mut q = d[lo as usize];
            if index > lo && d[hi as usize] != q {
                let h = index - lo;
                q = (1.0 - h) * q + h * d[hi as usize];
            }
            q
        })
        .collect()
}

/// `absSampleLmoments(sortedDataVector, numMoments, trim)`: trimmed sample L-moments, the
/// ratios from the third on, absolute values.
pub fn abs_sample_lmoments(sorted: &[f64], num_moments: usize, trim: usize) -> Vec<f64> {
    let num_obs = sorted.len();
    let n = num_obs as f64;
    let mut maxmom = num_moments + 2 * trim;
    let mut l = vec![0.0f64; maxmom];
    let sum: f64 = sorted.iter().sum();
    l[0] = sum / n;
    let mut lb = -n - 1.0;
    let mesh: Vec<f64> = (0..num_obs)
        .map(|_| {
            lb += 2.0;
            lb
        })
        .collect();
    let mut poly: Vec<f64> = mesh.iter().map(|v| v / (n - 1.0)).collect();
    let mut poly1 = vec![1.0f64; num_obs];
    let mut poly2 = vec![0.0f64; num_obs];
    let mut tmp: Vec<f64> = poly.iter().zip(sorted).map(|(p, x)| p * x).collect();
    l[1] = tmp.iter().sum::<f64>() / n;
    for j in 3..=maxmom {
        poly2.clone_from(&poly1);
        poly1.clone_from(&poly);
        let jf = j as f64;
        for i in 0..num_obs {
            poly[i] = ((2.0 * jf - 3.0) * mesh[i] * poly1[i]
                - (jf - 2.0) * (n + jf - 2.0) * poly2[i])
                / ((jf - 1.0) * (n - jf + 1.0));
            tmp[i] = poly[i] * sorted[i];
        }
        l[j - 1] = tmp.iter().sum::<f64>() / n;
    }
    if trim > 0 {
        for j in 1..=trim {
            maxmom -= 1;
            let jf = j as f64;
            for r in 1..(maxmom - j + 1) {
                let rf = r as f64;
                l[r - 1] = (l[r - 1] * (rf + jf) - l[r] * (rf + 1.0)) / (2.0 * rf + jf - 1.0);
            }
            l.truncate(maxmom - j);
        }
        for j in 1..=trim {
            // the C++ decrements maxmom here too, but never reads it again
            let jf = j as f64;
            let tf = trim as f64;
            let len = l.len();
            for r in 1..(len - 1) {
                let rf = r as f64;
                l[r - 1] = (l[r - 1] * (rf + jf + tf) + ((l[r] * (rf + 1.0) * (rf + tf)) / rf))
                    / (2.0 * rf + jf + tf - 1.0);
            }
            l.truncate(len - 1);
        }
    }
    for i in 2..l.len() {
        l[i] /= l[1];
    }
    l.iter_mut().for_each(|v| *v = v.abs());
    l
}

/// `scoreCandidateClusters` + `selectCandidateClusters`: score every candidate (dip and
/// L-moment shape of its margins), then greedily keep the best non-overlapping ones.
/// Returns a cluster number per cell (0 = unassigned).
pub fn select_candidate_clusters(candidates: &[Vec<usize>], data: &[Vec<f64>]) -> Vec<usize> {
    let n_cells = data.first().map_or(0, Vec::len);
    let n_cand = candidates.len();
    let col_num = data.len();
    let mut clustering = vec![0usize; n_cells];
    if n_cand == 0 {
        return clustering;
    }
    let mut dip_score = vec![0.0; n_cand];
    let mut linf = [vec![0.0; n_cand], vec![0.0; n_cand], vec![0.0; n_cand]];
    let mut tot = [vec![0.0; n_cand], vec![0.0; n_cand], vec![0.0; n_cand]];
    let mut tmp_dip = vec![0.0; col_num];
    let mut tmp_l = [vec![0.0; col_num], vec![0.0; col_num], vec![0.0; col_num]];
    let mut v = Vec::new();
    for (i, cand) in candidates.iter().enumerate() {
        for j in 0..col_num {
            v.clear();
            v.extend(cand.iter().map(|&k| data[j][k]));
            v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            tmp_dip[j] = single_dip(&v);
            let lm = abs_sample_lmoments(&v, 10, 2);
            tmp_l[0][j] = lm[1];
            tmp_l[1][j] = lm[2];
            tmp_l[2][j] = lm[3];
        }
        dip_score[i] = tmp_dip.iter().copied().fold(f64::INFINITY, f64::min);
        for m in 0..3 {
            linf[m][i] = tmp_l[m].iter().copied().fold(f64::NEG_INFINITY, f64::max);
            tot[m][i] = tmp_l[m].iter().sum();
        }
    }
    for m in 0..3 {
        let max_linf = linf[m].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let max_tot = tot[m].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        for i in 0..n_cand {
            linf[m][i] = 1.0 - linf[m][i] / max_linf;
            tot[m][i] = 1.0 - tot[m][i] / max_tot;
        }
    }
    let scores: Vec<f64> = (0..n_cand)
        .map(|i| {
            dip_score[i]
                + (linf[0][i] + tot[0][i]) / 2.0
                + (linf[1][i] + tot[1][i]) / 2.0
                + (linf[2][i] + tot[2][i]) / 2.0
        })
        .collect();
    // greedy selection; `anyIntersections` on sorted index lists
    let mut active = vec![true; n_cand];
    let mut current = 1usize;
    let mut num_active = n_cand;
    while num_active > 0 {
        let mut max_score = -1.0;
        let mut max_idx = 0;
        for i in 0..n_cand {
            if active[i] && scores[i] > max_score {
                max_score = scores[i];
                max_idx = i;
            }
        }
        let selected = &candidates[max_idx];
        active[max_idx] = false;
        for &k in selected {
            clustering[k] = current;
        }
        current += 1;
        for i in 0..n_cand {
            if active[i] && any_intersections(selected, &candidates[i]) {
                active[i] = false;
            }
        }
        num_active = active.iter().filter(|&&a| a).count();
    }
    clustering
}

fn any_intersections(v1: &[usize], v2: &[usize]) -> bool {
    if v1.is_empty() || v2.is_empty() {
        return false;
    }
    if v1[0] > v2[v2.len() - 1] || v1[v1.len() - 1] < v2[0] {
        return false;
    }
    // The C++ scans v2 from the start for every v1 value; sorted inputs make a merge exact.
    let (mut i, mut j) = (0, 0);
    while i < v1.len() && j < v2.len() {
        if v1[i] == v2[j] {
            return true;
        } else if v1[i] < v2[j] {
            i += 1;
        } else {
            j += 1;
        }
    }
    false
}

/// `annotateCluster`: per cell, `(score, boundaries)` per marker.
pub fn annotate_cluster(
    data: &[Vec<f64>],
    clustering: &[usize],
    boundaries: &[Vec<f64>],
    qs: &[f64; 3],
) -> Vec<Vec<(usize, usize)>> {
    let n_cols = data.len();
    let n_obs = clustering.len();
    let mut labels = vec![vec![(0usize, 0usize); n_cols]; n_obs];
    let mut uniq: Vec<usize> = clustering.to_vec();
    uniq.sort_unstable();
    uniq.dedup();
    for &c in &uniq {
        let lookup: Vec<usize> = (0..n_obs).filter(|&i| clustering[i] == c).collect();
        let mut lower = vec![0.0; n_cols];
        let mut med = vec![0.0; n_cols];
        let mut upper = vec![0.0; n_cols];
        let mut sub = Vec::with_capacity(lookup.len());
        for j in 0..n_cols {
            sub.clear();
            sub.extend(lookup.iter().map(|&i| data[j][i]));
            let q = r_quantile_cpp(&sub, qs);
            lower[j] = q[0];
            med[j] = q[1];
            upper[j] = q[2];
        }
        for &i in &lookup {
            for j in 0..n_cols {
                let b = &boundaries[j];
                let n_pts = b.len();
                let score = if n_pts == 1 {
                    let x = data[j][i];
                    let a = b[0];
                    if lower[j] >= a {
                        1
                    } else if upper[j] <= a {
                        0
                    } else if x >= a {
                        1
                    } else {
                        0
                    }
                } else {
                    b.iter().filter(|&&v| med[j] >= v).count()
                };
                labels[i][j] = (score, n_pts);
            }
        }
    }
    labels
}

/// `assignNumLabelToScoreVector`.
pub fn num_label(scored: &[(usize, usize)], names: &[String]) -> String {
    let mut s = String::new();
    for (j, &(score, n_pts)) in scored.iter().enumerate() {
        if score == 0 && n_pts == 0 {
            continue;
        }
        let set_score = score + 1;
        let set_size = n_pts + 1;
        if set_score <= set_size {
            s.push_str(&format!("{}~{}~{}~", names[j], set_score, set_size));
        } else {
            s.push_str(&format!("{}~AnnotationError~", names[j]));
        }
    }
    if s.is_empty() {
        "Uncertain".to_string()
    } else {
        s
    }
}

#[derive(Debug, Clone)]
pub struct ScampParams {
    pub gaussian_scale: f64,
    pub final_ann_qs: [f64; 3],
    pub seed: u64,
}

impl Default for ScampParams {
    fn default() -> Self {
        Self {
            gaussian_scale: 4.0,
            final_ann_qs: [0.4999, 0.5, 0.5001],
            seed: 123,
        }
    }
}

/// Outcome of one scamp run on one sample.
#[derive(Debug)]
pub struct ScampResult {
    /// Cluster number per cell (0 = never clustered).
    pub clustering: Vec<usize>,
    /// The label string per cell.
    pub labels: Vec<String>,
    pub n_candidates: usize,
    pub residual_rounds: usize,
}

/// scamp's `findCandidateClusters` on the single-thread random path: repeat
/// `candidateClusterSearch` with `seed += 1` until `maximumClusterNum` candidates are
/// collected. An iteration that finds nothing returns the whole set as one candidate and is
/// counted as a bad search when the data is larger than `10 * minClusterSize`; too many of
/// those abort the clustering in the C++ (which then loops forever in `cppNoisyScamp`), so it
/// is an unreachable case in practice and is left to the caller as an empty candidate list.
fn random_candidates(
    data: &[Vec<f64>],
    restriction: Option<&[Vec<u8>]>,
    sp: &SearchParams,
    seed: &mut u64,
) -> (Vec<Vec<usize>>, usize) {
    let n = data.first().map_or(0, Vec::len);
    let pathology_limit = sp.min_cluster_size * 10;
    let mut cands: Vec<Vec<usize>> = Vec::new();
    let mut bad = 0usize;
    while cands.len() < sp.max_clusters {
        *seed += 1;
        let r = search(data, restriction, sp, *seed);
        if r.vacuous && n > pathology_limit {
            bad += 1;
            if bad > 9 * (sp.max_clusters / 10) {
                return (Vec::new(), 0);
            }
        }
        cands.extend(r.candidates);
    }
    let n_candidates = cands.len();
    (cands, n_candidates)
}

/// `scamp(dataSet = raw selected markers, annForestVals = thresholds, randomSeed = seed, …)`
/// as FAUST calls it, through `cppNoisyScamp` with one iteration.
pub fn run_scamp(
    raw: &[Vec<f64>],
    names: &[String],
    thresholds: &[Vec<f64>],
    restriction: Option<&[Vec<u8>]>,
    params: &ScampParams,
) -> ScampResult {
    let n_cols = raw.len();
    let n_obs = raw.first().map_or(0, Vec::len);
    let scaled = r_scale(raw);
    let ann_gates = thresholds_to_scaled(raw, &scaled, thresholds);
    let sp = SearchParams::scamp(n_cols);
    // scamp(): noise, then the iterated random search; the seed keeps advancing through the
    // noise, every search iteration and every residual round (`randomSeed` is passed by
    // reference all the way down).
    let (noised, mut cur_seed) = add_noise_to_matrix(&scaled, params.gaussian_scale, params.seed);
    let (cands, n_candidates) = random_candidates(&noised, restriction, &sp, &mut cur_seed);
    let mut clustering = select_candidate_clusters(&cands, &noised);
    let mut c_max = clustering.iter().copied().max().unwrap_or(0);
    let mut residual_rounds = 0;
    loop {
        let zero_index: Vec<usize> = (0..n_obs).filter(|&i| clustering[i] == 0).collect();
        if zero_index.len() <= sp.min_cluster_size {
            break;
        }
        let subset: Vec<Vec<f64>> = (0..n_cols)
            .map(|j| zero_index.iter().map(|&i| noised[j][i]).collect())
            .collect();
        let still_viable = subset.iter().any(|col| {
            let mut s = col.clone();
            s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            single_dip(&s) < sp.dip_threshold
        });
        if !still_viable {
            break;
        }
        residual_rounds += 1;
        let res_subset: Option<Vec<Vec<u8>>> = restriction.map(|r| {
            (0..n_cols)
                .map(|j| zero_index.iter().map(|&i| r[j][i]).collect())
                .collect()
        });
        let (rcands, _) = random_candidates(&subset, res_subset.as_deref(), &sp, &mut cur_seed);
        // The C++ hands the scorer the FULL noised matrix while the candidates' indices are
        // positions in the subset, so the scores are computed on the wrong cells (the first
        // |subset| rows of every column). Reproduced as is: it is what the reference does.
        let rclust = select_candidate_clusters(&rcands, &noised);
        for (k, &c) in rclust.iter().enumerate() {
            if c > 0 {
                clustering[zero_index[k]] = c + c_max;
            }
        }
        c_max = clustering.iter().copied().max().unwrap_or(0);
    }
    let scored = annotate_cluster(&scaled, &clustering, &ann_gates, &params.final_ann_qs);
    let labels = scored.iter().map(|s| num_label(s, names)).collect();
    ScampResult {
        clustering,
        labels,
        n_candidates,
        residual_rounds,
    }
}
