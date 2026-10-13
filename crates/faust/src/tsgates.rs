//! Gate locations from the taut-string density and 1-D k-medians.
//!
//! Port of `src/cpPmden.cpp` + `src/tautstring.cpp` + `src/kkuiper.cpp` (Davies–Kovac taut
//! string density estimate with the Kuiper-based bandwidth search), `src/kMedDP.cpp` +
//! `src/fill_SMAWK.cpp` + `src/dynamic_prog.cpp` (exact 1-D k-medians by dynamic programming),
//! `src/tsGates.cpp` (`tsGates`: antimodes of the taut string, else k-medians boundaries) and
//! `src/tsModeEstimate.cpp`, `src/getTautStringApprox.cpp`.
//!
//! All inputs are **sorted ascending**. Outputs are the cut points **including the minimum and
//! the maximum** of the data, exactly as the C++ returns them.
//!
//! The arithmetic is kept in the C++ order so the results are bit-identical to R FAUST 0.5.5.
//! A few C++ quirks are reproduced on purpose and marked `QUIRK` below.
//!
//! Parity: `tests/parity_tsgates.rs` against `*_tsGates0.txt`, `*_tsGates2.txt`,
//! `*_tsGates3.txt` and `*_tautString.txt`.

/// What `cpPmden` returns (the C++ `stringInfo`).
///
/// `string` has `n - 1` values: the slope of the taut string on each gap between consecutive
/// observations, i.e. the density estimate. `knotst`/`knotsy` are the knot locations and heights
/// *after* the C++ post-processing (see `taut_string`), `nmax` the number of local maxima the
/// string sees (the mode count FAUST uses is `nmax - nmax / 2`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StringInfo {
    pub string: Vec<f64>,
    pub knotsind: Vec<usize>,
    pub knotst: Vec<f64>,
    pub knotsy: Vec<f64>,
    pub nknots: usize,
    pub nmax: usize,
}

// ------------------------------------------------------------------------------------------
// Kuiper-metric bandwidth bounds (pmden.R of ftnonpar, as truncated in cpPmden.cpp)
// ------------------------------------------------------------------------------------------

const MAX_KUIP_NR: usize = 19;

/// Safety cap on the bandwidth-search passes of [`cp_pmden`] (the C++ has none; real data
/// converges in tens of passes, the cap only stops degenerate inputs from spinning).
const MAX_PMDEN_ITERATIONS: usize = 20_000;

#[rustfmt::skip]
const KUIP_DIFF_BOUNDS: [[f64; 7]; MAX_KUIP_NR] = [
    [0.30750198, 0.22031415, 0.15749822, 0.10121863, 0.07178070, 0.05080985, 0.03215992],
    [0.26097572, 0.18867970, 0.13705738, 0.08810506, 0.06335109, 0.04492857, 0.02858946],
    [0.14897866, 0.11094044, 0.08092217, 0.05257830, 0.03748763, 0.02696284, 0.01707126],
    [0.12438262, 0.09529240, 0.06950244, 0.04524029, 0.03286277, 0.02348616, 0.01502098],
    [0.09319795, 0.07219891, 0.05339150, 0.03538433, 0.02598238, 0.01851747, 0.01174251],
    [0.08018187, 0.06280868, 0.04765458, 0.03164710, 0.02301769, 0.01648256, 0.01055136],
    [0.066165059, 0.052643314, 0.040499260, 0.027003915, 0.019671301, 0.014140314, 0.009111927],
    [0.05768921, 0.04729112, 0.03630651, 0.02443004, 0.01784079, 0.01289449, 0.00831476],
    [0.049727442, 0.041269518, 0.031941519, 0.021631778, 0.015981589, 0.011540437, 0.007469195],
    [0.044093021, 0.037163268, 0.029176417, 0.019971113, 0.014785146, 0.010818339, 0.006926491],
    [0.038637596, 0.033532556, 0.026499289, 0.018379105, 0.013643814, 0.009856368, 0.006386310],
    [0.035037285, 0.030666602, 0.024677670, 0.017132430, 0.012725232, 0.009228679, 0.005995483],
    [0.031642378, 0.028084749, 0.022794216, 0.015823040, 0.011929053, 0.008645996, 0.005649372],
    [0.029033131, 0.025996755, 0.021296584, 0.014971161, 0.011244505, 0.008220232, 0.005350641],
    [0.025883091, 0.024057281, 0.020032581, 0.014186956, 0.010545471, 0.007767515, 0.005061438],
    [0.023257645, 0.022536446, 0.018824915, 0.013486222, 0.010061760, 0.007427072, 0.004817209],
    [0.020675023, 0.020836156, 0.017644365, 0.012805928, 0.009588957, 0.007102941, 0.004595116],
    [0.019136349, 0.019431774, 0.016733216, 0.012146658, 0.009200609, 0.006771573, 0.004431340],
    [0.017891379, 0.018348230, 0.015944010, 0.011620645, 0.008794318, 0.006483762, 0.004265070],
];

const KUIP_DIFF_BOUNDS_X: [f64; 7] = [50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0];

// ------------------------------------------------------------------------------------------
// Small R-like helpers (cppApprox.cpp, rQuantile.cpp, medianAbsoluteDeviation.cpp)
// ------------------------------------------------------------------------------------------

/// `cppApprox`: linear interpolation of `(x, y)` at `v` (R `approx(..., rule = 2)`), `x`
/// ascending. Returns zeros when `x` is empty (the C++ would read out of bounds).
pub fn approx_linear(x: &[f64], y: &[f64], v: &[f64]) -> Vec<f64> {
    let nxy = x.len();
    let mut rvec = vec![0.0; v.len()];
    if nxy == 0 || y.len() < nxy {
        return rvec;
    }
    for (k, &cv) in v.iter().enumerate() {
        let mut i = 0usize;
        let mut j = nxy - 1;
        if cv < x[i] {
            rvec[k] = y[0];
            continue;
        }
        if cv > x[j] {
            rvec[k] = y[nxy - 1];
            continue;
        }
        while i + 1 < j {
            let ij = (i + j) / 2;
            if cv < x[ij] {
                j = ij;
            } else {
                i = ij;
            }
        }
        if cv == x[j] {
            rvec[k] = y[j];
        } else if cv == x[i] {
            rvec[k] = y[i];
        } else {
            rvec[k] = y[i] + (y[j] - y[i]) * ((cv - x[i]) / (x[j] - x[i]));
        }
    }
    rvec
}

/// `rQuantile`: R `quantile(type = 7)` of `data` (sorted internally) at probability `p`.
/// Returns `NaN` on empty input.
pub fn r_quantile(data: &[f64], p: f64) -> f64 {
    if data.is_empty() {
        return f64::NAN;
    }
    let mut d = data.to_vec();
    d.sort_by(f64::total_cmp);
    let n = d.len();
    let index = (n as f64 - 1.0) * p;
    let low = index.floor();
    let high = index.ceil();
    let li = (low.max(0.0) as usize).min(n - 1);
    let hi = (high.max(0.0) as usize).min(n - 1);
    let mut q = d[li];
    if index > low && d[hi] != q {
        let h_rem = index - low;
        q = ((1.0 - h_rem) * q) + (h_rem * d[hi]);
    }
    q
}

/// `medianAbsoluteDeviation`: R `mad(x)` (constant 1.4826).
pub fn median_absolute_deviation(data: &[f64]) -> f64 {
    let median = r_quantile(data, 0.5);
    let devs: Vec<f64> = data.iter().map(|v| (v - median).abs()).collect();
    r_quantile(&devs, 0.5) * 1.4826
}

// ------------------------------------------------------------------------------------------
// tautstring.cpp (strings.c of ftnonpar)
// ------------------------------------------------------------------------------------------

/// `tautString`: the taut string through the tube `[lower, upper]` around the empirical
/// distribution `fdist` at abscissae `t`, pinned to `y1` at the left and `yn` at the right.
/// `extrmean` is the `EXTRMEAN` flag (FAUST passes 1).
///
/// QUIRK (reproduced): the C++ strips **every** zero from `knotst` (meant for the unused tail)
/// and then truncates `knotsy` to that length, so a knot at `t == 0` shifts `knotsy` against
/// `knotst`. `knotsind` is returned truncated to `nknots`.
#[allow(clippy::too_many_arguments)]
pub fn taut_string(
    fdist: &[f64],
    t: &[f64],
    lower: &[f64],
    upper: &[f64],
    y1: f64,
    yn: f64,
    n: usize,
    extrmean: bool,
) -> StringInfo {
    if n < 2 || t.len() < n || fdist.len() < n || lower.len() < n || upper.len() < n {
        return StringInfo::default();
    }
    let mut i_taut_string = vec![0.0f64; n - 1];
    // One spare slot: the C++ sizes these `n` and relies on knot indices being strictly
    // increasing; a tie in `t` (0/0 = NaN derivatives) can break that, so we also stop adding
    // knots when the arrays are full instead of writing past the end.
    let mut knotsind = vec![0usize; n + 1];
    let mut knotst = vec![0.0f64; n + 1];
    let mut knotsy = vec![0.0f64; n + 1];
    let mut knotssign = vec![0i32; n + 1];
    let mut nknots = 1usize;
    let mut nmax = 0usize;

    let mut maxderiv = 1e+38f64;
    let mut minderiv = -maxderiv;
    let mut actind = 2usize;
    let mut maxind = 1usize;
    let mut minind = 1usize;
    let mut lastbound = 0i32;

    knotssign[0] = lastbound;
    knotsind[0] = 1;
    knotsy[0] = y1;
    knotst[0] = t[0];

    while actind <= n {
        let (newmaxderiv, newminderiv) = if actind < n {
            (
                (upper[actind - 1] - knotsy[nknots - 1]) / (t[actind - 1] - knotst[nknots - 1]),
                (lower[actind - 1] - knotsy[nknots - 1]) / (t[actind - 1] - knotst[nknots - 1]),
            )
        } else {
            let d = (yn - knotsy[nknots - 1]) / (t[actind - 1] - knotst[nknots - 1]);
            (d, d)
        };
        if newminderiv > maxderiv {
            if nknots >= n {
                break;
            }
            if lastbound == -1 {
                nmax += 1;
            }
            knotssign[nknots] = 1;
            knotsind[nknots] = maxind;
            knotsy[nknots] = upper[maxind - 1];
            knotst[nknots] = t[maxind - 1];
            nknots += 1;
            actind = maxind;
            lastbound = 1;
            maxderiv = 1e+38;
            minderiv = -1e+38;
        } else if newmaxderiv < minderiv {
            if nknots >= n {
                break;
            }
            if lastbound == 1 {
                nmax += 1;
            }
            knotsind[nknots] = minind;
            knotssign[nknots] = -1;
            knotsy[nknots] = lower[minind - 1];
            knotst[nknots] = t[minind - 1];
            nknots += 1;
            maxderiv = 1e+38;
            minderiv = -1e+38;
            lastbound = -1;
            actind = minind;
        } else {
            if newmaxderiv < maxderiv {
                maxderiv = newmaxderiv;
                maxind = actind;
            }
            if newminderiv > minderiv {
                minderiv = newminderiv;
                minind = actind;
            }
            if actind == n {
                if lastbound != 0 {
                    knotsind[nknots] = actind;
                    knotsy[nknots] = yn;
                    knotst[nknots] = t[actind - 1];
                    knotssign[nknots] = 0;
                    nknots += 1;
                } else {
                    lastbound = -1;
                    knotsind[0] = 1;
                    knotssign[0] = -1;
                    knotsy[0] = y1;
                    knotst[0] = t[0];
                    actind = 1;
                    minind = 1;
                    maxind = 1;
                }
            }
        }
        actind += 1;
    }

    for i in 0..nknots - 1 {
        let (lo, hi) = (knotsind[i], knotsind[i + 1]);
        for j in lo..hi {
            if j == 0 || j > n - 1 {
                continue;
            }
            i_taut_string[j - 1] = if knotssign[i] == knotssign[i + 1] || !extrmean {
                (knotsy[i + 1] - knotsy[i]) / (knotst[i + 1] - knotst[i])
            } else {
                (fdist[knotsind[i + 1] - 1] - fdist[knotsind[i] - 1]) / (knotst[i + 1] - knotst[i])
            };
        }
    }

    // QUIRK: remove *all* zeros (not just the unused tail), then truncate knotsy to match.
    knotst.retain(|&v| v != 0.0);
    knotsy.truncate(knotst.len());
    knotsind.truncate(nknots);
    StringInfo {
        string: i_taut_string,
        knotsind,
        knotst,
        knotsy,
        nknots,
        nmax,
    }
}

// ------------------------------------------------------------------------------------------
// kkuiper.cpp (kkuip.c of ftnonpar)
// ------------------------------------------------------------------------------------------

/// Scratch buffers for `difficultmax` (the C++ allocates four `n`-vectors per call).
struct KuiperScratch {
    mins: Vec<f64>,
    maxs: Vec<f64>,
    minis: Vec<usize>,
    maxis: Vec<usize>,
}

fn easymax(x: &[f64], left: usize, right: usize) -> (usize, usize, f64) {
    let mut min = 1e+38f64;
    let mut max = -1e+38f64;
    let mut mini = left;
    let mut maxi = left;
    for (i, &xi) in x.iter().enumerate().take(right + 1).skip(left) {
        if xi > max {
            max = xi;
            maxi = i;
        }
        if xi < min {
            min = xi;
            mini = i;
        }
    }
    if mini < maxi {
        (mini, maxi, max - min)
    } else {
        (maxi, mini, max - min)
    }
}

fn difficultmax(
    x: &[f64],
    left: usize,
    right: usize,
    s: &mut KuiperScratch,
) -> (usize, usize, f64) {
    let (mins, maxs, minis, maxis) = (&mut s.mins, &mut s.maxs, &mut s.minis, &mut s.maxis);
    let mut max = -1e+38f64;
    let mut maxi = left;
    let (erga, ergb);
    if x[left] < x[right] {
        maxs[left] = x[left];
        maxis[left] = left;
        for i in left + 1..=right {
            if x[i] > maxs[i - 1] {
                maxs[i] = x[i];
                maxis[i] = i;
            } else {
                maxs[i] = maxs[i - 1];
                maxis[i] = maxis[i - 1];
            }
        }
        mins[right] = x[right];
        minis[right] = right;
        for i in (left..right).rev() {
            if x[i] < mins[i + 1] {
                mins[i] = x[i];
                minis[i] = i;
            } else {
                mins[i] = mins[i + 1];
                minis[i] = minis[i + 1];
            }
        }
        for i in left..=right {
            if maxs[i] - mins[i] > max {
                max = maxs[i] - mins[i];
                maxi = i;
            }
        }
        erga = maxis[maxi];
        ergb = minis[maxi];
    } else {
        mins[left] = x[left];
        minis[left] = left;
        for i in left + 1..=right {
            if x[i] < mins[i - 1] {
                mins[i] = x[i];
                minis[i] = i;
            } else {
                mins[i] = mins[i - 1];
                minis[i] = minis[i - 1];
            }
        }
        maxs[right] = x[right];
        maxis[right] = right;
        for i in (left..right).rev() {
            if x[i] > maxs[i + 1] {
                maxs[i] = x[i];
                maxis[i] = i;
            } else {
                maxs[i] = maxs[i + 1];
                maxis[i] = maxis[i + 1];
            }
        }
        for i in left..=right {
            if maxs[i] - mins[i] > max {
                max = maxs[i] - mins[i];
                maxi = i;
            }
        }
        erga = minis[maxi];
        ergb = maxis[maxi];
    }
    let ergmax =
        (x[erga] - x[left]).abs() + (x[right] - x[ergb]).abs() - (x[right] - x[left]).abs();
    (erga, ergb, ergmax)
}

/// `kkuiper`: the first `kmax` Kuiper norms of `x` (the difference between the empirical
/// distribution and the taut string). Entries past an early exit stay 0.
pub fn kkuiper(x: &[f64], kmax: usize) -> Vec<f64> {
    let n = x.len();
    let mut norm = vec![0.0f64; kmax];
    if n == 0 || kmax == 0 {
        return norm;
    }
    let mut a = vec![0usize; kmax + 1];
    let mut b = vec![0usize; kmax + 1];
    let mut scratch = KuiperScratch {
        mins: vec![0.0; n],
        maxs: vec![0.0; n],
        minis: vec![0; n],
        maxis: vec![0; n],
    };

    let mut min = 1e+38f64;
    let mut max = -1e+38f64;
    let mut mini = 0usize;
    let mut maxi = 0usize;
    for (i, &xi) in x.iter().enumerate() {
        if xi > max {
            max = xi;
            maxi = i;
        }
        if xi < min {
            min = xi;
            mini = i;
        }
    }
    if mini < maxi {
        a[0] = mini;
        b[0] = maxi;
    } else {
        b[0] = mini;
        a[0] = maxi;
    }
    norm[0] = (x[b[0]] - x[a[0]]).abs();

    for k in 1..kmax {
        let mut abmax = -1e+38f64;
        let mut einordnen: isize = -1;
        let mut maxa = 0usize;
        let mut maxb = 0usize;
        if a[0] > 0 {
            let (ea, eb, em) = easymax(x, 0, a[0]);
            maxa = ea;
            maxb = eb;
            abmax = em;
            einordnen = 0;
        }
        for i in 0..k {
            let (newa, newb, newmax) = difficultmax(x, a[i], b[i], &mut scratch);
            if newmax > abmax {
                maxa = newa;
                maxb = newb;
                abmax = newmax;
                einordnen = i as isize;
            }
            let right = if i < k - 1 { a[i + 1] } else { n - 1 };
            let (newa, newb, newmax) = easymax(x, b[i], right);
            if newmax > abmax {
                maxa = newa;
                maxb = newb;
                abmax = newmax;
                einordnen = i as isize + 1;
            }
        }
        if einordnen < 0 {
            break;
        }
        let e = einordnen as usize;
        if e == k {
            a[k] = maxa;
            b[k] = maxb;
            norm[k] = norm[k - 1] + (x[maxb] - x[maxa]).abs();
        } else {
            for i in (e + 1..=k).rev() {
                a[i] = a[i - 1];
                b[i] = b[i - 1];
            }
            if maxa < a[e] {
                a[e] = maxa;
                b[e] = maxb;
                norm[k] = norm[k - 1] + (x[maxb] - x[maxa]).abs();
            } else {
                b[e] = maxa;
                a[e + 1] = maxb;
                norm[k] = norm[k - 1] - (x[a[e]] - x[b[e + 1]]).abs()
                    + (x[a[e]] - x[b[e]]).abs()
                    + (x[a[e + 1]] - x[b[e + 1]]).abs();
            }
        }
    }
    norm
}

// ------------------------------------------------------------------------------------------
// local_density.cpp (denlocal.f of ftnonpar)
// ------------------------------------------------------------------------------------------

/// `local_density`: marks in `kni` the observations where the taut string `cf` is locally too
/// dense (`icomax`) or too sparse (`irmax`) relative to the uniform grid `j / n`.
fn local_density(cf: &[f64], kni: &mut [u8], n: usize, icomax: usize, irmax: usize) {
    // 1-based I, J as in the Fortran.
    let mut i = 1usize;
    let mut j = 1usize;
    let mut icr = 0usize;
    loop {
        let mut ic = 0usize;
        loop {
            if i > n {
                break;
            }
            if cf[i - 1] <= (j as f64) / (n as f64) {
                if icr >= irmax {
                    kni[i - 1] = 1;
                    // C++ writes kni[I-2]; at I == 1 that is out of bounds (never reached:
                    // ICR is 0 the first time through).
                    if i >= 2 {
                        kni[i - 2] = 1;
                    }
                }
                icr = 0;
                ic += 1;
                if i < n {
                    i += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        if ic == 0 {
            icr += 1;
        }
        if ic >= icomax {
            for v in kni.iter_mut().take(i).skip(i - ic) {
                *v = 1;
            }
        }
        if i == n {
            break;
        }
        j += 1;
    }
}

// ------------------------------------------------------------------------------------------
// cpPmden.cpp (pmden of ftnonpar)
// ------------------------------------------------------------------------------------------

/// `cpPmden`: the taut-string density of sorted `x` with the Kuiper-metric bandwidth search
/// and (for `n < 10000`) local squeezing. Returns an empty [`StringInfo`] for `n < 2`.
pub fn cp_pmden(x_in: &[f64]) -> StringInfo {
    let nsamp = x_in.len();
    if nsamp < 2 {
        return StringInfo::default();
    }
    let mut x: Vec<f64> = x_in.to_vec();
    let locsq_factor = 0.9f64;
    let extrema_mean = true;

    let mut currbounds = [0.0f64; MAX_KUIP_NR];
    if nsamp > 5000 {
        let kuiper_multiplier = 5000.0f64.sqrt() / (nsamp as f64).sqrt();
        for (cb, row) in currbounds.iter_mut().zip(KUIP_DIFF_BOUNDS.iter()) {
            *cb = row[6] * kuiper_multiplier;
        }
    } else {
        let tmp_nsamp = [nsamp as f64];
        for (cb, row) in currbounds.iter_mut().zip(KUIP_DIFF_BOUNDS.iter()) {
            *cb = approx_linear(&KUIP_DIFF_BOUNDS_X, row, &tmp_nsamp)[0];
        }
    }

    let mut dataemp = vec![1.0 / ((nsamp - 1) as f64); nsamp];
    dataemp[0] = 0.0;

    let mut acc_min = 1.0f64;
    let max_diff = x[nsamp - 1] - x[0];
    for i in 0..nsamp - 1 {
        let current_ratio = (x[i + 1] - x[i]) / max_diff;
        if current_ratio < acc_min {
            acc_min = current_ratio;
        }
    }
    if acc_min < 1e-14 {
        // QUIRK: the C++ "accuracy" jitter loop computes `newaccx` but never stores it, so
        // the only effect is the round trip x -> x / accuracy -> * accuracy. The loop itself
        // (O(n^2)) has no side effects and is omitted.
        let accuracy = median_absolute_deviation(x_in) / 1000.0;
        for v in x.iter_mut() {
            *v = (*v / accuracy) * accuracy;
        }
    }

    let fdist_y: Vec<f64> = (0..nsamp)
        .map(|i| (i as f64) / ((nsamp - 1) as f64))
        .collect();
    let mut eps = vec![0.5f64; nsamp];
    eps[0] = 0.0;
    eps[nsamp - 1] = 0.0;

    let mut lower = vec![0.0f64; nsamp];
    let mut upper = vec![0.0f64; nsamp];
    let mut differences = vec![0.0f64; nsamp];
    let mut empirical_sum = vec![0.0f64; nsamp];
    let mut kni = vec![0u8; nsamp];
    let mut fts = StringInfo::default();

    // std::partial_sum: repeated addition (not i / (n - 1)).
    {
        let mut acc = dataemp[0];
        empirical_sum[0] = acc;
        for i in 1..nsamp {
            acc += dataemp[i];
            empirical_sum[i] = acc;
        }
    }

    let mut still_repeating = true;
    let mut iterations = 0usize;
    while still_repeating {
        for i in 0..nsamp {
            lower[i] = fdist_y[i] - eps[i];
            upper[i] = fdist_y[i] + eps[i];
        }
        fts = taut_string(
            &fdist_y,
            &x,
            &lower,
            &upper,
            upper[0],
            lower[nsamp - 1],
            nsamp,
            extrema_mean,
        );
        let x_string = &fts.string;

        let mut first_jump: Option<usize> = None;
        let mut last_jump: Option<usize> = None;
        for i in 0..nsamp.saturating_sub(2) {
            if x_string[i + 1] != x_string[i] {
                if first_jump.is_none() {
                    first_jump = Some(i);
                }
                last_jump = Some(i);
            }
        }
        if let (Some(fj), Some(lj)) = (first_jump, last_jump) {
            if x_string[fj] > x_string[fj + 1] {
                fts.nmax += 1;
            }
            if x_string[lj] < x_string[lj + 1] {
                fts.nmax += 1;
            }
        }

        let lastunif = approx_linear(&fts.knotst, &fts.knotsy, &x);
        for i in 0..nsamp {
            differences[i] = empirical_sum[i] - lastunif[i];
        }

        let currkkuip = kkuiper(&differences, MAX_KUIP_NR);
        let mut first_kuiper_ind: Option<usize> = None;
        if currkkuip[0] > (currbounds[0] + 1e-08) {
            first_kuiper_ind = Some(0);
        }
        for i in 1..MAX_KUIP_NR {
            if (currkkuip[i] - currkkuip[i - 1]) > (currbounds[i] + 1e-08)
                && first_kuiper_ind.is_none()
            {
                first_kuiper_ind = Some(i);
            }
        }

        iterations += 1;
        if let Some(fk) = first_kuiper_ind {
            let cur_kuiper_val = currbounds[fk] / 2.0;
            let mut changed = false;
            for e in eps.iter_mut() {
                if *e > 0.0 {
                    changed |= e.to_bits() != cur_kuiper_val.to_bits();
                    *e = cur_kuiper_val;
                }
            }
            // Guard (not in the C++): the loop is deterministic in `eps`, so an unchanged tube
            // means the C++ would spin forever (e.g. all observations tied, max_diff == 0).
            if !changed || iterations >= MAX_PMDEN_ITERATIONS {
                still_repeating = false;
            }
        } else if nsamp < 10000 {
            let irmax = ((nsamp as f64).ln().floor() as i64 + 3).max(0) as usize;
            let mut icomax = 1usize;
            let mut prp = (-1.0f64).exp();
            let mut currsum = 2.0 * (-1.0f64).exp();
            while currsum.ln() < (0.95f64.ln() / (nsamp as f64)) {
                icomax += 1;
                prp /= icomax as f64;
                currsum += prp;
            }
            kni.iter_mut().for_each(|v| *v = 0);
            local_density(&lastunif, &mut kni, nsamp, icomax, irmax);
            let mut kni_sum = 0usize;
            for (e, &k) in eps.iter_mut().zip(kni.iter()) {
                if k == 1 {
                    kni_sum += 1;
                    *e *= locsq_factor;
                }
            }
            if kni_sum == 0 || iterations >= MAX_PMDEN_ITERATIONS {
                still_repeating = false;
            }
        } else {
            still_repeating = false;
        }
    }
    fts
}

// ------------------------------------------------------------------------------------------
// kMedDP.cpp + fill_SMAWK.cpp + dynamic_prog.cpp (Ckmeans.1d.dp, L1 criterion)
// ------------------------------------------------------------------------------------------

/// `sabs`: the within-cluster L1 dissimilarity of `x[j..=i]` from prefix sums of shifted `x`.
#[inline]
fn sabs(j: usize, i: usize, sum_x: &[f64]) -> f64 {
    let sji = if j >= i {
        0.0
    } else if j > 0 {
        let l = (i + j) >> 1;
        if (i - j + 1) % 2 == 1 {
            -sum_x[l - 1] + sum_x[j - 1] + sum_x[i] - sum_x[l]
        } else {
            -sum_x[l] + sum_x[j - 1] + sum_x[i] - sum_x[l]
        }
    } else {
        let l = i >> 1;
        if (i + 1) % 2 == 1 {
            -sum_x[l - 1] + sum_x[i] - sum_x[l]
        } else {
            -sum_x[l] + sum_x[i] - sum_x[l]
        }
    };
    if sji < 0.0 { 0.0 } else { sji }
}

/// The two DP rows SMAWK works on: `q - 1` (read) and `q` (written).
struct DpRows<'a> {
    s_prev: &'a [f64],
    j_prev: &'a [usize],
    s_cur: &'a mut [f64],
    j_cur: &'a mut [usize],
    sum_x: &'a [f64],
}

fn reduce_in_place(imin: i64, imax: i64, istep: i64, js: &[usize], rows: &DpRows) -> Vec<usize> {
    let n_rows = (imax - imin) / istep + 1;
    let mut js_red = js.to_vec();
    if n_rows >= js.len() as i64 {
        return js_red;
    }
    let mut left: i64 = -1;
    let mut right: i64 = 0;
    let mut m = js_red.len();
    while (m as i64) > n_rows {
        let p = left + 1;
        let i = (imin + p * istep) as usize;
        let j = js_red[right as usize];
        let sl = rows.s_prev[j - 1] + sabs(j, i, rows.sum_x);
        let jplus1 = js_red[(right + 1) as usize];
        let slplus1 = rows.s_prev[jplus1 - 1] + sabs(jplus1, i, rows.sum_x);
        if sl < slplus1 && p < n_rows - 1 {
            left += 1;
            js_red[left as usize] = j;
            right += 1;
        } else if sl < slplus1 && p == n_rows - 1 {
            right += 1;
            js_red[right as usize] = j;
            m -= 1;
        } else {
            if p > 0 {
                js_red[right as usize] = js_red[left as usize];
                left -= 1;
            } else {
                right += 1;
            }
            m -= 1;
        }
    }
    let mut r = (left + 1) as usize;
    while r < m {
        js_red[r] = js_red[right as usize];
        right += 1;
        r += 1;
    }
    js_red.truncate(m);
    js_red
}

fn fill_even_positions(imin: i64, imax: i64, istep: i64, js: &[usize], rows: &mut DpRows) {
    let n = js.len();
    let istepx2 = istep << 1;
    let mut jl = js[0];
    let mut r = 0usize;
    let mut i = imin;
    while i <= imax {
        let iu = i as usize;
        while r < n && js[r] < jl {
            r += 1;
        }
        if r >= n {
            break;
        }
        rows.s_cur[iu] = rows.s_prev[js[r] - 1] + sabs(js[r], iu, rows.sum_x);
        rows.j_cur[iu] = js[r];

        let jh = if i + istep <= imax {
            rows.j_cur[(i + istep) as usize]
        } else {
            js[n - 1]
        };
        let jmax = jh.min(iu);
        let sjimin = sabs(jmax, iu, rows.sum_x);

        r += 1;
        while r < n && js[r] <= jmax {
            let jabs = js[r];
            if jabs > iu {
                break;
            }
            if jabs < rows.j_prev[iu] {
                r += 1;
                continue;
            }
            let s = sabs(jabs, iu, rows.sum_x);
            let sj = rows.s_prev[jabs - 1] + s;
            if sj <= rows.s_cur[iu] {
                rows.s_cur[iu] = sj;
                rows.j_cur[iu] = js[r];
            } else if rows.s_prev[jabs - 1] + sjimin > rows.s_cur[iu] {
                break;
            }
            r += 1;
        }
        r -= 1;
        jl = jh;
        i += istepx2;
    }
}

fn find_min_from_candidates(imin: i64, imax: i64, istep: i64, js: &[usize], rows: &mut DpRows) {
    let mut rmin_prev = 0usize;
    let mut i = imin;
    while i <= imax {
        let iu = i as usize;
        let rmin = rmin_prev;
        rows.s_cur[iu] = rows.s_prev[js[rmin] - 1] + sabs(js[rmin], iu, rows.sum_x);
        rows.j_cur[iu] = js[rmin];
        for (r, &j_abs) in js.iter().enumerate().skip(rmin + 1) {
            if j_abs < rows.j_prev[iu] {
                continue;
            }
            if j_abs > iu {
                break;
            }
            let sj = rows.s_prev[j_abs - 1] + sabs(j_abs, iu, rows.sum_x);
            if sj <= rows.s_cur[iu] {
                rows.s_cur[iu] = sj;
                rows.j_cur[iu] = j_abs;
                rmin_prev = r;
            }
        }
        i += istep;
    }
}

fn smawk(imin: i64, imax: i64, istep: i64, js: &[usize], rows: &mut DpRows) {
    if js.is_empty() {
        return;
    }
    if imax - imin <= 0 {
        find_min_from_candidates(imin, imax, istep, js, rows);
    } else {
        let js_odd = reduce_in_place(imin, imax, istep, js, rows);
        let istepx2 = istep << 1;
        let imin_odd = imin + istep;
        let imax_odd = imin_odd + (imax - imin_odd) / istepx2 * istepx2;
        smawk(imin_odd, imax_odd, istepx2, &js_odd, rows);
        fill_even_positions(imin, imax, istep, js, rows);
    }
}

/// `kMedDP`: exact 1-D k-medians of sorted `x` by dynamic programming (Ckmeans.1d.dp with the
/// L1 criterion, SMAWK row fill). Returns the 0-based cluster of each observation.
///
/// `O(k·n)` time, `O(k·n)` memory for the backtrack matrix. Returns an empty `Vec` when `x`
/// is empty or `k == 0`; `k` is clamped to `n` (the C++ requires `k <= n`).
pub fn k_med_dp(x: &[f64], k: usize) -> Vec<usize> {
    let n = x.len();
    if n == 0 || k == 0 {
        return Vec::new();
    }
    let k = k.min(n);

    // fill_dp_matrix (equally weighted, L1: sum_x_sq is never read).
    let shift = x[n / 2];
    let mut sum_x = vec![0.0f64; n];
    sum_x[0] = x[0] - shift;
    for i in 1..n {
        sum_x[i] = sum_x[i - 1] + x[i] - shift;
    }
    let mut s_prev = vec![0.0f64; n];
    let mut s_cur = vec![0.0f64; n];
    let mut j_rows: Vec<Vec<usize>> = Vec::with_capacity(k);
    j_rows.push(vec![0usize; n]);
    for (i, s) in s_prev.iter_mut().enumerate().skip(1) {
        *s = sabs(0, i, &sum_x);
    }
    for q in 1..k {
        let imin: i64 = if q < k - 1 {
            1.max(q) as i64
        } else {
            (n - 1) as i64
        };
        let imax: i64 = (n - 1) as i64;
        let mut j_cur = vec![0usize; n];
        {
            let js: Vec<usize> = (q..=n - 1).collect();
            let mut rows = DpRows {
                s_prev: &s_prev,
                j_prev: &j_rows[q - 1],
                s_cur: &mut s_cur,
                j_cur: &mut j_cur,
                sum_x: &sum_x,
            };
            smawk(imin, imax, 1, &js, &mut rows);
        }
        j_rows.push(j_cur);
        std::mem::swap(&mut s_prev, &mut s_cur);
    }

    // backtrack_L1 (centers / withinss are not part of the result).
    let mut cluster = vec![0usize; n];
    let mut cluster_right = n - 1;
    for q in (0..k).rev() {
        let cluster_left = j_rows[q][cluster_right].min(cluster_right);
        for c in cluster
            .iter_mut()
            .take(cluster_right + 1)
            .skip(cluster_left)
        {
            *c = q;
        }
        if q > 0 {
            if cluster_left == 0 {
                break;
            }
            cluster_right = cluster_left - 1;
        }
    }
    cluster
}

/// `findKmedGates`: `[min(x), max of class 0, ..., max of class k-1, max(x)]` (the C++ loops
/// over classes `0..k_estimate` where `k_estimate = max(class)`, so the last class's maximum
/// appears once, as `max(x)`).
pub fn find_kmed_gates(x: &[f64], class: &[usize], num_classes: usize) -> Vec<f64> {
    let mut gates = Vec::with_capacity(num_classes + 2);
    if x.is_empty() {
        return gates;
    }
    gates.push(min_element(x));
    let mut class_max: Vec<Option<f64>> = vec![None; num_classes];
    for (&v, &c) in x.iter().zip(class.iter()) {
        if c < num_classes {
            // std::max_element semantics: replace only when strictly greater.
            class_max[c] = Some(match class_max[c] {
                Some(m) if m < v => v,
                Some(m) => m,
                None => v,
            });
        }
    }
    for m in class_max.into_iter().flatten() {
        gates.push(m);
    }
    gates.push(max_element(x));
    gates
}

fn min_element(x: &[f64]) -> f64 {
    let mut m = x[0];
    for &v in &x[1..] {
        if v < m {
            m = v;
        }
    }
    m
}

fn max_element(x: &[f64]) -> f64 {
    let mut m = x[0];
    for &v in &x[1..] {
        if m < v {
            m = v;
        }
    }
    m
}

fn kmed_gates(x: &[f64], num_centers: usize) -> Vec<f64> {
    let k_clusters = k_med_dp(x, num_centers);
    let k_estimate = k_clusters.iter().copied().max().unwrap_or(0);
    find_kmed_gates(x, &k_clusters, k_estimate)
}

/// Collapse consecutive equal values (`std::unique`).
fn unique_consecutive(ys: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &v in ys {
        match out.last() {
            Some(&l) if l == v => {}
            _ => out.push(v),
        }
    }
    out
}

/// Indices of `ys` equal to `c`: `(count, sum of indices)`.
fn index_stats(ys: &[f64], c: f64) -> (usize, usize) {
    let mut count = 0usize;
    let mut sum = 0usize;
    for (i, &v) in ys.iter().enumerate() {
        if v == c {
            count += 1;
            sum += i;
        }
    }
    (count, sum)
}

// ------------------------------------------------------------------------------------------
// Public API: tsGates, tsModeEstimate, getTautStringApprox
// ------------------------------------------------------------------------------------------

/// The taut string of sorted `x` (`cpPmden(x).string`): the density estimate on each of the
/// `n - 1` gaps between consecutive observations. Empty for `n < 2`.
///
/// For the pairs the R `getTautStringApprox` returns, see [`taut_string_pairs`].
pub fn taut_string_approx(x: &[f64]) -> Vec<f64> {
    cp_pmden(x).string
}

/// `getTautStringApprox`: the taut string reduced to plotting pairs, `(stringX, stringY)`:
/// for every run of equal string values the data value where the run starts and where it
/// ends, each paired with the run's value.
///
/// QUIRK (reproduced): a run that starts at the last string position is dropped.
pub fn taut_string_pairs(x: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let full = cp_pmden(x).string;
    let mut out_x = Vec::new();
    let mut out_y = Vec::new();
    if full.is_empty() {
        return (out_x, out_y);
    }
    let epchk = full.len() - 1;
    let mut current_sp = 0usize;
    let mut current_sv = full[0];
    let mut account_for_end = true;
    for i in 1..full.len() {
        if full[i] != current_sv {
            if i == epchk {
                account_for_end = false;
            }
            out_x.push(x[current_sp]);
            out_y.push(current_sv);
            out_x.push(x[i]);
            out_y.push(current_sv);
            current_sv = full[i];
            current_sp = i;
        }
        if i == epchk && account_for_end {
            out_x.push(x[current_sp]);
            out_y.push(current_sv);
            out_x.push(x[i]);
            out_y.push(current_sv);
        }
    }
    (out_x, out_y)
}

/// Cut points for sorted `x`. `mode_prior` 0 lets the taut string decide (2 centres when it
/// sees one mode); `k > 0` forces `k` k-medians centres (FAUST passes 2/3/4 from the dip tests
/// when n ≤ 500).
///
/// Returns `[min, cuts..., max]`. Empty for `n < 2`.
pub fn ts_gates(x: &[f64], mode_prior: usize) -> Vec<f64> {
    let n = x.len();
    if n < 2 {
        return Vec::new();
    }
    let ts = cp_pmden(x);
    let num_modes = ts.nmax - ts.nmax / 2;
    let num_centers = if mode_prior == 0 { 2 } else { mode_prior };
    if num_modes <= 1 || mode_prior != 0 {
        return kmed_gates(x, num_centers);
    }
    let ys = &ts.string;
    let yvals = unique_consecutive(ys);
    let mut cut_values = Vec::new();
    if yvals.len() >= 3 {
        for i in 1..yvals.len() - 1 {
            if yvals[i] < yvals[i - 1] && yvals[i] < yvals[i + 1] {
                cut_values.push(yvals[i]);
            }
        }
    }
    let mut gate_indices = Vec::new();
    for &c_cut in &cut_values {
        let (count, index_sum) = index_stats(ys, c_cut);
        if count == 0 {
            continue;
        }
        let cut_index = index_sum / count;
        // reject thresholds that are point-masses arising from local squeezing in large samples
        if count > 5 || n < 5000 {
            gate_indices.push(cut_index);
        }
    }
    if gate_indices.is_empty() {
        return kmed_gates(x, num_centers);
    }
    let mut new_gates = Vec::with_capacity(gate_indices.len() + 2);
    new_gates.push(min_element(x));
    for ci in gate_indices {
        new_gates.push(x[ci.min(n - 1)]);
    }
    new_gates.push(max_element(x));
    new_gates.sort_by(f64::total_cmp);
    new_gates
}

/// Mode locations the taut string estimates (`tsModeEstimate`): `[min, mode, max]`, where the
/// mode is the data value under the first local maximum of the string; the sample median when
/// the string sees more than one mode or none. Empty for `n < 2`.
pub fn ts_mode_estimate(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n < 2 {
        return Vec::new();
    }
    let ts = cp_pmden(x);
    let num_modes = ts.nmax - ts.nmax / 2;
    let sample_median = r_quantile(x, 0.5);
    let finish = |mid: f64| {
        let mut g = vec![min_element(x), mid, max_element(x)];
        g.sort_by(f64::total_cmp);
        g
    };
    if num_modes > 1 {
        // C++: "Taut string disagrees with DIP; returning sample median."
        return finish(sample_median);
    }
    let ys = &ts.string;
    let yvals = unique_consecutive(ys);
    let mut first_cut: Option<f64> = None;
    if yvals.len() >= 3 {
        for i in 1..yvals.len() - 1 {
            if yvals[i] > yvals[i - 1] && yvals[i] > yvals[i + 1] {
                first_cut = Some(yvals[i]);
                break;
            }
        }
    }
    match first_cut {
        // C++: "Warning: tsModeEstimate string found no mode. Returning sample median."
        None => finish(sample_median),
        Some(c) => {
            let (count, index_sum) = index_stats(ys, c);
            if count == 0 {
                return finish(sample_median);
            }
            finish(x[(index_sum / count).min(n - 1)])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Brute-force 1-D k-medians: every contiguous partition of sorted `x` into `k` nonempty
    /// groups, minimum total absolute deviation from the group median (lower median, as the DP).
    fn brute_cost(x: &[f64], k: usize) -> f64 {
        fn cost(seg: &[f64]) -> f64 {
            let m = seg[(seg.len() - 1) / 2];
            seg.iter().map(|v| (v - m).abs()).sum()
        }
        fn rec(x: &[f64], k: usize) -> f64 {
            if k == 1 {
                return cost(x);
            }
            let mut best = f64::INFINITY;
            for cut in 1..=(x.len() - (k - 1)) {
                let c = cost(&x[..cut]) + rec(&x[cut..], k - 1);
                if c < best {
                    best = c;
                }
            }
            best
        }
        rec(x, k)
    }

    fn dp_cost(x: &[f64], cl: &[usize]) -> f64 {
        let k = cl.iter().copied().max().map_or(0, |m| m + 1);
        let mut total = 0.0;
        for c in 0..k {
            let seg: Vec<f64> = x
                .iter()
                .zip(cl)
                .filter(|&(_, &g)| g == c)
                .map(|(v, _)| *v)
                .collect();
            let m = seg[(seg.len() - 1) / 2];
            total += seg.iter().map(|v| (v - m).abs()).sum::<f64>();
        }
        total
    }

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / ((1u64 << 53) as f64)
    }

    #[test]
    fn k_med_dp_matches_brute_force() {
        let mut seed = 42u64;
        for n in 2..=14usize {
            for k in 1..=4usize.min(n) {
                for _ in 0..6 {
                    let mut x: Vec<f64> = (0..n)
                        .map(|_| (lcg(&mut seed) * 20.0).round() / 2.0)
                        .collect();
                    x.sort_by(f64::total_cmp);
                    let cl = k_med_dp(&x, k);
                    assert_eq!(cl.len(), n);
                    // contiguous, nondecreasing, exactly k groups
                    assert!(cl.windows(2).all(|w| w[1] == w[0] || w[1] == w[0] + 1));
                    assert_eq!(cl[0], 0);
                    assert_eq!(cl[n - 1], k - 1);
                    let got = dp_cost(&x, &cl);
                    let want = brute_cost(&x, k);
                    assert!(
                        (got - want).abs() <= 1e-9,
                        "n={n} k={k} x={x:?} got={got} want={want}"
                    );
                }
            }
        }
    }

    #[test]
    fn k_med_dp_small_cases() {
        assert_eq!(k_med_dp(&[], 2), Vec::<usize>::new());
        assert_eq!(k_med_dp(&[1.0], 2), vec![0]);
        assert_eq!(k_med_dp(&[1.0, 2.0], 0), Vec::<usize>::new());
        assert_eq!(
            k_med_dp(&[0.0, 0.1, 0.2, 10.0, 10.1], 2),
            vec![0, 0, 0, 1, 1]
        );
        assert_eq!(
            k_med_dp(&[0.0, 0.1, 5.0, 5.1, 10.0, 10.1], 3),
            vec![0, 0, 1, 1, 2, 2]
        );
    }

    #[test]
    fn find_kmed_gates_shape() {
        let x = [0.0, 0.1, 0.2, 10.0, 10.1];
        let g = find_kmed_gates(&x, &[0, 0, 0, 1, 1], 1);
        assert_eq!(g, vec![0.0, 0.2, 10.1]);
    }

    #[test]
    fn tiny_inputs_do_not_panic() {
        assert!(ts_gates(&[], 0).is_empty());
        assert!(ts_gates(&[1.0], 0).is_empty());
        assert!(ts_mode_estimate(&[]).is_empty());
        assert!(taut_string_approx(&[]).is_empty());
        assert_eq!(taut_string_pairs(&[1.0, 2.0]).0.len(), 0);
        let g = ts_gates(&[1.0, 2.0], 0);
        assert_eq!(g.len(), 3);
        let _ = ts_gates(&[1.0, 2.0, 3.0], 4);
        let _ = ts_mode_estimate(&[1.0, 2.0, 3.0]);
        let _ = ts_gates(&[1.0, 1.0, 1.0, 1.0], 0);
    }

    #[test]
    fn r_quantile_type7() {
        assert_eq!(r_quantile(&[1.0, 2.0, 3.0, 4.0], 0.5), 2.5);
        assert_eq!(r_quantile(&[3.0, 1.0, 2.0], 0.5), 2.0);
        assert!((median_absolute_deviation(&[1.0, 2.0, 3.0, 4.0, 5.0]) - 1.4826).abs() < 1e-12);
    }

    #[test]
    fn approx_rule2() {
        let x = [0.0, 1.0, 2.0];
        let y = [0.0, 10.0, 20.0];
        assert_eq!(
            approx_linear(&x, &y, &[-1.0, 0.5, 1.0, 3.0]),
            vec![0.0, 5.0, 10.0, 20.0]
        );
    }

    /// `cargo test -p faust --release --lib perf_one_million -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn perf_one_million_points() {
        let n = 1_000_000usize;
        let mut seed = 7u64;
        let mut x: Vec<f64> = (0..n)
            .map(|i| {
                let u = lcg(&mut seed).max(1e-12);
                let u2 = lcg(&mut seed);
                let g = (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
                g * 0.45 + if i % 3 == 0 { 0.0 } else { 2.5 }
            })
            .collect();
        x.sort_by(f64::total_cmp);
        let t = std::time::Instant::now();
        let g0 = ts_gates(&x, 0);
        let t0 = t.elapsed();
        let t = std::time::Instant::now();
        let g4 = ts_gates(&x, 4);
        let t4 = t.elapsed();
        let t = std::time::Instant::now();
        let m = ts_mode_estimate(&x);
        let tm = t.elapsed();
        eprintln!(
            "n={n}: tsGates(0)={g0:?} in {t0:?}; tsGates(4)={g4:?} in {t4:?}; mode={m:?} in {tm:?}"
        );
        assert!(g0.len() >= 3 && g4.len() == 5 && m.len() == 3);
    }
}
