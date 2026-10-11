# faust_operator — maintenance guide

Rust port of the R FAUST operator. Read this before touching anything; the numbers below are
measured, not designed.

## Architecture

| file | role |
|---|---|
| `src/main.rs`, `src/bin/dev.rs` | production (`--taskId`) and dev (`--workflowId/--stepId`) entry points |
| `src/lib.rs` | `execute()`: read the crosstab, split by sample, run FAUST, embed, write, upload |
| `src/input.rs`, `src/context.rs`, `src/upload.rs`, `src/tson.rs`, `src/progress.rs`, `src/pagecache.rs` | Tercen plumbing, copied from `umap_operator` (keep them in step) |
| `src/props.rs` | properties; the test asserts every one is in `operator.json` with the same default |
| `src/output.rs` | the `OperatorResult` TSON: per-cell table (`.ci`) + optional PNG table |
| `src/plots.rs` | the threshold histograms (R's `plotMarkerHistograms`) with `plotters` |
| `crates/faust/` | the algorithm, no Tercen types: see its module docs |

`crates/faust` modules, in pipeline order: `noise` (libstdc++-exact `rng`) → `forest` (search
over `dip` + `tsgates`) → `select` (depth scores, marker selection, threshold reconciliation,
gating, phenotype selection) → `scamp` → `embed`; `pipeline` composes them; `rrng` is R's
`set.seed`/`runif` for the one R draw FAUST makes.

## Crosstab contract

Rows = markers, columns = cells (+ `filename` or `sample_factor`), y = value. Per-cell output
keyed by `.ci`: `umap_1`, `umap_2` (double), `label` (string), `<marker>_gate` (double, 0/1/NA),
plus the PNG table joined on nothing. The matrix is gathered (cells × markers, f64): FAUST needs
every sample's whole matrix, so a stream makes no sense; `memory_model.json` declares the cost.

## Parity

Reference: R FAUST 0.5.5 + scamp 0.5.0 in `tercen/faust_operator:0.1.2` (R 4.0.4, g++ 9.4.0),
`threadNum = 1`. Fixtures: `dev/make_fixtures.R` run inside that image (`docker run --rm -v
$PWD/dev:/fix -v $PWD/tests:/operator/tests:ro --entrypoint Rscript tercen/faust_operator:0.1.2
/fix/make_fixtures.R /fix/out`), every number written with `%.17g`. Tests
`crates/faust/tests/parity_*.rs`:

| stage | fixture | claim |
|---|---|---|
| dip statistic, single/double/triple dip p-values | `*_singleDip.txt` … (46) | bit |
| taut string, `tsGates`, k-medians | `*_tsGates*.txt`, `*_tautString.txt` | bit |
| mt19937 + distributions, noise | `*_noise_s12[34].txt`, `ties_*` | bit |
| `set.seed(123); ceiling(runif(1,10,1e8))` | 28757760 | exact |
| annotation forest `gateData` | `*_forest_gateData/` (3 samples) | bit, same order |
| depth scores, selection, reconciliation, labels | `*_faust/` dumps, `*_annotationEmbedding.csv` | exact |

Three things in the C++ were reproduced rather than fixed, each marked in a doc comment: the
residual scamp search scores candidates on the wrong rows (full matrix indexed by subset
positions); `cppRnorm` builds a fresh distribution per draw; the gating keeps the node minimum
in the first child. Deliberate deviations: the embedding jitters and the histogram subsample
use a seeded generator (R's were unseeded), samples run in parallel, no `std::random_device`
path (seed 0 is refused).

**What is not reproducible even in R:** phase 2 with `threadNum > 1` (the operator ran 4),
because scamp merges thread results under a mutex and stops at a shared cap. Phase 1 is a
multiset of gates and is thread-order independent. So "same as the R operator" means the
R operator at one thread; against its 4-thread runs expect the same thresholds and a small
fraction of cells labelled differently.

## Seed / RNG

`seed` is FAUST's `seedValue`. The forest seed is `ceiling(runif(1, 10, 1e8))` after
`set.seed(seed)` (R's Mersenne-Twister, `rrng.rs`); the C++ engines are `std::mt19937` seeded
with it (+ per-column increments for the noise, +1 for the search), reproduced from GCC 9's
libstdc++ (`rng.rs`; uniform_int differs on GCC ≥ 11 and the tests pin the 9.4 numbers).

## Performance

Reference R (single thread, 10 markers, in the 0.1.2 image, this machine): 2 × 50k cells —
forest 364 s, scamp 470 s, UMAP 87 s, peak 1.0 GB, 761 phenotypes; 2 × 200k — 729 s / 180 s /
243 s, 2.1 GB, 984 phenotypes. The forest is linear in cells; scamp is capped at 50·d
candidates per sample.

Rust: see `STATUS.md` for the measured numbers and the memory model's two points. Where the
time went in the C++ and what this port does instead is in `crates/faust/src/forest.rs`.

## Tests

`cargo test` (unit + parity, ~1 min), `tests/test.json` (the platform test: goldens exported
from a Studio run of this operator on `tests/test_1_in.csv`, the R operator's own input).
Reproducibility is proved through the server: run the step twice and `cmp` the exports.

## Release

`operator.json` pins `ghcr.io/tercen/faust_operator:X.Y.Z`; bump the pin, CI green, tag, push,
in one go. `release.yml` builds the image (cache shared with `ci.yml`), pushes it, and runs the
install check on the org test server. The R operator lives on branch `r-legacy` (tag
`r-legacy-0.2.0`); its Docker Hub image `tercen/faust_operator:0.1.2` stays as the reference.
