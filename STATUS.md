# Status

## Measured (2026-10-12, 16-core laptop, release build)

Synthetic data shaped like `dev/profile_faust.R` (markers with 1–3 modes). R = FAUST 0.5.5 in
`tercen/faust_operator:0.1.2`, single thread; Rust = this crate.

| input | R forest | R scamp | R UMAP | R peak | Rust (1 thread) | Rust (1 thread/sample) | Rust peak |
|---|---|---|---|---|---|---|---|
| 2 × 50k cells × 10 markers | 364 s | 470 s | 87 s | 1.0 GB | 18.8 s (forest+scamp) | — | 62 MB |
| 2 × 200k cells × 10 markers | 729 s | 180 s | 243 s | 2.1 GB | 81 s | 40 s | 240 / 312 MB |

Forest alone (`examples/bench_forest.rs`): 50k × 10 in 13 s, 200k × 10 in 23 s (R: ~182 s and
~365 s per sample).

## Memory model

`memory_model.json`: 0.00015 MB per crosstab value + 1.5 × 200 MB. Points: 1 M values → 62 MB
(1 thread), 4 M values → 312 MB (2 samples in parallel); the UMAP on all cells adds ~70 B per
cell × selected markers (umap_operator's measurement). Refit from `stats_d_actual_ram_peak`
once there are production runs.

## Platform test

`tests/test.json`: the R operator's own `test_1_in.csv` (500 cells × 5 markers), goldens
exported from a Studio run at `threads = 1` (two runs byte-identical), compared with R2 ≥ 0.99.
The labels and all 2,500 gate values of that run are identical to the R operator's golden
(`r-legacy`), which came from its 4-thread run.

## Parity

Every stage is checked against values dumped from R FAUST 0.5.5 / scamp 0.5.0
(`crates/faust/tests/parity_*.rs`, 60+ tests): dip tests, taut-string gates, libstdc++
RNG/noise, annotation forests (unrestricted and with channel bounds), depth scores, marker
selection, threshold reconciliation, scamp labels, phenotype selection, final labels — all
bit-identical or string-identical on test_1 and two 20k-cell synthetic samples.

## Not done

- scamp's subsampled gates (`densitySubSampleThreshold`, nodes of ≥ 1e6 cells): `unimplemented!`;
  samples above 1 M cells per node are refused.
- `experimentalUnit`, `imputationHierarchy`, `supervisedList`, `channelBounds` inputs of
  `faust()` (the operator never set them).

## 1.1.3 (2026-10-14) — memory model refit on Tudor's prod runs

Measured peaks of 1.1.2 on prod (`stats_d_actual_ram_peak`): 500 cells 19 MB, 60k cells 170 MB,
75.5k cells 212 MB, 200k cells 518 MB (all 6 markers). Growth is ~2.7 MB per 1000 cells, i.e.
~0.00045 MB per crosstab value, dominated by the UMAP neighbour graph. The old coefficient
0.00015 booked 503 MB for the 200k run, below its 518 MB peak, so a run without a manual RAM
override would have been OOM-killed. New model: 0.0005 MB/value + 200 MB offset (×1.5 by the
server) → 301 MB (500 cells), 900 MB (200k cells), ~2 GB (584k cells, Tudor's largest R run).
No code change.

## 1.1.2 (2026-10-14) — one selected marker crashed the embedding

First prod run of 1.1.1 (Tudor Dimofte's workflow, 60k cells, 6 markers, no `sample_factor`
so one sample) selected a single marker (PAC1). The embedding matrix then has one column and
umaprs 0.2.0's PCA initialisation (used above 2000 points) indexed a second eigenvector that
does not exist: `index out of bounds: the len is 1 but the index is 1`, exit code 134.
Fixed in umaprs 0.2.1 (tercen/umaprs#3): components beyond the data's dimension start from
seeded unit-normal noise. Regression test `embedding_tests::umap_on_one_dimension_does_not_panic`.
No other change.

## 1.1.0 — threshold histograms

The histogram table is back: one PNG per selected marker (`hist_<marker>.png`, 800 × 600,
FAUST's recipe: ≤ 1000 cells per sample pooled, 1 %/99 % trim, FD/`pretty` bins, every
sample's thresholds as red dashed lines), as a second result table keyed by `.ri` with
`filename`, `mimetype`, `checksum`, `.content` — the UI lists them under Files. Drawn with
plotters into an RGB buffer, DejaVu Sans embedded (`fonts/`, Bitstream Vera licence), encoded
with `png`; the static image stays a few MB. `joinOperators` is left empty so the server joins
both tables itself.

Release note: tag `1.1.0` was cut with `operator.json` still pinning the 1.0.0 image, so the
release workflow's install check ran the old binary against the four-table golden and failed
(`task.test.operator.bad.nRelations`); the `1.1.0` image exists on ghcr but there is no release.
`1.1.1` is the same code with the container pin bumped — a release commit must bump
`operator.json` `container` and `Cargo.toml` `version` to the tag before tagging.
