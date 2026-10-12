# Overnight report — FAUST operator in Rust (2026-10-11 → 12)

Goal (Faris, 00:20): "create and test a rust faust operator ready in morning".

## Delivered

- **`tercen/faust_operator` PR #1** (`feature/rust` → `main`): the operator rewritten in Rust,
  version 1.0.0. The crate `crates/faust` is the port of RGLab/FAUST 0.5.5 + scamp 0.5.0; the
  shell is the umap_operator's Tercen plumbing. Same projection, same properties, same output
  columns as the R operator.
- **Parity with R, stage by stage, on fixtures dumped from the R packages** (60+ tests in
  `crates/faust/tests/parity_*.rs`, `cargo test` green):
  - dip tests, taut-string gates, k-medians, libstdc++ RNG and FAUST's noise: bit-identical;
  - annotation forests (with and without channel bounds), depth scores, marker selection,
    threshold reconciliation: bit-identical;
  - scamp labels, phenotype selection, final labels: identical on test_1 and two 20k-cell
    synthetic samples.
- **Through the server** (Studio, project `faust_rust_dev`): on `tests/test_1_in.csv` the
  Rust operator's 500 labels and 2,500 gate values equal the R operator's own golden; two runs
  at `threads = 1` export byte-identical. That export is the new platform test
  (`tests/test.json`).
- **Speed / memory** (STATUS.md), R single thread vs Rust:

  | input | R | Rust | peak R → Rust |
  |---|---|---|---|
  | 2 × 50k × 10 | 834 s | 18.8 s | 1.0 GB → 62 MB |
  | 2 × 200k × 10 | 909 s | 81 s (1 thread), 40 s (one thread per sample) | 2.1 GB → 312 MB |

  The image: 7.1 GB (R) → 9.9 MB, 3.0 MB compressed, static binary on `scratch`.
- **Production path exercised without publishing**: a platform-style `RunComputationTask` with a
  pre-allocated result file, run through the built image as uid 1000 under the memory the model
  books for test_1 (300 MB): exit 0, result uploaded into the pre-allocated file, peak 17 MB.

## What was learned about FAUST on the way (all reproduced, all in doc comments)

- scamp's candidate search is **iterated**: `findCandidateClusters` reruns the random search
  with `seed + 1` until `50 · markers` candidates exist; the seed keeps advancing through the
  residual rounds. The first port stopped after one search and 7 of 500 labels differed.
- The residual rounds score candidates on the **wrong rows** (the full matrix indexed by
  subset positions). Reproduced; it is what every R result contains.
- FAUST sees the expression **rounded to float32** (the GatingSet), and `channelBounds = ""`
  still **restricts** ~2 % of the cells per marker (empirical 1 %/99 % per sample, then
  5 %/95 % across samples); forest and scamp need the restriction matrix.
- R's `scale()` and `sum()` accumulate in 80-bit long double; the port emulates it.
- Phase 2 in R at `threadNum > 1` (what the R operator used) is not reproducible even in R.
  "Same as R" means the single-thread R run; the R operator's golden (4 threads) still matched
  on every cell here.

## Not done / decisions for you

1. **Released.** CI green on the branch, tag `1.0.0` cut at 00:06Z: the Release workflow
   pushed `ghcr.io/tercen/faust_operator:1.0.0` (package is public, anonymous pull works), the
   install check with the platform test passed on the org test server, and the GitHub release is
   published. **PR #1 is not merged** — that is your call; `operator.json` on the branch pins
   1.0.0, so installing the repo at the tag works now.
2. **The R operator's third table** (threshold histogram PNGs, joined on nothing) is not
   produced: declaring that empty-key join made the server drop the main table. If it matters,
   it should be a separate plot operator.
3. **Nodes of ≥ 1 M cells** (FAUST's `densitySubSampleThreshold`) use a subsampled gate in the
   C++ that is not ported; such samples are refused with a message.
4. After merge: move the R code to `r-legacy` (+ tag `r-legacy-0.2.0`), swap the library pin.
5. The Studio project `faust_rust_dev` can be deleted.
