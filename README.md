# FAUST operator

##### Description

FAUST — Full Annotation Using Shape-constrained Trees (Greene et al., *Patterns* 2021) —
discovers and annotates cell phenotypes per sample from a cytometry crosstab, then
standardises the annotation thresholds across samples so phenotypes are comparable.

Version 1.x is a **Rust implementation** (`crates/faust`): a port of the R packages
`RGLab/FAUST` 0.5.5 and `RGLab/scamp` 0.5.0 that the 0.x R operator ran. On the same seed it
returns the same marker selection, thresholds, per-cell phenotype labels and gates as the R
operator at `threadNum = 1`; the UMAP coordinates use the same landmark layout but come from
`umaprs`, so they differ from uwot's coordinate by coordinate. Since 1.2.0 (umaprs 0.4.0, with
uwot's spectral start) the maps match uwot's in structure: see umap_operator 2.2.0's validation.
The R code is kept on the `r-legacy` branch.

##### Usage

Input projection|.
---|---
`row`   | the variables (markers, channels)
`col`   | the observations (cells); a `filename` factor (or `sample_factor`) names the sample
`y-axis`| the measurement

Input parameters|.
---|---
`selectionQuantile`   | a marker is kept when this quantile of its depth score across samples is at least `depthScoreThreshold` (default 0.5)
`depthScoreThreshold` | depth score a marker must reach (default 0.01)
`seed`                | random seed; the same seed gives the same labels (default 123)
`sample_factor`       | column factor naming the sample; empty uses `filename` if projected, else one sample
`threads`             | samples processed in parallel; 0 = every core the task has (the result does not depend on it)


Output relations|.
---|---
`umap_1`, `umap_2`   | annotation embedding coordinates per cell
`label`              | FAUST phenotype, e.g. `CD3+CD4-CD8+`; cells of phenotypes below the occurrence threshold get `0_0_0`
`<marker>_gate`      | 1 for `+`, 0 for `-` on each selected marker (NA for markers with more than two levels)


##### Details

Per sample, FAUST adds a small amount of noise, grows the exhaustive depth-3 annotation forest
(one-dimensional gates from the dip test and the taut-string density), scores markers by how
consistently and how shallowly they separate the cells (depth score), keeps the markers whose
depth score is high in enough samples, and standardises the number and location of their
thresholds across samples. Then scamp clusters each sample on the selected markers and every
cell is labelled against the standardised thresholds. Phenotypes are kept by how often they
occur across samples (elbow rule above 4 samples, all phenotypes otherwise).

Parity with R is checked stage by stage in `crates/faust/tests/parity_*.rs` on fixtures dumped
from the R packages (`dev/make_fixtures.R`, run inside the 0.x operator's image). Differences
from the R operator, all deliberate: the embedding's jitters are drawn from a seeded generator
(R's were unseeded), samples run in parallel, and the threshold histograms are one PNG per
selected marker on the marker's row (`.ri`), all of its gates on the one plot (R wrote one
file per gate number, joined on nothing). Since 1.1.0.

##### References

- Greene E, Finak G, et al. New interpretable machine-learning method for single-cell data
  reveals correlates of clinical response to cancer immunotherapy. *Patterns* 2021.
- [RGLab/FAUST](https://github.com/RGLab/FAUST), [RGLab/scamp](https://github.com/RGLab/scamp)

##### See Also

[flowsom_operator](https://github.com/tercen/flowsom_operator),
[umap_operator](https://github.com/tercen/umap_operator)
