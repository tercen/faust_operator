//! The annotation embedding's input matrix (`R/makeAnnotationEmbedding.R`): per-marker
//! winsorising at the 1st/99th percentile, standardising inside each phenotype, and shifting
//! each phenotype to its own landmark on an equispaced grid. UMAP itself (uwot in R,
//! `umaprs` here) runs in the operator on the matrix this module returns.
