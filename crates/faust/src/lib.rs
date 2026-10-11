//! FAUST — Full Annotation Using Shape-constrained Trees — in Rust.
//!
//! A port of RGLab/FAUST 0.5.5 (R + C++) and RGLab/scamp 0.5.0, the packages the R
//! `tercen/faust_operator` runs. The reference is the **single-thread** R path
//! (`threadNum = 1`); every module here is checked against values dumped from that R, with
//! `%.17g`, in `tests/parity/`.
//!
//! Pipeline (FAUST's two phases, then the embedding):
//!
//! 1. **Annotation thresholds** — per sample: add noise ([`noise`]), grow the exhaustive
//!    depth-3 annotation forest ([`forest`]: one-dimensional gates from the dip test
//!    ([`dip`]) and the taut-string density ([`tsgates`])), score markers by depth
//!    ([`select`]), keep the markers whose depth score is high in enough samples, and
//!    standardise the number and location of the thresholds across samples ([`select`]).
//! 2. **Phenotype discovery** — per sample: scamp clustering ([`scamp`]) on the selected
//!    markers with the standardised thresholds as annotation boundaries; every cell gets a
//!    label like `CD3+CD4-CD8+`; phenotypes are kept by how often they occur across samples
//!    ([`select::phenotype_selection`]).
//! 3. **Annotation embedding** — ([`embed`]) the per-phenotype landmark layout the R
//!    `makeAnnotationEmbedding` feeds to UMAP (the UMAP itself lives in the operator).
//!
//! Module interfaces are fixed first so the ports can be written and checked independently.
pub mod dip;
pub mod embed;
pub mod forest;
pub mod noise;
pub mod rng;
pub mod scamp;
pub mod select;
pub mod tsgates;

/// A sample's expression matrix, column-major: `cols[j]` is marker `j` over all cells.
#[derive(Debug, Clone)]
pub struct Matrix {
    pub markers: Vec<String>,
    pub cols: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn n_cells(&self) -> usize {
        self.cols.first().map_or(0, Vec::len)
    }
    pub fn n_markers(&self) -> usize {
        self.cols.len()
    }
}

#[derive(Debug, thiserror::Error)]
pub enum FaustError {
    #[error("{0}")]
    Invalid(String),
}

pub type Result<T> = std::result::Result<T, FaustError>;
