//! scamp: selective clustering annotated using modes of projections (RGLab/scamp 0.5.0).
//!
//! Port of `src/scamp.cpp`, `src/selectCandidateClusters.cpp`, `src/absSampleLmoments.cpp`,
//! `src/annotateCluster.cpp`, `src/rQuantile.cpp`, `src/cppNoisyScamp.cpp` as FAUST calls them
//! (`numberIterations = 1`, `useAnnForest = TRUE`, `randomCandidateSearch = TRUE`, no restricted
//! values, no time budget). Implemented by the operator author on top of [`crate::forest`].
