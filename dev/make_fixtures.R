## Reference fixtures for the Rust port, from R FAUST 0.5.5 / scamp 0.5.0 (the image the
## operator runs). Everything numeric is written with %.17g so a Rust test can compare exactly.
## Datasets: the operator's own tests/test_1 (500 cells x 5 markers, 1 sample) and a synthetic
## 2-sample set (20k cells x 8 markers) with bimodal/trimodal markers.
suppressPackageStartupMessages({library(faust); library(scamp); library(flowCore); library(flowWorkspace)})
out <- commandArgs(trailingOnly = TRUE)[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
g17 <- function(x) sprintf("%.17g", x)
wvec <- function(x, name) writeLines(if (is.numeric(x)) g17(x) else as.character(x), file.path(out, name))
wmat <- function(m, name) {
  con <- file(file.path(out, name), "w")
  writeLines(paste(colnames(m), collapse = ","), con)
  for (i in seq_len(nrow(m))) writeLines(paste(if (is.numeric(m)) g17(m[i, ]) else m[i, ], collapse = ","), con)
  close(con)
}
## generic dump of an R object tree into files (lists -> directories)
dumpObj <- function(x, path) {
  if (is.list(x) && !is.data.frame(x)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    nms <- names(x); if (is.null(nms)) nms <- sprintf("i%03d", seq_along(x))
    nms[nms == ""] <- sprintf("i%03d", which(nms == ""))
    for (k in seq_along(x)) dumpObj(x[[k]], file.path(path, gsub("[^A-Za-z0-9_.+-]", "_", nms[k])))
  } else if (is.matrix(x) || is.data.frame(x)) {
    m <- as.matrix(x); if (is.null(colnames(m))) colnames(m) <- sprintf("c%d", seq_len(ncol(m)))
    if (is.null(rownames(m))) rownames(m) <- sprintf("r%d", seq_len(nrow(m)))
    con <- file(paste0(path, ".csv"), "w")
    writeLines(paste(c("rowname", colnames(m)), collapse = ","), con)
    for (i in seq_len(nrow(m))) writeLines(paste(c(rownames(m)[i], if (is.numeric(m)) g17(m[i, ]) else as.character(m[i, ])), collapse = ","), con)
    close(con)
  } else if (is.null(x)) {
    writeLines("NULL", paste0(path, ".txt"))
  } else {
    con <- file(paste0(path, ".txt"), "w")
    if (!is.null(names(x))) writeLines(paste(names(x), collapse = ","), con)
    writeLines(if (is.numeric(x)) g17(x) else as.character(x), con)
    close(con)
  }
}

## ---------- datasets ----------
t1 <- read.csv("/operator/tests/test_1_in.csv")
m1 <- matrix(t1$y_values, ncol = 5, dimnames = list(NULL, paste0("V", 1:5)))
# column j of the matrix = rows with gs0.variable == Vj, in rowId order
for (j in 1:5) m1[, j] <- t1$y_values[t1$gs0.variable == paste0("V", j)][order(t1$rowId[t1$gs0.variable == paste0("V", j)])]
wmat(m1, "test1_matrix.csv")

set.seed(7)
mkSample <- function(n, d, shift) {
  m <- matrix(0, n, d)
  for (j in seq_len(d)) {
    modes <- 1 + (j %% 3)          # 1, 2 or 3 modes
    comp <- sample.int(modes, n, replace = TRUE, prob = seq_len(modes) / sum(seq_len(modes)))
    m[, j] <- rnorm(n, mean = (comp - 1) * 2.5 + shift * 0.1, sd = 0.45)
  }
  colnames(m) <- paste0("M", seq_len(d)); m
}
s1 <- mkSample(20000, 8, 0); s2 <- mkSample(20000, 8, 1)
wmat(s1, "syn_s1.csv"); wmat(s2, "syn_s2.csv")

## ---------- primitives ----------
prim <- function(m, tag) {
  for (j in seq_len(ncol(m))) {
    x <- m[, j]; sx <- sort(x)
    wvec(faust:::singleDip(sx), sprintf("%s_%s_singleDip.txt", tag, colnames(m)[j]))
    wvec(faust:::tsGates(sx, 0L), sprintf("%s_%s_tsGates0.txt", tag, colnames(m)[j]))
    wvec(faust:::tsGates(sx, 2L), sprintf("%s_%s_tsGates2.txt", tag, colnames(m)[j]))
    wvec(faust:::tsGates(sx, 3L), sprintf("%s_%s_tsGates3.txt", tag, colnames(m)[j]))
    wvec(faust:::getTautStringApprox(sx), sprintf("%s_%s_tautString.txt", tag, colnames(m)[j]))
    if (length(sx) <= 500) {
      wvec(faust:::doubleDip(sx), sprintf("%s_%s_doubleDip.txt", tag, colnames(m)[j]))
      if (length(sx) <= 400) wvec(faust:::tripleDip(sx[1:400]), sprintf("%s_%s_tripleDip.txt", tag, colnames(m)[j]))
    }
    wvec(faust:::addNoiseToDataVector(x, 4, 123), sprintf("%s_%s_noise_s123.txt", tag, colnames(m)[j]))
    wvec(faust:::addNoiseToDataVector(x, 4, 124), sprintf("%s_%s_noise_s124.txt", tag, colnames(m)[j]))
  }
  # subsets of various sizes for the small-n dip branches
  sx <- sort(m[, 1])
  for (n in c(30, 100, 300, 400, 500)) {
    if (n <= length(sx)) {
      y <- sx[seq(1, length(sx), length.out = n)]
      wvec(y, sprintf("%s_sub%d_x.txt", tag, n))
      wvec(faust:::singleDip(y), sprintf("%s_sub%d_singleDip.txt", tag, n))
      wvec(faust:::doubleDip(y), sprintf("%s_sub%d_doubleDip.txt", tag, n))
      if (n <= 400) wvec(faust:::tripleDip(y), sprintf("%s_sub%d_tripleDip.txt", tag, n))
      wvec(faust:::tsGates(y, 0L), sprintf("%s_sub%d_tsGates0.txt", tag, n))
    }
  }
}
prim(m1, "test1"); prim(s1, "syn1")
cat("primitives done\n")

## ---------- annotation forest (phase 1 core), operator parameters ----------
forest <- function(m, tag) {
  af <- faust:::growAnnotationForest(dataSet = m, numberIterations = 1, pValueThreshold = 0.25,
    minimumClusterSize = 25, randomCandidateSearch = FALSE, maximumSearchDepth = 2L,
    numberOfThreads = 1, maximumGatingNum = 1e10, anyValueRestricted = FALSE,
    resValMatrix = matrix(0, nrow = 2, ncol = 2), cutPointUpperBound = 2, getDebugInfo = FALSE,
    randomSeed = 123, subSamplingThreshold = 1e6, subSampleSize = 1e6, subSampleIter = 1,
    recordCounts = FALSE, recordIndices = FALSE)
  dumpObj(af$gateData, file.path(out, sprintf("%s_forest_gateData", tag)))
  ePop <- apply(m, 2, function(x) 0); names(ePop) <- colnames(m)
  # effective sizes as .growForestForExpUnit computes them (no restricted values -> nrow)
  eSizes <- rep(nrow(m), ncol(m)); names(eSizes) <- colnames(m)
  pAnnF <- faust:::.parseAnnotationForest(af[["gateData"]], eSizes)
  dumpObj(pAnnF, file.path(out, sprintf("%s_forest_parsed", tag)))
}
forest(m1, "test1"); forest(s1, "syn1"); forest(s2, "syn2")
cat("forests done\n")

## ---------- full faust() runs: labels and thresholds ----------
runFaust <- function(mats, tag) {
  ffs <- lapply(mats, function(m) flowFrame(m)); fs <- as(ffs, "flowSet"); gs <- GatingSet(fs)
  proj <- file.path(tempdir(), paste0("FAUST_", tag)); unlink(proj, recursive = TRUE); dir.create(proj, recursive = TRUE)
  faust::faust(gatingSet = gs, startingCellPop = "root", projectPath = proj, annotationsApproved = TRUE,
    threadNum = 1, selectionQuantile = 0.5, seedValue = 123, depthScoreThreshold = 0.01,
    densitySubSampleThreshold = 1e6, densitySubSampleSize = 1e6, drawAnnotationHistograms = FALSE,
    plottingDevice = "png")
  # dump every rds under faustData, keeping the tree
  for (f in list.files(file.path(proj, "faustData"), pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)) {
    rel <- sub(paste0("^", proj, "/faustData/"), "", f); rel <- sub("\\.rds$", "", rel)
    dumpObj(readRDS(f), file.path(out, sprintf("%s_faust", tag), rel))
  }
  for (f in list.files(file.path(proj, "faustData"), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)) {
    rel <- sub(paste0("^", proj, "/faustData/"), "", f)
    dir.create(dirname(file.path(out, sprintf("%s_faust", tag), rel)), recursive = TRUE, showWarnings = FALSE)
    file.copy(f, file.path(out, sprintf("%s_faust", tag), rel), overwrite = TRUE)
  }
  snVec <- list.files(file.path(proj, "faustData", "sampleData"))
  ae <- faust::makeAnnotationEmbedding(projectPath = proj, sampleNameVec = snVec)
  write.csv(ae, file.path(out, sprintf("%s_annotationEmbedding.csv", tag)), row.names = FALSE)
  cat(sprintf("%s: %d cells, %d phenotypes\n", tag, nrow(ae), length(unique(ae$faustLabels))))
}
runFaust(list(sample = m1), "test1")
runFaust(list(s1 = s1, s2 = s2), "syn")
cat("fixtures done\n")
