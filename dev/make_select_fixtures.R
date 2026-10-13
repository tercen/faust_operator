## Extra R references for the Rust port of select.rs: synthetic reconciliation scenarios
## (imputation hierarchy, outliers, up/down conversion, preference + force lists), the
## conversion helpers on their own, and the GCM elbow. All numeric output is %.17g.
suppressPackageStartupMessages(library(faust))
out <- commandArgs(trailingOnly = TRUE)[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
g17 <- function(x) sprintf("%.17g", x)
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
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("NULL", paste0(path, ".txt"))
  } else {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    con <- file(paste0(path, ".txt"), "w")
    if (!is.null(names(x))) writeLines(paste(names(x), collapse = ","), con)
    writeLines(if (is.numeric(x)) g17(x) else as.character(x), con)
    close(con)
  }
}

## ---------- synthetic reconciliation ----------
units <- c("u2", "u1", "u3", "u4", "u6", "u5", "u7")   # analysisMap order (deliberately unsorted)
impH <- c(u1 = "A", u2 = "A", u3 = "A", u4 = "A", u5 = "B", u6 = "B", u7 = "C")
gates <- list(
  X = list(u1 = 1.0, u2 = 1.1, u3 = c(0.9, 2.5), u4 = 1.05, u5 = 1.2, u6 = NA, u7 = c(0.8, 1.3, 2.2)),
  Y = list(u1 = c(1, 3), u2 = c(1.1, 3.2), u3 = 1.05, u4 = c(0.95, 9.0), u5 = c(1.2, 3.3), u6 = NA, u7 = 3.0),
  Z = list(u1 = NA, u2 = NA, u3 = NA, u4 = NA, u5 = 2.0, u6 = NA, u7 = NA),
  W = list(u1 = 2.0, u2 = 2.2, u3 = c(1.9, 4.0), u4 = 2.1, u5 = c(2.05, 4.2), u6 = NA, u7 = c(1.8, 3.9))
)
runRecon <- function(tag, preferenceList, forceList) {
  proj <- file.path(tempdir(), paste0("recon_", tag)); unlink(proj, recursive = TRUE)
  dir.create(file.path(proj, "faustData", "metaData"), recursive = TRUE)
  dir.create(file.path(proj, "faustData", "expUnitData"), recursive = TRUE)
  am <- data.frame(sampleName = units, experimentalUnit = units, impH = impH[units], stringsAsFactors = FALSE)
  saveRDS(am, file.path(proj, "faustData", "metaData", "analysisMap.rds"))
  saveRDS("root", file.path(proj, "faustData", "metaData", "sanitizedCellPopStr.rds"))
  saveRDS(names(gates), file.path(proj, "faustData", "metaData", "initSelC.rds"))
  saveRDS(forceList, file.path(proj, "faustData", "metaData", "forceList.rds"))
  saveRDS(preferenceList, file.path(proj, "faustData", "metaData", "preferenceList.rds"))
  for (u in units) {
    dir.create(file.path(proj, "faustData", "expUnitData", u))
    pAnnF <- lapply(names(gates), function(ch) list(channelScore = 1, gates = gates[[ch]][[u]]))
    names(pAnnF) <- names(gates)
    saveRDS(pAnnF, file.path(proj, "faustData", "expUnitData", u, "root_pAnnF.rds"))
  }
  faust:::.reconcileAnnotationBoundaries(projectPath = proj, debugFlag = FALSE)
  dumpObj(readRDS(file.path(proj, "faustData", "gateData", "root_resListPrep.rds")), file.path(out, tag, "resListPrep"))
  dumpObj(readRDS(file.path(proj, "faustData", "gateData", "root_selectedChannels.rds")), file.path(out, tag, "selectedChannels"))
  dumpObj(readRDS(file.path(proj, "faustData", "gateData", "root_rawGateList.rds")), file.path(out, tag, "rawGateList"))
  dumpObj(readRDS(file.path(proj, "faustData", "metaData", "possibilityList.rds")), file.path(out, tag, "possibilityList"))
}
runRecon("recon_default", list(), list())
runRecon("recon_supervised", list(W = 2), list(X = 7.5, F = c(5, 6)))
cat("recon done\n")

## ---------- conversion helpers ----------
conv <- list(
  up1 = list(from = c(1.05), to = c(1.0, 3.1)),
  up2 = list(from = c(2.9, 3.0), to = c(1.0, 3.1, 5.0)),      # two from-gates map to one to-gate
  up3 = list(from = c(2.0), to = c(1.0, 3.0)),                  # exact tie: first (lower) wins
  up4 = list(from = c(0.8, 1.3, 2.2), to = c(1.1)),             # the inverted call in the fallback branch
  down1 = list(from = c(0.9, 2.5), to = c(1.05)),
  down2 = list(from = c(1.8, 3.9, 4.1), to = c(2.0, 4.0)),      # two from-gates map to the same to-gate
  down3 = list(from = c(2.0), to = c(1.0, 3.0)),                # fewer from than to (fallback branch)
  down4 = list(from = c(2.9, 3.1), to = c(3.0)),                # equal distances: last wins
  down5 = list(from = c(1.0, 2.0, 3.0), to = c(2.0, 4.0))
)
for (nm in names(conv)) {
  f <- if (startsWith(nm, "up")) faust:::.upConvert else faust:::.downConvert
  res <- tryCatch(f(conv[[nm]]$from, conv[[nm]]$to), error = function(e) "ERROR")
  dumpObj(conv[[nm]]$from, file.path(out, "convert", paste0(nm, "_from")))
  dumpObj(conv[[nm]]$to, file.path(out, "convert", paste0(nm, "_to")))
  dumpObj(res, file.path(out, "convert", paste0(nm, "_out")))
}
cat("convert done\n")

## ---------- elbow ----------
elbow <- list(
  e1 = c(50, 30, 20, 15, 12, 10, 8, 5, 3, 1),
  e2 = c(120, 60, 60, 40, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
  e3 = c(7, 7, 7, 7, 7),
  e4 = c(10, 9, 8, 1),
  e5 = c(100, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
  e6 = c(3, 2, 1)
)
for (nm in names(elbow)) {
  y <- elbow[[nm]]; df <- data.frame(x = seq_along(y), y = y)
  dumpObj(y, file.path(out, "elbow", paste0(nm, "_y")))
  dumpObj(faust:::.computeElbow(df), file.path(out, "elbow", paste0(nm, "_elbow")))
  if (length(y) >= 4) {
    lo <- max(2, ceiling(0.05 * nrow(df))); hi <- min(nrow(df) - 1, floor(0.95 * nrow(df)))
    g <- fdrtool::gcmlcm(x = df$x[lo:hi], y = df$y[lo:hi], type = "gcm")
    dumpObj(g$x.knots, file.path(out, "elbow", paste0(nm, "_xknots")))
    dumpObj(g$y.knots, file.path(out, "elbow", paste0(nm, "_yknots")))
  }
}
## isomean directly, with weights and ties
dumpObj(fdrtool:::pvt.isoMean(c(3, 1, 2, 2, 5, 4), c(1, 2, 1, 1, 3, 1)), file.path(out, "elbow", "isomean1"))
dumpObj(fdrtool:::pvt.isoMean(c(1.5, -2, -2, 0.25, 7, 7, 6.5), c(2, 1, 1, 1, 1, 2, 1)), file.path(out, "elbow", "isomean2"))
cat("elbow done\n")

## ---------- R numeric semantics spot checks ----------
num <- list(
  sum1 = c(0.1, 0.2, 0.3, 1e-17, 1e16, -1e16, 0.7),
  sum2 = c(0.79633578977392605, 0.83262122073205291, 1, 0.82915951718105707, 0.81744849519889695),
  med_even = c(3.7060300782780518, 2.4132228373732123, 1.1760337646717982, 3.7114253928275209),
  mad1 = c(1, 1.1, 0.95, 9.0, 1.2),
  q1 = c(0.3, 0.7, 0.1, 0.9, 0.55, 0.42)
)
dumpObj(num, file.path(out, "num", "in"))
dumpObj(list(sum1 = sum(num$sum1), sum2 = sum(num$sum2), med_even = median(num$med_even), mad1 = mad(num$mad1),
             q1_50 = as.numeric(quantile(num$q1, 0.5)), q1_30 = as.numeric(quantile(num$q1, 0.3)), q1_95 = as.numeric(quantile(num$q1, 0.95)),
             wm = weighted.mean(c(3.320695596908827, 3.7970430343892203, 3.8790180350111636), c(0.2, 0.3, 0.5)),
             qt = qt(0.975, df = 1), cumsum1 = cumsum(num$sum1)),
        file.path(out, "num", "out"))
cat("num done\n")
