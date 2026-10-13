## Per-sample 1%/99% quantiles of the float32-rounded expression matrices (what
## .processChannelBounds computes from the GatingSet data), and the empirical bounds rule.
out <- commandArgs(trailingOnly = TRUE)[1]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
g17 <- function(x) sprintf("%.17g", x)
f32 <- function(m) { r <- readBin(writeBin(as.numeric(m), raw(), size = 4), "double", size = 4, n = length(m)); matrix(r, nrow = nrow(m), dimnames = dimnames(m)) }
rd <- function(f) f32(as.matrix(read.csv(f, check.names = FALSE)))
sets <- list(test1 = list(sample = rd("/sp/faust/fix/out/test1_matrix.csv")),
             syn = list(s1 = rd("/sp/faust/fix/out/syn_s1.csv"), s2 = rd("/sp/faust/fix/out/syn_s2.csv")))
for (tag in names(sets)) {
  lows <- highs <- NULL
  for (sn in names(sets[[tag]])) {
    m <- sets[[tag]][[sn]]
    q01 <- apply(m, 2, function(x) quantile(x, probs = 0.01)); q99 <- apply(m, 2, function(x) quantile(x, probs = 0.99))
    writeLines(c(paste(colnames(m), collapse = ","), paste(g17(q01), collapse = ","), paste(g17(q99), collapse = ",")),
               file.path(out, sprintf("%s_%s_q01_q99.csv", tag, sn)))
    lows <- rbind(lows, q01); highs <- rbind(highs, q99)
  }
  low <- as.numeric(apply(lows, 2, function(x) quantile(x, probs = 0.05)))
  high <- as.numeric(apply(highs, 2, function(x) quantile(x, probs = 0.95)))
  writeLines(c(paste(colnames(lows), collapse = ","), paste(g17(low), collapse = ","), paste(g17(high), collapse = ",")),
             file.path(out, sprintf("%s_bounds.csv", tag)))
  ref <- read.csv(sprintf("/sp/faust/fix/out/%s_faust/metaData/channelBoundsUsedByFAUST/allSamples.csv", tag), check.names = FALSE)
  cat(tag, "bounds identical to faust dump:", identical(g17(low), g17(as.numeric(ref[1, -1]))), identical(g17(high), g17(as.numeric(ref[2, -1]))), "\n")
  for (sn in names(sets[[tag]])) {
    m <- sets[[tag]][[sn]]
    res <- matrix(0, nrow(m), ncol(m)); for (j in seq_len(ncol(m))) { res[m[, j] <= low[j], j] <- 1; res[m[, j] >= high[j], j] <- 2 }
    cat(tag, sn, "ePop:", apply(res, 2, function(x) sum(x == 0)), "\n")
  }
}
