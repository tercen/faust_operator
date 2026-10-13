cat("faust", as.character(packageVersion("faust")), "\n", file = "/o/ties_version.txt")
# A vector with ties at the bottom, inside and at the top of the unique range, plus a
# negative zero, so every branch of the uniform tie pass draws; row order is unsorted.
x <- c(1, 3, 1, 2, 5, 5, 3, 0, -0, 4, 5, 2.5, 1, 5, 3)
for (s in c(123, 7)) {
  y <- faust:::addNoiseToDataVector(x, 4, s)
  writeLines(sprintf("%.17g", y), sprintf("/o/ties_noise_s%d.txt", s))
}
# A single-column matrix via the matrix entry point is not exported; done.
writeLines(sprintf("%.17g", x), "/o/ties_x.txt")
