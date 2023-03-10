suppressPackageStartupMessages({
  library(tercen)
  library(tercenApi)
  library(dplyr)
  library(faust)
  library(tim)
})


ctx <- tercenCtx()

seed <- ctx$op.value("seed", as.double, 123)
depthScoreThreshold <- ctx$op.value("depthScoreThreshold", as.double, 0.01)
selectionQuantile <-  ctx$op.value("selectionQuantile", as.double, 0.5)

if("filename" %in% names(ctx$cnames)) {
  mat <- ctx %>% as.matrix() %>% t()
  colnames(mat) <- ctx$rselect()[[1]]
  df <- cbind(ctx$cselect("filename"), mat)
} else {
  mat <- ctx %>% as.matrix() %>% t()
  colnames(mat) <- ctx$rselect()[[1]]
  df <- as.data.frame(mat)
  df$filename <- "sample"
}

df_list <- split(df, df$filename)
ff_list <- lapply(df_list, function(x) {
  m <- x %>% dplyr::select(-filename) %>% as.matrix
  tim::matrix_to_flowFrame(m)
})

fs <- as(ff_list, "flowSet")
gs <- flowWorkspace::GatingSet(fs)

projPath <- file.path(tempdir(), "FAUST")
on.exit(unlink(projPath, recursive = TRUE))
dir.create(projPath, recursive = TRUE)

faust::faust(
  gatingSet           = gs,
  startingCellPop     = "root",
  projectPath         = projPath,
  annotationsApproved = TRUE,
  threadNum           = 4,
  selectionQuantile = selectionQuantile,
  seedValue = seed,
  depthScoreThreshold = depthScoreThreshold,
  densitySubSampleThreshold = 1e6,
  densitySubSampleSize = 1e6,
  drawAnnotationHistograms = FALSE,
  plottingDevice = "png"
)

snVec <- list.files(file.path(projPath, "faustData", "sampleData"))

annoEmbed <- faust::makeAnnotationEmbedding(	      
  projectPath=projPath,
  sampleNameVec=snVec
)

df_out <- annoEmbed %>%
  dplyr::select(umapX, umapY, faustLabels, contains("annotation")) %>%
  mutate(umapX = as.double(umapX)) %>%
  mutate(umapY = as.double(umapY)) %>%
  replace(. == "-", 0) %>%
  replace(. == "+", 1) %>%
  mutate_at(vars(contains("annotation")), as.numeric) %>%
  mutate(.i = seq_len(nrow(.)) - 1L) %>%
  rename(umap_1 = umapX, umap_2 = umapY, label = faustLabels)

colnames(df_out) <- gsub("_faust_annotation", "_gate", colnames(df_out))

# get plots and return table
diagnostic_plots <- list.files(path = projPath, pattern = ".png", recursive = TRUE, full.names = TRUE)
diagnostic_plots <- diagnostic_plots[grep("hist_", diagnostic_plots)]

df_out_png <- tim::png_to_df(diagnostic_plots)
# for( i in seq(1, length(diagnostic_plots))){
#   if(i == 1){
#     df_out_png <- tim::plot_file_to_df(diagnostic_plots[i])    
#   }else{
#     df_out_png <- rbind(df_out_png, tim::plot_file_to_df(diagnostic_plots[i]) )
#   }
# }
# df_out_png <- do.call(rbind, lapply(diagnostic_plots, tim::plot_file_to_df))

# output results
join_png = df_out_png %>% 
  ctx$addNamespace() %>%
  as_relation()

join_res = df_out %>%
  ctx$addNamespace() %>%
  as_relation() %>%
  left_join_relation(ctx$crelation, ".i", ctx$crelation$rids) %>%
  left_join_relation(join_png, list(), list()) %>%
  as_join_operator(ctx$cnames, ctx$cnames)

join_res %>%
  save_relation(ctx)
