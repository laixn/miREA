data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
     toy_background_MGI, toy_background_GGI, package = "miREA")

test_that("plotSummary(), plotHeatmap(), and plotHeatmapSankey() run without error", {
  input_data <- suppressWarnings(getAllInputData(
    methods = default_method_list, pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    GGI_source = "Omnipath", gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
    scoreFun = "rank"
  ))
  background <- list(
    TG = unique(toy_pathway$gene),
    MiR = unique(input_data$pathway$MiR$miRNA),
    Edge = unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))
  )
  result <- suppressWarnings(miREA(
    methods = default_method_list, input_data = input_data, background = background,
    minSize = 2, maxSize = 5000,
    pvalueType = list(TG = "std", MiR = "std", Edge = "neg"),
    pAdjMethod = "BH", pvalueCutoff = 1, iter = 20, ncores = 1
  ))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_no_error(plotSummary(result, penrichCutoff = 1))
  expect_no_error(plotHeatmap(method = "Edge_ORA", result = result, input_data = input_data,
                                penrichCutoff = 1))
  expect_no_error(plotHeatmapSankey(method = "Edge_ORA", result = result, input_data = input_data,
                                       penrichCutoff = 1))
})

test_that("plotHeatmap()/plotHeatmapSankey() inherit thresholds from input_data$threshold, and accept overrides", {
  input_data_strict <- getAllInputData(
    methods = c("Edge_ORA", "Edge_Score2D"), pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI,
    gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank",
    DElogFCCutoff = 1.5, DEpadjCutoff = 0.1
  )
  expect_equal(input_data_strict$threshold$DElogFCCutoff, 1.5)

  background <- list(
    Edge = unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))
  )
  result_strict <- miREA(
    methods = "Edge_ORA", input_data = input_data_strict, background = background,
    minSize = NULL, maxSize = NULL, pvalueCutoff = 1
  )

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  # no explicit threshold args: should silently pick up the strict thresholds baked into input_data
  expect_no_error(plotHeatmap(method = "Edge_ORA", result = result_strict, input_data = input_data_strict,
                                penrichCutoff = 1))
  expect_no_error(plotHeatmapSankey(method = "Edge_ORA", result = result_strict, input_data = input_data_strict,
                                       penrichCutoff = 1))
  # explicit args should still override input_data$threshold
  expect_no_error(plotHeatmap(method = "Edge_ORA", result = result_strict, input_data = input_data_strict,
                                penrichCutoff = 1, DElogFCCutoff = 1, DEpadjCutoff = 0.05))

  # input_data without a $threshold element (e.g. hand-built by a user) should fall back to the
  # historical hardcoded defaults rather than error
  input_data_notthr <- input_data_strict
  input_data_notthr$threshold <- NULL
  expect_no_error(plotHeatmap(method = "Edge_ORA", result = result_strict, input_data = input_data_notthr,
                                penrichCutoff = 1))
})
