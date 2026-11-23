data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
     toy_background_MGI, toy_background_GGI, package = "miREA")

test_that("getAllInputData() builds all expected components for the full method list", {
  input_data <- getAllInputData(
    methods = default_method_list,
    pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata,
    gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI,
    background_GGI = toy_background_GGI,
    GGI_source = "Omnipath",
    gene_mat = toy_gene_mat,
    mir_mat = toy_mir_mat,
    scoreFun = "rank"
  )

  expect_type(input_data, "list")
  expect_true(all(c("data", "pathway", "pathway_GGI", "GGI", "MGI", "dict_MGI") %in% names(input_data)))
  expect_true(all(default_method_list %in% names(input_data$data)))
  expect_s3_class(input_data$dict_MGI, "data.frame")
  expect_true(nrow(input_data$dict_MGI) > 0)
})

test_that("getAllInputData() falls back to the bundled background_MGI when NULL", {
  expect_message(
    input_data <- getAllInputData(
      methods = "TG_ORA", pathway = toy_pathway,
      mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata
    ),
    "Using the bundled default background_MGI"
  )
  data(background_MGI, package = "miREA")
  expect_equal(nrow(input_data$background_MGI), nrow(background_MGI))
})

test_that("getAllInputData() works for a single method without the full pipeline", {
  input_data <- getAllInputData(
    methods = "TG_ORA", pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI
  )
  expect_true(is.character(input_data$data$TG_ORA))
  expect_true(length(input_data$data$TG_ORA) > 0)
})

test_that("getAllInputData() accepts SummarizedExperiment and GeneSetCollection inputs", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("GSEABase")

  gene_se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = as.matrix(toy_gene_mat)))
  mir_se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = as.matrix(toy_mir_mat)))
  pathway_gsc <- GSEABase::GeneSetCollection(Map(function(genes, name) {
    GSEABase::GeneSet(genes, setName = name)
  }, split(toy_pathway$gene, toy_pathway$pathway), names(split(toy_pathway$gene, toy_pathway$pathway))))

  input_data_std <- getAllInputData(
    methods = default_method_list, pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
  )
  input_data_bioc <- getAllInputData(
    methods = default_method_list, pathway = pathway_gsc,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    gene_mat = gene_se, mir_mat = mir_se, scoreFun = "rank"
  )

  expect_equal(input_data_bioc$dict_MGI, input_data_std$dict_MGI)
  expect_equal(input_data_bioc$data$Edge_Topology, input_data_std$data$Edge_Topology)
})

test_that("getAllInputData() records the DE thresholds it used, and they actually affect the DE sets", {
  input_data_default <- getAllInputData(
    methods = c("TG_ORA", "Edge_ORA"), pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
  )
  expect_equal(
    input_data_default$threshold,
    list(DElogFCCutoff = 1, DEpadjCutoff = 0.05, corCutoff = 0, corPadjCutoff = 0.05)
  )

  input_data_strict <- getAllInputData(
    methods = c("TG_ORA", "Edge_ORA"), pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
    DElogFCCutoff = 10, DEpadjCutoff = 1e-6
  )
  expect_equal(input_data_strict$threshold$DElogFCCutoff, 10)
  expect_equal(input_data_strict$threshold$DEpadjCutoff, 1e-6)

  expect_true(length(input_data_strict$data$TG_ORA) <= length(input_data_default$data$TG_ORA))
  expect_true(length(input_data_strict$data$Edge_ORA) <= length(input_data_default$data$Edge_ORA))
})
