data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
     toy_background_MGI, toy_background_GGI, package = "miREA")

test_that("Edge_ORA, Edge_Score1D, Edge_Score2D return enrichment results for the signal pathway", {
  input_data <- getAllInputData(
    methods = c("Edge_ORA", "Edge_Score1D", "Edge_Score2D"), pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
    scoreFun = "rank"
  )
  background <- unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))

  ora <- Edge_ORA(DEMGI = input_data$data$Edge_ORA, background = background,
                   pathway = input_data$pathway$Edge, pvalueCutoff = 1)
  expect_s3_class(ora, "data.frame")
  expect_true("LINDGREN_BLADDER_CANCER_CLUSTER_2B" %in% ora$pathway)

  score1d <- Edge_Score1D(MGIList = input_data$data$Edge_Score1D, pathway = input_data$pathway$Edge,
                           pvalueType = "neg", pvalueCutoff = 1)
  expect_s3_class(score1d, "data.frame")
  expect_true("LINDGREN_BLADDER_CANCER_CLUSTER_2B" %in% score1d$pathway)

  score2d <- Edge_Score2D(MGI2D = input_data$data$Edge_Score2D, pathway = input_data$pathway$Edge,
                           background = background, iter = 20, ncores = 1, pvalueCutoff = 1)
  expect_s3_class(score2d, "data.frame")
  expect_true("LINDGREN_BLADDER_CANCER_CLUSTER_2B" %in% score2d$pathway)
})

test_that("Edge_ORA rejects an empty DEMGI vector", {
  expect_error(Edge_ORA(DEMGI = character(0), pathway = toy_pathway), "DEMGI is empty")
})

test_that("Edge_Topology returns enrichment results for the signal pathway", {
  input_data <- getAllInputData(
    methods = "Edge_Topology", pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    GGI_source = "Omnipath", gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
  )
  background <- unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))

  res <- Edge_Topology(
    MGI = input_data$data$Edge_Topology, background = background,
    pathway = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
    dict_MGI = input_data$dict_MGI, pvalueType = "neg",
    minSize = 2, maxSize = 5000, iter = 20, ncores = 1, pvalueCutoff = 1
  )
  expect_s3_class(res, "data.frame")
  expect_true("LINDGREN_BLADDER_CANCER_CLUSTER_2B" %in% res$pathway)
})

test_that("Edge_Network returns enrichment results for the signal pathway", {
  input_data <- getAllInputData(
    methods = "Edge_Network", pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    GGI_source = "Omnipath", gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
  )

  res <- Edge_Network(
    DEMGI = input_data$data$Edge_Network, pathway_gene = input_data$pathway$TG,
    pathway_MGI = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
    MGI = input_data$MGI, GGI = input_data$GGI,
    minSize = 2, maxSize = 5000, iter = 20, ncores = 1, pvalueCutoff = 1
  )
  expect_s3_class(res, "data.frame")
  expect_true("LINDGREN_BLADDER_CANCER_CLUSTER_2B" %in% res$pathway)
})

test_that("Edge_Score1D rejects invalid MGIList/pAdjMethod inputs", {
  expect_error(Edge_Score1D(MGIList = NULL, pathway = toy_pathway), "must be entered")
  expect_error(Edge_Score1D(MGIList = "not-numeric", pathway = toy_pathway), "must be a named number vector")
  expect_error(Edge_Score1D(MGIList = c(a = 1), pathway = toy_pathway, pAdjMethod = "bogus"), "Invalid pAdjMethod")
})

test_that("Edge_ORA rejects an invalid pAdjMethod", {
  expect_error(
    Edge_ORA(DEMGI = "a:b", background = "a:b", pathway = toy_pathway, pAdjMethod = "bogus"),
    "Invalid pAdjMethod"
  )
})

test_that("Edge_ORA errors when background is not supplied", {
  expect_error(
    Edge_ORA(DEMGI = "a:b", pathway = toy_pathway),
    "background must be supplied"
  )
})

test_that("Edge_ORA errors when all DEMGI are filtered out by the background", {
  expect_error(
    Edge_ORA(DEMGI = "nonexistent-mir:nonexistent-gene", background = "other-mir:other-gene",
             pathway = toy_pathway),
    "No valid DEMGI input"
  )
})

test_that("Edge_ORA filters DEMGI entries missing from the background with a message", {
  input_data <- getAllInputData(
    methods = "Edge_ORA", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
    gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI,
    gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
  )
  demgi <- input_data$data$Edge_ORA
  background_full <- unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))
  background_partial <- setdiff(background_full, demgi[1])
  expect_message(
    res <- Edge_ORA(DEMGI = demgi, background = background_partial,
                     pathway = input_data$pathway$Edge, pvalueCutoff = 1),
    "not included in the background"
  )
  expect_s3_class(res, "data.frame")
})

test_that("Edge_Topology rejects invalid pAdjMethod/pvalueType", {
  args <- list(MGI = data.frame(), background = character(),
               pathway = data.frame(), pathway_GGI = data.frame(), dict_MGI = data.frame())
  expect_error(do.call(Edge_Topology, c(args, list(pAdjMethod = "bogus"))), "Invalid pAdjMethod")
  expect_error(do.call(Edge_Topology, c(args, list(pvalueType = "bogus"))), "Invalid pvalueType")
})

test_that("Edge_Topology errors when no MGI pairs match the background", {
  expect_error(
    Edge_Topology(
      MGI = data.frame(MGI = "nonexistent:pair", miRNA = "x", gene = "y", strength = 1),
      background = "other:pair",
      pathway = data.frame(pathway = character(), miRNA = character(), gene = character()),
      pathway_GGI = data.frame(pathway = character(), from = character(), to = character()),
      dict_MGI = data.frame()
    ),
    "No valid MGI input"
  )
})

test_that("Edge_Topology handles n_pathway=0, the iter=1 degenerate branch, and ncores>1", {
  input_data <- getAllInputData(
    methods = "Edge_Topology", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
    gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI,
    background_GGI = toy_background_GGI, GGI_source = "Omnipath",
    gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
  )
  background <- unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))

  expect_warning(
    empty_res <- Edge_Topology(
      MGI = input_data$data$Edge_Topology, background = background,
      pathway = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
      dict_MGI = input_data$dict_MGI, minSize = 1e6, iter = 1
    ),
    "No pathways meet the size criteria"
  )
  expect_equal(nrow(empty_res), 0)

  res_iter1 <- suppressWarnings(Edge_Topology(
    MGI = input_data$data$Edge_Topology, background = background,
    pathway = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
    dict_MGI = input_data$dict_MGI, iter = 1, ncores = 1, pvalueCutoff = 1
  ))
  expect_s3_class(res_iter1, "data.frame")

  res_parallel <- Edge_Topology(
    MGI = input_data$data$Edge_Topology, background = background,
    pathway = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
    dict_MGI = input_data$dict_MGI, iter = 5, ncores = 2, pvalueCutoff = 1
  )
  expect_s3_class(res_parallel, "data.frame")
})

test_that("Edge_Topology_single warns on multi-pathway input and handles PW_GGI=NULL", {
  input_data <- getAllInputData(
    methods = "Edge_Topology", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
    gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI,
    background_GGI = toy_background_GGI, GGI_source = "Omnipath",
    gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
  )
  pathway <- as.data.frame(input_data$pathway$Edge)
  pathway$MGI <- paste0(pathway$miRNA, ":", pathway$gene)
  pathway_GGI <- as.data.frame(input_data$pathway_GGI)
  pw_names <- unique(pathway$pathway)[seq_len(2)]
  pw_two <- pathway[pathway$pathway %in% pw_names, ]
  ggi_two <- pathway_GGI[pathway_GGI$pathway %in% pw_names, ]
  expect_warning(
    miREA:::Edge_Topology_single(
      MGI = input_data$data$Edge_Topology, PW = pw_two, PW_GGI = ggi_two,
      dict_MGI = input_data$dict_MGI, pvalueType = "neg", iter = 5
    ),
    "not calculated one-by-one"
  )

  pw_name <- unique(pathway$pathway)[1]
  pw_one <- pathway[pathway$pathway == pw_name, ]
  res <- miREA:::Edge_Topology_single(
    MGI = input_data$data$Edge_Topology, PW = pw_one, PW_GGI = NULL,
    dict_MGI = input_data$dict_MGI, pvalueType = "neg", iter = 5
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
})

test_that("Edge_Score2D rejects invalid MGI2D/background inputs", {
  expect_error(
    Edge_Score2D(MGI2D = data.frame(MGI = character(), miRNA = character(), gene = character(),
                                     x = numeric(), y = numeric()),
                 pathway = toy_pathway, background = "a:b"),
    "No valid MGI pairs"
  )
  expect_error(
    Edge_Score2D(MGI2D = data.frame(MGI = "a:b", miRNA = "a", gene = "b", x = 0.1, y = 0.1),
                 pathway = toy_pathway, background = NULL),
    "background must be supplied"
  )
  expect_error(
    Edge_Score2D(MGI2D = data.frame(MGI = "a:b", miRNA = "a", gene = "b", x = 0.1, y = 0.1),
                 pathway = toy_pathway, background = "other:pair"),
    "No valid MGI input"
  )
})

test_that("Edge_Score2D handles n_pathway=0 early return and the ncores>1 parallel path", {
  input_data <- getAllInputData(
    methods = "Edge_Score2D", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
    gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI,
    gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
  )
  background <- unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))

  expect_warning(
    empty_res <- Edge_Score2D(
      MGI2D = input_data$data$Edge_Score2D, pathway = input_data$pathway$Edge,
      background = background, minSize = 1e6, iter = 1, pvalueCutoff = 1
    ),
    "No pathways meet the size criteria"
  )
  expect_equal(nrow(empty_res), 0)

  res_parallel <- Edge_Score2D(
    MGI2D = input_data$data$Edge_Score2D, pathway = input_data$pathway$Edge,
    background = background, iter = 5, ncores = 2, pvalueCutoff = 1
  )
  expect_s3_class(res_parallel, "data.frame")
})

test_that("Edge_Score2D_single warns on multi-pathway input and returns NULL for no overlap", {
  input_data <- getAllInputData(
    methods = "Edge_Score2D", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
    gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI,
    gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
  )
  mgi2d <- data.table::as.data.table(na.omit(input_data$data$Edge_Score2D))
  colnames(mgi2d) <- c("MGI", "miRNA", "gene", "x", "y")
  r_random <- sort(stats::runif(50, min = 0, max = 3))

  pw_names <- unique(input_data$pathway$Edge$pathway)[seq_len(2)]
  pw_two <- input_data$pathway$Edge[input_data$pathway$Edge$pathway %in% pw_names, ]
  expect_warning(
    miREA:::Edge_Score2D_single(MGI2D = mgi2d, PW = pw_two, r_random = r_random),
    "not calculated one-by-one"
  )

  pw_no_overlap <- data.frame(pathway = "fake_pw", miRNA = "fake_mir", gene = "fake_gene")
  expect_null(miREA:::Edge_Score2D_single(MGI2D = mgi2d, PW = pw_no_overlap, r_random = r_random))
})

test_that("Edge_Score2D/Edge_Topology/Edge_Network accept an explicit BPPARAM", {
  input_data <- getAllInputData(
    methods = c("Edge_Score2D", "Edge_Topology", "Edge_Network"), pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    GGI_source = "Omnipath", gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
  )
  background <- unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))

  score2d <- Edge_Score2D(
    MGI2D = input_data$data$Edge_Score2D, pathway = input_data$pathway$Edge,
    background = background, iter = 5, BPPARAM = BiocParallel::SerialParam(), pvalueCutoff = 1
  )
  expect_s3_class(score2d, "data.frame")

  topology <- Edge_Topology(
    MGI = input_data$data$Edge_Topology, background = background,
    pathway = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
    dict_MGI = input_data$dict_MGI, iter = 5, BPPARAM = BiocParallel::SerialParam(), pvalueCutoff = 1
  )
  expect_s3_class(topology, "data.frame")

  network <- Edge_Network(
    DEMGI = input_data$data$Edge_Network, pathway_gene = input_data$pathway$TG,
    pathway_MGI = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
    MGI = input_data$MGI, GGI = input_data$GGI,
    minSize = 2, maxSize = 5000, iter = 20, BPPARAM = BiocParallel::SerialParam(), pvalueCutoff = 1
  )
  expect_s3_class(network, "data.frame")
})

test_that("Edge_Network runs through the parallel RWR path with ncores > 1", {
  input_data <- getAllInputData(
    methods = "Edge_Network", pathway = toy_pathway,
    mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
    background_MGI = toy_background_MGI, background_GGI = toy_background_GGI,
    GGI_source = "Omnipath", gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
  )

  res <- Edge_Network(
    DEMGI = input_data$data$Edge_Network, pathway_gene = input_data$pathway$TG,
    pathway_MGI = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
    MGI = input_data$MGI, GGI = input_data$GGI,
    minSize = 2, maxSize = 5000, iter = 20, ncores = 2, pvalueCutoff = 1
  )
  expect_s3_class(res, "data.frame")
  expect_true("LINDGREN_BLADDER_CANCER_CLUSTER_2B" %in% res$pathway)
})
