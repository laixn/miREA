data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
    background_MGI,
    package = "miREA"
)

test_that("getAllInputData() does not modify its data.frame inputs", {
    pathway_df <- toy_pathway
    bg_df <- background_MGI
    gene_df <- toy_gene_DEdata
    mir_df <- toy_mir_DEdata

    pathway_snapshot <- pathway_df
    bg_snapshot <- bg_df
    gene_snapshot <- gene_df
    mir_snapshot <- mir_df

    invisible(getAllInputData(
        methods = "Edge_ORA", pathway = pathway_df, mir_DEdata = mir_df, gene_DEdata = gene_df,
        background_MGI = bg_df, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    ))

    expect_identical(pathway_df, pathway_snapshot)
    expect_identical(bg_df, bg_snapshot)
    expect_identical(gene_df, gene_snapshot)
    expect_identical(mir_df, mir_snapshot)
})

test_that("getAllInputData() does not modify its data.table inputs", {
    pathway_dt <- data.table::as.data.table(toy_pathway)
    bg_dt <- data.table::as.data.table(background_MGI)
    gene_dt <- data.table::as.data.table(toy_gene_DEdata)
    mir_dt <- data.table::as.data.table(toy_mir_DEdata)

    pathway_snapshot <- data.table::copy(pathway_dt)
    bg_snapshot <- data.table::copy(bg_dt)
    gene_snapshot <- data.table::copy(gene_dt)
    mir_snapshot <- data.table::copy(mir_dt)

    invisible(getAllInputData(
        methods = "Edge_ORA", pathway = pathway_dt, mir_DEdata = mir_dt, gene_DEdata = gene_dt,
        background_MGI = bg_dt, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    ))

    expect_identical(pathway_dt, pathway_snapshot)
    expect_identical(bg_dt, bg_snapshot)
    expect_identical(gene_dt, gene_snapshot)
    expect_identical(mir_dt, mir_snapshot)
})

test_that("Edge_Score2D() does not modify its data.table inputs", {
    input_data <- getAllInputData(
        methods = "Edge_Score2D", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
        scoreFun = "rank"
    )
    background <- unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))

    MGI2D_dt <- data.table::as.data.table(input_data$data$Edge_Score2D)
    pathway_dt <- data.table::as.data.table(input_data$pathway$Edge)
    MGI2D_snapshot <- data.table::copy(MGI2D_dt)
    pathway_snapshot <- data.table::copy(pathway_dt)

    invisible(Edge_Score2D(
        MGI2D = MGI2D_dt, pathway = pathway_dt, background = background,
        iter = 5, ncores = 1, pvalueCutoff = 1
    ))

    expect_identical(MGI2D_dt, MGI2D_snapshot)
    expect_identical(pathway_dt, pathway_snapshot)
})

test_that("Edge_Topology() does not modify its data.table inputs", {
    input_data <- getAllInputData(
        methods = "Edge_Topology", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    background <- unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))

    MGI_dt <- data.table::as.data.table(input_data$data$Edge_Topology)
    pathway_dt <- data.table::as.data.table(input_data$pathway$Edge)
    pathway_GGI_dt <- data.table::as.data.table(input_data$pathway_GGI)
    dict_MGI_dt <- data.table::as.data.table(input_data$dict_MGI)

    MGI_snapshot <- data.table::copy(MGI_dt)
    pathway_snapshot <- data.table::copy(pathway_dt)
    pathway_GGI_snapshot <- data.table::copy(pathway_GGI_dt)
    dict_MGI_snapshot <- data.table::copy(dict_MGI_dt)

    invisible(Edge_Topology(
        MGI = MGI_dt, background = background, pathway = pathway_dt, pathway_GGI = pathway_GGI_dt,
        dict_MGI = dict_MGI_dt, pvalueType = "neg", minSize = 2, maxSize = 5000,
        iter = 5, ncores = 1, pvalueCutoff = 1
    ))

    expect_identical(MGI_dt, MGI_snapshot)
    expect_identical(pathway_dt, pathway_snapshot)
    expect_identical(pathway_GGI_dt, pathway_GGI_snapshot)
    expect_identical(dict_MGI_dt, dict_MGI_snapshot)
})
