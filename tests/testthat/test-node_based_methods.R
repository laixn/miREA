data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
    background_MGI,
    package = "miREA"
)

test_that("TG_ORA and TG_Score return enrichment results for the signal pathway", {
    input_data <- getAllInputData(
        methods = c("TG_ORA", "TG_Score"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI
    )

    ora <- TG_ORA(
        TG = input_data$data$TG_ORA, background = unique(toy_pathway$gene),
        pathway = input_data$pathway$TG, pvalueCutoff = 1
    )
    expect_s3_class(ora, "data.frame")
    expect_true(nrow(ora) > 0)

    score <- TG_Score(
        geneList = input_data$data$TG_Score, pathway = input_data$pathway$TG,
        pvalueType = "std", pvalueCutoff = 1
    )
    expect_s3_class(score, "data.frame")
    expect_true(nrow(score) > 0)
})

test_that("TG_ORA rejects an empty TG vector", {
    expect_error(TG_ORA(TG = character(0), pathway = toy_pathway), "TG is empty")
})

test_that("MiR_ORA and MiR_Score return enrichment results for the signal pathway", {
    input_data <- getAllInputData(
        methods = c("MiR_ORA", "MiR_Score"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI
    )

    ora <- MiR_ORA(
        DEmiR = input_data$data$MiR_ORA, background = unique(input_data$pathway$MiR$miRNA),
        pathway = input_data$pathway$MiR, pvalueCutoff = 1
    )
    expect_s3_class(ora, "data.frame")
    expect_true(nrow(ora) > 0)

    score <- MiR_Score(
        miRList = input_data$data$MiR_Score, pathway = input_data$pathway$MiR,
        pvalueType = "std", pvalueCutoff = 1
    )
    expect_s3_class(score, "data.frame")
    expect_true(nrow(score) > 0)
})

test_that("MiR_ORA rejects an empty DEmiR vector", {
    expect_error(MiR_ORA(DEmiR = character(0), pathway = toy_pathway), "must be entered")
})

test_that("TG_ORA rejects an invalid pAdjMethod", {
    expect_error(
        TG_ORA(TG = "GENE1", pathway = toy_pathway, pAdjMethod = "bogus"),
        "Invalid pAdjMethod"
    )
})

test_that("TG_Score rejects invalid geneList/pAdjMethod inputs", {
    expect_error(TG_Score(geneList = NULL, pathway = toy_pathway), "must be entered")
    expect_error(TG_Score(geneList = "not-numeric", pathway = toy_pathway), "must be a named number vector")
    expect_error(TG_Score(geneList = c(a = 1), pathway = toy_pathway, pAdjMethod = "bogus"), "Invalid pAdjMethod")
})

test_that("MiR_ORA rejects an invalid pAdjMethod", {
    expect_error(
        MiR_ORA(DEmiR = "hsa-miR-1", pathway = toy_pathway, pAdjMethod = "bogus"),
        "Invalid pAdjMethod"
    )
})

test_that("MiR_ORA errors when all DEmiRs are filtered out by the background", {
    expect_error(
        MiR_ORA(DEmiR = "nonexistent-mir", background = "other-mir", pathway = toy_pathway),
        "No valid DEmiR input"
    )
})

test_that("MiR_ORA filters DEmiR entries missing from the background with a message", {
    input_data <- getAllInputData(
        methods = "MiR_ORA", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
        gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI
    )
    demir <- input_data$data$MiR_ORA
    background_full <- unique(input_data$pathway$MiR$miRNA)
    background_partial <- setdiff(background_full, demir[1])
    expect_message(
        res <- MiR_ORA(
            DEmiR = demir, background = background_partial,
            pathway = input_data$pathway$MiR, pvalueCutoff = 1
        ),
        "not included in the background"
    )
    expect_s3_class(res, "data.frame")
})

test_that("MiR_ORA falls back to the miRBase background when none is supplied", {
    input_data <- getAllInputData(
        methods = "MiR_ORA", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
        gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI
    )
    expect_message(
        res <- MiR_ORA(
            DEmiR = input_data$data$MiR_ORA, background = NULL,
            pathway = input_data$pathway$MiR, pvalueCutoff = 1
        ),
        "extract all mature miRNAs from miRBase"
    )
    expect_s3_class(res, "data.frame")
})

test_that("MiR_Score rejects invalid miRList/pAdjMethod inputs", {
    expect_error(MiR_Score(miRList = NULL, pathway = toy_pathway), "must be entered")
    expect_error(MiR_Score(miRList = "not-numeric", pathway = toy_pathway), "must be a named number vector")
    expect_error(MiR_Score(miRList = c(a = 1), pathway = toy_pathway, pAdjMethod = "bogus"), "Invalid pAdjMethod")
})
