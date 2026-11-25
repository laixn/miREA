data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
    background_MGI,
    package = "miREA"
)

test_that("getAllInputData() builds all expected components for the full method list", {
    input_data <- getAllInputData(
        methods = default_method_list,
        pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata,
        gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        GGI_source = "Reactome",
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
        background_MGI = background_MGI
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
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    input_data_bioc <- getAllInputData(
        methods = default_method_list, pathway = pathway_gsc,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = gene_se, mir_mat = mir_se, scoreFun = "rank"
    )

    expect_equal(input_data_bioc$dict_MGI, input_data_std$dict_MGI)
    expect_equal(input_data_bioc$data$Edge_Topology, input_data_std$data$Edge_Topology)
})

test_that("getAllInputData() records the DE thresholds it used, and they actually affect the DE sets", {
    input_data_default <- getAllInputData(
        methods = c("TG_ORA", "Edge_ORA"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
    )
    expect_equal(
        input_data_default$threshold,
        list(DElogFCCutoff = 1, DEpadjCutoff = 0.05, corCutoff = 0, corPadjCutoff = 0.05)
    )

    input_data_strict <- getAllInputData(
        methods = c("TG_ORA", "Edge_ORA"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
        DElogFCCutoff = 10, DEpadjCutoff = 1e-6
    )
    expect_equal(input_data_strict$threshold$DElogFCCutoff, 10)
    expect_equal(input_data_strict$threshold$DEpadjCutoff, 1e-6)

    expect_true(length(input_data_strict$data$TG_ORA) <= length(input_data_default$data$TG_ORA))
    expect_true(length(input_data_strict$data$Edge_ORA) <= length(input_data_default$data$Edge_ORA))
})

test_that("standardize_DEdata() errors (does not guess by position) when column names aren't recognized", {
    df <- data.frame(
        V1 = c("g1", "g2"), V2 = c(1, -2), V3 = c(3, -3), V4 = c(0.01, 0.2),
        stringsAsFactors = FALSE
    )
    expect_error(miREA:::standardize_DEdata(df, "gene"), "recognizable logFC/stat/padj columns")
})

test_that("standardize_DEdata() ignores extra columns and is idempotent on already-standardized input", {
    # extra columns beyond the 4 needed ones (as in full DESeq2::results()/limma::topTable() output)
    df <- data.frame(
        gene = c("g1", "g2"), log2FoldChange = c(1, -2), extra1 = c("x", "y"),
        stat = c(3, -3), extra2 = c(TRUE, FALSE), padj = c(0.01, 0.2),
        stringsAsFactors = FALSE
    )
    out <- miREA:::standardize_DEdata(df, "gene")
    expect_equal(colnames(out), c("gene", "geneLogFC", "geneStat", "genePadj"))
    expect_equal(out$gene, c("g1", "g2"))

    # re-standardizing an already-standardized table (as happens inside get_data_coerce_expr())
    # must be a no-op, not an error.
    out2 <- miREA:::standardize_DEdata(out, "gene")
    expect_equal(out2, out)
})

test_that("standardize_DEdata() matches DESeq2::results()-style columns by name, with ID from rownames", {
    df <- data.frame(
        baseMean = c(10, 20), log2FoldChange = c(1, -2), lfcSE = c(0.1, 0.2),
        stat = c(3, -3), pvalue = c(0.02, 0.3), padj = c(0.01, 0.2),
        row.names = c("g1", "g2")
    )
    out <- miREA:::standardize_DEdata(df, "gene")
    expect_equal(colnames(out), c("gene", "geneLogFC", "geneStat", "genePadj"))
    expect_equal(out$gene, c("g1", "g2"))
    expect_equal(out$geneLogFC, c(1, -2))
    expect_equal(out$geneStat, c(3, -3))
    expect_equal(out$genePadj, c(0.01, 0.2))
})

test_that("standardize_DEdata() matches limma::topTable()-style columns by name, with ID from rownames", {
    df <- data.frame(
        logFC = c(1, -2), AveExpr = c(5, 6), t = c(3, -3),
        P.Value = c(0.02, 0.3), adj.P.Val = c(0.01, 0.2), B = c(0.5, -0.5),
        row.names = c("g1", "g2")
    )
    out <- miREA:::standardize_DEdata(df, "gene")
    expect_equal(colnames(out), c("gene", "geneLogFC", "geneStat", "genePadj"))
    expect_equal(out$gene, c("g1", "g2"))
    expect_equal(out$geneLogFC, c(1, -2))
    expect_equal(out$geneStat, c(3, -3))
    expect_equal(out$genePadj, c(0.01, 0.2))
})

test_that("standardize_DEdata() uses the 'mir' prefix (not 'miRNA') for the miRNA entity's value columns", {
    df <- data.frame(log2FoldChange = 1, stat = 2, padj = 0.01, row.names = "mir1")
    out <- miREA:::standardize_DEdata(df, "miRNA")
    expect_equal(colnames(out), c("miRNA", "mirLogFC", "mirStat", "mirPadj"))
    expect_equal(out$miRNA, "mir1")
})

test_that("standardize_DEdata() errors clearly when columns can't be recognized and rownames are uninformative", {
    df <- data.frame(foo = c(1, 2), bar = c(3, 4), baz = c(0.1, 0.2))
    expect_error(miREA:::standardize_DEdata(df, "gene"), "recognizable logFC/stat/padj columns")

    df_recognized_no_id <- data.frame(log2FoldChange = c(1, 2), stat = c(3, 4), padj = c(0.1, 0.2))
    expect_error(miREA:::standardize_DEdata(df_recognized_no_id, "gene"), "has an ID column")
})

test_that("getAllInputData() accepts DESeq2::results()-style gene_DEdata and produces the same result as the positional form", {
    gene_DEdata_deseq2 <- data.frame(
        log2FoldChange = toy_gene_DEdata$log2FoldChange,
        stat = toy_gene_DEdata$stat,
        padj = toy_gene_DEdata$padj,
        row.names = toy_gene_DEdata$gene
    )

    input_data_positional <- getAllInputData(
        methods = "TG_ORA", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI
    )
    input_data_deseq2 <- getAllInputData(
        methods = "TG_ORA", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = gene_DEdata_deseq2,
        background_MGI = background_MGI
    )

    expect_equal(sort(input_data_deseq2$data$TG_ORA), sort(input_data_positional$data$TG_ORA))
})

test_that("getAllInputData() selects the requested assay via gene_assay/mir_assay", {
    skip_if_not_installed("SummarizedExperiment")

    # calc_MGI_corr() uses Spearman correlation, which is invariant to monotonic increasing
    # rescaling (e.g. x * 10); negating the matrix reverses the rank order instead, so the two
    # assays are guaranteed to produce different correlations.
    gene_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = as.matrix(toy_gene_mat), scaled = -as.matrix(toy_gene_mat))
    )

    input_data_assay1 <- getAllInputData(
        methods = "Edge_ORA", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = gene_se, mir_mat = toy_mir_mat, gene_assay = "counts"
    )
    input_data_assay2 <- getAllInputData(
        methods = "Edge_ORA", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = gene_se, mir_mat = toy_mir_mat, gene_assay = "scaled"
    )

    expect_false(isTRUE(all.equal(input_data_assay1$data$data, input_data_assay2$data$data)))

    input_data_direct <- getAllInputData(
        methods = "Edge_ORA", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
    )
    expect_equal(input_data_assay1$data$data, input_data_direct$data$data)
})

test_that("getAllInputData() falls back to bundled background_GGI_Reactome when background_GGI=NULL, GGI_source='Reactome'", {
    input_data <- getAllInputData(
        methods = "Edge_Topology", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    data(background_GGI_Reactome, package = "miREA")
    expected <- background_GGI_Reactome
    colnames(expected)[seq_len(2)] <- c("from", "to")
    expect_equal(input_data$background_GGI, expected)
})

test_that("getAllInputData() fetches Omnipath live when background_GGI=NULL and GGI_source is 'Omnipath' or NULL", {
    skip_if_not_installed("OmnipathR")
    input_data <- tryCatch(
        getAllInputData(
            methods = "Edge_Topology", pathway = toy_pathway,
            mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
            background_MGI = background_MGI, GGI_source = "Omnipath",
            gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
        ),
        error = function(e) {
            skip(paste("Omnipath/OmnipathR unavailable in this environment:", conditionMessage(e)))
        }
    )
    expect_equal(colnames(input_data$background_GGI), c("from", "to", "interact"))
    expect_true(nrow(input_data$background_GGI) > 0)

    # GGI_source = NULL is documented to behave the same as "Omnipath"
    input_data_null_source <- tryCatch(
        getAllInputData(
            methods = "Edge_Topology", pathway = toy_pathway,
            mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
            background_MGI = background_MGI, GGI_source = NULL,
            gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
        ),
        error = function(e) {
            skip(paste("Omnipath/OmnipathR unavailable in this environment:", conditionMessage(e)))
        }
    )
    expect_equal(colnames(input_data_null_source$background_GGI), c("from", "to", "interact"))
})

test_that("getAllInputData() uses a user-supplied background_GGI directly and ignores GGI_source", {
    custom_GGI <- data.frame(from = c("TP53", "MYC"), to = c("MDM2", "MAX"), stringsAsFactors = FALSE)
    expect_message(
        input_data <- getAllInputData(
            methods = "Edge_Topology", pathway = toy_pathway,
            mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
            background_MGI = background_MGI, background_GGI = custom_GGI, GGI_source = "Omnipath",
            gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
        ),
        "GGI_source is ignored"
    )
    expect_equal(input_data$background_GGI, custom_GGI)
})

test_that("fetch_background_GGI_omnipath() returns the expected structure", {
    skip_if_not_installed("OmnipathR")
    interactions <- tryCatch(
        miREA:::fetch_background_GGI_omnipath(),
        error = function(e) {
            skip(paste("Omnipath/OmnipathR unavailable in this environment:", conditionMessage(e)))
        }
    )
    expect_equal(colnames(interactions), c("source_genesymbol", "target_genesymbol", "interact"))
    expect_true(nrow(interactions) > 0)
})

test_that("getAllInputData() computes strength/MGI for Edge_ORA/Edge_Network alone, without any score method", {
    # Regression test: strength/MGI used to only be computed when methods
    # included Edge_Score1D/Edge_Score2D/Edge_Topology, which made
    # plotHeatmap()/plotHeatmapSankey() fail for input_data built from
    # Edge_ORA/Edge_Network alone.
    expect_message(
        input_data_ora <- getAllInputData(
            methods = "Edge_ORA", pathway = toy_pathway,
            mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
            background_MGI = background_MGI,
            gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
        ),
        "Use the default \\[rank\\] instead"
    )
    expect_true(all(c("strength", "MGI") %in% colnames(input_data_ora$data$data)))

    input_data_network <- getAllInputData(
        methods = "Edge_Network", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    expect_true(all(c("strength", "MGI") %in% colnames(input_data_network$data$data)))

    background <- list(Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene)))
    result_ora <- miREA(
        methods = "Edge_ORA", input_data = input_data_ora, background = background,
        minSize = NULL, maxSize = NULL, pvalueCutoff = 1
    )
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    expect_no_error(plotHeatmap(
        method = "Edge_ORA", result = result_ora, input_data = input_data_ora,
        penrichCutoff = 1
    ))
})
