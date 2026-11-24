data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
    background_MGI,
    package = "miREA"
)

test_that("plotSummary(), plotHeatmap(), and plotHeatmapSankey() run without error", {
    input_data <- suppressWarnings(cached_getAllInputData(
        methods = default_method_list, pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
        scoreFun = "rank"
    ))
    background <- list(
        TG = unique(toy_pathway$gene),
        MiR = unique(input_data$pathway$MiR$miRNA),
        Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))
    )
    result <- suppressWarnings(cached_miREA(
        methods = default_method_list, input_data = input_data, background = background,
        minSize = 2, maxSize = 5000,
        pvalueType = list(TG = "std", MiR = "std", Edge = "neg"),
        pAdjMethod = "BH", pvalueCutoff = 1, iter = 5, ncores = 1
    ))

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    expect_no_error(plotSummary(result, penrichCutoff = 1))
    expect_no_error(plotHeatmap(
        method = "Edge_ORA", result = result, input_data = input_data,
        penrichCutoff = 1
    ))
    expect_no_error(plotHeatmapSankey(
        method = "Edge_ORA", result = result, input_data = input_data,
        penrichCutoff = 1
    ))
})

test_that("plotHeatmap()/plotHeatmapSankey() inherit thresholds from input_data$threshold, and accept overrides", {
    input_data_strict <- getAllInputData(
        methods = c("Edge_ORA", "Edge_Score2D"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank",
        DElogFCCutoff = 1.5, DEpadjCutoff = 0.1
    )
    expect_equal(input_data_strict$threshold$DElogFCCutoff, 1.5)

    background <- list(
        Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))
    )
    result_strict <- miREA(
        methods = "Edge_ORA", input_data = input_data_strict, background = background,
        minSize = NULL, maxSize = NULL, pvalueCutoff = 1
    )

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    # no explicit threshold args: should silently pick up the strict thresholds baked into input_data
    expect_no_error(plotHeatmap(
        method = "Edge_ORA", result = result_strict, input_data = input_data_strict,
        penrichCutoff = 1
    ))
    expect_no_error(plotHeatmapSankey(
        method = "Edge_ORA", result = result_strict, input_data = input_data_strict,
        penrichCutoff = 1
    ))
    # explicit args should still override input_data$threshold
    expect_no_error(plotHeatmap(
        method = "Edge_ORA", result = result_strict, input_data = input_data_strict,
        penrichCutoff = 1, DElogFCCutoff = 1, DEpadjCutoff = 0.05
    ))

    # input_data without a $threshold element (e.g. hand-built by a user) should fall back to the
    # historical hardcoded defaults rather than error
    input_data_notthr <- input_data_strict
    input_data_notthr$threshold <- NULL
    expect_no_error(plotHeatmap(
        method = "Edge_ORA", result = result_strict, input_data = input_data_notthr,
        penrichCutoff = 1
    ))
})

test_that("plotHeatmap()/plotHeatmapSankey() always invisibly return list(ht_list, legend_list)", {
    input_data <- cached_getAllInputData(
        methods = c("Edge_ORA", "Edge_Score2D"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    background <- list(Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene)))
    result <- miREA(
        methods = "Edge_ORA", input_data = input_data, background = background,
        minSize = NULL, maxSize = NULL, pvalueCutoff = 1
    )

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    res_window <- plotHeatmap(method = "Edge_ORA", result = result, input_data = input_data, penrichCutoff = 1)
    expect_type(res_window, "list")
    expect_equal(names(res_window), c("ht_list", "legend_list"))

    plot_dir <- withr::local_tempdir()
    res_saved <- plotHeatmap(
        method = "Edge_ORA", result = result, input_data = input_data,
        penrichCutoff = 1, plot_path = plot_dir
    )
    expect_type(res_saved, "list")
    expect_equal(names(res_saved), c("ht_list", "legend_list"))
    expect_true(file.exists(file.path(plot_dir, "heatmap_Edge_ORA.pdf")))

    sankey_dir <- withr::local_tempdir()
    res_sankey <- plotHeatmapSankey(
        method = "Edge_ORA", result = result, input_data = input_data,
        penrichCutoff = 1, plot_path = sankey_dir
    )
    expect_type(res_sankey, "list")
    expect_equal(names(res_sankey), c("ht_list", "legend_list"))
    expect_true(file.exists(file.path(sankey_dir, "heatmap_sankey_Edge_ORA.pdf")))
})

test_that("plotHeatmap()/plotHeatmapSankey() close the pdf device they open (on.exit/dev.off)", {
    input_data <- cached_getAllInputData(
        methods = c("Edge_ORA", "Edge_Score2D"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    background <- list(Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene)))
    result <- miREA(
        methods = "Edge_ORA", input_data = input_data, background = background,
        minSize = NULL, maxSize = NULL, pvalueCutoff = 1
    )

    n_devices_before <- length(dev.list())
    plot_dir <- withr::local_tempdir()
    plotHeatmap(
        method = "Edge_ORA", result = result, input_data = input_data,
        penrichCutoff = 1, plot_path = plot_dir
    )
    expect_equal(length(dev.list()), n_devices_before)
})

test_that("plotHeatmap() derives legend titles from names(gene_annot)/names(mir_annot)", {
    input_data <- cached_getAllInputData(
        methods = c("Edge_ORA", "Edge_Score2D"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    background <- list(Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene)))
    result <- miREA(
        methods = "Edge_ORA", input_data = input_data, background = background,
        minSize = NULL, maxSize = NULL, pvalueCutoff = 1
    )

    some_genes <- unique(input_data$data$data$gene)[seq_len(3)]
    some_mirs <- unique(input_data$data$data$miRNA)[seq_len(3)]

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    res <- plotHeatmap(
        method = "Edge_ORA", result = result, input_data = input_data, penrichCutoff = 1,
        gene_annot = list(TP53_pathway = some_genes), mir_annot = list(let7_family = some_mirs)
    )

    lgd_titles <- vapply(res$legend_list, function(lgd) lgd@grob$children[[1]]$label, character(1))
    expect_true(any(grepl("TP53_pathway", lgd_titles)))
    expect_true(any(grepl("let7_family", lgd_titles)))
})

test_that("plotHeatmap() rejects an unsupported method and a non-edge-based method", {
    input_data <- getAllInputData(
        methods = c("Edge_ORA", "TG_ORA"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
    )
    background <- list(
        TG = unique(toy_pathway$gene),
        Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))
    )
    result <- miREA(
        methods = c("Edge_ORA", "TG_ORA"), input_data = input_data, background = background,
        minSize = NULL, maxSize = NULL, pvalueCutoff = 1
    )

    expect_error(
        plotHeatmap(method = "bogus_method", result = result, input_data = input_data),
        "No valid result for method"
    )
    expect_error(
        plotHeatmap(method = "TG_ORA", result = result, input_data = input_data),
        "Only useful when the methods include Edge-based method"
    )
    expect_error(
        plotHeatmap(method = c("Edge_ORA", "TG_ORA"), result = result, input_data = input_data),
        "Please specify only ONE"
    )
})

test_that("plotHeatmap() falls back to the top n_pathway smallest-padj pathways when none pass penrichCutoff", {
    input_data <- cached_getAllInputData(
        methods = c("Edge_ORA", "Edge_Score2D"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    background <- list(Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene)))
    result <- miREA(
        methods = "Edge_ORA", input_data = input_data, background = background,
        minSize = NULL, maxSize = NULL, pvalueCutoff = 1
    )

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    expect_warning(
        res <- plotHeatmap(
            method = "Edge_ORA", result = result, input_data = input_data,
            penrichCutoff = 1e-300, n_pathway = 3
        ),
        "No enriched items"
    )
    expect_type(res, "list")
    expect_equal(nrow(res$ht_list@ht_list$ht1@matrix), 3)
})

test_that("plotHeatmap() truncates to n_mir_heatmap miRNAs", {
    input_data <- cached_getAllInputData(
        methods = c("Edge_ORA", "Edge_Score2D"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI,
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    background <- list(Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene)))
    result <- miREA(
        methods = "Edge_ORA", input_data = input_data, background = background,
        minSize = NULL, maxSize = NULL, pvalueCutoff = 1
    )

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    res <- plotHeatmap(
        method = "Edge_ORA", result = result, input_data = input_data,
        penrichCutoff = 1, n_mir_heatmap = 2
    )
    expect_lte(nrow(res$ht_list@ht_list$ht2@matrix), 2)
})
