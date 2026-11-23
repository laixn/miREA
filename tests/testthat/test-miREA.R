data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
    background_MGI,
    package = "miREA"
)

test_that("miREA() runs all nine methods and returns a result per method", {
    input_data <- suppressWarnings(getAllInputData(
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

    result <- suppressWarnings(miREA(
        methods = default_method_list, input_data = input_data, background = background,
        minSize = 2, maxSize = 5000,
        pvalueType = list(TG = "std", MiR = "std", Edge = "neg"),
        pAdjMethod = "BH", pvalueCutoff = 1, iter = 20, ncores = 1
    ))

    expect_setequal(result$valid_methods, default_method_list)
    for (m in default_method_list) {
        expect_s3_class(result$result[[m]], "data.frame")
        expect_gt(nrow(result$result[[m]]), 0)
    }
})

test_that("miREA() records the DE thresholds used to build input_data in its parameter slot", {
    input_data <- getAllInputData(
        methods = "TG_ORA", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, DElogFCCutoff = 2, DEpadjCutoff = 0.01
    )
    result <- miREA(
        methods = "TG_ORA", input_data = input_data,
        background = list(TG = unique(toy_pathway$gene)), pvalueCutoff = 1
    )
    expect_equal(result$parameter$DE_threshold, input_data$threshold)
    expect_equal(result$parameter$DE_threshold$DElogFCCutoff, 2)
})

test_that("miREA() requires a background when none is supplied", {
    input_data <- suppressWarnings(getAllInputData(
        methods = "TG_ORA", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI
    ))
    expect_error(
        miREA(methods = "TG_ORA", input_data = input_data, background = NULL),
        "background must be supplied"
    )
})

test_that("miREA() validates input_data structure", {
    expect_error(
        miREA(methods = "TG_ORA", input_data = list(pathway = list()), background = list(TG = "a")),
        "missing the 'data' component"
    )
    expect_error(
        miREA(methods = "TG_ORA", input_data = list(data = list(TG_ORA = "a")), background = list(TG = "a")),
        "missing the 'pathway' component"
    )
    expect_error(
        miREA(
            methods = "TG_ORA",
            input_data = list(data = list(), pathway = list(TG = data.frame())),
            background = list(TG = "a")
        ),
        "missing in input_data\\$data"
    )
    expect_error(
        miREA(
            methods = "TG_ORA",
            input_data = list(data = list(TG_ORA = "a"), pathway = list()),
            background = list(TG = "a")
        ),
        "missing in input_data\\$pathway"
    )
})

test_that("miREA() warns and drops invalid method names, erroring if none remain", {
    input_data <- getAllInputData(
        methods = "TG_ORA", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
        gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI
    )
    expect_warning(
        res <- miREA(
            methods = c("TG_ORA", "bogus_method"), input_data = input_data,
            background = list(TG = unique(toy_pathway$gene))
        ),
        "cannot be specified"
    )
    expect_true("TG_ORA" %in% res$valid_methods)

    expect_error(
        miREA(methods = "bogus_method", input_data = input_data, background = list(TG = "a")),
        "No valid methods"
    )
})

test_that("miREA() requires pathway_GGI/dict_MGI/GGI/MGI for topology/network methods", {
    input_data <- getAllInputData(
        methods = "TG_ORA", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
        gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI
    )
    input_data$data$Edge_Topology <- data.frame()
    input_data$pathway$Edge <- data.frame()
    expect_error(
        miREA(methods = "Edge_Topology", input_data = input_data, background = list(Edge = "a")),
        "didn't specified pathway_GGI"
    )

    input_data$data$Edge_Network <- data.frame()
    expect_error(
        miREA(methods = "Edge_Network", input_data = input_data, background = list(Edge = "a")),
        "didn't specified pathway_GGI"
    )
})

test_that("miREA() validates pvalueType and ncores", {
    input_data <- getAllInputData(
        methods = "TG_ORA", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
        gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI
    )
    background <- list(TG = unique(toy_pathway$gene))

    expect_error(
        miREA(
            methods = "TG_ORA", input_data = input_data, background = background,
            pvalueType = list(MiR = "std")
        ),
        "missing in pvalueType"
    )

    expect_error(
        miREA(methods = "TG_ORA", input_data = input_data, background = background, ncores = -1),
        "ncores must be a positive integer"
    )
    expect_error(
        miREA(methods = "TG_ORA", input_data = input_data, background = background, ncores = 1.5),
        "ncores must be a positive integer"
    )

    expect_message(
        res <- miREA(methods = "TG_ORA", input_data = input_data, background = background, ncores = NULL),
        "Use unparallel computing"
    )
    expect_true("TG_ORA" %in% res$valid_methods)
})

test_that("miREA() accepts a plain character vector background for a single target", {
    input_data <- getAllInputData(
        methods = "TG_ORA", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
        gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI
    )
    result <- miREA(
        methods = "TG_ORA", input_data = input_data,
        background = unique(toy_pathway$gene)
    )
    expect_true("TG_ORA" %in% result$valid_methods)
})

test_that("miREA() rejects a plain character vector background for multiple targets", {
    input_data <- getAllInputData(
        methods = c("TG_ORA", "MiR_ORA"), pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
        gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI
    )
    expect_error(
        miREA(
            methods = c("TG_ORA", "MiR_ORA"), input_data = input_data,
            background = unique(toy_pathway$gene)
        ),
        "please specify background for all the targets"
    )
})

test_that("miREA() continues running other methods after one method errors, and fills missing ncores entries", {
    input_data <- getAllInputData(
        methods = c("TG_ORA", "Edge_Score2D", "Edge_Topology"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    # deliberately break TG_ORA's input so it errors inside miREA()'s tryCatch
    input_data$data$TG_ORA <- character(0)
    background <- list(
        TG = unique(toy_pathway$gene),
        Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))
    )

    expect_warning(
        result <- miREA(
            methods = c("TG_ORA", "Edge_Score2D", "Edge_Topology"), input_data = input_data,
            background = background, pvalueType = list(TG = "std", Edge = "neg"),
            minSize = NULL, maxSize = NULL, iter = 5,
            ncores = list(Edge_Score2D = 2), pvalueCutoff = 1
        ),
        "didn't customized the ncores"
    )

    expect_false("TG_ORA" %in% result$valid_methods)
    expect_true(all(c("Edge_Score2D", "Edge_Topology") %in% result$valid_methods))

    expect_true("TG_ORA" %in% names(result$failed_methods))
    expect_type(result$failed_methods[["TG_ORA"]], "character")
    expect_false("TG_ORA" %in% names(result$result))
})
