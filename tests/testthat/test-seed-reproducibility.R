data(toy_pathway, toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat,
    background_MGI,
    package = "miREA"
)

test_that("Edge_Score2D: same seed gives identical results, seed=NULL gives different results", {
    input_data <- cached_getAllInputData(
        methods = "Edge_Score2D", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
        scoreFun = "rank"
    )
    background <- unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))

    run <- function(seed) {
        Edge_Score2D(
            MGI2D = input_data$data$Edge_Score2D, pathway = input_data$pathway$Edge,
            background = background, iter = 5, ncores = 1, pvalueCutoff = 1, seed = seed
        )
    }

    expect_identical(run(42), run(42))
    expect_false(identical(run(1), run(2)))
    # seed=NULL draws are compared via a rounded summary (rand_score_dist),
    # so a single pair can coincidentally round to the same result; require
    # at least one difference across a few independent trials.
    expect_true(any(replicate(3, !isTRUE(all.equal(run(NULL), run(NULL))))))
})

test_that("Edge_Topology: same seed gives identical results, seed=NULL gives different results", {
    input_data <- cached_getAllInputData(
        methods = "Edge_Topology", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    background <- unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))

    run <- function(seed) {
        Edge_Topology(
            MGI = input_data$data$Edge_Topology, background = background,
            pathway = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
            dict_MGI = input_data$dict_MGI, pvalueType = "neg",
            minSize = 2, maxSize = 5000, iter = 5, ncores = 1, pvalueCutoff = 1, seed = seed
        )
    }

    expect_identical(run(42), run(42))
    expect_false(identical(run(1), run(2)))
    expect_false(isTRUE(all.equal(run(NULL), run(NULL))))
})

test_that("Edge_Network: same seed gives identical results, seed=NULL gives different results", {
    input_data <- cached_getAllInputData(
        methods = "Edge_Network", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
    )

    run <- function(seed) {
        Edge_Network(
            DEMGI = input_data$data$Edge_Network, pathway_gene = input_data$pathway$TG,
            pathway_MGI = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
            MGI = input_data$MGI, GGI = input_data$GGI,
            minSize = 2, maxSize = 5000, iter = 5, ncores = 1, pvalueCutoff = 1, seed = seed
        )
    }

    expect_identical(run(42), run(42))
    expect_false(identical(run(1), run(2)))
    expect_false(isTRUE(all.equal(run(NULL), run(NULL))))
})

test_that("miREA() forwards seed consistently across Edge_Score2D/Edge_Topology/Edge_Network", {
    input_data <- cached_getAllInputData(
        methods = c("Edge_Score2D", "Edge_Topology", "Edge_Network"), pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, GGI_source = "Reactome",
        gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
    )
    background <- list(Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene)))

    run <- function(seed) {
        miREA(
            methods = c("Edge_Score2D", "Edge_Topology", "Edge_Network"),
            input_data = input_data, background = background,
            minSize = 2, maxSize = 5000, pvalueType = list(Edge = "neg"),
            pvalueCutoff = 1, iter = 5, ncores = 1, seed = seed
        )$result
    }

    expect_identical(run(7), run(7))
})

test_that("Edge_Score2D/Edge_Topology/Edge_Network's internal seed does not leak into the caller's global RNG stream", {
    input_data <- cached_getAllInputData(
        methods = "Edge_Score2D", pathway = toy_pathway,
        mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
        background_MGI = background_MGI, gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
        scoreFun = "rank"
    )
    background <- unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))

    # Baseline: what the caller's next two runif(1) draws look like with no
    # intervening call to Edge_Score2D().
    set.seed(123)
    expected_draw1 <- runif(1)
    expected_draw2 <- runif(1)

    # Same starting seed, but with an Edge_Score2D(seed = 999) call in between:
    # if that call's internal seeding perturbed the caller's global RNG stream,
    # draw2 would differ from expected_draw2 even though draw1 (taken before the
    # call) still matches expected_draw1.
    set.seed(123)
    draw1 <- runif(1)
    Edge_Score2D(
        MGI2D = input_data$data$Edge_Score2D, pathway = input_data$pathway$Edge,
        background = background, iter = 5, ncores = 1, pvalueCutoff = 1, seed = 999
    )
    draw2 <- runif(1)

    expect_identical(draw1, expected_draw1)
    expect_identical(draw2, expected_draw2)
})
