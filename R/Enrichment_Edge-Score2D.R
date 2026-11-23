#' Edge-Score2D (MiRNA-Gene Interaction 2D Edge Score-based Analysis)
#'
#' Perform miRNA-gene interaction (MGI) enrichment analysis based on a two-dimensional
#' edge score (expression correlation and differential-expression-informed ratio),
#' using a permutation test against randomized correlation values.
#'
#' @param MGI2D A five-column dataframe, including MGI miRNA gene cor(x) ranknormratio(y).
#' @param pathway A three-column dataframe contains pathway miRNA gene.
#' @param background A character vector contains background MGIs in the format of [miRNA]:[gene].
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default_pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If no adjustment is desired, set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values above the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs less than the minSize threshold.
#' @param maxSize The maximum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs more than the maxSize threshold.
#' @param iter Number of iterations for permutation test
#' @param ncores Number of cores for calculating randomized r during the permutation test.
#'        Ignored if `BPPARAM` is supplied.
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] instance controlling how the
#'        permutation test is parallelized. If NULL (the default), one is constructed
#'        automatically from `ncores` and `seed` (serial if `ncores` is NULL/1, otherwise
#'        [BiocParallel::SnowParam()]).
#' @param seed Optional integer seed for the permutation test. If NULL (the
#'        default), no internal seed is set and the permutation draws depend
#'        on the caller's current random state. If supplied, the permutation
#'        is made reproducible via [withr::with_seed()] when running serially,
#'        or via `BPPARAM`'s `RNGseed` when `ncores > 1`/a parallel `BPPARAM`
#'        is used; neither alters the caller's global random state.
#'
#' @return a dataframe contains enrichment result for Edge-Score2D, including
#'        `nMoreExtreme` (the number of permutation scores at least as
#'        extreme as the observed score, out of `iter`), `p_value` (a
#'        normal-approximation p-value from the permutation null
#'        distribution), and `p_empirical` (the empirical permutation
#'        p-value, `(nMoreExtreme + 1) / (iter + 1)`).
#' @examples
#' data(toy_pathway)
#' data(toy_mir_DEdata)
#' data(toy_gene_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
#' data(background_MGI)
#' input_data <- getAllInputData(
#'     methods = "Edge_Score2D", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
#'     gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI,
#'     gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
#' )
#' background <- unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))
#' res <- Edge_Score2D(
#'     MGI2D = input_data$data$Edge_Score2D, pathway = input_data$pathway$Edge,
#'     background = background, iter = 100, ncores = 1, pvalueCutoff = 1
#' )
#' head(res)
#' @export
Edge_Score2D <- function(MGI2D, pathway, background = NULL, pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL, iter = 1000, ncores = NULL, BPPARAM = NULL, seed = NULL) {
    message("\n  Start Edge_Score2D analysis ...")

    MGI2D <- data.table::copy(data.table::as.data.table(MGI2D))
    MGI2D <- na.omit(MGI2D)
    if (is.null(MGI2D) || !is.data.frame(MGI2D)) {
        stop("Please make sure the input MGI2D is a dataframe!")
    }
    if (nrow(MGI2D) == 0) {
        stop("No valid MGI pairs in MGI2D!")
    }

    if (is.null(background)) {
        stop(
            "background must be supplied: a character vector of MGIs in '[miRNA]:[gene]' ",
            "format. miREA does not ship a bundled default background; see the package ",
            "vignette for how to build one from background_MGI, or download the full ",
            "background used in the miREA manuscript from its Zenodo deposit."
        )
    }

    # filtering MGI in the background
    if (all(MGI2D$MGI %in% background) == FALSE) {
        MGI2D <- MGI2D[MGI %in% background]
        if (nrow(MGI2D) == 0) {
            stop("No valid MGI input from background list!")
        }
    }
    colnames(MGI2D) <- c("MGI", "miRNA", "gene", "x", "y")

    r_all <- sqrt(((1 - MGI2D$x) / 2)^2 + (MGI2D$y)^2)


    # filtering pathway in the background
    pathway <- data.table::copy(data.table::as.data.table(pathway))
    pathway[, MGI := paste0(miRNA, ":", gene)]
    pathway <- pathway[MGI %in% MGI2D$MGI]
    pathway <- pathway[MGI2D, on = .(MGI), nomatch = 0L]


    ### filtering pathway size
    pathway <- filter_pathway_by_size(pathway, minSize, maxSize)

    n_pathway <- length(unique(pathway$pathway))
    if (n_pathway == 0) {
        warning("No pathways meet the size criteria.")
        return(data.frame(
            pathway = character(), MGISize = integer(),
            score = numeric(), rand_score_dist = character(),
            p_value = numeric(), padj = numeric(),
            nMoreExtreme = integer(), p_empirical = numeric(),
            TopMGI = character(), stringsAsFactors = FALSE
        ))
    }

    bpparam <- resolve_BPPARAM(ncores, seed = seed, BPPARAM = BPPARAM)
    serial <- is.null(BPPARAM) && (is.null(ncores) || ncores <= 1)

    # generate random pathways
    message("    Generate random MGI set radius on ", BiocParallel::bpnworkers(bpparam), " cores...")

    pathway_sizes <- table(pathway$pathway)
    unique_sizes <- sort(unique(pathway_sizes))

    nchunks <- min(BiocParallel::bpnworkers(bpparam), length(unique_sizes))
    size_groups <- split(
        unique_sizes,
        rep(seq_len(nchunks), length.out = length(unique_sizes))
    )

    draw_random_radius <- function(size_group) {
        sub_result <- list()
        for (size in size_group) {
            sub_result[[as.character(size)]] <- vapply(seq_len(iter), function(i) sum(sample(r_all, size)), numeric(1))
        }
        sub_result
    }

    # Serial execution runs in-process, so plain R-level seeding via
    # withr::with_seed() is sufficient and avoids depending on BiocParallel's
    # SerialParam RNG handling (observed to repeat draws across separate
    # calls on some platform/BiocParallel combinations even with no RNGseed
    # set). True parallel execution still needs BiocParallel's own RNGseed,
    # since workers don't inherit this process's RNG stream.
    r_random_list <- if (serial) {
        if (is.null(seed)) {
            lapply(size_groups, draw_random_radius)
        } else {
            withr::with_seed(seed, lapply(size_groups, draw_random_radius))
        }
    } else {
        BiocParallel::bplapply(size_groups, draw_random_radius, BPPARAM = bpparam)
    }
    r_random_list <- unlist(r_random_list, recursive = FALSE)
    names(r_random_list) <- sub("^\\d+\\.", "", names(r_random_list))
    r_random_list <- lapply(r_random_list, sort)

    message("    Calculate pathway radius and significant level...")
    result <- Edge_Score2D_simple(MGI2D = MGI2D, r_random_list = r_random_list, pathway = pathway, iter = iter)
    result <- result %>%
        dplyr::mutate(padj = stats::p.adjust(p_value, method = pAdjMethod)) %>%
        dplyr::relocate(padj, .after = p_value) %>%
        dplyr::filter(
            p_value <= pvalueCutoff,
            padj <= pvalueCutoff
        ) %>%
        arrange(padj)

    message("  Edge-Score2D analysis have finished!")
    return(result)
}


Edge_Score2D_simple <- function(MGI2D, r_random_list, pathway, iter) {
    pathway_names <- unique(pathway$pathway)
    pathway_list <- split(pathway, pathway$pathway)
    compute_pathway <- function(i) {
        pathway_name <- pathway_names[i]
        pw_single <- pathway_list[[pathway_name]]
        pw_size <- nrow(pw_single)
        r_random <- r_random_list[[as.character(pw_size)]]

        result <- Edge_Score2D_single(MGI2D = MGI2D, PW = pw_single, r_random = r_random)

        message("      ", i, " / ", length(pathway_names), " , ", pathway_name, " : ", round(result$p_value, 4))

        return(as.data.frame(result))
    }

    pathway_result <- bind_rows(lapply(seq_len(length(pathway_names)), compute_pathway))
    return(pathway_result)
}


Edge_Score2D_single <- function(MGI2D, PW, r_random) {
    MGI2D <- na.omit(MGI2D)
    if (is.null(MGI2D) || !is.data.frame(MGI2D)) {
        stop("Please make sure the input MGI2D is a dataframe!")
    }
    if (nrow(MGI2D) == 0) {
        stop("No valid MGI pairs in MGI2D!")
    }
    if (length(unique(PW$pathway)) != 1) {
        warning("The pathway is not calculated one-by-one, please check the PW. Calculate the first one instead.")
        name <- unique(PW$pathway)[1]
        PW <- PW %>% dplyr::filter(pathway == name)
    }

    points_data <- semi_join(MGI2D, PW, by = c("miRNA", "gene"))

    if (nrow(points_data) == 0) {
        return(NULL)
    }


    points_data <- points_data %>%
        mutate(r = sqrt(((1 - x) / 2)^2 + y^2))

    r <- sum(points_data$r)
    count_leq <- findInterval(r, r_random) # count <= r
    nMoreExtreme <- length(r_random) - count_leq
    p_empirical <- (nMoreExtreme + 1) / (length(r_random) + 1)

    sd_r_random <- sd(r_random)
    if (is.na(sd_r_random) || sd_r_random == 0) {
        p_value <- 1.0000
    } else {
        z <- (r - mean(r_random)) / sd_r_random
        p_value <- pnorm(z, lower.tail = FALSE)
    }


    threshold <- quantile(points_data$r, 0.9)
    top_pairs <- points_data %>%
        filter(r >= threshold) %>%
        mutate(MGI = paste0(miRNA, ":", gene)) %>%
        pull(MGI)

    ## permutation test
    PW_name <- unique(PW$pathway)
    MGI_list <- PW$MGI
    size <- nrow(PW)

    random_ds <- round(quantile(r_random, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 4)
    random_mean <- round(mean(r_random, na.rm = TRUE), 4)
    random_ds <- c(
        min = random_ds[1],
        q25 = random_ds[2],
        median = random_ds[3],
        mean = random_mean,
        q75 = random_ds[4],
        max = random_ds[5]
    )

    result <- data.frame(
        "pathway" = PW_name,
        "MGISize" = size,
        "score" = r,
        "rand_score_dist" = paste(random_ds, collapse = ","),
        "p_value" = p_value,
        "nMoreExtreme" = nMoreExtreme,
        "p_empirical" = p_empirical,
        "TopMGI" = paste(top_pairs, collapse = "/"),
        # "rand_score_dist" = paste(round(r_random, 4), collapse = ","),
        stringsAsFactors = FALSE
    )
    return(result)
}
