#' Edge-Network (MGI Edge Network Enrichment Analysis)
#'
#' Perform miRNA-gene interaction (MGI) network enrichment analysis, based on
#' Random Walk with Restart (RWR), by evaluating how relevant the
#' differentially expressed MGI set is to the pathway MGI set.
#'
#' @param DEMGI A character vector contains a list of interested MGIs in the format of [miRNA]:[gene].
#' @param pathway_gene A two-column dataframe, including pathway and gene. Equals to input_data$pathway$TG.
#' @param pathway_MGI A three-column dataframe, including pathway, miRNA, gene. Equals to input_data$pathway$Edge.
#' @param pathway_GGI A three column dataframe, including pathway, from, to. Equals to input_data$pathway_GGI.
#' @param MGI A two-column dataframe, including MGIs across all DEMGIs and pathways built from background.
#' @param GGI A two-column dataframe, including GGIs across all DEMGIs and pathways, including internal GGIs and GGIs linking different pathways and DEMGIs.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default_pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If no adjustment is desired, set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values above the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs less than the minSize threshold.
#' @param maxSize The maximum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs more than the maxSize threshold.
#' @param iter Number of iterations for permutation test
#' @param ncores Number of cores for permutation test (parallel computing random RWR procedure).
#'        Ignored if `BPPARAM` is supplied.
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] instance controlling how the
#'        permutation test is parallelized. If NULL (the default), one is constructed
#'        automatically from `ncores` (serial if `ncores` is NULL/1, otherwise
#'        [BiocParallel::SnowParam()]).
#' @param seed Optional integer seed for the permutation test. If NULL (the
#'        default), no internal seed is set and the permutation draws depend
#'        on the caller's current random state. If supplied, the permutation
#'        is made reproducible via [withr::with_seed()], which does not alter
#'        the caller's global random state.
#' @return A dataframe containing the enrichment result for Edge-Network,
#'        including `nMoreExtreme` (the number of permutation scores at
#'        least as extreme as the observed score, out of `iter`), `p_value`
#'        (a normal-approximation p-value from the permutation null
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
#'     methods = "Edge_Network", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
#'     gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI,
#'     GGI_source = "Reactome",
#'     gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
#' )
#' res <- Edge_Network(
#'     DEMGI = input_data$data$Edge_Network, pathway_gene = input_data$pathway$TG,
#'     pathway_MGI = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
#'     MGI = input_data$MGI, GGI = input_data$GGI,
#'     minSize = 2, maxSize = 5000, iter = 10, ncores = 1, pvalueCutoff = 1
#' )
#' head(res)
#' @export
# pathway_GGI: only contain GGIs inside each pathway, pathway, from_gene, to_gene
# GGI: both GGIs inside each pathway, among different pathways, and for the DEMGI
# MGI: all MGIs in the network
Edge_Network <- function(DEMGI,
                         pathway_gene, pathway_MGI, pathway_GGI,
                         MGI, GGI,
                         pAdjMethod = "BH", pvalueCutoff = 0.05,
                         minSize = NULL, maxSize = NULL,
                         iter = 1000, ncores = NULL, BPPARAM = NULL, seed = NULL) {
    message("\n  Start Edge-Network analysis ...")
    if (is.null(DEMGI)) {
        stop("A valid character DEMGI vector contains the differential expressed MGI must be entered!")
    }

    if (!pAdjMethod %in% default_pAdjMethod) {
        stop(
            "Invalid pAdjMethod: ", pAdjMethod,
            "\n Please choose one of: ", paste(default_pAdjMethod, collapse = ",")
        )
    }

    prepped <- edge_network_prepare_pathway(pathway_gene, pathway_MGI, pathway_GGI, MGI, GGI, minSize, maxSize)
    pathway <- prepped$pathway
    pathway_gene <- prepped$pathway_gene
    pathway_GGI <- prepped$pathway_GGI
    MGI <- prepped$MGI
    GGI <- prepped$GGI
    MGISize <- prepped$MGISize

    n_pathway <- length(unique(pathway$pathway))
    if (n_pathway == 0) {
        warning("No pathways meet the size criteria.")
        return(data.frame(
            pathway = character(), n_MGI = integer(), score = numeric(), rand_score_dist = character(),
            z = numeric(), p_value = numeric(),
            nMoreExtreme = integer(), p_empirical = numeric(), padj = numeric(),
            contribute_edge = character(), n_contribute_edge = integer(),
            stringsAsFactors = FALSE
        ))
    }

    built <- edge_network_build_graph(pathway = pathway, pathway_gene = pathway_gene, MGI = MGI, GGI = GGI, DEMGI = DEMGI)

    restart_matrix <- edge_network_build_restart_matrix(
        edge_ids = built$edge_ids, all_edges = built$all_edges, pathway = pathway, DEMGI = DEMGI,
        DEMGI_out = built$DEMGI_out, DEMGI_in = built$DEMGI_in, DEMGI_connect = built$DEMGI_connect,
        iter = iter, seed = seed
    )

    if (is.null(BPPARAM) && (is.null(ncores) || ncores == 1)) {
        ncores <- 1
        message("You are use unparallel computing. It could cause huge amount of time.
            \nFor parallel computing, please specify [ncores] parameter to be larger than 1 (We highly recommend that).\n")
    }
    message("    Build transition matrix...")
    transition_matrix <- edge_network_build_transition_matrix(graph = built$graph, all_edges = built$all_edges)

    bpparam <- resolve_BPPARAM(ncores, BPPARAM = BPPARAM)
    message("    Execute RWR with ", iter, " iterations on ", BiocParallel::bpnworkers(bpparam), " cores...")
    rwr <- edge_network_run_rwr(transition_matrix = transition_matrix, restart_matrix = restart_matrix, BPPARAM = bpparam)

    message("    Summarise by pathways...")
    scored <- edge_network_pathway_scores(pathway = pathway, pathway_GGI = pathway_GGI, obs_score = rwr$obs_score, perm_scores = rwr$perm_scores)

    final <- edge_network_finalize_results(
        obs_result = scored$obs_result, perm_result = scored$perm_result, contributing_mgi = scored$contributing_mgi,
        MGISize = MGISize, pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, n_iter = scored$n_iter
    )

    message("  Edge-Network analysis have finished!")
    return(final)
}

#' @noRd
edge_network_prepare_pathway <- function(pathway_gene, pathway_MGI, pathway_GGI, MGI, GGI, minSize, maxSize) {
    # filtering GGI
    colnames(GGI)[seq_len(2)] <- c("from", "to")
    colnames(pathway_GGI)[seq_len(3)] <- c("pathway", "from", "to")
    colnames(pathway_gene)[seq_len(2)] <- c("pathway", "gene")

    # filtering pathway
    colnames(pathway_MGI)[seq_len(3)] <- c("pathway", "miRNA", "gene")
    colnames(MGI)[seq_len(2)] <- c("miRNA", "gene")
    pathway <- pathway_MGI %>%
        mutate(MGI = paste0(miRNA, ":", gene))
    MGISize <- pathway_MGI %>%
        group_by(pathway) %>%
        summarise(n_MGI = n()) %>%
        ungroup()

    pathway <- filter_pathway_by_size(pathway, minSize, maxSize)

    list(
        pathway = pathway, pathway_gene = pathway_gene, pathway_GGI = pathway_GGI,
        MGI = MGI, GGI = GGI, MGISize = MGISize
    )
}

#' @noRd
edge_network_build_graph <- function(pathway, pathway_gene, MGI, GGI, DEMGI) {
    ### build global graph
    MGI$type <- "MGI"
    colnames(MGI)[seq_len(3)] <- c("from", "to", "type")
    GGI$type <- "GGI"
    all_edges <- dplyr::distinct(rbind(MGI, GGI)) %>%
        dplyr::mutate(edge = paste0(from, ":", to))
    all_edges$weight <- 1 ### set weight as 1

    graph <- igraph::graph_from_data_frame(all_edges, directed = TRUE)
    edge_ids <- paste0(all_edges$from, ":", all_edges$to)
    names(edge_ids) <- seq_along(edge_ids)
    orig_edges_df <- igraph::as_data_frame(graph, what = "edges")
    E(graph)$name <- paste0(orig_edges_df$from, ":", orig_edges_df$to)

    # Divide DEMGI into two groups: isolated/non-isolated
    isolated_edges <- base::setdiff(DEMGI, pathway$MGI)
    isolated_df <- data.frame(
        miRNA = sub(":.*", "", isolated_edges),
        gene = sub(".*:", "", isolated_edges),
        stringsAsFactors = FALSE
    )
    pathway_genes <- unique(pathway_gene$gene)
    connect_genes <- unique(orig_edges_df$from[
        orig_edges_df$type == "GGI" & orig_edges_df$to %in% pathway_genes
    ])

    DEMGI_out <- isolated_df %>%
        filter(!(gene %in% connect_genes)) %>%
        mutate(MGI = paste0(miRNA, ":", gene)) %>%
        pull(MGI)
    DEMGI_connect <- setdiff(isolated_edges, DEMGI_out)
    DEMGI_in <- setdiff(DEMGI, isolated_edges)
    if (length(DEMGI_in) + length(DEMGI_connect) == 0) {
        stop("All DEMGIs are not included in any pathway. Edge_Network doesn't work in this situation. Please try other methods.\n")
    }

    list(
        graph = graph, all_edges = all_edges, edge_ids = edge_ids,
        DEMGI_out = DEMGI_out, DEMGI_in = DEMGI_in, DEMGI_connect = DEMGI_connect
    )
}

#' @noRd
edge_network_build_restart_matrix <- function(edge_ids, all_edges, pathway, DEMGI,
                                              DEMGI_out, DEMGI_in, DEMGI_connect,
                                              iter, seed) {
    message("    Build restart list...")
    restart_list <- as.integer(edge_ids %in% DEMGI)
    names(restart_list) <- edge_ids

    num_DE_in <- length(DEMGI_in)
    num_DE_connect <- length(DEMGI_connect)
    MGI_out_indices <- which(all_edges$type == "MGI" & edge_ids %in% DEMGI_out) # fixed as 1
    MGI_pw_indices <- which(all_edges$type == "MGI" & edge_ids %in% pathway$MGI) # all MGIs that in the randomization range
    MGI_connect_indices <- setdiff(which(all_edges$type == "MGI"), c(MGI_out_indices, MGI_pw_indices))

    draw_random_restart <- function(i) {
        new_restart <- rep(0, length(edge_ids))
        new_restart[MGI_out_indices] <- 1 # fixed DEMGI that outside the pathways as 1
        selected_pw_indices <- sample(MGI_pw_indices, num_DE_in, replace = FALSE)
        selected_connect_indices <- sample(MGI_connect_indices, num_DE_connect, replace = FALSE)
        selected_indices <- c(selected_pw_indices, selected_connect_indices)
        new_restart[selected_indices] <- 1
        new_restart
    }
    random_restarts <- if (is.null(seed)) {
        vapply(seq_len(iter), draw_random_restart, numeric(length(edge_ids)))
    } else {
        withr::with_seed(seed, vapply(seq_len(iter), draw_random_restart, numeric(length(edge_ids))))
    }

    restart_matrix <- cbind(restart_list, random_restarts)
    colnames(restart_matrix) <- c("OBS", paste0("P", seq(iter)))

    restart_matrix <- as(Matrix::Matrix(restart_matrix, sparse = TRUE), "dgCMatrix")
    message("      Size of sparse restart_matrix: ", object.size(restart_matrix) / (1024 * 1024), " MB.")
    restart_matrix
}

#' @noRd
edge_network_build_transition_matrix <- function(graph, all_edges) {
    lg <- igraph::make_line_graph(graph)
    V(lg)$name <- E(graph)$name
    edge_adjM <- as_adjacency_matrix(lg, sparse = TRUE)
    row_sums <- Matrix::rowSums(edge_adjM)
    row_sums[row_sums == 0] <- 1
    D_inv <- Matrix::Diagonal(x = 1 / row_sums)
    transition_matrix <- D_inv %*% edge_adjM
    transition_matrix <- Matrix::t(transition_matrix) # each row is a target edge, each column is a source edge, so each colsum should be 1 or 0
    colnames(transition_matrix) <- rownames(transition_matrix)

    # add self-loop for GGIs who don't have downstream edge
    ggi_nodes <- all_edges$edge[all_edges$type == "GGI"]
    cols_to_fix <- which(colnames(transition_matrix) %in% ggi_nodes & Matrix::colSums(transition_matrix) == 0)

    if (length(cols_to_fix) > 0) {
        diag_fix <- Matrix::sparseMatrix(
            i = cols_to_fix,
            j = cols_to_fix,
            x = rep(1, length(cols_to_fix)),
            dims = dim(transition_matrix),
            dimnames = dimnames(transition_matrix)
        )
        transition_matrix <- transition_matrix + diag_fix
    }
    message("      Size of transition_matrix: ", object.size(transition_matrix) / (1024 * 1024), " MB.")
    gc()
    transition_matrix
}

#' @noRd
edge_network_make_rwr <- function(transition_matrix) {
    function(seed.vec) {
        alpha <- 0.5
        delta <- 1e-6
        max_steps <- 100
        p <- seed.vec
        new_p <- p # pre-allocate, to avoid repeatedly creating the object
        diff_val <- Inf
        step <- 0

        repeat {
            tmp <- transition_matrix %*% p
            tmp@x <- tmp@x * (1 - alpha)
            new_p <- tmp
            new_p <- new_p + alpha * seed.vec

            # check difference
            diff_val <- Matrix::colSums(abs(new_p - p))
            step <- step + 1

            if (diff_val < delta || step >= max_steps) break

            # swap references instead of copying (in-place reference swap avoids holding two large objects)
            tmp <- p
            p <- new_p
            new_p <- tmp
        }
        p
    }
}

#' @noRd
edge_network_run_rwr <- function(transition_matrix, restart_matrix, BPPARAM) {
    run_rwr <- edge_network_make_rwr(transition_matrix)

    obs_vec <- restart_matrix[, 1, drop = FALSE] # drop = FALSE
    perm_mat <- restart_matrix[, -1, drop = FALSE]
    n_cols <- ncol(perm_mat)

    obs_score <- run_rwr(obs_vec) # calculate the observation scores

    # divide the permutation columns into one block per worker, where each
    # column of a block represents a single permutation
    ncores <- BiocParallel::bpnworkers(BPPARAM)
    avg <- ceiling(n_cols / ncores)
    block_indices <- split(seq_len(n_cols), rep(seq_len(ncores), each = avg, length.out = n_cols))

    perm_scores_list <- BiocParallel::bplapply(block_indices, function(idxs) {
        do.call(cbind, lapply(idxs, function(i) run_rwr(perm_mat[, i, drop = FALSE])))
    }, BPPARAM = BPPARAM)
    perm_scores <- do.call(cbind, perm_scores_list)

    list(obs_score = obs_score, perm_scores = perm_scores)
}

#' @noRd
edge_network_pathway_scores <- function(pathway, pathway_GGI, obs_score, perm_scores) {
    pathway <- pathway %>% mutate(edge = paste0(miRNA, ":", gene))
    pathway_GGI <- pathway_GGI %>% mutate(edge = paste0(from, ":", to))
    pw_all <- rbind(pathway[, c("pathway", "edge")], pathway_GGI[, c("pathway", "edge")])
    pw_dt <- as.data.table(pw_all)

    obs_dt <- data.table(edge = rownames(obs_score), score = as.numeric(obs_score))
    setkey(obs_dt, edge)
    obs_result <- pw_dt[, .(score = sum(obs_dt[edge, score])), by = pathway]

    # get top contributed MGIs
    edge_score <- merge(pw_dt, obs_dt, by = "edge", all.x = TRUE)
    contributing_mgi <- edge_score %>%
        group_by(pathway) %>%
        arrange(desc(score)) %>%
        summarise(
            contribute_edge = paste(edge[score > 0], collapse = "/"),
            n_contribute_edge = sum(score > 0),
            .groups = "drop"
        )

    # get permutation score for each pathway under each iteration
    edges <- unique(pw_dt$edge)
    pathways <- unique(pw_dt$pathway)
    P <- sparseMatrix(
        i = match(pw_dt$pathway, pathways),
        j = match(pw_dt$edge, edges),
        x = 1,
        dims = c(length(pathways), length(edges)),
        dimnames = list(pathways, edges)
    )

    # matrix-multiply across all permutation columns at once
    perm_scores <- perm_scores[match(colnames(P), rownames(perm_scores)), , drop = FALSE]
    perm_result_mat <- P %*% perm_scores # pathways x iterations

    # convert back to data.table
    # Matrix::summary() only returns the stored (non-zero) entries of the sparse
    # product; pathway/iteration combinations whose propagated score is exactly
    # zero are legitimately missing from perm_coo, not numerically negligible.
    # edge_network_finalize_results() accounts for this explicitly via n_iter,
    # rather than assuming perm_result contains every (pathway, iter) pair.
    perm_coo <- Matrix::summary(perm_result_mat) # i,j,x
    perm_result <- data.table(
        pathway = rownames(perm_result_mat)[perm_coo$i],
        iter = colnames(perm_result_mat)[perm_coo$j],
        score = perm_coo$x
    )

    list(
        obs_result = obs_result, contributing_mgi = contributing_mgi, perm_result = perm_result,
        n_iter = ncol(perm_result_mat)
    )
}

#' @noRd
edge_network_finalize_results <- function(obs_result, perm_result, contributing_mgi,
                                          MGISize, pAdjMethod, pvalueCutoff, n_iter) {
    # perm_result only stores the non-zero permutation scores (see
    # edge_network_pathway_scores()); reconstruct each pathway's true null
    # distribution (length n_iter, zero-padded) one pathway at a time so the
    # mean/sd/quantiles and nMoreExtreme count the missing zeros correctly,
    # without ever materializing the full pathway-by-iteration matrix.
    all_pathways <- obs_result$pathway
    score_list <- split(perm_result$score, perm_result$pathway)
    obs_score_by_pathway <- setNames(obs_result$score, obs_result$pathway)

    pathway_stats <- rbindlist(lapply(all_pathways, function(pw) {
        stored <- score_list[[pw]]
        if (is.null(stored)) stored <- numeric(0)
        full <- c(stored, rep(0, n_iter - length(stored)))
        obs <- obs_score_by_pathway[[pw]]

        sd_perm <- sd(full)
        mu_perm <- mean(full)
        if (is.na(sd_perm) || sd_perm == 0) {
            z <- NA_real_
            p_value <- 1.0000
        } else {
            z <- (obs - mu_perm) / sd_perm
            p_value <- pnorm(z, lower.tail = FALSE)
        }
        nMoreExtreme <- sum(full >= obs)
        p_empirical <- (nMoreExtreme + 1) / (n_iter + 1)

        rand_score_dist <- paste(
            sprintf("%.4f", min(full)),
            sprintf("%.4f", quantile(full, 0.25)),
            sprintf("%.4f", median(full)),
            sprintf("%.4f", mu_perm),
            sprintf("%.4f", quantile(full, 0.75)),
            sprintf("%.4f", max(full)),
            sep = ","
        )
        data.table(
            pathway = pw, rand_score_dist = rand_score_dist,
            z = z, p_value = p_value, nMoreExtreme = nMoreExtreme, p_empirical = p_empirical
        )
    }))

    final <- merge(obs_result, pathway_stats, by = "pathway")
    final[, padj := p.adjust(p_value, method = pAdjMethod)]
    final <- merge(final, contributing_mgi, by = "pathway", all.x = TRUE)
    final <- merge(MGISize, final, by = "pathway")
    final <- final %>%
        dplyr::filter(p_value <= pvalueCutoff, padj <= pvalueCutoff) %>%
        dplyr::arrange(padj) %>%
        dplyr::relocate(rand_score_dist, .after = score)
    final
}
