#' Edge-Topology (MGI Edge Pathway Topology Enrichment Analysis)
#'
#' Perform miRNA-gene interaction (MGI) network enrichment analysis, based on
#' internal pathway topological importance and edge score importance, by
#' evaluating whether the MGIs with stronger regulatory strength also have
#' more topological importance inside the pathway.
#'
#' @param MGI A four-column dataframe, including MGI miRNA gene strength.
#' @param pvalueType How to calculate p value. See default_pvalueType for all possible options. Default pvalueType is "pos".
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negative ("neg") enrichment).
#' The "std" is used when interested in both positive and negative scores.
#' @param background A character vector contains background MGI in the format of [miRNA]:[gene].
#' @param pathway A three-column dataframe, including pathway, miRNA, gene. Equals to input_data$pathway$Edge.
#' @param pathway_GGI A three column dataframe, including pathway, from, to. Equals to input_data$pathway_GGI.
#' @param dict_MGI A four-column dataframe, including MGI miRNA gene strength. Containing all possible MGIs that could be used in permutation test. Equals to input_data$dict_MGI.
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values above the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs less than the minSize threshold.
#' @param maxSize The maximum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs more than the maxSize threshold.
#' @param iter Maximum number of permutations for the permutation test. The
#'        actual number of permutations run for a given pathway can be lower:
#'        internally, permutations are drawn in up to three sequential
#'        stages (roughly 10%, 50%, and 100% of `iter`), stopping early once
#'        enough extreme permutation scores have been observed to settle the
#'        p-value. The number of permutations actually used for a pathway is
#'        reported in the `nPerm` column of the result.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default_pAdjMethod to find all possible values.
#' @param ncores Number of cores for permutation test (parallel computing different pathways).
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
#' @return A dataframe containing the enrichment result for Edge-Topology,
#'        including `nPerm` (the actual number of permutations run for that
#'        pathway, see `iter`), `nMoreExtreme` (the number of those
#'        permutations at least as extreme as the observed score), `p_value`
#'        (a normal-approximation p-value from the permutation null
#'        distribution), and `p_empirical` (the empirical permutation
#'        p-value, `(nMoreExtreme + 1) / (nPerm + 1)`).
#' @examples
#' data(toy_pathway); data(toy_mir_DEdata); data(toy_gene_DEdata)
#' data(toy_gene_mat); data(toy_mir_mat)
#' data(toy_background_MGI); data(toy_background_GGI)
#' input_data <- getAllInputData(
#'   methods = "Edge_Topology", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
#'   gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI,
#'   background_GGI = toy_background_GGI, GGI_source = "Omnipath",
#'   gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
#' )
#' background <- unique(paste0(toy_background_MGI$miRNA, ":", toy_background_MGI$gene))
#' res <- Edge_Topology(
#'   MGI = input_data$data$Edge_Topology, background = background,
#'   pathway = input_data$pathway$Edge, pathway_GGI = input_data$pathway_GGI,
#'   dict_MGI = input_data$dict_MGI, pvalueType = "neg",
#'   minSize = 2, maxSize = 5000, iter = 20, ncores = 1, pvalueCutoff = 1
#' )
#' head(res)
#' @export
Edge_Topology <- function(MGI, pvalueType = "neg",
                           background, pathway, pathway_GGI, dict_MGI,
                           iter = 1000, pAdjMethod = "BH", pvalueCutoff = 0.05,
                           minSize = NULL, maxSize = NULL, ncores = NULL, BPPARAM = NULL, seed = NULL){
  message("\n  Start Edge-Topology analysis ...")
  # check pAdjMethod
  if (!pAdjMethod %in% default_pAdjMethod){
    stop("Invalid pAdjMethod: ", pAdjMethod,
         "\n Please choose one of: ", paste(default_pAdjMethod, collapse = ", "))
  }
  # check pvalueType
  if (!pvalueType %in% default_pvalueType){
    stop("Invalid pvalueType: ", pvalueType,
         "\n Please choose one of: ", paste(default_pvalueType, collapse = ", "))
  }

  # filtering MGI in the background
  if (all(MGI$MGI %in% background) == FALSE) {
    MGI <- MGI %>% dplyr::filter(MGI %in% background)
    if (nrow(MGI) == 0) {
      stop("No valid MGI input from background list!")
    }
  }

  # pathway
  colnames(pathway)[seq_len(3)] <- c("pathway", "miRNA", "gene")
  colnames(pathway_GGI)[seq_len(3)] <- c("pathway", "from", "to")
  ## filtering pathway in the background
  pathway <- pathway %>%
    mutate(MGI = paste0(miRNA,":",gene)) %>%
    filter(MGI %in% background)
  ## filtering pathway size
  pathway_counts <- table(pathway$pathway)

  # NOTE: `pathway` (the data.table) has a column also named "pathway", so
  # `pathway$pathway` inside `pathway[...]` gets NSE-resolved against the
  # *column* (an atomic vector) first, raising "$ operator is invalid for
  # atomic vectors". Use the bare column reference instead.
  if (!is.null(minSize)){
    pathways_to_remove_min <- names(pathway_counts[pathway_counts < minSize])
    pathway <- pathway[!pathway %in% pathways_to_remove_min, ]
  }
  if (!is.null(maxSize)){
    pathways_to_remove_max <- names(pathway_counts[pathway_counts > maxSize])
    pathway <- pathway[!pathway %in% pathways_to_remove_max, ]
  }

  n_pathway <- length(unique(pathway$pathway))
  if (n_pathway == 0) {
    warning("No pathways meet the size criteria.")
    return(data.frame(
      pathway = character(), n_MGI = integer(),
      score = numeric(), rand_score_dist = character(),
      p_value = numeric(), padj = numeric(),
      nPerm = integer(), nMoreExtreme = integer(), p_empirical = numeric(),
      TopMGI = character(),
      stringsAsFactors = FALSE
    ))}


  # check dict
  dict_MGI <- data.table::as.data.table(dict_MGI)
  dict_MGI <- dict_MGI[, .(MGI, miRNA, gene, strength)]

  message("    Start analysis on ", n_pathway, " pathways...")
  serial <- is.null(BPPARAM) && (is.null(ncores) || ncores == 1)
  if (serial){
    ncores <- 1
    message("You are use unparallel computing. It could cause huge amount of time.
            \nFor parallel computing, please specify [ncores] parameter to be larger than 1 (We highly recommend that).\n")
  }
  message("    Calculate pathway scores and significant level on ",
          if (serial) 1 else BiocParallel::bpnworkers(resolve_BPPARAM(ncores, BPPARAM = BPPARAM)), " cores...")
  if (serial){
    result <- Edge_Topology_simple(MGI = MGI, dict_MGI = dict_MGI, pathway = pathway, pathway_GGI = pathway_GGI, pvalueType = pvalueType, iter = iter, seed = seed) # background = background,
  } else {
    result <- Edge_Topology_parallel(MGI = MGI, dict_MGI = dict_MGI, pathway = pathway, pathway_GGI = pathway_GGI, pvalueType = pvalueType, iter = iter, ncores = ncores, BPPARAM = BPPARAM, seed = seed) # background = background,
  }
  result <- result %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = pAdjMethod)) %>%
    dplyr::relocate(padj, .after = p_value) %>%
    dplyr::filter(
      p_value <= pvalueCutoff,
      padj <= pvalueCutoff
    ) %>%
    arrange(padj)

  message("  Edge-Topology analysis have finished!")
  return(result)
}


Edge_Topology_simple <- function(MGI, dict_MGI, pathway, pathway_GGI, pvalueType, iter, seed = NULL){ # background,


  pathway_list <- split(pathway, pathway$pathway)
  pathway_GGI_list <- split(pathway_GGI, pathway_GGI$pathway)
  pathway_name <- unique(pathway$pathway)
  n_pathway <- length(pathway_name)

  compute_pathway <- function(i){
    name <- pathway_name[i]
    pw_single <- pathway_list[[name]]
    pw_GGI_single <- pathway_GGI_list[[name]]

    # Derive a pathway-specific seed from the pathway's own name (rather than
    # from i, which is only a chunk-local index under Edge_Topology_parallel)
    # so permutations are independent across pathways while remaining fully
    # reproducible given the top-level seed.
    pathway_seed <- if (is.null(seed)) NULL else seed + sum(utf8ToInt(name) * seq_along(utf8ToInt(name)))

    result <- Edge_Topology_single(MGI = MGI, PW = pw_single, PW_GGI = pw_GGI_single, dict_MGI = dict_MGI, pvalueType = pvalueType, iter = iter, seed = pathway_seed) # background = background,

    message("      ", i, " / ", n_pathway, " , ", name, " : ", round(result$p_value, 4))

    return(result)
  }
  pathway_result <- dplyr::bind_rows(lapply(seq_len(n_pathway), compute_pathway))
  return(pathway_result)
}

## Edge_Topology_parallel
Edge_Topology_parallel <- function(MGI, dict_MGI, pathway, pathway_GGI, pvalueType, iter, ncores, BPPARAM = NULL, seed = NULL){ # background,
  bpparam <- resolve_BPPARAM(ncores, BPPARAM = BPPARAM)
  nchunks <- BiocParallel::bpnworkers(bpparam)

  pathway_counts <- table(pathway$pathway)
  pathway_counts <- sort(pathway_counts, decreasing = TRUE, na.last = NA)
  pathway_names_sorted <- names(pathway_counts)

  # initialize the task list of each chunk
  if (length(pathway_counts) < nchunks) {
    nchunks <- length(pathway_counts)
    warning("Since your task can only be divided into ", nchunks, " parts, modify the number of cores to ", nchunks, " \n")
  }
  task_list <- vector("list", nchunks)
  task_sizes <- rep(0, nchunks)  # record total pathway size of each chunk

  # greedy algorithm：distribute the largest pathway to the smallest load chunk
  for (i in seq_along(pathway_names_sorted)) {
    min_core <- which.min(task_sizes)
    task_list[[min_core]] <- c(task_list[[min_core]], pathway_names_sorted[i])
    task_sizes[min_core] <- task_sizes[min_core] + pathway_counts[i]
  }
  message(paste("      Total MGI size of each core is:", paste(task_sizes, collapse = " ")))

  results <- BiocParallel::bplapply(seq_len(nchunks), function(i) {
    pw <- pathway[pathway %in% task_list[[i]]]
    ggi <- pathway_GGI[pathway %in% task_list[[i]]]
    Edge_Topology_simple(MGI = MGI, dict_MGI = dict_MGI,
                          pathway = pw, pathway_GGI = ggi, # background,
                          pvalueType = pvalueType,
                          iter = iter, seed = seed)
  }, BPPARAM = bpparam)

  pathway_result <- combine_parallel(results)
  return(pathway_result)
}



Edge_Topology_single <- function(MGI, PW, PW_GGI, dict_MGI, pvalueType, iter, seed = NULL){
  if (is.null(PW_GGI)){
    PW_GGI <- data.table(pathway = character(), from = character(), to = character(), stringsAsFactors = FALSE)
  }
  if (length(unique(PW$pathway)) != 1){
    warning("The pathway is not calculated one-by-one, please check the PW. Calculate the first one instead.\n")
    name <- unique(PW$pathway)[1]
    PW <- PW %>% dplyr::filter(pathway == name)
    PW_GGI <- PW_GGI %>% dplyr::filter(pathway == name)
  }
  ## filter MGI and dict_MGI
  PW <- data.table::as.data.table(PW)
  MGI <- data.table::as.data.table(MGI)
  dict_MGI <- data.table::as.data.table(dict_MGI)

  data.table::setkey(PW, miRNA, gene, MGI)
  data.table::setkey(MGI, miRNA, gene, MGI)
  MGI <- MGI[PW, nomatch = 0]

  dict_MGI <- dict_MGI[miRNA %in% unique(PW$miRNA) & gene %in% unique(PW$gene)]
  dict_MGI <- dict_MGI[, setNames(strength, MGI)]

  PW_name <- unique(PW$pathway)
  mirnas <- unique(PW$miRNA)

  PW_edge <- rbind(
    data.table(from = PW$miRNA, to = PW$gene),
    data.table(from = PW_GGI$from, to = PW_GGI$to)
  )
  PW_edge <- unique(PW_edge)

  g <- igraph::graph_from_data_frame(d = PW_edge, directed = TRUE)
  edge_eb <- igraph::edge_betweenness(g)
  edge_dt <- data.table(igraph::as_data_frame(g, what = "edges"))[, EB := edge_eb]
  edge_dt[, EB := fifelse(is.na(EB), 0, EB)]
  edge_dt[, centrality := min_max(EB) + 0.01]

  # === score for original pathway === #
  PW_dt <- data.table(PW)
  setnames(PW_dt, c("miRNA", "gene"), c("from", "to"))
  PW_dt <- merge(PW_dt, edge_dt, by = c("from", "to"), all.x = TRUE)
  # PW_dt[, EB := fifelse(is.na(EB), 0, EB)]
  # PW_dt[, centrality := min_max(EB)]
  PW_dt[MGI, strength := i.strength, on = "MGI"]
  PW_dt[is.na(strength), strength := 0]
  PW_dt[, score := centrality * strength]
  PW_dt <- unique(PW_dt)
  PW_score <- sum(PW_dt$score, na.rm = TRUE)

  # === Random Permutation === #
  PW_GGI <- data.table(PW_GGI)
  g_MGI <- igraph::graph_from_data_frame(PW_dt[, .(from, to)], directed = TRUE) # original MGI network
  all_vertices <- unique(c(PW_GGI$from, PW_GGI$to, PW_dt$from, PW_dt$to))
  # GGI network with both MGI and GGI vertices (remain the same during permutation)
  basic_net <- igraph::graph_from_data_frame(PW_GGI[, .(from, to)], directed = TRUE, vertices = all_vertices)

  compute_random_score <- function(i){
    g_rand <- igraph::rewire(g_MGI, with = keeping_degseq(niter = ecount(g_MGI) * 50))
    g_combined <- igraph::union(g_rand, basic_net)

    rand_edge_eb <- igraph::edge_betweenness(g_combined)

    rand_MGI <- data.table(igraph::as_data_frame(g_rand, what = "edges"))
    rand_edge_dt <- data.table(igraph::as_data_frame(g_combined, what = "edges"))[, EB := rand_edge_eb]
    rand_edge_dt[, EB := fifelse(is.na(EB), 0, EB)]
    rand_edge_dt[, centrality := min_max(EB)+0.01]
    rand_PW_dt <- merge(rand_MGI, rand_edge_dt, by = c("from", "to"), all.x = TRUE)
    rand_PW_dt[, MGI := paste0(from, ":", to)]
    rand_PW_dt[, strength := dict_MGI[MGI]]
    rand_PW_dt[is.na(strength), strength := 0]
    rand_PW_dt[, score := centrality * strength]
    rand_PW_dt <- unique(rand_PW_dt)
    sum(rand_PW_dt$score, na.rm = TRUE)
  }

  # === Stepwise Permutation Process === #
  iter1 <- min(iter, max(10, ceiling(0.1 * iter)))
  iter2 <- min(iter, max(50, ceiling(0.5 * iter)))

  run_stepwise_permutation <- function(){
    s_random1 <- vapply(seq_len(iter1), compute_random_score, numeric(1))
    # n_extreme1 <- sum(s_random1 >= PW_score)
    n_up1 <- sum(s_random1 >= PW_score)
    n_down1 <- sum(s_random1 <= PW_score)

    n_extreme1 <- switch(
      pvalueType,
      "pos" = n_up1,
      "neg" = n_down1,
      "std" = min(n_up1, n_down1)
    )
    if (n_extreme1 > 0.5 * iter1){
      s_random <- s_random1
      n_up <- n_up1
      n_down <- n_down1
      nMoreExtreme <- n_extreme1
      nPerm <- iter1
    } else if (iter2 <= iter1) {
      s_random <- s_random1
      n_up <- n_up1
      n_down <- n_down1
      nMoreExtreme <- n_extreme1
      nPerm <- iter1
    } else {
      s_random2 <- vapply((iter1 + 1):iter2, compute_random_score, numeric(1))
      s_random <- c(s_random1, s_random2)
      n_up2 <- sum(s_random >= PW_score)
      n_down2 <- sum(s_random <= PW_score)
      n_extreme2 <- switch(
        pvalueType,
        "pos" = n_up2,
        "neg" = n_down2,
        "std" = min(n_up2, n_down2)
      )
      if (n_extreme2 > 0.1 * iter2 || iter <= iter2){
        n_up <- n_up2
        n_down <- n_down2
        nMoreExtreme <- n_extreme2
        nPerm <- iter2
      } else {
        s_random3 <- vapply((iter2 + 1):iter, compute_random_score, numeric(1))
        s_random <- c(s_random, s_random3)
        n_up <- sum(s_random >= PW_score)
        n_down <- sum(s_random <= PW_score)
        nMoreExtreme <- switch(
          pvalueType,
          "pos" = n_up,
          "neg" = n_down,
          "std" = min(n_up, n_down)
        )
        nPerm <- iter
      }
    }
    list(s_random = s_random, n_up = n_up, n_down = n_down, nMoreExtreme = nMoreExtreme, nPerm = nPerm)
  }

  stepwise_result <- if (is.null(seed)) {
    run_stepwise_permutation()
  } else {
    withr::with_seed(seed, run_stepwise_permutation())
  }
  s_random <- stepwise_result$s_random
  n_up <- stepwise_result$n_up
  n_down <- stepwise_result$n_down
  nMoreExtreme <- stepwise_result$nMoreExtreme
  nPerm <- stepwise_result$nPerm

  # z-score
  sd_s_random <- sd(s_random)
  if (is.na(sd_s_random) || sd_s_random == 0){
    p_value <- 1.0000
  } else {
    z <- (PW_score - mean(s_random)) / sd_s_random
    p_value <- switch(
      pvalueType,
      "pos" = pnorm(z, lower.tail = FALSE),
      "neg" = pnorm(z),
      "std" = 2 * pnorm(-abs(z))
    )
  }

  random_ds <- round(quantile(s_random, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 4)
  random_mean <- round(mean(s_random, na.rm = TRUE), 4)
  random_ds <- c(
    min = random_ds[1],
    q25 = random_ds[2],
    median = random_ds[3],
    mean = random_mean,
    q75 = random_ds[4],
    max = random_ds[5]
  )

  # calculate top contributed MGIs
  n_top <- max(1, floor(nrow(PW_dt) * 0.1))
  if (pvalueType == "pos" || (pvalueType == "std" && n_up < n_down)) {
    PW_leading_edge <- PW_dt[order(-score)][seq_len(n_top)][, paste0(from, ":", to)]
  } else {
    PW_leading_edge <- PW_dt[order(score)][seq_len(n_top)][, paste0(from, ":", to)]
  }

  p_empirical <- (nMoreExtreme + 1) / (nPerm + 1)

  result <- data.frame(
    pathway = PW_name,
    n_MGI = nrow(PW_dt),
    score = PW_score,
    rand_score_dist = paste(random_ds, collapse = ","),
    p_value = p_value,
    nPerm = nPerm,
    nMoreExtreme = nMoreExtreme,
    p_empirical = p_empirical,
    TopMGI = paste(PW_leading_edge, collapse = "/"),
    #"MGI_list" = paste(MGI_list, collapse = "/"),
    stringsAsFactors = FALSE
  )
  return(result)
}
