#' Title: Edge-Network (MGI Edge Network Enrichment Analysis that based on Random Walk with Restart (RWR))
#' Description: Perform MGI network enrichment analysis by evaluating how relevant the differential expressed MGI set is to the pathway MGI set based on RWR.
#' @param DEMGI A character vector contains a list of interested MGIs in the format of [miRNA]:[gene].
#' @param pathway_gene A two-column dataframe, including pathway and gene. Equals to input_data$pathway$TG.
#' @param pathway_MGI A three-column dataframe, including pathway, miRNA, gene. Equals to input_data$pathway$Edge.
#' @param pathway_GGI A three column dataframe, including pathway, from, to. Equals to input_data$pathway_GGI.
#' @param MGI A two-column dataframe, including MGIs across all DEMGIs and pathways built from background.
#' @param GGI A two-column dataframe, including GGIs across all DEMGIs and pathways, including internal GGIs and GGIs linking different pathways and DEMGIs.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default.pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If you don't want any adjustment, please set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values below the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs less than the minSize threshold.
#' @param maxSize The maximum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs more than the maxSize threshold.
#' @param iter Number of iterations for permutation test
#' @param ncores Number of cores for permutation test (parallel computing random RWR procedure).

#' @return A dataframe containing the enrichment result for Edge-Network.


# pathway_GGI: only contain GGIs inside each pathway, pathway, from_gene, to_gene
# GGI: both GGIs inside each pathway, among different pathways, and for the DEMGI
# MGI: all MGIs in the network
Edge_Network <- function(DEMGI,
                         pathway_gene, pathway_MGI, pathway_GGI,
                         MGI, GGI,
                         pAdjMethod = "BH", pvalueCutoff = 0.05,
                         minSize = NULL, maxSize = NULL,
                         iter = 1000, ncores = NULL){
  cat("\n  Start Edge-Network analysis ...\n")
  if (is.null(DEMGI)) {
    stop("A valid character DEMGI vector contains the differential expressed MGI must be entered!")
  }

  if (!pAdjMethod %in% default.pAdjMethod){
    stop(paste("Invalid pAdjMethod:", pAdjMethod,
               "\n Please choose one of:", paste(default.pAdjMethod, collapse = ",")))
  }


  # filtering GGI
  colnames(GGI)[1:2] <- c("from", "to")
  colnames(pathway_GGI)[1:3] <- c("pathway", "from", "to")
  colnames(pathway_gene)[1:2] <- c("pathway", "gene")

  # filtering pathway
  colnames(pathway_MGI)[1:3] <- c("pathway", "miRNA", "gene")
  colnames(MGI)[1:2] <- c("miRNA", "gene")
  pathway <- pathway_MGI %>%
    mutate(MGI = paste0(miRNA,":",gene)) #%>%
    #filter(MGI %in% background)
  MGISize <- pathway_MGI %>% group_by(pathway) %>%
    summarise(n_MGI = n()) %>% ungroup()
  rm(pathway_MGI)

  pathway_counts <- table(pathway$pathway)
  if (!is.null(minSize)){
    pathways_to_remove_min <- names(pathway_counts[pathway_counts < minSize])
    pathway <- pathway[!pathway$pathway %in% pathways_to_remove_min, ]
  }
  if (!is.null(maxSize)){
    pathways_to_remove_max <- names(pathway_counts[pathway_counts > maxSize])
    pathway <- pathway[!pathway$pathway %in% pathways_to_remove_max, ]
  }

  n_pathway <- length(unique(pathway$pathway))
  if (n_pathway == 0) {
    warning("No pathways meet the size criteria.")
    #return(NULL)
    return(data.frame(
      pathway = character(), n_MGI = integer(), score = numeric(), rand_score_dist = character(), p_value = numeric(),  padj = numeric(),
      contribute_MGI = character(), n_contribute_MGI = integer(),
      stringsAsFactors = FALSE
    ))}


  ### build global graph
  MGI$type <- "MGI"
  colnames(MGI)[1:3] <- c("from", "to", "type")
  GGI$type <- "GGI"
  all_edges <- dplyr::distinct(rbind(MGI, GGI)) %>%
    dplyr::mutate(edge = paste0(from, ":", to))
  all_edges$weight <- 1  ### set weight as 1

  graph <- igraph::graph_from_data_frame(all_edges, directed = TRUE)
  edge_ids <- paste0(all_edges$from, ":", all_edges$to)
  names(edge_ids) <- seq_along(edge_ids)
  orig_edges_df <- igraph::as_data_frame(graph, what = "edges")
  edge_names <- paste0(orig_edges_df$from, ":", orig_edges_df$to)
  E(graph)$name <- edge_names

  # Divide DEMGI into two groups: isolated/non-isolated
  isolated_edges <- base::setdiff(DEMGI, pathway$MGI)
  isolated_df <- data.frame(
    miRNA = sub(":.*", "", isolated_edges),
    gene  = sub(".*:", "", isolated_edges),
    stringsAsFactors = FALSE
  )
  pathway_genes <- unique(pathway_gene$gene)
  orig_edges_df <- igraph::as_data_frame(graph, what = "edges")
  connect_genes <- unique(orig_edges_df$from[
    orig_edges_df$type == "GGI" & orig_edges_df$to %in% pathway_genes])

  DEMGI_out <- isolated_df %>%
    filter(!(gene %in% connect_genes)) %>%
    mutate(MGI = paste0(miRNA, ":", gene)) %>%
    pull(MGI)
  DEMGI_connect <- setdiff(isolated_edges, DEMGI_out)
  DEMGI_in <- setdiff(DEMGI, isolated_edges)
  # DEMGI_in <- setdiff(DEMGI, DEMGI_out)
  if (length(DEMGI_in) + length(DEMGI_connect) == 0){ #
    stop("All DEMGIs are not included in any pathway. Edge_Network doesn't work in this situation. Please try other methods.\n")
  }

  ### Step 1. build restart list
  cat("    Build restart list...\n")
  restart_list <- as.integer(edge_ids %in% DEMGI)
  names(restart_list) = edge_ids

  # num_DE = sum(restart_list)
  num_DE_in <- length(DEMGI_in)
  num_DE_connect <- length(DEMGI_connect)
  MGI_out_indices <- which(all_edges$type == "MGI" & edge_ids %in% DEMGI_out) # fixed as 1
  MGI_pw_indices <- which(all_edges$type == "MGI" & edge_ids %in% pathway$MGI) # all MGIs that in the randomization range
  MGI_connect_indices <- setdiff(which(all_edges$type == "MGI"), c(MGI_out_indices, MGI_pw_indices))
  # MGI_pw_indices <- setdiff(which(all_edges$type == "MGI"), MGI_out_indices)

  set.seed(12345)
  seeds_null = sample(10000, iter, replace = FALSE)

  restart_matrix = cbind(
    restart_list,
    sapply(1:iter, function(i){
      set.seed(seeds_null[i])
      new_restart <- rep(0, length(edge_ids))
      new_restart[MGI_out_indices] <-1 # fixed DEMGI that outside the pathways as 1
      # selected_indices <- sample(MGI_indices, num_DE, replace = FALSE)
      selected_pw_indices <- sample(MGI_pw_indices, num_DE_in, replace = FALSE)
      selected_connect_indices <- sample(MGI_connect_indices, num_DE_connect, replace = FALSE)
      selected_indices <- c(selected_pw_indices,selected_connect_indices)
      new_restart[selected_indices] <- 1
      return(new_restart)
    })
  )
  colnames(restart_matrix) = c('OBS', paste0("P", seq(iter)))

  # cat("      Size of restart_matrix:", object.size(restart_matrix)/(1024*1024), "MB.\n")
  restart_matrix <- as(Matrix::Matrix(restart_matrix, sparse = TRUE), "dgCMatrix")
  cat("      Size of sparse restart_matrix:", object.size(restart_matrix)/(1024*1024), "MB.\n")
  rm(orig_edges_df, isolated_df, MGI, GGI)

  # ### Step 2. Sparse Edge-based Transition Matrix
  if (is.null(ncores) || ncores == 1){
    ncores = 1
    warning("You are use unparallel computing. It could cause huge amount of time!!!
            \nFor parallel computing, please specify [ncores] parameter to be larger than 1 (We highly recommend that)!!!\n")
  }
  cat("    Build transition matrix...\n")
  # cat("      Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

  lg <- igraph::make_line_graph(graph)
  V(lg)$name <- E(graph)$name
  edge_adjM <- as_adjacency_matrix(lg, sparse = TRUE)
  row_sums <- Matrix::rowSums(edge_adjM)
  row_sums[row_sums == 0] <- 1
  D_inv <- Matrix::Diagonal(x = 1 / row_sums)
  transition_matrix <- D_inv %*% edge_adjM
  transition_matrix <- Matrix::t(transition_matrix) # each row is a target edge, each column is a source edge, so each colsum should be 1 or 0
  colnames(transition_matrix) = rownames(transition_matrix)

  # add self-loop for GGIs who don't have downstream edge
  ggi_nodes <- all_edges$edge[all_edges$type == "GGI"]
  cols_to_fix <- which(colnames(transition_matrix) %in% ggi_nodes & Matrix::colSums(transition_matrix) == 0)

  if (length(cols_to_fix) > 0) {
    diag_fix <- Matrix::sparseMatrix(i = cols_to_fix,
                                     j = cols_to_fix,
                                     x = rep(1, length(cols_to_fix)),
                                     dims = dim(transition_matrix),
                                     dimnames = dimnames(transition_matrix))
    transition_matrix <- transition_matrix + diag_fix
  }
  cat("      Size of transition_matrix:", object.size(transition_matrix)/(1024*1024), "MB.\n")
  rm(diag_fix, D_inv, edge_adjM, graph)
  gc()

  ### Step 3. RWR Execution
  run_rwr <- function(seed.vec){
    alpha <- 0.5
    delta <- 1e-6
    max_steps <- 100
    p <- seed.vec
    new_p <- p # 预分配，避免重复创建对象
    diff_val <- Inf
    step <- 0

    repeat {
      tmp <- transition_matrix %*% p
      tmp@x <- tmp@x * (1 - alpha) 
      new_p <- tmp
      new_p <- new_p + alpha * seed.vec

      # check difference
      # diff_val <- sum(abs(new_p - p))
      diff_val <- Matrix::colSums(abs(new_p - p))
      step <- step + 1

      if (diff_val < delta || step >= max_steps) break

      # swap 引用而不是复制（就地交换引用避免两份大对象）
      tmp <- p
      p <- new_p
      new_p <- tmp

    }
    return(p)
  }


  obs_vec <- restart_matrix[, 1, drop = FALSE] # drop = FALSE
  perm_mat <- restart_matrix[, -1, drop = FALSE]
  n_cols <- ncol(perm_mat)



  cat("    Execute RWR with", iter, "iterations on", ncores, "cores...\n")

  obs_score <- run_rwr(obs_vec)  # calculate the observation scores

  if (ncores > 1){
    # divided restart_matrix to ncores blocks (sub-matrix), where each column represent a single permutation
    avg <- ceiling(n_cols / ncores)
    block_indices <- split(seq_len(n_cols), rep(seq_len(ncores), each = avg, length.out = n_cols))



    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl, varlist = c("transition_matrix", "run_rwr", "perm_mat"), envir = environment())  # , "edge_ids", "restart_matrix"
    parallel::clusterEvalQ(cl, library(Matrix))

    perm_scores_list <- parallel::parLapply(cl, block_indices, function(idxs) {
      do.call(cbind, lapply(idxs, function(i) run_rwr(perm_mat[, i, drop = FALSE])))
    })
    parallel::stopCluster(cl)
    perm_scores <- do.call(cbind, perm_scores_list)
  } else {
    # perm_scores <- apply(perm_mat, 2, run_rwr)
    perm_scores <- do.call(cbind, lapply(seq_len(ncol(perm_mat)), function(i) {
      run_rwr(perm_mat[, i, drop = FALSE])
    }))
  }


  ### Step 4. pathway aggregation
  cat("    Summarise by pathways...\n")
  pathway <- pathway %>% mutate(edge = paste0(miRNA, ":", gene))
  pathway_GGI <- pathway_GGI %>% mutate(edge = paste0(from, ":", to))
  pw_all <- rbind(pathway[, c("pathway", "edge")], pathway_GGI[, c("pathway", "edge")])
  pw_dt <- as.data.table(pw_all)
  rm(pathway, pathway_GGI, pw_all)

  # obs_dt <- data.table(edge = edge_ids, score = obs_score)
  obs_dt <- data.table(edge = rownames(obs_score), score = as.numeric(obs_score))
  setkey(obs_dt, edge)
  obs_result <- pw_dt[,.(score = sum(obs_dt[edge,score])), by = pathway]
  # obs_result <- merge(obs_dt, pw_dt, by = "edge", allow.cartesian = TRUE)[, .(score = sum(score)), by = pathway]

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
  rm(obs_dt, edge_score)
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

  # 对每列 permutation 做矩阵乘法
  perm_scores <- perm_scores[match(colnames(P), rownames(perm_scores)), , drop = FALSE]
  perm_result_mat <- P %*% perm_scores  # pathways x iterations

  # 转回 data.table
  perm_coo <- Matrix::summary(perm_result_mat)  # i,j,x
  perm_result <- data.table(
    pathway = rownames(perm_result_mat)[perm_coo$i],
    iter = colnames(perm_result_mat)[perm_coo$j],
    score = perm_coo$x
  )

  rm(perm_scores, perm_coo)

  random_stats <- perm_result[, .(
    rand_score_dist = paste(
      sprintf("%.4f", min(score)),
      sprintf("%.4f", quantile(score, 0.25)),
      sprintf("%.4f", median(score)),
      sprintf("%.4f", mean(score)),
      sprintf("%.4f", quantile(score, 0.75)),
      sprintf("%.4f", max(score)),
      sep = ","
    )
  ), by = pathway]

  result <- merge(perm_result, obs_result, by = "pathway", suffixes = c("_perm", "_obs"))
  rm(perm_result)

  # calculate p-values
  ctx <- V8::v8()
  decimal_code <- paste(readLines("Javascript/decimal.js"), collapse = "\n")
  getP_code     <- paste(readLines("Javascript/GetpValueFromZ.js"), collapse = "\n")
  ctx$eval(decimal_code)
  ctx$eval(getP_code)

  pval_dt <- result[, {
    sd_perm <- sd(score_perm)
    mu_perm <- mean(score_perm)
    obs <- unique(score_obs)
    if (is.na(sd_perm) || sd_perm == 0) {
      .(z = NA_real_,
        p_value = 1.0000)
    } else {
      z <- (obs - mu_perm) / sd_perm
      p <- ctx$call("GetpValueFromZ", z, "right")
      .(z = as.numeric(z),
        p_value = as.numeric(p))
    }
  }, by = pathway]
  final <- merge(obs_result, pval_dt, by = "pathway")
  final <- merge(final, random_stats, by = "pathway")
  final[, padj := p.adjust(p_value, method = pAdjMethod)]
  final <- merge(final, contributing_mgi, by = "pathway", all.x = TRUE)
  final <- merge(MGISize, final, by = "pathway")
  final <- final %>% dplyr::filter(p_value <= pvalueCutoff, padj <= pvalueCutoff) %>%
    dplyr::arrange(padj) %>% dplyr::relocate(rand_score_dist, .after = score)


  cat("  Edge-Network analysis have finished!\n")
  return(final)
}


