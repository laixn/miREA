#' Title: TG-ORA (empirical) (Target Gene Empirical Over-Representation Analysis)
#' Discription: Permutation-based (empirical) alternative to TG_ORA(). Instead of a plain
#'        hypergeometric test on the target-gene set, the query is defined at the miRNA
#'        level (DEmiR) and translated to target genes via a miRNA-target gene interaction
#'        (MGI) table; the same-size resampling of miRNA identities (translated to their
#'        pooled target genes via the SAME MGI table) is used to build an empirical null
#'        distribution for the pathway overlap. This is a faithful, minimal-change port of
#'        miRNetR's algo=='emp' branch (Chang et al. 2020, GPL-3.0; https://github.com/xia-lab/miRNetR/tree/master;
#'        utils_mir_target_enrich.R -> my.mir.target.enrich() / GetRandomMirTargetGenes()),
#'        already bit-exact validated against a live miRNetR run this session (see
#'        major2/miRNet/utils_custom_mir_enrich.R, PerformCustomMirTargetORA() -- same core
#'        algorithm, this function only adapts the argument/return shapes to miREA's own
#'        data conventions). The p-value uses miRNet's own formula unmodified,
#'        nMoreExtreme / iter (no pseudocount), so it can be exactly 0 when the observed hit
#'        count exceeds every one of the `iter` permutation draws.
#' @param DEmiR A character vector of differentially expressed (query) miRNAs. Unlike
#'        TG_ORA()'s `TG` (an already-derived gene list), this is the function's defining
#'        input: the empirical null needs miRNA-level identity to resample, so the target
#'        gene set is derived *internally* from DEmiR + background_MGI (mirroring miRNet's
#'        own design), not supplied directly.
#' @param background_MGI A two-column dataframe (or matrix) of the miRNA-target gene
#'        interaction (MGI) database: column 1 = miRNA ID, column 2 = target gene ID. This
#'        is the "MGI database" -- swap it for any custom miRNA-target resource.
#' @param pathway A two-column dataframe contains pathway and their member genes.
#' @param background_gene A character vector of background/universe genes. If NULL
#'        (default), the universe is all genes appearing in `pathway` (unique(pathway[[2]])),
#'        matching PerformCustomMirTargetORA()'s `universe` and TG_ORA()'s implicit
#'        clusterProfiler behaviour.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and
#'        control the false discovery rate. Use default.pAdjMethod to find all possible
#'        values. If not specified, Benjamini & Hochberg (BH) method will be automatically
#'        used. If you don't want any adjustment, please set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways
#'        with padj values below the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#'        (Kept for interface parity with TG_ORA()/TG_Score(); this function itself does not
#'        drop any pathway based on pvalueCutoff -- all tested pathways are returned, and
#'        pvalueCutoff is expected to be applied by the caller, exactly like the
#'        n_enrich = sum(padj < 0.05) convention used throughout this project.)
#' @param minSize The minimum gene set sizes allowed for analysis, filtering out pathways
#'        with number of member genes less than the minSize threshold.
#' @param maxSize The maximum gene set sizes allowed for analysis, filtering out pathways
#'        with number of member genes more than the maxSize threshold.
#' @param iter Number of permutations for the empirical null. Default iter is 1000 -- same
#'        name/default as the `iter` parameter already used by Edge_2Ddist/Edge_Topology/
#'        Edge_Network in miREA() (code/github/R/miREA.R), not a new naming convention.
#' @param replace Whether miRNA resampling is done with replacement. Default TRUE, faithful
#'        to miRNet's shipped GetRandomMirTargetGenes() (`sample(mirs, qSize, replace=TRUE)`)
#'        -- this is what was bit-exact validated against a live miRNetR run earlier this
#'        session. major2/01_run_bias_corrected_tg_ora.R's TG_ORA_emp() instead deliberately
#'        uses replace = FALSE (a documented, justified deviation, found to reduce
#'        anti-conservative bias on that analysis's own FPR benchmark) -- set replace = FALSE
#'        here to reproduce that variant's resampling scheme with this function instead.
#' @return A dataframe containing the enrichment result for TG-ORA (empirical): pathway,
#'        n_hit, n_pathway_size, p_value, padj. Same shape as the existing TG_ORA_emp() in
#'        major2/01_run_bias_corrected_tg_ora.R and consumable by the same
#'        n_enrich = sum(padj < 0.05) convention used everywhere downstream in this project.
#'
TG_ORA_emp <- function(DEmiR, background_MGI, pathway, background_gene = NULL,
                        pAdjMethod = "BH", pvalueCutoff = 0.05,
                        minSize = NULL, maxSize = NULL, iter = 1000, replace = TRUE) {
  cat("\n  Start TG-ORA (empirical) analysis ...\n")

  if (is.null(DEmiR)) {
    stop("A valid character vector contains the differentially expressed miRNAs must be entered!")
  }

  if (!pAdjMethod %in% default.pAdjMethod) {
    stop(paste("Invalid pAdjMethod:", pAdjMethod,
               "\n Please choose one of:", paste(default.pAdjMethod, collapse = ",")))
  }

  colnames(pathway)[1:2] <- c("Term", "Gene")
  colnames(background_MGI)[1:2] <- c("miRNA", "Gene")

  if (is.null(background_gene)) {
    background_gene <- unique(pathway$Gene)
  }

  path_ids <- unique(pathway$Term)
  members_list <- split(pathway$Gene, pathway$Term)[path_ids]
  set_size <- lengths(members_list)
  size_ok <- rep(TRUE, length(set_size))
  if (!is.null(minSize)) size_ok <- size_ok & (set_size >= minSize)
  if (!is.null(maxSize)) size_ok <- size_ok & (set_size <= maxSize)
  path_ids <- path_ids[size_ok]
  members_list <- members_list[path_ids]
  set_size <- set_size[path_ids]

  DEmiR <- unique(DEmiR)
  all_mirs <- unique(background_MGI$miRNA)
  mir2gene <- split(background_MGI$Gene, background_MGI$miRNA)

  q_mirs <- DEmiR[DEmiR %in% all_mirs]
  if (length(q_mirs) == 0) {
    stop("None of the query miRNAs (DEmiR) were found in background_MGI (no known targets).")
  }

  get_targets <- function(mirs) {
    unique(unlist(mir2gene[mirs], use.names = FALSE))
  }

  TG <- get_targets(q_mirs)
  TG <- TG[TG %in% background_gene]

  hit_obs <- vapply(members_list, function(g) length(intersect(TG, g)), numeric(1))

  perm.out <- matrix(0L, nrow = length(path_ids), ncol = iter)
  for (i in seq_len(iter)) {
    rand_mir <- sample(all_mirs, length(q_mirs), replace = replace)
    rand_TG <- get_targets(rand_mir)
    rand_TG <- rand_TG[rand_TG %in% background_gene]
    perm.out[, i] <- vapply(members_list, function(g) length(intersect(rand_TG, g)), numeric(1))
  }

  # empirical p from permutation: nMoreExtreme / iter -- miRNet's own shipped formula
  # (`perm.out - hit.num > 0` / iter), unmodified, no pseudocount. Bit-exact validated against
  # a live miRNetR run this session (major2/08_validate_TG_ORA_emp_vs_live_miRNet.R). Can be
  # exactly 0 when hit_obs exceeds every one of the `iter` permutation draws.
  nMoreExtreme <- rowSums(perm.out > hit_obs)
  p_value <- nMoreExtreme / iter
  padj <- p.adjust(p_value, method = pAdjMethod)

  ora_result <- data.frame(pathway = path_ids, n_hit = hit_obs, n_pathway_size = set_size,
                            p_value = p_value, padj = padj, stringsAsFactors = FALSE,
                            row.names = NULL)

  cat("  TG-ORA (empirical) analysis have finished! \n")

  # class(ora_result) = "TG_ORA_emp"
  return(ora_result) # a dataframe
}
