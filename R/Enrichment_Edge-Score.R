#' Title: Edge-Score (MiRNA-Gene Interaction Edge Score-based Analysis)
#' Description: Perform miRNA-gene interaction (MGI) enrichment analysis based on whether the pathway MGI are located unexpectedly forward towards the ranked MGI list.
#' @param MGIList A named vector shows the ranked score of MGIs.
#' @param pathway A three-column dataframe contains pathway, miRNA, gene.
#' @param pvalueType How to calculate p value. See default.pvalueType for all possible options. Default pvalueType is "pos".
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negative ("neg") enrichment).
#' The "std" is used when interested in both positive and negative scores.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default.pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If you don't want any adjustment, please set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values below the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs less than the minSize threshold.
#' @param maxSize The maximum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for Edge-Score.


Edge_Score <- function(MGIList,  pathway, pvalueType = "neg", pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL){
  cat("\n  Start Edge-Score analysis ...\n")

  if (is.null(MGIList)) {
    stop("A valid named vector MGIList contains the MGIs with their scores must be entered!")
  }
  if (!is.vector(MGIList) || !is.numeric(MGIList)) {
    stop("MGIList must be a named number vector, where the name is MGI, and the value is the ranked score!")
  }

  if (!pAdjMethod %in% default.pAdjMethod){
    stop(paste("Invalid pAdjMethod:", pAdjMethod,
               "\n Please choose one of:", paste(default.pAdjMethod, collapse = ",")))
  }

  colnames(pathway)[1:3] <- c("pathway", "miRNA", "gene")
  pathway <- pathway %>%
    dplyr::mutate(MGI = paste0(miRNA,":",gene)) %>%
    dplyr::select(pathway, MGI)

  MGIList <- sort(MGIList, decreasing = TRUE, na.last = NA)

  colnames(pathway)[1:2] <- c("Term", "MGI")
  gsea <- clusterProfiler::GSEA(geneList = MGIList, TERM2GENE = pathway, scoreType = pvalueType, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minGSSize = minSize, maxGSSize = maxSize)
  gsea_result <- as.data.frame(gsea) %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)

  cat("  Edge-Score analysis have finished!\n")

  # class(gsea_result) = "Edge_Score"
  return(gsea_result)
}
