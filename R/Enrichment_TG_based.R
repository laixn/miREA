#' Title: TG-ORA (Target Gene Over-Representation Analysis)
#' Discription: Perform target gene enrichment analysis based on the overlap between differentially expressed gene set and pathway gene set.
#' @param TG A character vector contains selected target genes of DEmiRs.
#' @param background A character vector contains background genes.
#' @param pathway A two-column dataframe contains pathway and their member genes.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default.pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If you don't want any adjustment, please set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values below the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum gene set sizes allowed for analysis, filtering out pathways with number of member genes less than the minSize threshold.
#' @param maxSize The maximum gene set sizes allowed for analysis, filtering out pathways with number of member genes more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for TG-Score.
#'

# if background is null, will use background provided in clusterProfiler
TG_ORA <- function(TG, background = NULL, pathway, pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL){
  cat("\n  Start TG-ORA analysis ...\n")

  if (is.null(TG)) {
    stop("A valid character vector contains the differential expressed genes must be entered!")
  }

  if (!pAdjMethod %in% default.pAdjMethod){
    stop(paste("Invalid pAdjMethod:", pAdjMethod,
               "\n Please choose one of:", paste(default.pAdjMethod, collapse = ",")))
  }

  colnames(pathway)[1:2] <- c("Term", "Gene")
  ora <- clusterProfiler::enricher(gene = TG, universe = background, TERM2GENE = pathway, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = 1, minGSSize = minSize, maxGSSize = maxSize)
  ora_result <- as.data.frame(ora)  %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)
  ## if use ora@result, all pathways with valid DEG will remain, no matter whether the pvalue pass the Cutoff threshold.
  cat("  TG-ORA analysis have finished! \n")

  # class(ora_result) = "TG_ORA"
  return(ora_result) # a dataframe

}


#' Title: TG-Score (Target Gene Score-based Analysis)
#' Description: Perform target gene enrichment analysis based on whether the pathway genes are located unexpectedly forward towards the ranked gene list.
#' @param geneList A named vector shows the ranked score of target genes.
#' @param pathway A two-column dataframe contains pathway and their member genes.
#' @param pvalueType ow to calculate p value. See default.pvalueType for all possible options. Default pvalueType is "pos".
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negative ("neg") enrichment).
#' The "std" is used when interested in both positive and negative scores.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default.pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If you don't want any adjustment, please set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values below the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum gene set sizes allowed for analysis, filtering out pathways with number of member genes less than the minSize threshold.
#' @param maxSize The maximum gene set sizes allowed for analysis, filtering out pathways with number of member genes more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for TG-Score.
#'
TG_Score <- function(geneList, pathway, pvalueType = "std", pAdjMethod = "BH",  pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL){ # pvalueType = "std", background = NULL
  cat("\n  Start TG-Score analysis ...\n")

  if (is.null(geneList)) {
    stop("A valid named vector geneList contains the genes with their scores must be entered!")
  }
  if (!is.vector(geneList) || !is.numeric(geneList)) {
    stop("geneList must be a named number vector, where the name is gene, and the value is the ranked score!")
  }

  if (!pAdjMethod %in% default.pAdjMethod){
    stop(paste("Invalid pAdjMethod:", pAdjMethod,
               "\n Please choose one of:", paste(default.pAdjMethod, collapse = ",")))
  }

  geneList <- sort(geneList, decreasing = TRUE, na.last = NA)

  colnames(pathway)[1:2] <- c("Term", "Gene") # scoreType = pvalueType,
  gsea <- clusterProfiler::GSEA(geneList = geneList, TERM2GENE = pathway,  scoreType = pvalueType, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minGSSize = minSize, maxGSSize = maxSize)
  gsea_result <- as.data.frame(gsea) %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)
  ## gsea@result will only remain pathways that pass the cutoff threshold. But use as.data.frame() as well for harmonization with ORA analysis
  # gsea_result <- do.call(data.frame, gsea_result)

  cat("  TG-Score analysis have finished!\n")

  # class(gsea_result) = "TG_Score"
  return(gsea_result)
}
