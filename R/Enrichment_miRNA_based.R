#' Title: MiR-ORA (MiRNA Over-Representation Analysis)
#' Discription: Perform miRNA enrichment analysis based on the overlap between differentially expressed miRNA set and pathway miRNA set. Pathway miRNA set was mapped from pathway gene set.
#' @param DEMiR A character vector contains selected differential expressed miRNAs.
#' @param background A character vector contains background miRNAs.
#' @param pathway A two-column dataframe contains pathway and their member miRNAs.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default.pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If you don't want any adjustment, please set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values below the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum miRNA set sizes allowed for analysis, filtering out pathways with number of member miRNAs less than the minSize threshold.
#' @param maxSize The maximum miRNA set sizes allowed for analysis, filtering out pathways with number of member miRNAs more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for MiR-ORA.
#'

MiR_ORA <- function(DEmiR, background = NULL, pathway, pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL){
  cat("\n  Start MiR-ORA analysis ...\n")

  if (is.null(DEmiR)) {
    stop("A valid character vector contains the differential expressed miRNAs must be entered!")
  }

  if (!pAdjMethod %in% default.pAdjMethod){
    stop(paste("Invalid pAdjMethod:", pAdjMethod,
               "\n Please choose one of:", paste(default.pAdjMethod, collapse = ",")))
  }

  if (is.null(background)){
    cat("    No background data is given, extract all mature miRNAs from miRBase version 22 instead. \n")
    # background <- get_background_data("MiR_ORA")
    background <- unique(miRBaseConverter::getAllMiRNAs(version="v22", type="mature", species="hsa")$Name)
  }

  # filtering DEmiRs in the background
  if (all(DEmiR %in% background) == FALSE) {
    cat(length(DEmiR[!DEmiR %in% background]), "of", length(DEmiR),
        "miRNAs are not included in the background. Those are: \n",
        paste(DEmiR[!DEmiR %in% background], collapse = " "),
        "\n Those miRNAs are deleted. \n", sep = " ")
    DEMiR <- DEmiR[DEmiR %in% background]
    if (is.null(DEmiR)) {
      stop("No valid DEmiR input from background list!")
    }
  }

  # filtering pathway in the background
  pathway <- pathway %>%
    filter(.[[2]] %in% background)

  colnames(pathway)[1:2] <- c("Term", "MiRNA")
  ora <- clusterProfiler::enricher(gene = DEmiR, universe = background, TERM2GENE = pathway, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = 1, minGSSize = minSize, maxGSSize = maxSize)
  ora_result <- as.data.frame(ora) %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)
  cat("  MiR-ORA analysis have finished! \n")

  # class(ora_result) = "MiR_ORA"
  return(ora_result) # a dataframe

}


#' Title: MiR-Score (MiRNA Score-based Analysis)
#' Description: Perform miRNA enrichment analysis based on whether the pathway miRNA are located unexpectedly forward towards the ranked miRNA list.
#' @param miRList A named vector shows the ranked score of miRNAs.
#' @param pathway A two-column dataframe contains pathway and their member miRNAs.
#' @param pvalueType ow to calculate p value. See default.pvalueType for all possible options. Default pvalueType is "pos".
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negative ("neg") enrichment).
#' The "std" is used when interested in both positive and negative scores.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default.pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If you don't want any adjustment, please set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values below the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum miRNA set sizes allowed for analysis, filtering out pathways with number of member miRNAs less than the minSize threshold.
#' @param maxSize The maximum miRNA set sizes allowed for analysis, filtering out pathways with number of member miRNAs more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for MiR-Score.
#'
#'   background
MiR_Score <- function(miRList, pathway, pvalueType = "std", pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL){ # pvalueType = "std",
  cat("\n  Start MiR-Score analysis ...\n")

  if (is.null(miRList)) {
    stop("A valid named vector miRList contains the miRNAs with their scores must be entered!")
  }
  if (!is.vector(miRList) || !is.numeric(miRList)) {
    stop("miRList must be a named number vector, where the name is miRNA, and the value is the ranked score!")
  }

  if (!pAdjMethod %in% default.pAdjMethod){
    stop(paste("Invalid pAdjMethod:", pAdjMethod,
               "\n Please choose one of:", paste(default.pAdjMethod, collapse = ",")))
  }

  miRList <- sort(miRList, decreasing = TRUE, na.last = NA)

  colnames(pathway)[1:2] <- c("Term", "MiRNA")
  gsea <- clusterProfiler::GSEA(geneList = miRList, TERM2GENE = pathway, scoreType = pvalueType, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff,  minGSSize = minSize, maxGSSize = maxSize) #  scoreType = pvalueType,
  gsea_result <- as.data.frame(gsea) %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)
  # gsea_result <- do.call(data.frame, gsea_result)
  cat("  MiR-Score analysis have finished!\n")

  #class(gsea_result) = "MiR_Score"
  return(gsea_result)
}
