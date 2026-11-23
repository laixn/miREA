#' MiR-ORA (MiRNA Over-Representation Analysis)
#'
#' Perform miRNA enrichment analysis based on the overlap between differentially expressed miRNA set and pathway miRNA set. Pathway miRNA set was mapped from pathway gene set.
#'
#' @param DEmiR A character vector contains selected differential expressed miRNAs.
#' @param background A character vector contains background miRNAs.
#' @param pathway A two-column dataframe contains pathway and their member miRNAs.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default_pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If no adjustment is desired, set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values above the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum miRNA set sizes allowed for analysis, filtering out pathways with number of member miRNAs less than the minSize threshold.
#' @param maxSize The maximum miRNA set sizes allowed for analysis, filtering out pathways with number of member miRNAs more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for MiR-ORA.
#' @examples
#' data(toy_pathway); data(toy_mir_DEdata); data(toy_gene_DEdata); data(toy_background_MGI)
#' input_data <- getAllInputData(
#'   methods = "MiR_ORA", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
#'   gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI
#' )
#' res <- MiR_ORA(
#'   DEmiR = input_data$data$MiR_ORA, background = unique(input_data$pathway$MiR$miRNA),
#'   pathway = input_data$pathway$MiR, pvalueCutoff = 1
#' )
#' head(res)
#' @export
MiR_ORA <- function(DEmiR, background = NULL, pathway, pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL){
  message("\n  Start MiR-ORA analysis ...")

  if (is.null(DEmiR) || length(DEmiR) == 0) {
    stop("A valid character vector contains the differential expressed miRNAs must be entered!")
  }

  if (!pAdjMethod %in% default_pAdjMethod){
    stop(paste("Invalid pAdjMethod:", pAdjMethod,
               "\n Please choose one of:", paste(default_pAdjMethod, collapse = ",")))
  }

  if (is.null(background)){
    message("    No background data is given, extract all mature miRNAs from miRBase version 22 instead. ")
    background <- unique(miRBaseConverter::getAllMiRNAs(version="v22", type="mature", species="hsa")$Name)
  }

  # filtering DEmiRs in the background
  if (all(DEmiR %in% background) == FALSE) {
    message(paste(length(DEmiR[!DEmiR %in% background]), "of", length(DEmiR),
        "miRNAs are not included in the background. Those are: \n",
        paste(DEmiR[!DEmiR %in% background], collapse = " "),
        "\n Those miRNAs are deleted. "))
    DEmiR <- DEmiR[DEmiR %in% background]
    if (length(DEmiR) == 0) {
      stop("No valid DEmiR input from background list!")
    }
  }

  # filtering pathway in the background
  pathway <- pathway %>%
    filter(.[[2]] %in% background)

  colnames(pathway)[seq_len(2)] <- c("Term", "MiRNA")
  ora <- clusterProfiler::enricher(gene = DEmiR, universe = background, TERM2GENE = pathway, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = 1, minGSSize = minSize, maxGSSize = maxSize)
  ora_result <- as.data.frame(ora) %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)
  message("  MiR-ORA analysis have finished! ")

  return(ora_result) # a dataframe

}


#' MiR-Score (MiRNA Score-based Analysis)
#'
#' Perform miRNA enrichment analysis based on whether the pathway miRNA are located unexpectedly forward towards the ranked miRNA list.
#'
#' @param miRList A named vector shows the ranked score of miRNAs.
#' @param pathway A two-column dataframe contains pathway and their member miRNAs.
#' @param pvalueType How to calculate p value. See default_pvalueType for all possible options. Default pvalueType is "pos".
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negative ("neg") enrichment).
#' The "std" is used when interested in both positive and negative scores.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default_pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If no adjustment is desired, set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values above the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum miRNA set sizes allowed for analysis, filtering out pathways with number of member miRNAs less than the minSize threshold.
#' @param maxSize The maximum miRNA set sizes allowed for analysis, filtering out pathways with number of member miRNAs more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for MiR-Score.
#' @examples
#' data(toy_pathway); data(toy_mir_DEdata); data(toy_gene_DEdata); data(toy_background_MGI)
#' input_data <- getAllInputData(
#'   methods = "MiR_Score", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
#'   gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI
#' )
#' res <- MiR_Score(
#'   miRList = input_data$data$MiR_Score, pathway = input_data$pathway$MiR,
#'   pvalueType = "std", pvalueCutoff = 1
#' )
#' head(res)
#' @export
MiR_Score <- function(miRList, pathway, pvalueType = "std", pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL){ # pvalueType = "std",
  message("\n  Start MiR-Score analysis ...")

  if (is.null(miRList)) {
    stop("A valid named vector miRList contains the miRNAs with their scores must be entered!")
  }
  if (!is.vector(miRList) || !is.numeric(miRList)) {
    stop("miRList must be a named number vector, where the name is miRNA, and the value is the ranked score!")
  }

  if (!pAdjMethod %in% default_pAdjMethod){
    stop(paste("Invalid pAdjMethod:", pAdjMethod,
               "\n Please choose one of:", paste(default_pAdjMethod, collapse = ",")))
  }

  miRList <- sort(miRList, decreasing = TRUE, na.last = NA)

  colnames(pathway)[seq_len(2)] <- c("Term", "MiRNA")
  gsea <- clusterProfiler::GSEA(geneList = miRList, TERM2GENE = pathway, scoreType = pvalueType, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff,  minGSSize = minSize, maxGSSize = maxSize) #  scoreType = pvalueType,
  gsea_result <- as.data.frame(gsea) %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)
  message("  MiR-Score analysis have finished!")

  return(gsea_result)
}
