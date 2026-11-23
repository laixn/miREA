#' TG-ORA (Target Gene Over-Representation Analysis)
#'
#' Perform target gene enrichment analysis based on the overlap between differentially expressed gene set and pathway gene set.
#'
#' @param TG A character vector contains selected target genes of DEmiRs.
#' @param background A character vector contains background genes.
#' @param pathway A two-column dataframe contains pathway and their member genes.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default_pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If no adjustment is desired, set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values above the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum gene set sizes allowed for analysis, filtering out pathways with number of member genes less than the minSize threshold.
#' @param maxSize The maximum gene set sizes allowed for analysis, filtering out pathways with number of member genes more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for TG-ORA.
#' @examples
#' data(toy_pathway); data(toy_mir_DEdata); data(toy_gene_DEdata); data(toy_background_MGI)
#' input_data <- getAllInputData(
#'   methods = "TG_ORA", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
#'   gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI
#' )
#' res <- TG_ORA(
#'   TG = input_data$data$TG_ORA, background = unique(toy_pathway$gene),
#'   pathway = input_data$pathway$TG, pvalueCutoff = 1
#' )
#' head(res)
#' @export
#'
# if background is null, will use background provided in clusterProfiler
TG_ORA <- function(TG, background = NULL, pathway, pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL){
  message("\n  Start TG-ORA analysis ...")

  if (is.null(TG) || length(TG) == 0) {
    stop("A valid character vector contains the differential expressed genes must be entered! ",
         "TG is empty: no genes passed the upstream differential-expression filters.")
  }

  if (!pAdjMethod %in% default_pAdjMethod){
    stop("Invalid pAdjMethod: ", pAdjMethod,
         "\n Please choose one of: ", paste(default_pAdjMethod, collapse = ","))
  }

  colnames(pathway)[seq_len(2)] <- c("Term", "Gene")
  ora <- clusterProfiler::enricher(gene = TG, universe = background, TERM2GENE = pathway, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = 1, minGSSize = minSize, maxGSSize = maxSize)
  ora_result <- as.data.frame(ora)  %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)
  ## if use ora@result, all pathways with valid DEG will remain, no matter whether the pvalue pass the Cutoff threshold.
  message("  TG-ORA analysis have finished! ")

  return(ora_result) # a dataframe

}


#' TG-Score (Target Gene Score-based Analysis)
#'
#' Perform target gene enrichment analysis based on whether the pathway genes are located unexpectedly forward towards the ranked gene list.
#'
#' @param geneList A named vector shows the ranked score of target genes.
#' @param pathway A two-column dataframe contains pathway and their member genes.
#' @param pvalueType How to calculate p value. See default_pvalueType for all possible options. Default pvalueType is "pos".
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negative ("neg") enrichment).
#' The "std" is used when interested in both positive and negative scores.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default_pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If no adjustment is desired, set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values above the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum gene set sizes allowed for analysis, filtering out pathways with number of member genes less than the minSize threshold.
#' @param maxSize The maximum gene set sizes allowed for analysis, filtering out pathways with number of member genes more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for TG-Score.
#' @examples
#' data(toy_pathway); data(toy_mir_DEdata); data(toy_gene_DEdata); data(toy_background_MGI)
#' input_data <- getAllInputData(
#'   methods = "TG_Score", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
#'   gene_DEdata = toy_gene_DEdata, background_MGI = toy_background_MGI
#' )
#' res <- TG_Score(
#'   geneList = input_data$data$TG_Score, pathway = input_data$pathway$TG,
#'   pvalueType = "std", pvalueCutoff = 1
#' )
#' head(res)
#' @export
TG_Score <- function(geneList, pathway, pvalueType = "std", pAdjMethod = "BH",  pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL){ # pvalueType = "std", background = NULL
  message("\n  Start TG-Score analysis ...")

  if (is.null(geneList)) {
    stop("A valid named vector geneList contains the genes with their scores must be entered!")
  }
  if (!is.vector(geneList) || !is.numeric(geneList)) {
    stop("geneList must be a named number vector, where the name is gene, and the value is the ranked score!")
  }

  if (!pAdjMethod %in% default_pAdjMethod){
    stop("Invalid pAdjMethod: ", pAdjMethod,
         "\n Please choose one of: ", paste(default_pAdjMethod, collapse = ","))
  }

  geneList <- sort(geneList, decreasing = TRUE, na.last = NA)

  colnames(pathway)[seq_len(2)] <- c("Term", "Gene") # scoreType = pvalueType,
  gsea <- clusterProfiler::GSEA(geneList = geneList, TERM2GENE = pathway,  scoreType = pvalueType, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minGSSize = minSize, maxGSSize = maxSize)
  gsea_result <- as.data.frame(gsea) %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)
  ## gsea@result will only remain pathways that pass the cutoff threshold. But use as.data.frame() as well for harmonization with ORA analysis

  message("  TG-Score analysis have finished!")

  return(gsea_result)
}
