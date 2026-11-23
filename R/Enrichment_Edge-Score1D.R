#' Edge-Score1D (MiRNA-Gene Interaction Edge Score-based Analysis)
#'
#' Perform miRNA-gene interaction (MGI) enrichment analysis based on whether the pathway MGI are located unexpectedly forward towards the ranked MGI list.
#'
#' @param MGIList A named vector shows the ranked score of MGIs.
#' @param pathway A three-column dataframe contains pathway, miRNA, gene.
#' @param pvalueType How to calculate p value. See default_pvalueType for all possible options. Default pvalueType is "pos".
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negative ("neg") enrichment).
#' The "std" is used when interested in both positive and negative scores.
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default_pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If no adjustment is desired, set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values above the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs less than the minSize threshold.
#' @param maxSize The maximum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs more than the maxSize threshold.
#' @return A dataframe containing the enrichment result for Edge-Score1D.
#' @examples
#' data(toy_pathway)
#' data(toy_mir_DEdata)
#' data(toy_gene_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
#' data(background_MGI)
#' input_data <- getAllInputData(
#'     methods = "Edge_Score1D", pathway = toy_pathway, mir_DEdata = toy_mir_DEdata,
#'     gene_DEdata = toy_gene_DEdata, background_MGI = background_MGI,
#'     gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
#' )
#' res <- Edge_Score1D(
#'     MGIList = input_data$data$Edge_Score1D, pathway = input_data$pathway$Edge,
#'     pvalueType = "neg", pvalueCutoff = 1, minSize = 2, maxSize = 5000
#' )
#' head(res)
#' @export
Edge_Score1D <- function(MGIList, pathway, pvalueType = "neg", pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL) {
    message("\n  Start Edge-Score1D analysis ...")

    if (is.null(MGIList)) {
        stop("A valid named vector MGIList contains the MGIs with their scores must be entered!")
    }
    if (!is.vector(MGIList) || !is.numeric(MGIList)) {
        stop("MGIList must be a named number vector, where the name is MGI, and the value is the ranked score!")
    }

    if (!pAdjMethod %in% default_pAdjMethod) {
        stop(
            "Invalid pAdjMethod: ", pAdjMethod,
            "\n Please choose one of: ", paste(default_pAdjMethod, collapse = ",")
        )
    }

    colnames(pathway)[seq_len(3)] <- c("pathway", "miRNA", "gene")
    pathway <- pathway %>%
        dplyr::mutate(MGI = paste0(miRNA, ":", gene)) %>%
        dplyr::select(pathway, MGI)

    MGIList <- sort(MGIList, decreasing = TRUE, na.last = NA)

    colnames(pathway)[seq_len(2)] <- c("Term", "MGI")
    gsea <- clusterProfiler::GSEA(geneList = MGIList, TERM2GENE = pathway, scoreType = pvalueType, pAdjustMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minGSSize = minSize, maxGSSize = maxSize)
    gsea_result <- as.data.frame(gsea) %>% rename(pathway = ID, padj = p.adjust, p_value = pvalue)

    message("  Edge-Score1D analysis have finished!")

    return(gsea_result)
}
