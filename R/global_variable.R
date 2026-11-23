#' Supported enrichment methods
#'
#' Character vector of the nine node-based and edge-based enrichment methods
#' supported by [miREA()] and [getAllInputData()].
#'
#' @format A character vector of length 9.
#' @return A character vector of length 9.
#' @examples
#' default_method_list
#' @export
default_method_list <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score1D", "Edge_Score2D", "Edge_Topology", "Edge_Network")

#' Supported score transformations
#'
#' Character vector of the score transformations supported by the
#' `scoreFun` argument of [getAllInputData()].
#'
#' @format A character vector of length 4.
#' @return A character vector of length 4.
#' @examples
#' default_scoreFun
#' @export
default_scoreFun <- c("none", "rank", "sigmoid", "normalize")

#' Supported p-value adjustment methods
#'
#' Character vector of the multiple-testing correction methods supported by
#' the `pAdjMethod` argument of [miREA()] and the individual enrichment
#' method functions (see [stats::p.adjust()]).
#'
#' @format A character vector of length 8.
#' @return A character vector of length 8.
#' @examples
#' default_pAdjMethod
#' @export
default_pAdjMethod <- c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none")

#' Supported predefined pathway sets
#'
#' Character vector naming the pathway sets shipped with the miREA
#' manuscript's companion data deposit.
#'
#' @format A character vector of length 3.
#' @return A character vector of length 3.
#' @examples
#' default_pathway_list
#' @export
default_pathway_list <- c("Reactome", "hallmark", "TP")

#' Cancer types used for the TP pathway set
#'
#' Character vector of TCGA cancer-type abbreviations used to build the
#' `"TP"` (cancer-specific) pathway set in [default_pathway_list].
#'
#' @format A character vector of length 16.
#' @return A character vector of length 16.
#' @examples
#' default_TP_cancer
#' @export
default_TP_cancer <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "KIRC",
                      "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD",
                      "READ", "STAD", "THCA", "UCEC")

#' Supported GSEA score types
#'
#' Character vector of the score types supported by the `pvalueType`
#' argument of the scoring-based enrichment methods (`TG_Score`, `MiR_Score`,
#' `Edge_Score1D`, `Edge_Score2D`, `Edge_Topology`).
#'
#' @format A character vector of length 3.
#' @return A character vector of length 3.
#' @examples
#' default_pvalueType
#' @export
default_pvalueType <- c("pos", "std", "neg")

#' Supported gene-gene interaction sources
#'
#' Character vector of the background gene-gene interaction (GGI) sources
#' supported by the `GGI_source` argument of [getAllInputData()].
#'
#' @format A character vector of length 2.
#' @return A character vector of length 2.
#' @examples
#' default_GGI_source
#' @export
default_GGI_source <- c("Omnipath", "Reactome")

#' @noRd
combine_parallel <- function(res) {
  dplyr::bind_rows(lapply(res, as.data.frame))
}

#' @noRd
coerce_expr_mat <- function(x) {
  if (is.null(x) || !methods::is(x, "SummarizedExperiment")) return(x)
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required to pass a SummarizedExperiment as gene_mat/mir_mat.")
  }
  as.data.frame(SummarizedExperiment::assay(x))
}

#' @noRd
coerce_pathway <- function(pathway) {
  if (!methods::is(pathway, "GeneSetCollection")) return(pathway)
  if (!requireNamespace("GSEABase", quietly = TRUE)) {
    stop("Package 'GSEABase' is required to pass a GeneSetCollection as pathway.")
  }
  do.call(rbind, lapply(pathway, function(gs) {
    data.frame(pathway = GSEABase::setName(gs), gene = GSEABase::geneIds(gs), stringsAsFactors = FALSE)
  }))
}
