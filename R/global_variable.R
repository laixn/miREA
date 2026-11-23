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
default_TP_cancer <- c(
  "BLCA", "BRCA", "CESC", "COAD", "ESCA", "KIRC",
  "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD",
  "READ", "STAD", "THCA", "UCEC"
)

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
fetch_background_GGI_omnipath <- function() {
  if (!requireNamespace("OmnipathR", quietly = TRUE)) {
    stop(
      "Package 'OmnipathR' is required to fetch a default background_GGI from Omnipath. ",
      "Please install OmnipathR, or supply background_GGI directly."
    )
  }
  # Returns Omnipath's full interaction set, with no license/source filtering:
  # Omnipath aggregates many underlying databases under different licenses,
  # some of which restrict commercial use (see README "Third-party data
  # licensing"). Users who need a license-restricted subset should filter it
  # themselves (e.g. via OmnipathR::resource_info()) and supply the result
  # via background_GGI.
  interactions <- tryCatch(
    OmnipathR::all_interactions(),
    error = function(e) {
      stop(
        "Failed to fetch a default background_GGI from Omnipath (", conditionMessage(e), "). ",
        "Please supply background_GGI directly instead."
      )
    }
  )

  interactions <- interactions[, c(
    "source_genesymbol", "target_genesymbol", "is_directed",
    "is_stimulation", "is_inhibition", "n_references"
  )]
  interactions <- interactions[interactions$n_references != 0, ]
  interactions <- interactions[!grepl(",", interactions$source_genesymbol) &
    !grepl("_", interactions$source_genesymbol) &
    !grepl("_", interactions$target_genesymbol), ]
  undirected <- interactions[interactions$is_directed == 0, c(2, 1, 3:6)]
  names(undirected) <- names(interactions)
  interactions <- rbind(undirected, interactions)
  interactions$interact <- interactions$is_stimulation - interactions$is_inhibition
  interactions <- interactions[, c("source_genesymbol", "target_genesymbol", "interact")]
  interactions[!grepl("^hsa-", interactions$source_genesymbol) &
    !grepl("^hsa-", interactions$target_genesymbol), ]
}

#' @noRd
fetch_background_GGI_reactome <- function() {
  # Pinned to the exact dated Reactome FI release used to build the
  # manuscript's background_GGI_Reactome (reactome.org hosts past releases
  # permanently at stable URLs, so this snapshot will not drift over time,
  # unlike the live Omnipath fetch above). Reactome's data is CC0, so no
  # license filtering is needed here.
  url <- "https://reactome.org/download/tools/ReactomeFIs/FIsInGene_04142025_with_annotations.txt.zip"
  tmp_zip <- tempfile(fileext = ".zip")
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(c(tmp_zip, tmp_dir), recursive = TRUE), add = TRUE)

  tryCatch(
    {
      utils::download.file(url, destfile = tmp_zip, quiet = TRUE, mode = "wb")
      utils::unzip(tmp_zip, exdir = tmp_dir)
    },
    error = function(e) {
      stop(
        "Failed to fetch a default background_GGI from Reactome (", conditionMessage(e), "). ",
        "Please supply background_GGI directly instead."
      )
    }
  )

  txt_file <- list.files(tmp_dir, pattern = "_with_annotations\\.txt$", full.names = TRUE)
  if (length(txt_file) == 0) {
    stop(
      "Failed to fetch a default background_GGI from Reactome: downloaded archive did not ",
      "contain the expected annotation file. Please supply background_GGI directly instead."
    )
  }
  reactome_ggi <- utils::read.table(txt_file[1],
    sep = "\t", header = TRUE, quote = "",
    stringsAsFactors = FALSE
  )

  # Undirected/mixed-sign relations are split into two directed rows, exactly
  # as in the original pathway_generate.Rmd expand_direction() logic.
  expand_direction <- function(gene1, gene2, dir) {
    switch(dir,
      "-"   = list(c(gene1, gene2, "-"), c(gene2, gene1, "-")),
      "<->" = list(c(gene1, gene2, "<-"), c(gene2, gene1, "->")),
      "|-|" = list(c(gene1, gene2, "|-"), c(gene2, gene1, "-|")),
      "|->" = list(c(gene1, gene2, "|-"), c(gene2, gene1, "->")),
      "<-|" = list(c(gene1, gene2, "<-"), c(gene2, gene1, "-|")),
      list(c(gene1, gene2, dir))
    )
  }
  expanded <- do.call(rbind, unlist(
    Map(expand_direction, reactome_ggi$Gene1, reactome_ggi$Gene2, reactome_ggi$Direction),
    recursive = FALSE
  ))
  result <- as.data.frame(expanded, stringsAsFactors = FALSE, row.names = seq_len(nrow(expanded)))
  colnames(result) <- c("source_gene", "target_gene", "direction")
  result
}

#' @noRd
combine_parallel <- function(res) {
  dplyr::bind_rows(lapply(res, as.data.frame))
}

#' @noRd
resolve_BPPARAM <- function(ncores, seed = NULL, BPPARAM = NULL) {
  if (!is.null(BPPARAM)) {
    return(BPPARAM)
  }
  if (is.null(ncores) || ncores <= 1) {
    BiocParallel::SerialParam(RNGseed = seed)
  } else {
    BiocParallel::SnowParam(workers = ncores, RNGseed = seed)
  }
}

#' @noRd
coerce_expr_mat <- function(x, assay = 1) {
  if (is.null(x) || !methods::is(x, "SummarizedExperiment")) {
    return(x)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required to pass a SummarizedExperiment as gene_mat/mir_mat.")
  }
  as.data.frame(SummarizedExperiment::assay(x, i = assay))
}

#' @noRd
filter_pathway_by_size <- function(pathway, minSize = NULL, maxSize = NULL) {
  pathway_counts <- table(pathway$pathway)
  # NOTE: `pathway` has a column also named "pathway", so `pathway$pathway`
  # inside `pathway[...]` gets NSE-resolved against the *column* (an atomic
  # vector) first, raising "$ operator is invalid for atomic vectors". Use
  # the bare column reference instead.
  if (!is.null(minSize)) {
    pathways_to_remove_min <- names(pathway_counts[pathway_counts < minSize])
    pathway <- pathway[!pathway %in% pathways_to_remove_min, ]
  }
  if (!is.null(maxSize)) {
    pathways_to_remove_max <- names(pathway_counts[pathway_counts > maxSize])
    pathway <- pathway[!pathway %in% pathways_to_remove_max, ]
  }
  pathway
}

#' @noRd
coerce_pathway <- function(pathway) {
  if (!methods::is(pathway, "GeneSetCollection")) {
    return(pathway)
  }
  if (!requireNamespace("GSEABase", quietly = TRUE)) {
    stop("Package 'GSEABase' is required to pass a GeneSetCollection as pathway.")
  }
  do.call(rbind, lapply(pathway, function(gs) {
    data.frame(pathway = GSEABase::setName(gs), gene = GSEABase::geneIds(gs), stringsAsFactors = FALSE)
  }))
}
