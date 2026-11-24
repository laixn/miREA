#' MGI edge-based microRNA-oriented enrichment analysis
#'
#' The main entry point for miREA, implementing four node-based methods and
#' five edge-based methods for miRNA-oriented enrichment analysis.
#'
#' @param methods Method used for enrichment analysis, see default_method_list for all supportive options. Methods could be several, like c("Edge_ORA", "Edge_Score1D").
#' @param input_data Input data list, containing the input data for all methods, including data, pathway, pathway_GGI, GGI, and dict_MGI.
#'        Use getAllInputData() function to generate. You could also generate it by yourself.
#'          input_data$data, each sub-element named by the method, contains the input data in the
#'            following formats.
#'            TG_ORA uses DEG, MiR_ORA uses DEmiR, Edge_ORA uses DEMGI, and Edge_Network uses DEMGI.
#'              Each is a character vector, containing all selected differential expressed items (genes/miRNAs/MGI).
#'            TG_Score uses geneList, MiR_Score uses miRList, and Edge_Score1D uses MGIList.
#'              Each is a named numeric vector, where the names are items (genes/miRNAs/MGI), and the corresponding value is their score, by which the vector is ranked.
#'            Edge_Topology uses MGIdata,
#'              a five-column dataframe, containing MGI miRNA gene strength.
#'            Edge_Score2D uses MGI2D,
#'              a five-column dataframe, including MGI miRNA gene cor(x) ranknormratio(y).
#'          input_data$pathway, each sub-element named by the object, which could be TG, MiR, and Edge.
#'          input_data$pathway_GGI, a dataframe containing pathway and their directed GGI involved. Only useful when methods contain Edge_Topology and Edge_Network.
#'          input_data$MGI, a dataframe containing MGIs across all the pathways and DEMGIs. Only useful for Edge_Network.
#'          input_data$GGI, a dataframe containing directed GGIs across all the pathways and DEMGIs. Only useful for Edge_Network.
#'          input_data$dict_MGI, a dataframe containing all possible MGIs and their expression profiles. Only useful for Edge_Topology.
#'
#' @param background Contains background items (genes for TG_ method, miRNAs for MiR_ method, MGIs for Edge_ method).
#'                    If methods contains more than one targets, it should be a list where each sub-element is a character vector named by the target (TG, MiR, or Edge).
#'                    If methods contains only one target, it could also be a character vector.
#'                    Required. miREA does not ship a bundled default background, so background = NULL raises an error.
#' @param minSize The minimum pathway set sizes allowed for analysis, filtering out pathways with number of members less than the minSize threshold.
#' @param maxSize The maximum pathway set sizes allowed for analysis, filtering out pathways with number of members more than the maxSize threshold.
#' @param pvalueType How to calculate p value. See default_pvalueType for all possible options. Default pvalueType is list(TG = "std", MiR = "std", Edge = "neg").
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negative ("neg") enrichment).
#' The "std" is used when interested in both positive and negative scores, which is the same as the original GSEA algorithm.
#' Only used in Score-based method (Edge_Score2D, Edge_Topology, TG_Score, MiR_Score).
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default_pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If no adjustment is desired, set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values above the specified pvalueCutoff. Default pvalueCutoff is 1.
#' @param iter Number of iterations when perform permutation test. Default iter is 1000.
#' @param ncores Number of cores for parallel computing. Only useful if methods contain Edge_Score2D, Edge_Topology, and Edge_Network.
#'        It accepts either a integer that is universal for all methods, or a list specify the number of cores for each method separately.
#'        Our recommendation is to specify different nodes for each method. E.g., list(Edge_Score2D = 4, Edge_Topology = 32, Edge_Network = 8).
#'        If NULL, unparallel computing will be used.
#' @param seed Optional integer seed for the permutation-based methods
#'        (Edge_Score2D, Edge_Topology, Edge_Network). If NULL (the default),
#'        no internal seed is set and the permutation draws depend on the
#'        caller's current random state. If supplied, it is forwarded to all
#'        three methods. Edge_Topology and Edge_Network make it reproducible
#'        via [withr::with_seed()]. Edge_Score2D does the same when running
#'        serially, or uses it as `BPPARAM`'s `RNGseed` when run in parallel
#'        (`ncores > 1`). None of these alter the caller's global random
#'        state.


#' @return A [miREAResult-class] object containing the per-method result tables
#' (`result`), per-method elapsed time (`time`), the methods that completed
#' without error (`valid_methods`), the methods that errored along with their
#' error message (`failed_methods`), and the call parameters (`parameter`).
#' `parameter$DE_threshold` records the `DElogFCCutoff`/`DEpadjCutoff`/`corCutoff`/
#' `corPadjCutoff` values used to build `input_data` (see [getAllInputData()]),
#' or `NULL` if `input_data` does not carry a `threshold` element.
#' Supports `$` access to these fields for backward compatibility.
#' @references Zhang Z, Lai X (2026). miREA: a network-based tool for
#' microRNA-oriented enrichment analysis. bioRxiv. \doi{10.64898/2026.02.27.708509}
#' @examples
#' data(toy_pathway)
#' data(toy_gene_DEdata)
#' data(toy_mir_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
#' data(background_MGI)
#' input_data <- getAllInputData(
#'   methods = default_method_list, pathway = toy_pathway,
#'   mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
#'   background_MGI = background_MGI, GGI_source = "Reactome",
#'   gene_mat = toy_gene_mat, mir_mat = toy_mir_mat,
#'   scoreFun = "rank"
#' )
#' background <- list(
#'   TG = unique(toy_pathway$gene),
#'   MiR = unique(input_data$pathway$MiR$miRNA),
#'   Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))
#' )
#' result <- miREA(
#'   methods = default_method_list, input_data = input_data, background = background,
#'   minSize = NULL, maxSize = NULL,
#'   pvalueType = list(TG = "std", MiR = "std", Edge = "neg"),
#'   pAdjMethod = "BH", pvalueCutoff = 1, iter = 10, ncores = 1
#' )
#' result$valid_methods
#' @export
miREA <- function(methods, input_data, background = NULL,
                  minSize = NULL, maxSize = NULL, pvalueType = list(TG = "std", MiR = "std", Edge = "neg"),
                  pAdjMethod = "BH", pvalueCutoff = 1, iter = 1000, ncores = NULL, seed = NULL) {
  if (!is.list(input_data)) stop("input_data must be a list!")
  methods <- miREA_validate_methods(methods)
  targets <- unique(vapply(strsplit(methods, "_"), `[`, character(1), 1))
  miREA_validate_input_data(input_data, methods, targets)
  background <- miREA_validate_background(background, targets)
  miREA_validate_pvalueType(pvalueType, targets)
  ncores <- miREA_validate_ncores(ncores, methods)

  message("Start Enrichment analysis...")
  message("=== Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ===")
  all_result <- list()
  failed_methods <- list()
  time <- data.frame(method = character(), time = numeric(), stringsAsFactors = FALSE)
  n_pathway <- length(unique(input_data$pathway[[1]]$pathway))

  method_calls <- miREA_method_calls(input_data, background, pvalueType, pAdjMethod, pvalueCutoff, minSize, maxSize, iter, ncores, seed)
  for (m in names(method_calls)) {
    if (!m %in% methods) next
    run <- run_miREA_method(m, method_calls[[m]])
    if (!is.null(run$error)) {
      failed_methods[[m]] <- run$error
    } else {
      all_result[[m]] <- run$result
      time <- rbind(time, data.frame(method = m, time = run$elapsed, stringsAsFactors = FALSE))
    }
  }

  message("\nmiREA analysis have finished!...")
  message("=== End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ===")

  miREAResult(
    result = all_result, time = time, valid_methods = names(all_result),
    failed_methods = failed_methods,
    parameter = list(method = methods, n_pathway = n_pathway, minSize = minSize, maxSize = maxSize, pvalueType = pvalueType, pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, iter = iter, ncores = ncores, seed = seed, DE_threshold = input_data[["threshold"]])
  )
}

#' @noRd
miREA_validate_methods <- function(methods) {
  invalid_methods <- setdiff(methods, default_method_list)
  if (length(invalid_methods) > 0) {
    warning(
      "Please check the methods you have input. These methods cannot be specified: ",
      paste(invalid_methods, collapse = ", "), ". They will be removed.\n"
    )
    methods <- intersect(methods, default_method_list)
    if (length(methods) == 0) {
      stop("No valid methods! Please select one or more from default_method_list!")
    }
  }
  methods
}

#' @noRd
miREA_validate_input_data <- function(input_data, methods, targets) {
  # ensure input_data have all input prepared.
  if (!"data" %in% names(input_data)) {
    stop("input_data is missing the 'data' component.")
  }
  missing_data <- methods[!methods %in% names(input_data$data)]
  if (length(missing_data) > 0) {
    stop(
      "The following methods are missing in input_data$data: ",
      paste(missing_data, collapse = ", ")
    )
  }

  if (!"pathway" %in% names(input_data)) {
    stop("input_data is missing the 'pathway' component.")
  }
  missing_pathway <- targets[!targets %in% names(input_data$pathway)]
  if (length(missing_pathway) > 0) {
    stop(
      "The following method targets are missing in input_data$pathway: ",
      paste(missing_pathway, collapse = ", ")
    )
  }

  if ("Edge_Topology" %in% methods) {
    if (!"pathway_GGI" %in% names(input_data)) {
      stop("You have used topology method but didn't specified pathway_GGI in input_data!")
    } else if (!"dict_MGI" %in% names(input_data)) {
      stop("You have used topology method but didn't specified dict_MGI in input_data!")
    }
  }

  if ("Edge_Network" %in% methods) {
    if (!"pathway_GGI" %in% names(input_data)) {
      stop("You have used network method but didn't specified pathway_GGI in input_data!")
    } else if (!"GGI" %in% names(input_data)) {
      stop("You have used network method but didn't specified GGI in input_data!")
    } else if (!"MGI" %in% names(input_data)) {
      stop("You have used network method but didn't specified MGI in input_data!")
    }
  }
  invisible(NULL)
}

#' @noRd
miREA_validate_background <- function(background, targets) {
  if (!is.null(background)) {
    if (!is.list(background)) {
      if (length(targets) > 1) {
        stop("You have selected methods contain multiple targets, please specify background for all the targets.")
      }
      background <- setNames(list(background), targets)
    } else if (any(!targets %in% names(background))) {
      stop(
        "The following background targets are missing in background: ",
        paste(targets[!targets %in% names(background)], collapse = ", ")
      )
    }
  } else {
    stop(
      "background must be supplied for the selected methods. miREA does not ship ",
      "a bundled default background; see the package vignette for how to build ",
      "one (e.g. from background_MGI for Edge_ methods, or miRBaseConverter::",
      "getAllMiRNAs() for MiR_ methods), or download the full background used in ",
      "the miREA manuscript from its Zenodo deposit."
    )
  }
  background
}

#' @noRd
miREA_validate_pvalueType <- function(pvalueType, targets) {
  missing_pvalueType <- targets[!targets %in% names(pvalueType)]
  if (length(missing_pvalueType) > 0) {
    stop(
      "The following method targets are missing in pvalueType: ",
      paste(missing_pvalueType, collapse = ", ")
    )
  }
  invisible(NULL)
}

#' @noRd
miREA_validate_ncores <- function(ncores, methods) {
  methods_ncores <- c("Edge_Score2D", "Edge_Topology", "Edge_Network")
  methods_ncores <- intersect(methods_ncores, methods)
  if (!is.null(ncores) && is.list(ncores)) {
    missing_ncores <- setdiff(methods_ncores, names(ncores))
    exist_ncores <- intersect(methods_ncores, names(ncores))
    if (length(missing_ncores) > 0) {
      if (length(exist_ncores) > 0) {
        warning("You have selected ", paste(missing_ncores, collapse = ", "), " methods. But you didn't customized the ncores used when parallel computing. Use the same as other methods.\n")
        max_core <- max(unlist(ncores[exist_ncores]))
        for (m in missing_ncores) {
          ncores[[m]] <- max_core
        }
      } else {
        warning("You have selected ", paste(missing_ncores, collapse = ", "), " methods. But you didn't customized the ncores used when parallel computing. Use unparallel computing instead.\n")
      }
    }
  } else if (is.integer(ncores) || is.numeric(ncores)) {
    if (length(ncores) != 1 ||
      is.na(ncores) ||
      !is.finite(ncores) ||
      ncores < 1 ||
      ncores != as.integer(ncores)) {
      stop("ncores must be a positive integer.")
    }
    ncores <- as.integer(ncores)
    ncores <- setNames(as.list(rep(ncores, length(methods_ncores))), methods_ncores)
  } else if (is.null(ncores)) {
    message("You didn't specify number of cores [ncores]. Use unparallel computing instead.\n")
    ncores <- setNames(as.list(rep(1L, length(methods_ncores))), methods_ncores)
  }
  ncores
}

#' @noRd
miREA_method_calls <- function(input_data, background, pvalueType, pAdjMethod, pvalueCutoff, minSize, maxSize, iter, ncores, seed) {
  list(
    TG_ORA = function() {
      TG_ORA(
        TG = input_data$data$TG_ORA, background = background$TG, pathway = input_data$pathway$TG,
        pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize
      )
    },
    TG_Score = function() {
      TG_Score(
        geneList = input_data$data$TG_Score, pathway = input_data$pathway$TG, pvalueType = pvalueType$TG,
        pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize
      )
    },
    MiR_ORA = function() {
      MiR_ORA(
        DEmiR = input_data$data$MiR_ORA, background = background$MiR, pathway = input_data$pathway$MiR,
        pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize
      )
    },
    MiR_Score = function() {
      MiR_Score(
        miRList = input_data$data$MiR_Score, pathway = input_data$pathway$MiR, pvalueType = pvalueType$MiR,
        pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize
      )
    },
    Edge_ORA = function() {
      Edge_ORA(
        DEMGI = input_data$data$Edge_ORA, background = background$Edge, pathway = input_data$pathway$Edge,
        pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize
      )
    },
    Edge_Score1D = function() {
      Edge_Score1D(
        MGIList = input_data$data$Edge_Score1D, pathway = input_data$pathway$Edge, pvalueType = pvalueType$Edge,
        pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize
      )
    },
    Edge_Score2D = function() {
      Edge_Score2D(
        MGI2D = input_data$data$Edge_Score2D, pathway = input_data$pathway$Edge, background = background$Edge,
        iter = iter, ncores = ncores$Edge_Score2D, pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize, seed = seed
      )
    },
    Edge_Topology = function() {
      Edge_Topology(
        MGI = input_data$data$Edge_Topology, background = background$Edge, pathway = input_data$pathway$Edge,
        pathway_GGI = input_data$pathway_GGI, dict_MGI = input_data$dict_MGI, ncores = ncores$Edge_Topology, pvalueType = pvalueType$Edge,
        iter = iter, pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize, seed = seed
      )
    },
    Edge_Network = function() {
      Edge_Network(
        DEMGI = input_data$data$Edge_Network,
        pathway_gene = input_data$pathway$TG, pathway_MGI = input_data$pathway$Edge,
        pathway_GGI = input_data$pathway_GGI,
        MGI = input_data$MGI, GGI = input_data$GGI, ncores = ncores$Edge_Network,
        pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize, iter = iter, seed = seed
      )
    }
  )
}

#' @noRd
run_miREA_method <- function(name, call_fn) {
  tryCatch(
    {
      elapsed <- system.time({
        result <- call_fn()
      })[["elapsed"]]
      message("    - ", name, " finished, time cost: ", elapsed, " s. ")
      list(result = result, elapsed = elapsed, error = NULL)
    },
    error = function(e) {
      message("    - ", name, " failed: ", e$message)
      list(result = NULL, elapsed = NA_real_, error = e$message)
    }
  )
}
