#' Title: MGI edge-based MicroRNA-oriented Enrichment Analysis (miREA)
#' Description: The main function to realize miREA which contains four node-based methods and five edge-based methods.
#' @param methods Method used for enrichment analysis, see default.method_list for all supportive options. Methods could be several, like c("Edge_ORA", "Edge_Score").
#' @param input_data Input data list. Contains all kind of input data for all methods, including data, pathway, pathway_GGI, GGI, dict_MGI. 
#'        Use get_all_input_data() function to generate. You could also generate it by yourself.
#'          input_data$data: each sub-element is named by the method, containing the input data format.
#'            TG_ORA: DEG, MiR_ORA: DEmiR, Edge_ORA: DEMGI, Edge_Network: DEMGI
#'              All of them is a character vector, containing all selected differential expressed items (genes/miRNAs/MGI)
#'            TG_Score: geneList, MiR_Score: miRList, Edge_Score: MGIList
#'              All of them is a named numeric vector, where the names are items (genes/miRNAs/MGI), corresponding value is their score, which is ranked by the score.
#'            Edge_Topology: MGIdata
#'              A five-column dataframe, containing MGI miRNA gene strength.
#'            Edge_2Ddist: MGI2D
#'              A five-column dataframe, including MGI miRNA gene cor(x) ranknormratio(y).
#'          input_data$pathway: each sub-element is named by the object, which could be TG, MiR, and Edge.
#'          input_data$pathway_GGI: a dataframe containing pathway and their directed GGI involved. Only useful when methods contain Edge_Topology and Edge_Network.
#'          input_data$MGI: a dataframe containing MGIs across all the pathways and DEMGIs. Only useful for Edge_Network.
#'          input_data$GGI: a dataframe containing directed GGIs across all the pathways and DEMGIs. Only useful for Edge_Network.
#'          input_data$dict_MGI: a dataframe containing all posible MGIs and their expression profiles. Only useful for Edge_Topology.
#'
#' @param background Contains background items (genes for TG_ method, miRNAs for MiR_ method, MGIs for Edge_ method).
#'                    If methods contains more than one targets, it should be a list where each sub-element is a character vector named by the target (TG, MiR, or Edge).
#'                    If methods contains only one target, it could also be a character vector.
#'                    If NULL, use the global background instead.
#' @param minSize The minimum pathway set sizes allowed for analysis, filtering out pathways with number of members less than the minSize threshold.
#' @param maxSize The maximum pathway set sizes allowed for analysis, filtering out pathways with number of members more than the maxSize threshold.
#' @param pvalueType How to calculate p value. See default.pvalueType for all possible options. Default pvalueType is list(TG = "std", MiR = "std", Edge = "neg").
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negative ("neg") enrichment).
#' The "std" is used when interseted both positive and negative scores, which is the same with the original GSEA algorithm.
#' Only used in Score-based method (Edge_2Ddist, Edge_TopoScore, TG_Score, MiR_Score).
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default.pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If you don't want any adjustment, please set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values below the specified pvalueCutoff. Default pvalueCutoff is 1.
#' @param iter Number of iteractions when perform permutation test. Default iter is 1000.
#' @param ncores Number of cores for parallel computing. Only useful if methods contain Edge_2Ddist, Edge_Topology, and Edge_Network.
#'        It accepts either a integer that is universal for all methods, or a list specify the number of cores for each method separately.
#'        Our recommendation is to specify different nodes for each method. E.g., list(Edge_2Ddist = 4, Edge_Topology = 32, Edge_Network = 8).
#'        If NULL, unparallel computing will be used.



#' @return A list for the enrichment result, contains result, time, valid methods, and parameter.
#' @examples

miREA <- function(methods, input_data, background = NULL,
                  minSize = NULL, maxSize = NULL, pvalueType = list(TG = "std", MiR = "std", Edge = "neg"),
                  pAdjMethod = "BH", pvalueCutoff = 1, iter = 1000, ncores = NULL){
  if (!is.list(input_data)) stop("input_data must be a list!")
  invalid_methods <- setdiff(methods, default.method_list)
  if (length(invalid_methods) > 0){
    warning("Please check the methods you have input. These methods cannot be specified: ",
            paste(invalid_methods, collapse = ", "), ". They will be removed.\n")
    methods <- intersect(methods, default.method_list)
    if (length(methods) == 0){
      stop("No valid methods! Please select one or more from default.method_list!")
    }
  }

  targets <- unique(sapply(strsplit(methods, "_"), `[`, 1))
  # ensure input_data have all input prepared.
  if (!"data" %in% names(input_data)) {
    stop("input_data is missing the 'data' component.")
  }
  missing_data <- methods[!methods %in% names(input_data$data)]
  if (length(missing_data) > 1){ #####################################需要改成>0
    stop("The following methods are missing in input_data$data: ",
         paste(missing_data, collapse = ", "))
  }


  if (!"pathway" %in% names(input_data)) {
    stop("input_data is missing the 'pathway' component.")
  }
  missing_pathway <- targets[!targets %in% names(input_data$pathway)]
  if (length(missing_pathway) > 0){
    stop("The following method targets are missing in input_data$pathway: ",
         paste(missing_pathway, collapse = ", "))
  }

  if ("Edge_Topology" %in% methods){
    if (! "pathway_GGI" %in% names(input_data)){
      stop("You have used topology method but didn't specified pathway_GGI in input_data!")
    } else if (! "dict_MGI" %in% names(input_data)){
      stop("You have used topology method but didn't specified dict_MGI in input_data!")
    }
  }

  if ("Edge_Network" %in% methods){
    if (! "pathway_GGI" %in% names(input_data)){
      stop("You have used network method but didn't specified pathway_GGI in input_data!")
    } else if (! "GGI" %in% names(input_data)){
      stop("You have used network method but didn't specified GGI in input_data!")
    } else if (! "MGI" %in% names(input_data)){
      stop("You have used network method but didn't specified MGI in input_data!")
    }
  }

  # check background
  if (!is.null(background)){
    if(!is.list(background) && length(targets) > 1){
      stop("You have selected methods contain multiple targets, please specify background for all the targets.")
    } else {
      if (!is.list(background) || any(!targets %in% names(background))){
        stop("The following background targets are missing in background: ",
             paste(targets[!targets %in% names(background)], collapse = ", "))
      }
    }
  } else {
    cat("\nNo background has specified. Use the global one instead!\n")
    # background_edge <- read.csv(paste0("/scratch/project_2011179/raw_data/background/old/",cancer, "_background.csv"), header = TRUE)
    load("data/raw_data/background/background_MGI.RData")
    background_edge <- paste0(background_MGI$miRNA, ":", background_MGI$gene)

    background <- list(TG = NULL,
                       MiR = unique(miRBaseConverter::getAllMiRNAs(version="v22", type="mature", species="hsa")$Name),
                       Edge = background_edge)
  }

  # check pvalueType
  missing_pvalueType <- targets[!targets %in% names(pvalueType)]
  if (length(missing_pvalueType) > 0){
    stop("The following method targets are missing in pvalueType: ",
         paste(missing_pvalueType, collapse = ", "))
  }

  # check ncores
  methods_ncores <- c("Edge_2Ddist", "Edge_Topology", "Edge_Network")
  methods_ncores <- intersect(methods_ncores, methods)
  if (!is.null(ncores) && is.list(ncores)){
    missing_ncores <- setdiff(methods_ncores, names(ncores))
    exist_ncores <- intersect(methods_ncores, names(ncores))
    if (length(missing_ncores) > 0){
      if (length(exist_ncores) > 0){
        warning("You have selected ", paste(missing_ncores, collapse = ", "), " methods. But you didn't customized the ncores used when parallel computing. Use the same as other methods.\n")
        max_core <- max(unlist(ncores[exist_ncores]))
        for (m in missing_ncores) {
          ncores[[m]] <- max_core
        }
      } else {
        warning("You have selected ", paste(missing_ncores, collapse = ", "), " methods. But you didn't customized the ncores used when parallel computing. Use unparallel computing instead.\n")
      }
    }
  } else if (is.integer(ncores) || is.numeric(ncores)){
    ncores <- as.integer(ncores)
    ncores <- setNames(as.list(rep(ncores, length(methods_ncores))), methods_ncores)
  } else if (is.null(ncores)){
    ncores <- setNames(as.list(rep(ncores, length(methods_ncores))), methods_ncores)
  }


  cat("Start Enrichment analysis...\n")
  cat("=== Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
  all_result <- list()
  time <- data.frame(method = character(), time = numeric(), stringsAsFactors = FALSE)
  n_pathway <- length(unique(input_data$pathway[[1]]$pathway))
  tryCatch({
    if ("TG_ORA" %in% methods){
      time_TG_ORA <- system.time({
        result_TG_ORA <- TG_ORA(TG = input_data$data$TG_ORA, background = background$TG, pathway = input_data$pathway$TG,
                                pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize)
      })
      #result_TG_ORA <- as.data.frame(result_TG_ORA)
      all_result[["TG_ORA"]] <- result_TG_ORA
      time <- rbind(time, data.frame(method = "TG_ORA", time = time_TG_ORA[["elapsed"]], stringsAsFactors = FALSE))
      cat("    - TG_ORA finished, time cost:", time_TG_ORA[["elapsed"]],"s. \n")
    }}, error = function(e) {
      message("TG_ORA error: ", e$message)
    })
  tryCatch({
    if ("TG_Score" %in% methods){
      time_TG_Score <- system.time({
        result_TG_Score <- TG_Score(
          geneList = input_data$data$TG_Score, pathway = input_data$pathway$TG, pvalueType = pvalueType$TG,
          pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize)
      })
      #result_TG_Score <- as.data.frame(result_TG_Score)
      all_result[["TG_Score"]] <- result_TG_Score
      time <- rbind(time, data.frame(method = "TG_Score", time = time_TG_Score[["elapsed"]], stringsAsFactors = FALSE))
      cat("    - TG_Score finished, time cost:", time_TG_Score[["elapsed"]],"s. \n")
    }
  }, error = function(e) {
    message("TG_Score error: ", e$message)
  })
  tryCatch({
    if ("MiR_ORA" %in% methods){
      time_MiR_ORA <- system.time({
        time_MiR_ORA <- system.time({
        result_MiR_ORA <- MiR_ORA(DEmiR = input_data$data$MiR_ORA, background = background$MiR, pathway = input_data$pathway$MiR,
                                  pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize)
      })
      })
      time <- rbind(time, data.frame(method = "MiR_ORA", time = time_MiR_ORA[["elapsed"]], stringsAsFactors = FALSE))
      #result_MiR_ORA <- as.data.frame(result_MiR_ORA)
      all_result[["MiR_ORA"]] <- result_MiR_ORA
      cat("    - MiR_ORA finished, time cost:", time_MiR_ORA[["elapsed"]],"s. \n")
    }
  }, error = function(e) {
    message("MiR_ORA error: ", e$message)
  })
  tryCatch({
    if ("MiR_Score" %in% methods){
      time_MiR_Score <- system.time({
        result_MiR_Score <- MiR_Score(miRList = input_data$data$MiR_Score, pathway = input_data$pathway$MiR, pvalueType = pvalueType$MiR,
                                      pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize)
      })
      #result_MiR_Score <- as.data.frame(result_MiR_Score)
      time <- rbind(time, data.frame(method = "MiR_Score", time = time_MiR_Score[["elapsed"]], stringsAsFactors = FALSE))
      all_result[["MiR_Score"]] <- result_MiR_Score
      cat("    - MiR_Score finished, time cost:", time_MiR_Score[["elapsed"]],"s. \n")
    }
  }, error = function(e) {
    message("MiR_Score error: ", e$message)
  })
  tryCatch({
    if ("Edge_ORA" %in% methods){
      time_Edge_ORA <- system.time({
        result_Edge_ORA <- Edge_ORA(DEMGI = input_data$data$Edge_ORA, background = background$Edge, pathway = input_data$pathway$Edge,
                                    pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize)
      })
      time <- rbind(time, data.frame(method = "Edge_ORA", time = time_Edge_ORA[["elapsed"]], stringsAsFactors = FALSE))
      all_result[["Edge_ORA"]] <- result_Edge_ORA
      cat("    - Edge_ORA finished, time cost:", time_Edge_ORA[["elapsed"]],"s. \n")
    }
  }, error = function(e) {
    message("Edge_ORA error: ", e$message)
  })
  tryCatch({
    if ("Edge_Score" %in% methods){
      time_Edge_Score <- system.time({
        result_Edge_Score <- Edge_Score(MGIList = input_data$data$Edge_Score, pathway = input_data$pathway$Edge, pvalueType = pvalueType$Edge,
                                        pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize)
      })
      #result_Edge_Score <- as.data.frame(result_Edge_Score)
      time <- rbind(time, data.frame(method = "Edge_Score", time = time_Edge_Score[["elapsed"]], stringsAsFactors = FALSE))
      all_result[["Edge_Score"]] <- result_Edge_Score
      cat("    - Edge_Score finished, time cost:", time_Edge_Score[["elapsed"]],"s. \n")
    }
  }, error = function(e) {
    message("Edge_Score error: ", e$message)
  })
  tryCatch({
    if ("Edge_2Ddist" %in% methods){
      time_Edge_2Ddist <- system.time({
        result_Edge_2Ddist <- Edge_2Ddist(MGI2D = input_data$data$Edge_2Ddist, pathway = input_data$pathway$Edge, background = background$Edge,
                                          iter = iter, ncores = ncores$Edge_2Ddist, pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize)
      })
      #result_Edge_2Ddist <- as.data.frame(result_Edge_2Ddist)
      time <- rbind(time, data.frame(method = "Edge_2Ddist", time = time_Edge_2Ddist[["elapsed"]], stringsAsFactors = FALSE))
      all_result[["Edge_2Ddist"]] <- result_Edge_2Ddist
      cat("    - Edge_2Ddist finished, time cost:", time_Edge_2Ddist[["elapsed"]],"s. \n")
    }
  }, error = function(e) {
    message("Edge_2Ddist error: ", e$message)
  })
  tryCatch({
    if ("Edge_manova" %in% methods){
      input <- input_data$data$Edge_2Ddist %>% dplyr::mutate(normcor = (1-cor)/2) %>% dplyr::select(MGI, miRNA, gene, normcor, normratio)
      time_Edge_manova <- system.time({
        result_Edge_manova <- Edge_manova(MGI2D = input, pathway = input_data$pathway$Edge, pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize)
      })
      #result_Edge_2Ddist <- as.data.frame(result_Edge_2Ddist)
      time <- rbind(time, data.frame(method = "Edge_manova", time = time_Edge_manova[["elapsed"]], stringsAsFactors = FALSE))
      all_result[["Edge_manova"]] <- result_Edge_manova
      cat("    - Edge_manova finished, time cost:", time_Edge_manova[["elapsed"]],"s. \n")
    }
  }, error = function(e) {
    message("Edge_manova error: ", e$message)
  })
  tryCatch({
    if ("Edge_Topology" %in% methods){
      time_Edge_Topology <- system.time({
        result_Edge_Topology <- Edge_Topology(MGI = input_data$data$Edge_Topology, background = background$Edge, pathway = input_data$pathway$Edge,
                                              pathway_GGI = input_data$pathway_GGI, dict_MGI = input_data$dict_MGI, ncores = ncores$Edge_Topology, pvalueType = pvalueType$Edge,
                                              iter = iter, pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize)
      })
      # result_Edge_Topology <- as.data.frame(result_Edge_Topology)
      time <- rbind(time, data.frame(method = "Edge_Topology", time = time_Edge_Topology[["elapsed"]], stringsAsFactors = FALSE))
      all_result[["Edge_Topology"]] <- result_Edge_Topology
      cat("    - Edge_Topology finished, time cost:", time_Edge_Topology[["elapsed"]],"s. \n")
    }
  }, error = function(e) {
    message("Edge_Topology error: ", e$message)
  })
  tryCatch({
    if ("Edge_Network" %in% methods){
      time_Edge_Network <- system.time({
        result_Edge_Network <- Edge_Network(DEMGI = input_data$data$Edge_Network,
                                            pathway_gene = input_data$pathway$TG, pathway_MGI = input_data$pathway$Edge,
                                            pathway_GGI = input_data$pathway_GGI,
                                            MGI = input_data$MGI, GGI = input_data$GGI, ncores = ncores$Edge_Network,
                                            pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, minSize = minSize, maxSize = maxSize, iter = iter)
      })
      #result_Edge_Network <- as.data.frame(result_Edge_Network)
      time <- rbind(time, data.frame(method = "Edge_Network", time = time_Edge_Network[["elapsed"]], stringsAsFactors = FALSE))
      all_result[["Edge_Network"]] <- result_Edge_Network
      cat("    - Edge_Network finished, time cost:", time_Edge_Network[["elapsed"]],"s. \n")
    }
  }, error = function(e) {
    message("Edge_Network error: ", e$message)
  })

  cat("\nmiREA analysis have finished!...\n")
  cat("=== End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")

  return(list(result = all_result, time = time, valid_methods = names(all_result),
              parameter = list(method = methods, n_pathway = n_pathway, minSize = minSize, maxSize = maxSize, pvalueType = pvalueType, pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff, iter = iter, ncores = ncores)))
}
