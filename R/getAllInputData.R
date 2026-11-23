# if TG_method, initially use the differential expressed miRNAs (interested miRNAs), and link to genes (as interested genes, or geneList (need gene_data)).
# if MiR_method, use mir_DEdata directly
# if Edge_method, need both mir_DEdata and gene_DEdata

#' Generate input data for miREA
#'
#' Functions to generate input data for miREA, based on the original expression and differential expression results for miRNAs and genes.
#'
#' @param methods compulsory, see default_method_list for all supportive methods. You could either input one method, or a vector contains several methods you want to try.
#' @param pathway compulsory. It should be a two-column pathway, containing pathway and gene, or a GSEABase GeneSetCollection object.
#' @param mir_DEdata compulsory, should be a four-column dataframe, containing miRNA, log2FC, stat, padj
#' @param gene_DEdata optional, except for TG-ORA, MiR-ORA, and MiR-Score methods, should be a four-column dataframe, containing gene, log2FC, stat, padj
#' @param background_MGI optional. It should be either a character vector in the format of [miRNA]:[gene], or a two-column dataframe, containing miRNA gene. If NULL, miREA falls back to its bundled `background_MGI` dataset (the full manuscript background).
#' @param background_GGI optional, it should be a two-column data frame that contains source gene and target gene. Combined with *GGI_source*. If NULL, miREA falls back to its bundled `background_GGI_Omnipath`/`background_GGI_Reactome` dataset (the full manuscript background for that source), selected via GGI_source. These are frozen snapshots; to fetch a fresh/current background instead, see `inst/scripts/build_background_GGI.R`.
#' @param GGI_source optional, select from "Omnipath" or "Reactome". Selects which bundled background_GGI snapshot to fall back to when background_GGI is not provided. The default is Omnipath.
#' @param gene_mat optional, needed for Edge_ method, used for calculating correlations between miRNAs and genes. It should be a raw gene expression matrix (each row a gene, each column a sample), either as a dataframe/matrix or as a SummarizedExperiment object.
#' @param mir_mat optional, needed for Edge_ method, used for calculating correlations between miRNAs and genes. It should be a raw miRNA expression matrix (each row a miRNA, each column a sample), either as a dataframe/matrix or as a SummarizedExperiment object.
#' @param scoreFun optional, needed for Edge_Score1D, Edge_Score2D, Edge_Topology. See default_scoreFun for all supportive methods.
#' @note gene_mat and mir_mat should have same sample name format. If the number of shared samples below 50 % of the minimum of gene samples and miRNA samples, a warning would happen to remind you recheck them.
#' @return A list that contains all the data formats needed for miREA.
#' @references Zhang Z, Lai X (2026). miREA: a network-based tool for
#' microRNA-oriented enrichment analysis. bioRxiv. \doi{10.64898/2026.02.27.708509}
#' @examples
#' data(toy_pathway)
#' data(toy_gene_DEdata)
#' data(toy_mir_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
#' data(toy_background_MGI)
#' data(toy_background_GGI)
#' input_data <- getAllInputData(
#'   methods = default_method_list,
#'   pathway = toy_pathway,
#'   mir_DEdata = toy_mir_DEdata,
#'   gene_DEdata = toy_gene_DEdata,
#'   background_MGI = toy_background_MGI,
#'   background_GGI = toy_background_GGI,
#'   GGI_source = "Omnipath",
#'   gene_mat = toy_gene_mat,
#'   mir_mat = toy_mir_mat,
#'   scoreFun = "rank"
#' )
#' names(input_data)
#' @export
getAllInputData <- function(methods, pathway,
                               mir_DEdata, gene_DEdata = NULL,
                               background_MGI = NULL, background_GGI = NULL, GGI_source = "Omnipath",
                               gene_mat = NULL, mir_mat = NULL,
                               scoreFun = NULL){
  message("Generate all input data...")
  pathway <- coerce_pathway(pathway)
  gene_mat <- coerce_expr_mat(gene_mat)
  mir_mat <- coerce_expr_mat(mir_mat)
  sum_TG <- sum(grepl("TG_", methods))
  sum_MiR <- sum(grepl("MiR_", methods))
  sum_Edge <- sum(grepl("Edge_", methods))
  target_method <- c("Edge_Topology", "Edge_Network")
  score_method <- c("Edge_Score1D", "Edge_Score2D", "Edge_Topology")
  colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
  # check scoreFun
  if (length(intersect(score_method, methods)) > 0 && is.null(scoreFun)){
    warning("No [scoreFun] parameter has been specified. Use the default [rank] instead!\n")
    scoreFun <- "rank"
  }
  if (length(intersect(score_method, methods)) > 0 && (! scoreFun %in% default_scoreFun)){
    stop("The [scoreFun] you input is not valid. Please select one from the default_scoreFun!")
  }

  result <- list()

  # check background_MGI
  # NOTE: nearly every method (all but MiR_Score) derives its input via
  # background_MGI somewhere below (e.g. TG_ORA maps significant miRNAs to
  # target genes through it), so this is required unconditionally rather
  # than only for Edge_/MiR_ methods.
  if (is.null(background_MGI)){
    message("No background_MGI is given. Using the bundled default background_MGI instead.\n")
    utils::data("background_MGI", package = "miREA", envir = environment())
  }
  colnames(background_MGI)[seq_len(2)] <- c("miRNA", "gene")

  # check background_GGI
  if (length(intersect(target_method, methods)) != 0){
    if (is.null(background_GGI)){
      if (! GGI_source %in% default_GGI_source){
        stop("You tend to extract default background_GGI for ", GGI_source, ". But [GGI_source] you have specified are not involved. Please use default_GGI_source to see all supporting options!")
      }
      bundled_env <- new.env()
      if (GGI_source == "Omnipath"){
        message("No background_GGI is given. Using the bundled default background_GGI_Omnipath instead.\n")
        utils::data("background_GGI_Omnipath", package = "miREA", envir = bundled_env)
        background_GGI <- bundled_env[["background_GGI_Omnipath"]]
      } else {
        message("No background_GGI is given. Using the bundled default background_GGI_Reactome instead.\n")
        utils::data("background_GGI_Reactome", package = "miREA", envir = bundled_env)
        background_GGI <- bundled_env[["background_GGI_Reactome"]]
      }
    }
    colnames(background_GGI)[seq_len(2)] <- c("from", "to")
  }

  # check expression data
  if (!is.null(gene_DEdata)){
    if (ncol(gene_DEdata) != 4 || !is.data.frame(gene_DEdata)){
      stop("Please make sure the input gene_DEdata is a four-column dataframe contains gene, log2FC, stat, and padj!")
    }
    colnames(gene_DEdata) <- c("gene", "geneLogFC", "geneStat", "genePadj")
    gene_DEdata <- gene_DEdata %>% distinct(gene, .keep_all = TRUE)
    gene_DEdata <- data.table::as.data.table(gene_DEdata)
  }

  if (ncol(mir_DEdata) != 4 || !is.data.frame(mir_DEdata)){
    stop("Please make sure the input mir_DEdata is a four-column dataframe contains miRNA, log2FC, stat, and padj!")
  }
  colnames(mir_DEdata) <- c("miRNA", "mirLogFC", "mirStat", "mirPadj")
  mir_DEdata <- mir_DEdata %>% distinct(miRNA, .keep_all = TRUE)
  mir_DEdata <- data.table::as.data.table(mir_DEdata)

  # ensure that mir_mat and gene_mat have their names as rownames
  if (!is.null(mir_mat) && !is.null(gene_mat)){
    if (!is.numeric(mir_mat[[1]]) && !is.integer(mir_mat[[1]])) {
      # we suppose the first column is miRNA name in that case.
      colnames(mir_mat)[1] <- "miRNA"
      mir_mat <- mir_mat %>% distinct(.[1], .keep_all = TRUE) # only keep the first appearance if duplicated values.
      rownames(mir_mat) <- mir_mat[[1]]
      mir_mat <- mir_mat[, -1, drop = FALSE]
    }
    if (!is.numeric(gene_mat[[1]]) && !is.integer(gene_mat[[1]])) {
      colnames(gene_mat)[1] <- "gene"
      gene_mat <- gene_mat %>% distinct(.[1], .keep_all = TRUE)
      rownames(gene_mat) <- gene_mat[[1]]
      gene_mat <- gene_mat[, -1, drop = FALSE]
    }
  }
  

  raw_data <- list(
    methods = methods,
    mir_DEdata = mir_DEdata, gene_DEdata = gene_DEdata,
    gene_mat = gene_mat, mir_mat = mir_mat,
    scoreFun = scoreFun
    # check_mir_names = check_mir_names, mir_version = mir_version,
    # check_gene_names = check_gene_names
  )
  result[["raw_data"]] <- raw_data
  result[["background_MGI"]] <- background_MGI
  result[["background_GGI"]] <- background_GGI


  message("  Start generate expression data...")
  data <- get_data(methods = methods, mir_DEdata = mir_DEdata, gene_DEdata = gene_DEdata,
                   background_MGI = background_MGI, gene_mat = gene_mat, mir_mat = mir_mat, scoreFun = scoreFun)
  result[["data"]] <- data

  message("  Start generate pathway data...")
  pathway_data <- get_pathway_data(methods = methods, pathway = pathway, background_MGI = background_MGI)
  result[["pathway"]] <- pathway_data

  if (length(intersect(target_method, methods)) != 0){
    message("  Start generate pathway_GGI data...")
    pathway_GGI <- get_pathway_GGI_data(methods = methods, pathway = pathway, background_GGI = background_GGI)
    result[["pathway_GGI"]] <- pathway_GGI

    if ("Edge_Network" %in% methods){
      message("  Start generate GGI data for whole network...")
      GGI <- get_GGI_data(methods = methods, DEMGI = result[["data"]][["Edge_Network"]], pathway = pathway, background_GGI = background_GGI)
      result[["GGI"]] <- GGI
      message("  Start generate MGI data for whole network...")
      MGI <- get_MGI_data(
        methods = methods, DEMGI = result[["data"]][["Edge_Network"]], pathway = pathway,
        pathway_MGI = result[["pathway"]][["Edge"]], background_MGI = background_MGI)
      result[["MGI"]] <- MGI
    }

    if ("Edge_Topology" %in% methods){
      message("  Start generate dict_MGI data for permutation...")
      dict_MGI <- get_dict_MGI(methods = methods, mir_mat = mir_mat, gene_mat = gene_mat, data = data$data,
                               pathway_MGI = pathway_data$Edge, scoreFun = scoreFun)
      result[["dict_MGI"]] <- dict_MGI
    }

  }

  message("All the input data is ready!")
  return(result)
}

## get all input data, including:
### data, pathway, pathway_GGI, dict_MGI, dict_GGI
## all the original input include:
### method, gene_DEdata, mir_DEdata, pathway, gene_mat, mir_mat, background_MGI, background_GGI, scoreFun
### Notice: gene_mat, mir_mat should only conclude tumor samples


# gene_mat: Edge_ORA
# mir_mat: Edge_ORA
# gene_DEdata: TG_ORA, TG_Score
# mir_DEdata: MiR_ORA, MiR_Score
# background_MGI: could be either two-column dataframe, or a character vector with format [miRNA]:[gene]
#                 Edge_ORA, Edge_Score1D, Edge_Topology, Edge_Network
# background_GGI
# pathway
# scoreFun: Edge_Score1D

# gene_DEdata: except for TG_ORA, MiR_ORA and MiR_Score, should be four-column dataframe, containing gene, log2FC, stat, padj
# mir_DEdata: compulsory, should be four-column dataframe, containing miRNA, log2FC, stat, padj
# gene_mat: only needed for Edge_ method
# mir_mat: only needed for Edge_ method
# scoreFun: only needed for Edge_Score1D, Edge_Score2D, Edge_Topology, choose one from rank, sigmoid, normalize
# background_MGI: only needed for Edge_ method
# will return a list, contains the input of each method
get_data <- function(methods, mir_DEdata, gene_DEdata = NULL,
                     background_MGI, gene_mat = NULL, mir_mat = NULL, scoreFun){
  scoreFun <- get_data_validate_inputs(methods, gene_DEdata, gene_mat, mir_mat, background_MGI, scoreFun)

  coerced <- get_data_coerce_expr(gene_DEdata, mir_DEdata)
  gene_DEdata <- coerced$gene_DEdata
  mir_DEdata <- coerced$mir_DEdata
  background_MGI <- get_data_coerce_background_mgi(background_MGI)

  result <- get_data_node_based(methods, mir_DEdata, gene_DEdata, background_MGI)

  if (sum(grepl("Edge_", methods)) > 0){
    result <- c(result, get_data_edge_based(methods, mir_DEdata, gene_DEdata, background_MGI, mir_mat, gene_mat, scoreFun))
  }
  result
}

#' @noRd
get_data_validate_inputs <- function(methods, gene_DEdata, gene_mat, mir_mat, background_MGI, scoreFun) {
  invalid_method <- base::setdiff(methods, default_method_list)
  if (length(invalid_method) > 0){
    stop("The method(s):", invalid_method, "you are input is(are) not one of the supportive methods. Please revise the method(s)!")
  }
  sum_Edge <- sum(grepl("Edge_", methods))
  no_gene_DEdata <- c("TG_ORA", "MiR_ORA", "MiR_Score")
  score_method <- c("Edge_Score1D", "Edge_Score2D", "Edge_Topology")
  if (length(setdiff(methods, no_gene_DEdata)) > 0 && is.null(gene_DEdata)) {
    stop("Please make sure you have input gene_DEdata!")
  }
  if (sum_Edge != 0 && (is.null(gene_mat) || is.null(mir_mat))){
    stop("Please make sure you have input gene_mat and mir_mat")
  }
  if (any(score_method %in% methods) && is.null(scoreFun)){
    warning("You didn't specify scoreFun for transform expression score, set 'rank' automatically.\n")
    scoreFun <- "rank"
  }
  # NOTE: unreachable in practice — getAllInputData() resolves background_MGI
  # (falling back to the bundled default) before calling get_data(), which is
  # the only caller of this function. Kept as a defensive safety net.
  if (is.null(background_MGI)){
    stop("background_MGI must be supplied: a two-column data.frame of miRNA, gene ",
         "pairs, or use the bundled default background_MGI via getAllInputData().")
  }
  scoreFun
}

#' @noRd
get_data_coerce_expr <- function(gene_DEdata, mir_DEdata) {
  if (!is.null(gene_DEdata)){
    if (ncol(gene_DEdata) != 4 || !is.data.frame(gene_DEdata)){
      stop("Please make sure the input gene_DEdata is a four-column dataframe contains gene, log2FC, stat, and padj!")
    }
    colnames(gene_DEdata) <- c("gene", "geneLogFC", "geneStat", "genePadj")
    gene_DEdata <- data.table::as.data.table(gene_DEdata)
  }

  if (ncol(mir_DEdata) != 4 || !is.data.frame(mir_DEdata)){
    stop("Please make sure the input mir_DEdata is a four-column dataframe contains miRNA, log2FC, stat, and padj!")
  }
  colnames(mir_DEdata) <- c("miRNA", "mirLogFC", "mirStat", "mirPadj")
  mir_DEdata <- data.table::as.data.table(mir_DEdata)

  list(gene_DEdata = gene_DEdata, mir_DEdata = mir_DEdata)
}

#' @noRd
get_data_coerce_background_mgi <- function(background_MGI) {
  if (is.data.frame(background_MGI)){
    if (ncol(background_MGI) != 2){
      stop("Please make sure the input background_MGI is a two-column dataframe contains miRNA and gene")
    }
    colnames(background_MGI) <- c("miRNA", "gene")
  } else if (is.character(background_MGI) && is.vector(background_MGI)){
    background_MGI <- data.frame(
      miRNA = sub(":.*", "", background_MGI),
      gene = sub(".*:", "", background_MGI),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Please make sure the input background_MGI is either a two-column dataframe contains miRNA and gene,
         or a character vector with each element in the format of [miRNA]:[gene]!")
  }
  colnames(background_MGI) <- c("miRNA", "gene")
  data.table::as.data.table(background_MGI)
}

#' @noRd
get_data_node_based <- function(methods, mir_DEdata, gene_DEdata, background_MGI) {
  result <- list()
  if (sum(methods %in% c("TG_ORA", "TG_Score", "MiR_ORA")) > 0){
    interested_miRNAs <- mir_DEdata[abs(mirLogFC) > 1 & mirPadj < 0.05, unique(miRNA)]
  }

  if ("TG_ORA" %in% methods){
    result[["TG_ORA"]] <- background_MGI[miRNA %in% interested_miRNAs, unique(gene)]
  }
  if ("TG_Score" %in% methods){
    genes <- background_MGI[miRNA %in% interested_miRNAs, unique(gene)]
    data <- gene_DEdata[gene %in% genes, setNames(geneStat, gene)]
    result[["TG_Score"]] <- sort(data, decreasing = TRUE)
  }
  if ("MiR_ORA" %in% methods){
    result[["MiR_ORA"]] <- interested_miRNAs
  }
  if ("MiR_Score" %in% methods){
    data <- mir_DEdata[, setNames(mirStat, miRNA)]
    result[["MiR_Score"]] <- sort(data, decreasing = TRUE)
  }
  result
}

#' @noRd
get_data_edge_based <- function(methods, mir_DEdata, gene_DEdata, background_MGI, mir_mat, gene_mat, scoreFun) {
  DE_method <- c("Edge_ORA", "Edge_Network")
  score_method <- c("Edge_Score1D", "Edge_Score2D", "Edge_Topology")
  result <- list()

  dict <- merge(background_MGI, mir_DEdata[, .(miRNA, mirLogFC, mirPadj)], by = "miRNA", all = FALSE)
  dict <- merge(dict, gene_DEdata[, .(gene, geneLogFC, genePadj)], by = "gene", all = FALSE)
  cor <- data.table::as.data.table(calc_MGI_corr(dict = dict, mir_mat = mir_mat, gene_mat = gene_mat))
  dict <- merge(dict, cor, by = c("miRNA", "gene"), all = FALSE)
  dict <- na.omit(dict)
  data.table::setcolorder(dict, c("miRNA", "gene", "cor", "corPadj", "mirLogFC", "mirPadj", "geneLogFC", "genePadj"))
  result[["data"]] <- dict

  if (any(DE_method %in% methods)){
    sel <- dict[
      cor < 0 & corPadj < 0.05 &
        abs(mirLogFC) > 1 & mirPadj < 0.05 &
        abs(geneLogFC) > 1 & genePadj < 0.05
    ]
    data <- if (nrow(sel) > 0) paste0(sel$miRNA, ":", sel$gene) else character(0)
    for (sub_method in intersect(DE_method, methods)){
      result[[sub_method]] <- data
    }
  }
  if (any(score_method %in% methods)){
    dict[, ratio := (geneLogFC / mirLogFC) * ((1-genePadj)/(1+mirPadj))]
    if (scoreFun == "rank") {
      dict[, normratio := ranknorm(ratio)]
    } else if (scoreFun == "sigmoid") {
      dict[, normratio := sigmoid(ratio)]
    } else if (scoreFun == "normalize") {
      dict[, normratio := max_min(ratio)]
    } else {
      dict[, normratio := NA_real_]
    }
    dict[, `:=`( # `:=`() form modifies multiple columns; infix `:=` modifies a single column
      strength = cor * normratio,
      MGI = paste0(miRNA, ":", gene)
    )]
    data.table::setorder(dict, strength)
    dict <- unique(dict, by = "MGI")

    if ("Edge_Score1D" %in% methods) {
      result[["Edge_Score1D"]] <- setNames(as.numeric(dict$strength), dict$MGI)
    }
    if ("Edge_Topology" %in% methods) {
      result[["Edge_Topology"]] <- as.data.frame(dict[, .(MGI, miRNA, gene, strength)])
    }
    if ("Edge_Score2D" %in% methods) {
      result[["Edge_Score2D"]] <- as.data.frame(dict[, .(MGI, miRNA, gene, cor, normratio)])
    }
  }
  result
}


get_pathway_data <- function(methods, pathway, background_MGI){
  colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
  sum_TG <- sum(grepl("TG_|^Edge_Network$", methods)) # Edge-Network needs pathway gene set to map GGI.
  sum_MiR <- sum(grepl("MiR_", methods))
  sum_Edge <- sum(grepl("Edge_", methods))

  result <- list()

  if (sum_TG > 0){
    result[["TG"]] <- pathway
  }
  if (sum_MiR + sum_Edge > 0){
    colnames(background_MGI)[seq_len(2)] <- c("miRNA", "gene")
    background_MGI <- data.table::as.data.table(background_MGI)
    pathway <- data.table::as.data.table(pathway)

    pathway_genes <- pathway[, .(genes = list(gene)), by = pathway]
    result_list <- lapply(seq_len(nrow(pathway_genes)), function(i) {
      pw <- pathway_genes$pathway[i]
      genes <- pathway_genes$genes[[i]]
      background_MGI[gene %in% genes, .(pathway = pw, miRNA, gene)]
    })
    pathway_MGI <- unique(rbindlist(result_list))
    if (sum_MiR > 0){
      result[["MiR"]] <- unique(pathway_MGI[, .(pathway, miRNA)])
    }
    if (sum_Edge > 0){
      result[["Edge"]] <- pathway_MGI
    }
  }
  return(result)
}


# only useful for Edge_Topology, Edge_Network
get_pathway_GGI_data <- function(methods, pathway, background_GGI){
  target_method <- c("Edge_Network", "Edge_Topology")
  if (length(intersect(target_method, methods)) == 0){
    warning("The methods you have selected don't need gene-gene interactions data. Skip the pathway_GGI process.\n")
    return(NULL)
  } else {
    colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
    colnames(background_GGI)[seq_len(2)] <- c("from", "to")
    background_GGI <- data.table::as.data.table(background_GGI)
    pathway <- data.table::as.data.table(pathway)

    pathway_genes <- pathway[, .(genes = list(gene)), by = pathway]
    result_list <- vector("list", nrow(pathway_genes))

    result_list <- lapply(seq_len(nrow(pathway_genes)), function(i) {
      pw <- pathway_genes$pathway[i]
      genes <- pathway_genes$genes[[i]]

      background_GGI[from %in% genes & to %in% genes,
                     .(pathway = pw, from, to)]  # , type = "GGI"
    })

    pathway_GGI <- unique(rbindlist(result_list))
    return(pathway_GGI)
  }
}

# DEMGI: result[["data"]][["Edge_Network"]]
get_GGI_data <- function(methods, DEMGI, pathway, background_GGI){ # gene_mat, gene_DEdata,
  if (! "Edge_Network" %in% methods){
    warning("The methods you have selected don't need global GGI data. Skip the get_GGI_data process.\n")
    return(NULL)
  } else{
    DEMGI <- data.frame(
      miRNA = sub(":.*", "", DEMGI),
      gene = sub(".*:", "", DEMGI),
      stringsAsFactors = FALSE)

    colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
    colnames(background_GGI)[seq_len(2)] <- c("from", "to")
    background_GGI <- data.table::as.data.table(background_GGI)
    pathway <- data.table::as.data.table(pathway)
    DEMGI <- data.table::as.data.table(DEMGI)

    pathway_genes <- unique(pathway$gene)
    DEMGI_genes <- unique(DEMGI$gene)
    global_genes <- unique(c(pathway_genes, DEMGI_genes))
    global_GGI <- background_GGI[
      from %in% global_genes & to %in% global_genes,
      .(from = from, to = to)
    ]
    global_GGI <- as.data.frame(global_GGI)
    return(global_GGI)
  }
}

get_MGI_data <- function(methods, DEMGI, pathway, pathway_MGI, background_MGI){
  if (! "Edge_Network" %in% methods){
    warning("The methods you have selected don't need global GGI data. Skip the get_GGI_data process.\n")
    return(NULL)
  } else{
    DEMGI <- data.frame(
      miRNA = sub(":.*", "", DEMGI),
      gene = sub(".*:", "", DEMGI),
      stringsAsFactors = FALSE)

    colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
    colnames(background_MGI)[seq_len(2)] <- c("miRNA", "gene")
    colnames(pathway_MGI)[seq_len(3)] <- c("pathway", "miRNA", "gene")
    background_MGI <- data.table::as.data.table(background_MGI)
    pathway <- data.table::as.data.table(pathway)
    pathway_MGI <- data.table::as.data.table(pathway_MGI)
    DEMGI <- data.table::as.data.table(DEMGI)

    pathway_genes <- unique(pathway$gene)
    pathway_mirs <- unique(pathway_MGI$miRNA)
    DEMGI_genes <- unique(DEMGI$gene)
    DEMGI_mirs <- unique(DEMGI$miRNA)
    global_genes <- unique(c(pathway_genes, DEMGI_genes))
    global_mirs <- unique(c(pathway_mirs, DEMGI_mirs))

    global_MGI <- background_MGI[
      gene %in% global_genes & miRNA %in% global_mirs,
      .(miRNA = miRNA, gene = gene)
    ]
    global_MGI <- as.data.frame(global_MGI)
    return(global_MGI)
  }
}


# useful in Edge_Topology and Edge_Network
# data is data$data get in the get_data() function
get_dict_MGI <- function(methods, mir_mat, gene_mat, data, pathway_MGI, scoreFun){
  if (! "Edge_Topology" %in% methods){
    warning("The methods you have selected don't need dict_MGI data. Skip the get_dict_MGI process.\n")
    return(NULL)
  } else{
    # ensure that mir_mat and gene_mat have their names as rownames
    if (!is.numeric(mir_mat[[1]]) && !is.integer(mir_mat[[1]])) {
      rownames(mir_mat) <- mir_mat[[1]]
      mir_mat <- mir_mat[, -1, drop = FALSE]
    }
    if (!is.numeric(gene_mat[[1]]) && !is.integer(gene_mat[[1]])) {
      rownames(gene_mat) <- gene_mat[[1]]
      gene_mat <- gene_mat[, -1, drop = FALSE]
    }

    # check shared samples
    common_samples <- intersect(colnames(mir_mat), colnames(gene_mat))
    if (length(common_samples) == 0) {
      stop("No matched samples between miRNA expression and gene expression profiles,
         please ensure mir_mat and gene_mat have shared colnames as sample names!")
    }
    if (length(common_samples) < min(ncol(mir_mat), ncol(gene_mat)) * 0.5) {
      warning("The number of matched samples between mir_mat and gene_mat is:", length(common_samples),
              ", which is below 50% of the total number of samples.
            We suggest to double check the colnames of miRNA and gene expression profiles have shared colnames as sample names!")
    }

    colnames(pathway_MGI)[seq_len(3)] <- c("pathway", "miRNA", "gene")
    pathway_MGI <- data.table::as.data.table(pathway_MGI)

    data <- data.table::as.data.table(data)

    # build mat that have both miRNA and gene with only shared samples.
    mir_sub <- mir_mat[, common_samples, drop = FALSE]
    gene_sub <- gene_mat[, common_samples, drop = FALSE]

    # build dict dataframe: all possible combinations between miRNAs and genes,
    # who are existed in pathway and meanwhile have valid expresison profiles.
    unique_miRNAs <- intersect(unique(pathway_MGI$miRNA), unique(data$miRNA))
    unique_genes <- intersect(unique(pathway_MGI$gene), unique(data$gene))
    dict <- data.table::CJ(miRNA = unique_miRNAs, gene = unique_genes)

    # merge with gene and miRNA differential expressed data
    mir_info <- unique(data[, .(miRNA, mirLogFC, mirPadj)])
    gene_info <- unique(data[, .(gene, geneLogFC, genePadj)])

    dict <- merge(dict, mir_info, by = "miRNA", all.x = TRUE)
    dict <- merge(dict, gene_info, by = "gene", all.x = TRUE)

    results_mat <- calc_MGI_corr(dict = dict, mir_mat = mir_sub, gene_mat = gene_sub)
    results_mat <- data.table::as.data.table(results_mat)

    dict <- merge(dict, results_mat[, .(miRNA, gene, cor, corPadj)], by = c("miRNA", "gene"), all.x = TRUE)
    dict <- dict[complete.cases(dict)]

    # dict[, corPadj := p.adjust(pvalue, method = "BH")]

    dict[, ratio := (geneLogFC/mirLogFC) * ((1-genePadj)/(1+mirPadj))]
    if (scoreFun == "rank") {
      dict[, normratio := ranknorm(ratio)]
    } else if (scoreFun == "sigmoid") {
      dict[, normratio := sigmoid(ratio)]
    } else if (scoreFun == "normalize") {
      dict[, normratio := max_min(ratio)]
    } else {
      dict[, normratio := NA_real_]
    }

    dict[, strength := cor * normratio]
    dict[, MGI := paste0(miRNA, ":", gene)]
    data.table::setorder(dict, -strength)
    dict <- unique(dict)
    dict <- dict[, .(MGI, miRNA, gene, cor, corPadj, mirLogFC, mirPadj, geneLogFC, genePadj,
                     ratio, normratio, strength)]
    dict <- as.data.frame(dict)
    return(dict)
  }
}



#' @param dict data.frame/data.table with at least `miRNA` and `gene` columns naming the pairs to correlate.
#' @param mir_mat miRNA-by-sample expression matrix (or data.frame with miRNA names in the first column).
#' @param gene_mat gene-by-sample expression matrix (or data.frame with gene names in the first column).
#' @return A data.frame with `miRNA`, `gene`, `cor` (Spearman rho), and `corPadj` (BH-adjusted p-value).
#' @details Computes Spearman correlation for every (miRNA, gene) pair in `dict` by rank-transforming
#' each expression matrix once, then evaluating all pairs with matched-index matrix algebra
#' (no per-pair loop). Pairs are processed in chunks of 2e5 to bound peak memory on very large
#' `dict` inputs (e.g. whole-pathway background correlation for `Edge_Topology`).
#' @noRd
calc_MGI_corr <- function(dict, mir_mat, gene_mat) {
  # ensure that mir_mat and gene_mat have their names as rownames
  if (!is.numeric(mir_mat[[1]]) && !is.integer(mir_mat[[1]])) {
    rownames(mir_mat) <- mir_mat[[1]]
    mir_mat <- mir_mat[, -1, drop = FALSE]
  }
  if (!is.numeric(gene_mat[[1]]) && !is.integer(gene_mat[[1]])) {
    rownames(gene_mat) <- gene_mat[[1]]
    gene_mat <- gene_mat[, -1, drop = FALSE]
  }

  # check shared samples
  common_samples <- intersect(colnames(mir_mat), colnames(gene_mat))
  if (length(common_samples) == 0) {
    stop("No matched samples between miRNA expression and gene expression profiles,
         please ensure mir_mat and gene_mat have shared colnames as sample names!")
  }
  if (length(common_samples) < min(ncol(mir_mat), ncol(gene_mat)) * 0.5) {
    warning("The number of matched samples between mir_mat and gene_mat is:", length(common_samples),
            ", which is below 50% of the total number of samples.
            We suggest to double check the colnames of miRNA and gene expression profiles have shared colnames as sample names!")
  }

  # build mat that have both miRNA and gene with only shared samples.
  mir_sub <- mir_mat[, common_samples, drop = FALSE]
  gene_sub <- gene_mat[, common_samples, drop = FALSE]

  # rank-transform each row once; Spearman correlation is Pearson correlation on ranks.
  rank_mir  <- t(apply(as.matrix(mir_sub), 1, rank))
  rank_gene <- t(apply(as.matrix(gene_sub), 1, rank))

  idx_mir_all  <- match(dict$miRNA, rownames(mir_sub))
  idx_gene_all <- match(dict$gene, rownames(gene_sub))
  n_samples <- ncol(rank_mir)
  n_pairs   <- nrow(dict)

  chunk_size <- 2e5
  rho_all  <- numeric(n_pairs)
  pval_all <- numeric(n_pairs)
  for (start in seq(1, n_pairs, by = chunk_size)) {
    idx <- start:min(start + chunk_size - 1, n_pairs)
    Rm <- rank_mir[idx_mir_all[idx], , drop = FALSE]
    Rg <- rank_gene[idx_gene_all[idx], , drop = FALSE]
    mean_m <- rowMeans(Rm)
    mean_g <- rowMeans(Rg)
    Rm_c <- Rm - mean_m
    Rg_c <- Rg - mean_g
    cov_mg <- rowSums(Rm_c * Rg_c)
    sd_m <- sqrt(rowSums(Rm_c^2))
    sd_g <- sqrt(rowSums(Rg_c^2))
    rho <- cov_mg / (sd_m * sd_g)
    tstat <- rho * sqrt((n_samples - 2) / (1 - rho^2))
    rho_all[idx]  <- rho
    pval_all[idx] <- 2 * stats::pt(-abs(tstat), df = n_samples - 2)
  }

  data.frame(
    miRNA = dict$miRNA,
    gene = dict$gene,
    cor = rho_all,
    corPadj = stats::p.adjust(pval_all, method = "BH")
  )
}


min_max <- function(x) {
  diff <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (diff != 0){
    scaled <- (x - min(x, na.rm = TRUE)) / diff
  } else {
    scaled <- 0
  }
}

ranknorm <- function(x){
  ranks <- rank(x)
  result <- (max(ranks) + 1 - ranks)/length(ranks)
  return(result)
}

ranknorm_pos <- function(x){
  ranks <- rank(x)
  result <- ranks/length(ranks)
  return(result)
}

sigmoid <- function(x){
  1 / (1 + exp(x))
}

sigmoid_pos <- function(x){
  1 / (1 + exp(-x))
}

max_min <- function(x){
  diff <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (diff != 0){
    scaled <- 1 - (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    # scaled * 0.98 + 0.01
  } else {
    scaled <- 0
  }

}
