# if TG_method, initially use the differential expressed miRNAs (interested miRNAs), and link to genes (as interested genes, or geneList (need gene_data)).
# if MiR_method, use mir_DEdata directly
# if Edge_method, need both mir_DEdata and gene_DEdata

#' Description: Functions to generate input data for miREA, based on the original expression and differential expression results for miRNAs and genes.
#' 
#' @param method compulsory, see default.method_list for all supportive methods. You could either input one method, or a vector contains several methods you want to try.
#' @param pathway compulsory. It should be a two-column pathway, containing pathway and gene.
#' @param mir_DEdata: compulsory, should be a four-column dataframe, containing miRNA, log2FC, stat, padj
#' @param gene_DEdata: optional, except for MiR_ based methods, should be a four-column dataframe, containing gene, log2FC, stat, padj
#' @param background_MGI: optional, but if doesn't include, we would load our predefined background. It should be either a character vector in the format of [miRNA]:[gene], or a two-column dataframe, containing miRNA gene
#' @param background_GGI: optional, it should be a two-column data frame that contains source gene and target gene. Combined with *GGI_source*. If NULL, we would use GGI_source to determine the background GGIs.
#' @param GGI_source: optional, select from "Omnipath" or "Reactome".if you don't have specified background_GGI and neither GGI_source. The default Omnipath will be used.
#' @param gene_mat: optional, needed for Edge_ method, used for calculating correlations between miRNAs and genes. It should be raw gene expression matrix, where each row is a gene, and each column represents a sample.
#' @param mir_mat: optional, needed for Edge_ method, used for calculating correlations between miRNAs and genes. It should be raw miRNA expression matrix, where each row is a miRNA, and each column represents a sample.
#' @param scoreFun: optional, needed for Edge_Score, Edge_2Ddist, Edge_Topology. See default.scoreFun for all supportive methods.
#' @note gene_mat and mir_mat should have same sample name format. If the number of shared samples below 50 % of the minimum of gene samples and miRNA samples, a warning would happen to remind you recheck them.

#' @return A list that contains all the data formats needed for miREA.

get_all_input_data <- function(methods, pathway,
                               mir_DEdata, gene_DEdata = NULL,
                               background_MGI = NULL, background_GGI = NULL, GGI_source = "Omnipath",
                               gene_mat = NULL, mir_mat = NULL,
                               scoreFun = NULL){
  cat("Generate all input data...\n")
  sum_TG <- sum(grepl("TG_", methods))
  sum_MiR <- sum(grepl("MiR_", methods))
  sum_Edge <- sum(grepl("Edge_", methods))
  target_method <- c("Edge_Topology", "Edge_Network")
  score_method <- c("Edge_Score", "Edge_2Ddist", "Edge_Topology")
  colnames(pathway)[1:2] <- c("pathway", "gene")
  # check scoreFun
  if (length(intersect(score_method, methods)) > 0 && is.null(scoreFun)){
    warning("No [scoreFun] parameter has been specified. Use the default [rank] instead!\n")
    scoreFun = "rank"
  }
  if (length(intersect(score_method, methods)) > 0 && (! scoreFun %in% default.scoreFun)){
    stop("The [scoreFun] you input is not valid. Please select one from the default.scoreFun!")
  }

  result <- list()

  # check background_MGI
  if (sum_MiR + sum_Edge > 0){
    if (is.null(background_MGI)){
      warning("You didn't specify background_MGI. Read the default one instead.\n")
      load("data/raw_data/background/background_MGI.RData")
    }
    colnames(background_MGI)[1:2] <- c("miRNA", "gene")

  }

  # check background_GGI
  if (length(intersect(target_method, methods)) != 0){
    if (is.null(background_GGI)){
      if (! GGI_source %in% default.GGI_source){
        stop("You tend to extract default background_GGI for", GGI_source,". But [GGI_source] you have specified are not involved. Please use default.GGI_source() to see all supporting options!")
      }
      warning("No background_GGI is given. Use global GGI list from ", GGI_source, " instead.\n")
      load(paste0("data/raw_data/background/background_GGI_", GGI_source, ".RData"))
    }
    colnames(background_GGI)[1:2] <- c("from", "to")
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


  cat("  Start generate expression data...\n")
  data <- get_data(methods = methods, mir_DEdata = mir_DEdata, gene_DEdata = gene_DEdata,
                   background_MGI = background_MGI, gene_mat = gene_mat, mir_mat = mir_mat, scoreFun = scoreFun)
  result[["data"]] <- data

  cat("  Start generate pathway data...\n")
  pathway_data <- get_pathway_data(methods = methods, pathway = pathway, background_MGI = background_MGI)
  result[["pathway"]] <- pathway_data

  if (length(intersect(target_method, methods)) != 0){
    cat("  Start generate pathway_GGI data...\n")
    pathway_GGI <- get_pathway_GGI_data(methods = methods, pathway = pathway, background_GGI = background_GGI)
    result[["pathway_GGI"]] <- pathway_GGI

    if ("Edge_Network" %in% methods){
      cat("  Start generate GGI data for whole network...\n")
      GGI <- get_GGI_data(methods = methods, DEMGI = result[["data"]][["Edge_Network"]], pathway = pathway, background_GGI = background_GGI)
      result[["GGI"]] <- GGI
      cat("  Start generate MGI data for whole network...\n")
      MGI <- get_MGI_data(
        methods = methods, DEMGI = result[["data"]][["Edge_Network"]], pathway = pathway,
        pathway_MGI = result[["pathway"]][["Edge"]], background_MGI = background_MGI)
      result[["MGI"]] <- MGI
    }

    if ("Edge_Topology" %in% methods){
      cat("  Start generate dict_MGI data for permutation...\n")
      dict_MGI <- get_dict_MGI(methods = methods, mir_mat = mir_mat, gene_mat = gene_mat, data = data$data,
                               pathway_MGI = pathway_data$Edge, scoreFun = scoreFun)
      result[["dict_MGI"]] <- dict_MGI
    }

  }

  cat("All the input data is ready!\n")
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
#                 Edge_ORA, Edge_Score, Edge_TopoDE, Edge_network
# background_GGI
# pathway
# scoreFun: Edge_Score

# gene_DEdata: except for MiR_ORA and MiR_Score, should be four-column dataframe, containing gene, log2FC, stat, padj
# mir_DEdata: compulsory, should be four-column dataframe, containing miRNA, log2FC, stat, padj
# gene_mat: only needed for Edge_ method
# mir_mat: only needed for Edge_ method
# scoreFun: only needed for Edge_Score, Edge_2Ddist, Edge_TopoScore, Edge_TopoScore_weight, choose one from rank, sigmoid, normalize
# background_MGI: only needed for Edge_ method
# will return a list, contains the input of each method
get_data <- function(methods, mir_DEdata, gene_DEdata = NULL,
                     background_MGI, gene_mat = NULL, mir_mat = NULL, scoreFun){
  invalid_method <- base::setdiff(methods, default.method_list)
  if (length(invalid_method) > 0){
    stop("The method(s):", invalid_method, "you are input is(are) not one of the supportive methods. Please revise the method(s)!")
  }
  sum_TG <- sum(grepl("TG_", methods))
  sum_MiR <- sum(grepl("MiR_", methods))
  sum_Edge <- sum(grepl("Edge_", methods))
  DE_method <- c("Edge_ORA", "Edge_Network")
  score_method <- c("Edge_Score", "Edge_2Ddist", "Edge_Topology")
  if ((sum_TG != 0 || sum_Edge != 0) && is.null(gene_DEdata)){
    stop("Please make sure you have input gene_DEdata!")
  }
  if (sum_Edge != 0 && (is.null(gene_mat) || is.null(mir_mat))){
    stop("Please make sure you have input gene_mat and mir_mat")
  }
  if (any(score_method %in% methods) && is.null(scoreFun)){
    warning("You didn't specify scoreFun for transform expression score, set 'rank' automatically.\n")
    scoreFun = "rank"
  }
  if (sum_Edge != 0 && is.null(background_MGI)){
    warning("You didn't specify background_MGI. Read the default one instead.\n")
    load("data/raw_data/background/background_mirv22_geneHGNC.RData")
  }

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


  if (!is.null(background_MGI)){
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
    background_MGI <- data.table::as.data.table(background_MGI)
  }

  result <- list()

  if (sum(methods %in% c("TG_ORA", "TG_Score", "MiR_ORA")) > 0){
    interested_miRNAs <- mir_DEdata[abs(mirLogFC) > 1 & mirPadj < 0.05, unique(miRNA)]
  }


  if ("TG_ORA" %in% methods){
    data <- background_MGI[miRNA %in% interested_miRNAs, unique(gene)]
    result[["TG_ORA"]] <- data
  }
  if ("TG_Score" %in% methods){
    genes <- background_MGI[miRNA %in% interested_miRNAs, unique(gene)]
    data <- gene_DEdata[gene %in% genes,
                               setNames(geneStat, gene)]
    data <- sort(data, decreasing = TRUE)
    result[["TG_Score"]] <- data
  }
  if ("MiR_ORA" %in% methods){
    result[["MiR_ORA"]] <- interested_miRNAs
  }
  if ("MiR_Score" %in% methods){
    data <- mir_DEdata[, setNames(mirStat, miRNA)]
    data <- sort(data, decreasing = TRUE)
    result[["MiR_Score"]] <- data

  }
  if (sum(grepl("Edge_", methods)) > 0){
    dict <- merge(background_MGI, mir_DEdata[, .(miRNA, mirLogFC, mirPadj)], by = "miRNA", all = FALSE)
    dict <- merge(dict, gene_DEdata[, .(gene, geneLogFC, genePadj)], by = "gene", all = FALSE)
    cor <- data.table::as.data.table(calc_MGI_corr(dict = dict, mir_mat = mir_mat, gene_mat = gene_mat))
    dict <- merge(dict, cor, by = c("miRNA", "gene"), all = FALSE)
    dict <- na.omit(dict)
    data.table::setcolorder(dict, c("miRNA", "gene", "cor", "corPadj", "mirLogFC", "mirPadj", "geneLogFC", "genePadj"))
    result[["data"]] <- dict

    if (any(DE_method %in% methods)){
      data <- dict[
        cor < 0 & corPadj < 0.05 &
          abs(mirLogFC) > 1 & mirPadj < 0.05 &
          abs(geneLogFC) > 1 & genePadj < 0.05,
        paste0(miRNA, ":", gene)
      ]
      # data <- paste0(data$miRNA, ":", data$gene)
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
      #dict[, normcor := (1-cor)/2]

      dict[, `:=`( #`:=`表示修改多个列，:=是只修改一个列
        strength = cor * normratio,
        MGI = paste0(miRNA, ":", gene)
      )]
      setorder(dict, strength)
      dict <- unique(dict, by = "MGI")

      if ("Edge_Score" %in% methods) {
        result[["Edge_Score"]] <- setNames(as.numeric(dict$strength), dict$MGI)
      }

      if ("Edge_Topology" %in% methods) {
        result[["Edge_Topology"]] <- as.data.frame(dict[, .(MGI, miRNA, gene, strength)])
      }

      # if ("Edge_manova" %in% methods) {
      #   result[["Edge_manova"]] <- as.data.frame(dict[, .(MGI, miRNA, gene, normcor, normratio, cor, ratio)])
      # }
      if ("Edge_2Ddist" %in% methods) {
        result[["Edge_2Ddist"]] <- as.data.frame(dict[, .(MGI, miRNA, gene, cor, normratio)])
      }
    }
  }
  return(result)
}


get_pathway_data <- function(methods, pathway, background_MGI){
  colnames(pathway)[1:2] <- c("pathway", "gene")
  sum_TG <- sum(grepl("TG_", methods))
  sum_MiR <- sum(grepl("MiR_", methods))
  sum_Edge <- sum(grepl("Edge_", methods))

  result <- list()

  if (sum_TG > 0){
    result[["TG"]] <- pathway
  }
  if (sum_MiR + sum_Edge > 0){
    # if (is.null(background_MGI)){
    #   cat("No background_MGI is given. Use global MGI list instead.\n")
    #   load("data/raw_data/background/background_MGI.RData")
    # }
    colnames(background_MGI)[1:2] <- c("miRNA", "gene")
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


# only useful for Edge_TopoDE, Edge_TopoScore, Edge_network, Edge_TopoScore_weight
get_pathway_GGI_data <- function(methods, pathway, background_GGI){
  target_method <- c("Edge_Network", "Edge_Topology")
  if (length(intersect(target_method, methods)) == 0){
    warning("The methods you have selected don't need gene-gene interactions data. Skip the pathway_GGI process.\n")
    return(NULL)
  } else {
    # if (is.null(background_GGI)){
    #   cat("No background_GGI is given. Use global GGI list instead.\n")
    #   load("data/raw_data/background/background_GGI_HGNC.RData") # variable name is background_GGI
    # }
    colnames(pathway)[1:2] <- c("pathway", "gene")
    colnames(background_GGI)[1:2] <- c("from", "to")
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

# DEMGI: result[["data"]][["Edge_network"]]
get_GGI_data <- function(methods, DEMGI, pathway, background_GGI){ # gene_mat, gene_DEdata,
  if (! "Edge_Network" %in% methods){
    warning("The methods you have selected don't need global GGI data. Skip the get_GGI_data process.\n")
    return(NULL)
  } else{
    # if (is.null(background_GGI)){
    #   cat("No background_GGI is given. Use global GGI list instead.\n")
    #   load("data/raw_data/background/background_GGI.RData") # variable name is background_GGI
    # }
    DEMGI <- data.frame(
      miRNA = sub(":.*", "", DEMGI),
      gene = sub(".*:", "", DEMGI),
      stringsAsFactors = FALSE)

    colnames(pathway)[1:2] <- c("pathway", "gene")
    colnames(background_GGI)[1:2] <- c("from", "to")
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
    # if (is.null(background_MGI)){
    #   cat("No background_MGI is given. Use global MGI list instead.\n")
    #   load("data/raw_data/background/background_mirv22_geneHGNC.RData")
    # }
    DEMGI <- data.frame(
      miRNA = sub(":.*", "", DEMGI),
      gene = sub(".*:", "", DEMGI),
      stringsAsFactors = FALSE)

    colnames(pathway)[1:2] <- c("pathway", "gene")
    colnames(background_MGI)[1:2] <- c("miRNA", "gene")
    colnames(pathway_MGI)[1:3] <- c("pathway", "miRNA", "gene")
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


# useful in TopoDE, TopoScore, TopoScore_weight, and network
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

    colnames(pathway_MGI)[1:3] <- c("pathway", "miRNA", "gene")
    pathway_MGI <- data.table::as.data.table(pathway_MGI)

    data <- data.table::as.data.table(data)

    # build mat that have both miRNA and gene with only shared samples.
    mir_sub <- mir_mat[, common_samples, drop = FALSE]
    gene_sub <- gene_mat[, common_samples, drop = FALSE]
    # mat <- rbind(mir_sub, gene_sub)

    # build dict dataframe: all possible combinations between miRNAs and genes,
    # who are existed in pathway and meanwhile have valid expresison profiles.
    unique_miRNAs <- intersect(unique(pathway_MGI$miRNA), unique(data$miRNA))
    unique_genes <- intersect(unique(pathway_MGI$gene), unique(data$gene))
    dict <- data.table::CJ(miRNA = unique_miRNAs, gene = unique_genes)
    # dict <- expand.grid(miRNA = unique_miRNAs, gene = unique_genes, stringsAsFactors = FALSE)

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

  # 取出子数据
  pairs <- as.matrix(dict[, c("miRNA", "gene")])

  # 准备mir和gene数据（注意转置）
  mir_list <- as.list(as.data.frame(t(mir_sub)))
  gene_list <- as.list(as.data.frame(t(gene_sub)))

  # 转到Python
  py_mir_mat <- r_to_py(mir_list)
  py_gene_mat <- r_to_py(gene_list)
  py_pairs <- r_to_py(pairs)

  # 预先定义Python函数，一次性生效
  py_run_string("
import numpy as np
from scipy.stats import spearmanr

def batch_spearman(mir_mat, gene_mat, pairs):
    mirnas = sorted(set(miR for miR, _ in pairs))
    genes = sorted(set(gene for _, gene in pairs))

    mir_matrix = {miR: np.array(mir_mat[miR]) for miR in mirnas}
    gene_matrix = {gene: np.array(gene_mat[gene]) for gene in genes}

    results = []
    for miR, gene in pairs:
        mir_vec = mir_matrix[miR]
        gene_vec = gene_matrix[gene]
        cor, pval = spearmanr(mir_vec, gene_vec)
        results.append((miR, gene, cor, pval))
    return results
")

  # 调用Python计算
  py_results <- py$batch_spearman(py_mir_mat, py_gene_mat, py_pairs)

  # 结果拉回R
  result_df <- as.data.frame(do.call(rbind, py_results))
  colnames(result_df) <- c("miRNA", "gene", "cor", "pvalue")
  result_df$cor <- as.numeric(result_df$cor)
  result_df$pvalue <- as.numeric(result_df$pvalue)
  result_df <- data.frame(
    miRNA = unlist(result_df$miRNA),
    gene = unlist(result_df$gene),
    cor = unlist(result_df$cor),
    pvalue = unlist(result_df$pvalue)
  )
  result_df <- result_df %>%
    dplyr::mutate(corPadj = p.adjust(pvalue, method = "BH")) %>%
    dplyr::select(miRNA, gene, cor, corPadj)
  return(result_df)
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