setwd("/scratch/project_2011179/code/miREA/")
.libPaths(c("/projappl/project_2011179/rpackages_440", .libPaths()))  ## only use when on CSC
base_dir <- "/scratch/project_2011179/code/miREA/analysis/2.5_robustness/threshold_result/"

source("R/global_variable.R")
source("R/lib.R")

if (!dir.exists(base_dir)){dir.create(base_dir)}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_enrichment.R <pathway>") #  <method_list> <pAdjMethod>
}
#
path_name <- args[1]

cancer_list <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "READ", "STAD", "THCA", "UCEC")
scoreFun = "rank"
pvalueCutoff = 1
minSize = NULL
maxSize = NULL
pAdjMethod = "BH"
penrichCutoff = 0.05
iter = 1000

logFC_thresholds <- c(1, log2(1.5), log2(1))
padj_thresholds <- c(0.05, 0.1, 0.2)


# functions for getting data ----
get_thresholds_test_data <- function(methods = c("Edge_ORA", "Edge_Network"),
                                     pathway, logFC_thresholds, padj_thresholds,
                               mir_DEdata, gene_DEdata = NULL,
                               background_MGI = NULL, background_GGI = NULL, GGI_source = "Omnipath",
                               gene_mat = NULL, mir_mat = NULL){
  cat("Generate all input data...\n")
  sum_TG <- sum(grepl("TG_", methods))
  sum_MiR <- sum(grepl("MiR_", methods))
  sum_Edge <- sum(grepl("Edge_", methods))
  target_method <- c("Edge_Topology", "Edge_Network")
  score_method <- c("Edge_Score", "Edge_2Ddist", "Edge_Topology")
  colnames(pathway)[1:2] <- c("pathway", "gene")

  result <- list()

  # check background_MGI
  if (is.null(background_MGI)){
    warning("You didn't specify background_MGI. Read the default one instead.\n")
    load("data/raw_data/background/background_MGI.RData")
  }
  colnames(background_MGI)[1:2] <- c("miRNA", "gene")


  # check background_GGI
  if (is.null(background_GGI)){
    if (! GGI_source %in% default.GGI_source){
      stop("You tend to extract default background_GGI for", GGI_source,". But [GGI_source] you have specified are not involved. Please use default.GGI_source() to see all supporting options!")
    }
    warning("No background_GGI is given. Use global GGI list from ", GGI_source, " instead.\n")
    load(paste0("data/raw_data/background/background_GGI_", GGI_source, ".RData"))
  }
  colnames(background_GGI)[1:2] <- c("from", "to")

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

  # raw_data <- list(
  #   methods = methods,
  #   mir_DEdata = mir_DEdata, gene_DEdata = gene_DEdata,
  #   gene_mat = gene_mat, mir_mat = mir_mat,
  #   scoreFun = scoreFun
  #   # check_mir_names = check_mir_names, mir_version = mir_version,
  #   # check_gene_names = check_gene_names
  # )
  # result[["raw_data"]] <- raw_data
  result[["background_MGI"]] <- background_MGI
  result[["background_GGI"]] <- background_GGI


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

  }

  cat("  Start generate expression data...\n")
  # data <- get_data(methods = methods, mir_DEdata = mir_DEdata, gene_DEdata = gene_DEdata,
  #                  background_MGI = background_MGI, gene_mat = gene_mat, mir_mat = mir_mat, scoreFun = scoreFun)
  dict <- merge(background_MGI, mir_DEdata[, .(miRNA, mirLogFC, mirPadj)], by = "miRNA", all = FALSE)
  dict <- merge(dict, gene_DEdata[, .(gene, geneLogFC, genePadj)], by = "gene", all = FALSE)
  cor <- data.table::as.data.table(calc_MGI_corr(dict = dict, mir_mat = mir_mat, gene_mat = gene_mat))
  dict <- merge(dict, cor, by = c("miRNA", "gene"), all = FALSE)
  dict <- na.omit(dict)
  data.table::setcolorder(dict, c("miRNA", "gene", "cor", "corPadj", "mirLogFC", "mirPadj", "geneLogFC", "genePadj"))
  # result[["data"]][["data"]] <- dict
  dict <- as.data.table(dict)

  test_data_list <- list()
  for (padj_threshold in padj_thresholds){
    for (logFC_threshold in logFC_thresholds){
      data <- dict[
        cor < 0 & corPadj < padj_threshold &
          abs(mirLogFC) > logFC_threshold & mirPadj < padj_threshold &
          abs(geneLogFC) > logFC_threshold & genePadj < padj_threshold,
        paste0(miRNA, ":", gene)
      ]
      name <- paste0("padj", padj_threshold, "_logFC", logFC_threshold)
      test_data_list[[name]] <- list(padj_threshold = padj_threshold, logFC_threshold = logFC_threshold, DEMGI = data)
    }
  }

  # result[["data"]][["Edge_Network"]] <- data
  # result[["data"]][["Edge_ORA"]] <- data
  # data <- paste0(data$miRNA, ":", data$gene)


  output <- list(result, test_data_list)
  cat("All the input data is ready!\n")
  return(output)
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
  scaled <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  # scaled * 0.98 + 0.01
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
  scaled <- 1 - (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  scaled * 0.98 + 0.01
}





# get input data ----
if (path_name == "Reactome"){
  pathway <- read.csv("data/raw_data/pathway/Reactome/Reactome_gene.csv", header = TRUE)
  GGI_source = "Reactome"
} else if (path_name == "hallmark"){
  pathway <- read.csv("data/raw_data/pathway/hallmark/hallmark_gene.csv", header = TRUE)
  GGI_source = "Omnipath"
}  else if (path_name == "TN"){
  pathway <- read.csv("data/raw_data/pathway/TN/TN_gene.csv", header = TRUE)
  GGI_source = "Omnipath"
} else if (path_name == "TP"){
  pathway <- read.csv(paste0("data/raw_data/pathway/TP/cancer_specific/TP_", cancer, "_gene.csv"), header = TRUE)
  GGI_source <- "Omnipath"
}


methods <- c("Edge_ORA", "Edge_Network")
for (cancer in cancer_list){
  cat("Start processing", cancer, "...\n")
  cancer_path <- paste0(base_dir, cancer, "/")
  if (!dir.exists(cancer_path)) {dir.create(cancer_path)}

  gene_DEdata <- read.csv(paste0("/scratch/project_2011179/processed_data/DE/gene/", cancer, "_DEG.csv"))
  mir_DEdata <- read.csv(paste0("/scratch/project_2011179/processed_data/DE/miRNA/", cancer, "_DEmiR.csv"))

  gene_mat <- read.csv(paste0("/scratch/project_2011179/raw_data/TCGA/paired/", cancer, "_gene_TP.csv"),
                       header = TRUE, check.names = FALSE)
  mir_mat <- read.csv(paste0("/scratch/project_2011179/raw_data/TCGA/paired/", cancer, "_mir_TP.csv"),
                      header = TRUE, check.names = FALSE)
  gene_DEdata <- gene_DEdata %>% select(gene, log2FoldChange, stat, padj)
  mir_DEdata <- mir_DEdata %>% select(miRNA, log2FoldChange, stat, padj)

  input_data <- get_thresholds_test_data(methods = methods, pathway = pathway, mir_DEdata = mir_DEdata, gene_DEdata = gene_DEdata,
                                         background_MGI = NULL, background_GGI = NULL, GGI_source = GGI_source,
                                         gene_mat = gene_mat, mir_mat = mir_mat, padj_thresholds = padj_thresholds, logFC_thresholds = logFC_thresholds)
  input_loc <- paste0(cancer_path, cancer, "_", path_name, "_input_data.RData")
  save(input_data, file = input_loc)
}


