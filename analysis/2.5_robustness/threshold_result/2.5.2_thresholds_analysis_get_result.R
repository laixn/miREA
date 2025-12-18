args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript run_enrichment.R <cancer> <pathway> <ncores>") #  <method_list> <pAdjMethod>
}
#
cancer <- args[1]
path_name <- args[2]
padj_threshold <- args[3]
logFC_threshold <- args[4]

if (logFC_threshold == 0.58){
  logFC_threshold = 0.584962500721156
}

# cancer_list <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "READ", "STAD", "THCA", "UCEC")
ncores = 8
# scoreFun = "rank"
pvalueCutoff = 1
minSize = NULL
maxSize = NULL
pAdjMethod = "BH"
# penrichCutoff = 0.05
iter = 1000
methods = c("Edge_ORA", "Edge_Network")

# logFC_thresholds <- c(1, log2(1.5), log2(1))
# padj_thresholds <- c(0.05, 0.1, 0.2)

.libPaths(c("/projappl/project_2011179/rpackages_440", .libPaths()))  ## only use when on CSC


setwd("/scratch/project_2011179/code/miREA")
source("R/lib.R")
load("function.RData")
base_dir <- "/scratch/project_2011179/code/miREA/analysis/2.5_robustness/threshold_result/"
cancer_path <- paste0(base_dir, cancer, "/")
rdata_path <- paste0(cancer_path, "RData/")
result_path <- paste0(cancer_path,"result/")

# if (!dir.exists(pathway_path)) {dir.create(pathway_path)}
if (!dir.exists(cancer_path)) {dir.create(cancer_path)}
# if (!dir.exists(plot_path)) {dir.create(plot_path)}
if (!dir.exists(rdata_path)) {dir.create(rdata_path)}
if (!dir.exists(result_path)) {dir.create(result_path)}

if (path_name == "Reactome"){
  pathway <- read.csv("data/raw_data/pathway/Reactome/Reactome_gene.csv", header = TRUE)
} else if (path_name == "hallmark"){
  pathway <- read.csv("data/raw_data/pathway/hallmark/hallmark_gene.csv", header = TRUE)
}  else if (path_name == "TN"){
  pathway <- read.csv("data/raw_data/pathway/TN/TN_gene.csv", header = TRUE)
} else if (path_name == "TP"){
  pathway <- read.csv(paste0("data/raw_data/pathway/TP/cancer_specific/TP_", cancer, "_gene.csv"), header = TRUE)
}

output_file  <- paste0(cancer_path, "log.txt")
sink(output_file, split = TRUE, append = TRUE)

# for (cancer in cancer_list){
input_loc <- paste0(cancer_path, cancer, "_", path_name, "_input_data.RData")
load(input_loc)
cat("Start processing", cancer, "for", path_name, ".\n")
basic_data <- input_data[[1]]
#for (i in 1:length(input_data[[2]])){
name <- paste0("padj", padj_threshold, "_logFC", logFC_threshold)
data <- input_data[[2]][[name]]
padj_threshold <- data$padj_threshold
logFC_threshold <- data$logFC_threshold
# cat("Start processing padj_threshold:", padj_threshold, "logFC_threshold:", logFC_threshold, "...\n")
input <- basic_data
input[["pathway"]][["TG"]] <- pathway
input[["data"]][["Edge_ORA"]] <- data$DEMGI
input[["data"]][["Edge_Network"]] <- data$DEMGI
rm(input_data)
gc()
cat("Enrichment analysis on padj_threshold:", padj_threshold, ",logFC_threshold:", logFC_threshold, ". In total", length(data$DEMGI), "DEMGIs", "...\n")
cat("\n  === Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
# setwd("/scratch/project_2011179/code/miREA/")
# background <- list(TG = NULL, MiR = unique(input_data$background_MGI$miRNA), Edge = input_data$background_MGI)
result <- miREA(methods = methods, input_data = input, background = NULL, minSize = NULL, maxSize = NULL,
                pvalueType = list(TG = "std", MiR = "std", Edge = "neg"), pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff,
                iter = iter, ncores = ncores)
gc()
save(result, file = paste0(rdata_path, path_name, "_padj", padj_threshold, "_logFC", logFC_threshold, "_result.RData"))
ORA_result <- result$result$Edge_ORA
ORA_result$n_DEMGI = length(data$DEMGI)
net_result <- result$result$Edge_Network
if (!is.null(net_result)){
  net_result$n_DEMGI = length(data$DEMGI)
  write.csv(net_result, file = paste0(result_path, path_name, "_padj", padj_threshold, "_logFC", logFC_threshold, "_Network.csv"), row.names = FALSE)
}

write.csv(ORA_result, file = paste0(result_path, path_name, "_padj", padj_threshold, "_logFC", logFC_threshold, "_ORA.csv"), row.names = FALSE)

cat("\n  === End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
#}
# }
