# This is the code for 2.4. Comparisons between integrating gene data or not
# Test Edge_ORA and Edge_Network results for only having miRNA data
# methods <- c("Edge_ORA", "Edge_Network")
# paths <- c("TP", "TN")

#### script for parallel computing based on (cancer, pathway)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_enrichment.R <cancer>") #  <method_list> <pAdjMethod>
}
#
cancer <- args[1]

path_name = "TP"
ncore <- 16
ncores <- list("Edge_Topology" = ncore, "Edge_Network" = min(8, ncore), "Edge_2Ddist" = 1)


setwd("/scratch/project_2011179/code/miREA")
load("function.RData")
.libPaths(c("/projappl/project_2011179/rpackages_440", .libPaths()))  ## only use when on CSC
source("R/lib.R")

scoreFun = "rank"
pvalueCutoff = 1
minSize = NULL
maxSize = NULL
pAdjMethod = "BH"
penrichCutoff = 0.05
iter = 1000

base_dir <- "analysis/2.6_data_source/miR_only_result/"
if (!dir.exists(base_dir)) {dir.create(base_dir)}
pathway_path <- paste0(base_dir, path_name, "/")
if (!dir.exists(pathway_path)) {dir.create(pathway_path)}

cancer_path <- paste0(pathway_path, cancer,"/")
plot_path <- paste0(cancer_path, "plot/")
rdata_path <- paste0(cancer_path, "RData/")
result_path <- paste0(cancer_path,"result/")

if (!dir.exists(cancer_path)) {dir.create(cancer_path)}
if (!dir.exists(plot_path)) {dir.create(plot_path)}
if (!dir.exists(rdata_path)) {dir.create(rdata_path)}
if (!dir.exists(result_path)) {dir.create(result_path)}

output_file  <- paste0(cancer_path, "log.txt")
sink(output_file, split = TRUE, append = TRUE)

##### 1. Data preparation----

methods <- c("Edge_ORA", "Edge_Network")

cat("Enrichment Analysis Based on cancer:", cancer, ", pathway:", path_name,".\n")
cat("Step 1. Get all input data prepared.\n")

input_loc <- paste0("data/input_data/", path_name, "/", cancer, "_", path_name, "_input_data.RData")
load(input_loc)
DEMGI <- input_data$data$data %>% filter(abs(mirLogFC) > 1, mirPadj < 0.05) %>% pull(MGI) %>% unique()

input_data$data$Edge_ORA <- DEMGI
input_data$data$Edge_Network <- DEMGI



##### 2. Enrichment analysis ----
cat("Step 2. Enrichment analysis...\n")
cat("\n  === Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
setwd("/scratch/project_2011179/code/miREA/")
background <- list()
if ("TG" %in% names(input_data$pathway)){
  background[["TG"]] <- unique(input_data$pathway$TG$gene)
}
if ("MiR" %in% names(input_data$pathway)){
  background[["MiR"]] <- unique(input_data$background_MGI$miRNA)
}
if ("background_MGI" %in% names(input_data)){
  background[["Edge"]] <- unique(paste0(input_data$background_MGI$miRNA, ":", input_data$background_MGI$gene))
}
result <- miREA(methods = methods, input_data = input_data, background = background, minSize = NULL, maxSize = NULL,
                pvalueType = list(TG = "std", MiR = "std", Edge = "neg"), pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff,
                iter = iter, ncores = ncores)
save(result, file = paste0(rdata_path, "result.RData"))


##### 3. Summary ----
cat("Step 3. summary...\n")
cat("\n  === Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
enrich <- list()
n_gene_path <- length(unique(input_data$pathway$TG$pathway))
n_mir_path <- length(unique(input_data$pathway$MiR$pathway))
n_edge_path <- length(unique(input_data$pathway$Edge$pathway))
summary <- data.frame(pathway = character(), method = character(), n_enrich = integer(), n_pathway = integer(), time = double(), stringsAsFactors = FALSE)
for (method in names(result$result)){
  result_df <- result$result[[method]]
  time_df <- result$time
  enrich_item <- result$result[[method]] %>% filter(padj < penrichCutoff) %>% pull(pathway)
  write.csv(result_df, file = paste0(result_path, cancer, "_", method, ".csv"), row.names = FALSE)
  enrich[[method]] <- enrich_item
  if(grepl("TG_", method)){
    n_path <- n_gene_path
  } else if (grepl("MiR_", method)){
    n_path = n_mir_path
  } else if (grepl("Edge_", method)){
    n_path = n_edge_path
  } else {
    n_path = NA
  }
  summary <- rbind(summary,
                   data.frame(
                     pathway = path_name, method = method, n_enrich = length(enrich_item), n_pathway = n_path, time = time_df[time_df$method == method, "time"],
                     stringsAsFactors = FALSE))
}
summary <- summary %>% mutate(positive_rate = n_enrich / n_pathway)

write.csv(summary, file = paste0(cancer_path, cancer, "_", path_name, "_result.csv"))
save(enrich, file = paste0(rdata_path,  "enrich.RData"))

##### 4. Visualization ----
cat("Step 4. Visualization...\n")
cat("\n  === Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
summary_plot <- plot_summary(result = result, penrichCutoff = penrichCutoff, plot_path,
                             fill_col = c("TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
                                          "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_2Ddist" = "#626FB3",
                                          "Edge_Topology" = "#EAA5C2","Edge_Network" = "#FA8072"))


cat("\nMission completed! Please check your folder:\n    ", cancer_path, "\n")
cat("\n  === End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
while (sink.number() > 0) {
  sink(NULL)
}
