# This is an example to illustrate how miREA works
# BLCA enrichment on 10 cancer hallmarks, for all supported methods.
# stage 1. get input data
# stage 2. enrichment result
# stage 3. visualization and summary

setwd("/scratch/project_2011179/code/miREA/") # change your own directory here
.libPaths(c("/projappl/project_2011179/rpackages_440", .libPaths()))  ## Our special case, please remove this line when using.

# 0. Load miREA functions and packages ----
load("function.RData")
source("R/lib.R")

# 1. Setting all parameters ----
cancer = "BLCA"
path_name = "hallmark"
ncores = 8

# read pathway
# (if your have other pathways, please specify [pathway] (a dataframe contains pathway and gene),
# and GGI_source (support Omnipath, and Reactome)/background_GGI (your own background, a dataframe includes from and to) by yourself)

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

methods <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network")
scoreFun = "rank"
pvalueCutoff = 1
minSize = NULL
maxSize = NULL
pAdjMethod = "BH"
penrichCutoff = 0.05
iter = 1000

pathway_path <- paste0("results/", path_name, "/")
cancer_path <- paste0(pathway_path, cancer,"/")
plot_path <- paste0(cancer_path, "plot/")
rdata_path <- paste0(cancer_path, "RData/")
result_path <- paste0(cancer_path,"result/")

if (!dir.exists(pathway_path)) {dir.create(pathway_path)}
if (!dir.exists(cancer_path)) {dir.create(cancer_path)}
if (!dir.exists(plot_path)) {dir.create(plot_path)}
if (!dir.exists(rdata_path)) {dir.create(rdata_path)}
if (!dir.exists(result_path)) {dir.create(result_path)}

output_file  <- paste0(cancer_path, "log.txt")
sink(output_file, split = TRUE, append = TRUE)

# 2. Data preparation ----
if (exists("cancer") & exists("path_name")){
  cat("Enrichment Analysis Based on cancer:", cancer, ", pathway:", path_name,".\n")
}

cat("Step 1. Get all input data prepared.\n")
if (!dir.exists("data/input_data/")){
  dir.exists("data/input_data/")
}
input_loc <- paste0("data/input_data/", path_name, "/", cancer, "_", path_name, "_input_data.RData")
if (file.exists(input_loc)){
  load(input_loc)
} else {
  gene_DEdata <- read.csv(paste0("data/raw_data/cancer_data/DEG/", cancer, "_DEG.csv"))
  mir_DEdata <- read.csv(paste0("data/raw_data/cancer_data/DEmiR/", cancer, "_DEmiR.csv"))

  gene_mat <- read.csv(paste0("data/raw_data/cancer_data/paired/", cancer, "_gene_TP.csv"),
                       header = TRUE, check.names = FALSE)
  mir_mat <- read.csv(paste0("data/raw_data/cancer_data/paired/", cancer, "_mir_TP.csv"),
                      header = TRUE, check.names = FALSE)
  gene_DEdata <- gene_DEdata %>% select(gene, log2FoldChange, stat, padj)
  mir_DEdata <- mir_DEdata %>% select(miRNA, log2FoldChange, stat, padj)

  load("data/raw_data/background/background_MGI.RData")
  time_get_data <- system.time({
    input_data <- get_all_input_data(
      methods = methods, pathway = pathway, mir_DEdata = mir_DEdata, gene_DEdata = gene_DEdata,
      background_MGI = background_MGI, background_GGI = NULL, GGI_source = GGI_source,
      gene_mat = gene_mat, mir_mat = mir_mat, scoreFun = scoreFun) # check_mir_names = TRUE, mir_version = "v22", check_gene_names = TRUE
  })
  save(input_data, file = input_loc)
  cat("  === Total time cost: ", time_get_data[["elapsed"]], " s.\n")
}

# 3. Enrichment ----
cat("Step 2. Enrichment analysis...\n")
cat("\n  === Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
background <- list()
if ("TG" %in% names(input_data$pathway)){
  background[["TG"]] <- unique(input_data$pathway$TG$gene)
}
if ("MiR" %in% names(input_data$pathway)){
  background[["MiR"]] <- unique(input_data$pathway$MiR$miRNA)
}
if ("background_MGI" %in% names(input_data)){
  background[["Edge"]] <- unique(paste0(input_data$background_MGI$miRNA, ":", input_data$background_MGI$gene))
}
result <- miREA(methods = methods, input_data = input_data, background = background, minSize = NULL, maxSize = NULL,
                pvalueType = list(TG = "std", MiR = "std", Edge = "neg"), pAdjMethod = pAdjMethod, pvalueCutoff = pvalueCutoff,
                iter = iter, ncores = ncores)
save(result, file = paste0(rdata_path, "result.RData"))

# 4. Summary----
cat("Step 3. summary...\n")
cat("\n  === Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
enrich <- list()
if ("TG" %in% names(input_data$pathway)){
  n_gene_path <- length(unique(input_data$pathway$TG$pathway))
}
if ("MiR" %in% names(input_data$pathway)){
  n_mir_path <- length(unique(input_data$pathway$MiR$pathway))
}
if ("Edge" %in% names(input_data$pathway)){
  n_edge_path <- length(unique(input_data$pathway$Edge$pathway))
}

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

write.csv(summary, file = paste0(cancer_path, cancer, "_", path_name, "_result.csv"), row.names = FALSE)
save(enrich, file = paste0(rdata_path,  "enrich.RData"))

##### 5. Visualization ----
cat("Step 4. Visualization...\n")
cat("\n  === Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
summary_plot <- plot_summary(result = result, penrichCutoff = 0.05, plot_path,
                             fill_col = c("TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
                                          "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_2Ddist" = "#626FB3",
                                          "Edge_Topology" = "#EAA5C2","Edge_Network" = "#FA8072"))


# background <- unique(paste0(input_data$background_MGI$miRNA, ":", input_data$background_MGI$gene))
for (method in methods[grepl("Edge", methods)]){
  cat("    Heatmap plot for", method, ":\n")
  temp_result <- result
  temp_result$result <- result$result[method]
  if (length(enrich[[method]]) > 1){

    heatmap_sankey_plot <- plot_heatmap_sankey(method = method, result = result, input_data = input_data, n_mir_heatmap = 20, n_mir_sankey = 10, n_pathway = 10,
                                               # background = background,
                                               plot_path = plot_path,  penrichCutoff = penrichCutoff,
                                               height = NULL, width = NULL, sankey_prop = 1.5)
  }
}



cat("\nMission completed! Please check your folder:\n    ", cancer_path, "\n")
cat("\n  === End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
while (sink.number() > 0) {
  sink(NULL)
}
