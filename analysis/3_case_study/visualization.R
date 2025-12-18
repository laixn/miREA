top_miRNA <- c("hsa-miR-4665-5p", "hsa-miR-103a-2-5p", "hsa-miR-301a-3p", "hsa-miR-301b-3p", "hsa-miR-4778-3p", "hsa-miR-93-5p", "hsa-miR-17-5p",
          "hsa-miR-183-5p", "hsa-miR-20a-5p", "hsa-miR-7g-3p", "hsa-miR-130b-3p", "hsa-miR-19a-3p", "hsa-miR-93-3p", "hsa-miR-106b-5p",
          "hsa-miR-205-3p",
          "hsa-miR-4728-5p", "hsa-miR-149-3p", "hsa-miR-19b-3p", "hsa-miR-203a-3p", "hsa-miR-429")
# sum(miRs %in% cancer_list$CMC_miR)
#
# miRs_sankey <- c("hsa-miR-301b-3p", "hsa-miR-183-5p", "hsa-miR-17-5p", "hsa-miR-4665-5p", "hsa-miR-20a-5p", "hsa-miR-301a-3p",
#                  "hsa-miR-93-5p", "hsa-let-7g-3p", "hsa-miR-103a-3-5p", "hsa-miR-4778-3p")
# sum(miRs_sankey %in% cancer_list$CMC_miR)
#
# miR_ref <- cancer_list$CMC_miR
# miR_ref <- data.frame(miR = miR_ref, stringsAsFactors =  FALSE)
# oncoKB_ref <- data.frame(oncoKB = cancer_list$oncoKB_gene, stringsAsFactors = FALSE)
# cosmic_ref <- data.frame(cosmic = cancer_list$COSMIC_gene, stringsAsFactors = FALSE)
#
# genes <- c("CXCR5", "PCDHGC5", "PCDHGB3", "PCDHGA12", "PCDHGA5", "PCDHGA9", "PCDHGB7", "ACTN2", "GNAL", "PTGIS", "PPP1R1A",
#            "MAP1B", "PDE2A", "ADARB1", "PRKG1", "ITPR1", "SLIT2", "RUNX1T1", "DCHS1", "EBF1", "EPHA6", "TGFB1I1", "PCDHGC3",
#            "AKAP12", "PER1", "GADD45B", "LIFR", "PDE4D", "JCAD", "ITGA1", "RECK", "TLN1", "COLEC12", "SOCS3", "LEPR", "PCDH7",
#            "DIRAS1", "PARVA", "MEF2C", "PTPRD", "MRAS", "TTBK2", "RPS6KA5", "IL33", "DCUN1D3","SCN1B", "NFATC3", "C1S",
#            "F13A1", "IMPA2", "ARHGAP28", "MYO9A", "CNTFR", "LDLRAD4", "CSNK1A1", "NPR3", "TUBB6", "PPARA", "GPC5", "CNTNAP3",
#            "BCL9L", "SIGLEC11", "AFF4", "SMARCA2", "DUSP22", "OSMR", "S1PR2", "SHB", "GNE", "RSIP1", "RNH1", "ITGAL", "PLCXD3",
#            "TNFRSF10C", "CBL", "MOB3A", "CD3D", "FPR1", "METRN", "RICTOR", "BMP6", "DAB2", "PAX5", "ATP6V0E1", "HLA-DMB",
#            "APLNR", "DEK", "HPSE", "CDC25B", "PDCD1LG2")
# sum(genes %in% oncoKB_ref$oncoKB)
# sum(genes %in% cosmic_ref$cosmic)

setwd("/scratch/project_2011179/code/miREA/") # change your own directory here
library(dplyr)
library(ggalluvial)
library(ComplexHeatmap)
library(circlize) # colorRamp2()
library(ggnewscale) # new_scale_fill()

result_dir <- "analysis/3_case_study/"
plot_path <- "analysis/3_case_study/plot/"
cancer = "BLCA"
path_name = "hallmark"
penrichCutoff = 0.05
methods <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network")
# load result
enrich_dir <- paste0("results/", path_name, "/", cancer, "/")
load(paste0(enrich_dir, "RData/result.RData"))
invalid_methods <- setdiff(names(result$result), methods)
result$result[[invalid_methods]] <- NULL

# load input data
input_dir <- paste0("data/input_data/", path_name, "/")
load(paste0(input_dir, cancer, "_", path_name, "_input_data.RData"))

# load gene and miRNA annotation
load("data/cancer_list/cancer_list.RData")

# load plot function
source("R/plot_heatmap_sankey.R")
source("R/plot_summary.R")


gene_annot <- list("COSMIC (v103, 18 Nov 2025)" = cancer_list$COSMIC_gene,
                   "OncoKB (v24, Nov 2025)" = cancer_list$oncoKB_gene)
mir_annot <- list("CMC (Suszynska,M. et al. 2024)" = cancer_list$CMC_miR)

ht_sankey_plot <- plot_heatmap_sankey(method = "Edge_Network", result = result, input_data = input_data, n_mir_heatmap = 20, n_mir_sankey = 10,  n_pathway = 10, plot_path = plot_path, penrichCutoff = penrichCutoff,
                                height = 17, width = 25, sankey_prop = 1.5, gene_annot = gene_annot, mir_annot = mir_annot, annot_color = "palegreen3") # #A6DDBA"


summary_plot <- plot_summary(result = result, penrichCutoff = penrichCutoff, plot_path = plot_path,
                             fill_col = c("TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
                                          "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_2Ddist" = "#626FB3",
                                          "Edge_Topology" = "#EAA5C2","Edge_Network" = "#FA8072"))
