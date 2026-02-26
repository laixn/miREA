# This is the script to evaluate abilities to identify cancer-related miRNAs and genes
# generates Figure 2D: summary.pdf

## Start directly from 2. draw if previous data has been already prepared.
# CMC: 139/165 cancer-related miRNAs from cancer miRNAs census, https://academic.oup.com/nar/article/52/4/1628/7585677#441564665
# oncoKB: 1200/1209 genes, last updated 11/24/2025, https://www.oncokb.org/cancer-genes
# COSMIC: 586 genes, cancer gene census only in Tier 1, https://cancer.sanger.ac.uk/cosmic/census?tier=1#cl_overview
# enriched pathways in cancer hallmarks for different cancers
library(dplyr)
library(biomaRt)
library(readr) # cols()
library(tidyr) # separate_rows()

.libPaths(c("/projappl/project_2011179/rpackages_440", .libPaths()))  ## Our special case, please remove this line when using.
library(miRBaseConverter)

library(ggplot2)

library(patchwork)
setwd("/scratch/project_2011179/code/miREA/") # change your own directory here

result_dir <- "analysis/2.4_cancer_identify_ability/"
if (!dir.exists(result_dir)){
  dir.create(result_dir)
}

fill_col = c("TG_Score" = "#07AEE3",
             "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_Network" = "#FA8072")

# 0. transform gene and miRNA names ----
cosmic_gene <- read.csv("data/cancer_list/raw/cancerGeneList_COSMIC.csv")
oncoKB_gene <- readr::read_tsv("data/cancer_list/raw/cancerGeneList_OncoKB.tsv", col_types = cols())
miR <- read.csv("data/cancer_list/raw/cancerMiRList.csv", header = TRUE, check.names = FALSE)

# transform gene names
ensembl <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)
oncoKB_gene_map <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters    = "entrezgene_id",
  values     = unique(oncoKB_gene$`Entrez Gene ID`),
  mart       = ensembl
)
oncoKB_gene <- oncoKB_gene %>%
  left_join(oncoKB_gene_map, by = c("Entrez Gene ID" = "entrezgene_id")) %>%
  dplyr::select(hgnc_symbol, dplyr::everything())
cat("Number of consistent gene names in oncoKB:", sum(oncoKB_gene$`Hugo Symbol` == oncoKB_gene$hgnc_symbol, na.rm = TRUE), "/", nrow(oncoKB_gene), "\n")

cosmic_gene_map <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters    = "entrezgene_id",
  values     = unique(cosmic_gene$Entrez.GeneId),
  mart       = ensembl
)
cosmic_gene <- cosmic_gene %>%
  left_join(cosmic_gene_map, by = c("Entrez.GeneId" = "entrezgene_id")) %>%
  dplyr::select(hgnc_symbol, dplyr::everything())
cat("Number of consistent gene names in COSMIC:", sum(cosmic_gene$`Gene.Symbol` == cosmic_gene$hgnc_symbol, na.rm = TRUE), "/", nrow(cosmic_gene), "\n")
# transform mature miRNAs
miR <- miR %>%
  mutate(strand = `predominantly expressed miRNA (miRNA-strand balance)`)
miR_long <- miR %>%
  separate_rows(strand, sep = "/")
miR_long <- miR_long %>%
  mutate(
    mature_miRNA = paste0("hsa-", `miRNA ID (miRBase)`, "-", strand)
  ) %>%
  dplyr::select(mature_miRNA, dplyr::everything())
cat("Total mature miRNA:", length(unique(miR_long$mature_miRNA)), "\n")

precursor_list <- unique(miR_long$`miRNA precursor/locus ID (miRBase)`)

mat_map <- miRNA_PrecursorToMature(
  precursor_list,
  "v22"
)
# mat_map_df <- data.frame(ID = precursor_list, precursor = mat_map$OriginalName, mature1 = mat_map$Mature1, mature2 = mat_map$Mature2)
mat_map_df <- na.omit(rbind(mat_map %>% dplyr::select(OriginalName, mature = Mature1), mat_map %>% dplyr::select(OriginalName, mature = Mature2)))

miR_long_df <- miR_long %>%
  left_join(mat_map_df, by = c("miRNA precursor/locus ID (miRBase)" = "OriginalName")) %>%
  dplyr::select(mature, strand, dplyr::everything()) %>%
  dplyr::filter(
    (strand == "5p" & grepl("-5p$", mature)) |
      (strand == "3p" & grepl("-3p$", mature))
  ) %>% distinct()

# mir_list <- unique(miR_long$mature_miRNA)
# acc <- miRNA_NameToAccession(mir_list, version = "v22")
# valid_miR <- acc$Target[!is.na(acc$Accession)]
# miR_long_valid <- miR_long %>%
#   filter(mature_miRNA %in% valid_miR)
cat("Valid in miRBase:", sum(miR_long_df$mature == miR_long_df$mature_miRNA), "\n")

cancer_list <- list(oncoKB_gene = unique(oncoKB_gene$hgnc_symbol),
                    COSMIC_gene = unique(cosmic_gene$hgnc_symbol),
                    CMC_miR = unique(miR_long_df$mature))
write.csv(oncoKB_gene, file = "data/cancer_list/cancerGeneList_oncoKB.csv", row.names = FALSE)
write.csv(cosmic_gene, file = "data/cancer_list/cancerGeneList_COSMIC.csv", row.names = FALSE)
write.csv(miR_long_df, file = "data/cancer_list/cancerMiRList.csv", row.names = FALSE)
save(cancer_list, file = "data/cancer_list/cancer_list.RData")

# 1. generate summary dataframe ----
paths <- c("hallmark")
cancers <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "KICH",  "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "READ", "STAD", "THCA", "UCEC")
load("data/cancer_list/cancer_list.RData")

summary <- data.frame(pathway = character(), cancer = character(), method = character(), n_enrich_pw = integer(),
                      n_core_miR = integer(), n_CMC = integer(), n_CMC_core = integer(), prop_inter_core = numeric(), prop_inter_CMC = numeric(),
                      n_core_gene = integer(), n_oncoKB = integer(), n_oncoKB_core = integer(), prop_oncoKB_inter_core = numeric(), prop_inter_oncoKB = numeric(),
                      n_COSMIC = integer(), n_COSMIC_core = integer(), prop_COSMIC_inter_core = numeric(), prop_inter_COSMIC = numeric(), stringsAsFactors = FALSE
)
for (path_name in paths){
  cat("Start processing pathway:", path_name, "\n")
  for (cancer in cancers){
    cat("--", cancer, "\n")

    load(paste0("data/input_data/hallmark/", cancer, "_", path_name, "_input_data.RData"))
    load(paste0("results/", path_name, "/", cancer, "/RData/result.RData"))

    pathway_MGI = input_data$pathway$Edge

    enrich_TG_Score <- result$result$TG_Score %>% filter(padj < 0.05)
    core_TG_Score <- unique(unlist(strsplit(enrich_TG_Score$core_enrichment, "/")))

    enrich_Edge_Score <- result$result$Edge_Score %>% filter(padj < 0.05)
    core_Edge_Score <- unique(unlist(strsplit(enrich_Edge_Score$core_enrichment, "/")))

    enrich_Edge_ORA <- result$result$Edge_ORA %>% filter(padj < 0.05)
    DEMGI = input_data$data$Edge_ORA
    core_Edge_ORA <- unique(pathway_MGI %>% mutate(MGI = paste0(miRNA, ":", gene)) %>% filter(pathway %in% enrich_Edge_ORA$pathway) %>% filter(MGI %in% DEMGI) %>% pull(MGI))

    enrich_Edge_Network <- result$result$Edge_Network %>% filter(padj < 0.05)
    core_Edge_Network <- unique(unlist(strsplit(enrich_Edge_Network$contribute_MGI, "/")))

    # check those cancer miRs in the certain enriched pathways
    # check those oncoKB genes in the certain enriched pathways
    # check those COSMIC genes in the certain enriched pathways

    ## TG_Score
    if (nrow(enrich_TG_Score) > 0){
      split_list <- strsplit(core_TG_Score, ":")
      all_elements <- unlist(split_list)
      DEmiR_TG <- input_data$data$data %>% filter(mirPadj < 0.05, abs(mirLogFC) > 1)
      core_gene_TG_Score <- unique(all_elements[!grepl("^hsa-", all_elements)])
      core_miR_TG_Score <- unique(DEmiR_TG %>% filter(gene %in% core_gene_TG_Score) %>% pull(miRNA))

      CMC_miR <- unique(pathway_MGI %>% filter(miRNA %in% cancer_list$CMC_miR, pathway %in% enrich_TG_Score$pathway) %>% pull(miRNA))
      n_core_miR = length(core_miR_TG_Score)
      n_CMC = length(CMC_miR)
      n_CMC_core = length(intersect(CMC_miR, core_miR_TG_Score))
      prop_inter_core = length(intersect(CMC_miR, core_miR_TG_Score))/length(core_miR_TG_Score)
      prop_inter_CMC = length(intersect(CMC_miR, core_miR_TG_Score))/length(CMC_miR)

      gene <- unique(pathway_MGI %>% filter(gene %in% cancer_list$oncoKB_gene, pathway %in% enrich_TG_Score$pathway) %>% pull(gene))
      n_core_gene = length(core_gene_TG_Score)
      n_oncoKB = length(gene)
      n_oncoKB_core = length(intersect(gene, core_gene_TG_Score))
      prop_oncoKB_inter_core = n_oncoKB_core / n_core_gene
      prop_inter_oncoKB = n_oncoKB_core / n_oncoKB

      gene <- unique(pathway_MGI %>% filter(gene %in% cancer_list$COSMIC_gene, pathway %in% enrich_TG_Score$pathway) %>% pull(gene))
      n_COSMIC = length(gene)
      n_COSMIC_core = length(intersect(gene, core_gene_TG_Score))
      prop_COSMIC_inter_core = n_COSMIC_core / n_core_gene
      prop_inter_COSMIC = n_COSMIC_core / n_COSMIC

      TG_Score_row <- data.frame(
        pathway = path_name, cancer = cancer, method = "TG_Score", n_enrich_pw = nrow(enrich_TG_Score),
        n_core_miR = n_core_miR, n_CMC = n_CMC, n_CMC_core = n_CMC_core,
        prop_inter_core = prop_inter_core, prop_inter_CMC = prop_inter_CMC,
        n_core_gene = n_core_gene, n_oncoKB = n_oncoKB, n_oncoKB_core = n_oncoKB_core,
        prop_oncoKB_inter_core = prop_oncoKB_inter_core, prop_inter_oncoKB = prop_inter_oncoKB,
        n_COSMIC = n_COSMIC, n_COSMIC_core = n_COSMIC_core,
        prop_COSMIC_inter_core = prop_COSMIC_inter_core, prop_inter_COSMIC = prop_inter_COSMIC, stringsAsFactors = FALSE
      )
    } else {
      TG_Score_row <- data.frame(
        pathway = path_name, cancer = cancer, method = "TG_Score", n_enrich_pw = nrow(enrich_TG_Score),
        n_core_miR = NA, n_CMC = NA, n_CMC_core = NA,
        prop_inter_core = NA, prop_inter_CMC = NA,
        n_core_gene = NA, n_oncoKB = NA, n_oncoKB_core = NA,
        prop_oncoKB_inter_core = NA, prop_inter_oncoKB = NA,
        n_COSMIC = NA, n_COSMIC_core = NA,
        prop_COSMIC_inter_core = NA, prop_inter_COSMIC = NA, stringsAsFactors = FALSE
      )
    }


    # Edge_ORA
    if (nrow(enrich_Edge_ORA) > 0){
      split_list <- strsplit(core_Edge_ORA, ":")
      all_elements <- unlist(split_list)
      core_miR_Edge_ORA  <- unique(all_elements[grepl("^hsa-", all_elements)])
      core_gene_Edge_ORA <- unique(all_elements[!grepl("^hsa-", all_elements)])

      CMC_miR <- unique(pathway_MGI %>% filter(miRNA %in% cancer_list$CMC_miR, pathway %in% enrich_Edge_ORA$pathway) %>% pull(miRNA))
      n_core_miR = length(core_miR_Edge_ORA)
      n_CMC = length(CMC_miR)
      n_CMC_core = length(intersect(CMC_miR, core_miR_Edge_ORA))
      prop_inter_core = length(intersect(CMC_miR, core_miR_Edge_ORA))/length(core_miR_Edge_ORA)
      prop_inter_CMC = length(intersect(CMC_miR, core_miR_Edge_ORA))/length(CMC_miR)

      gene <- unique(pathway_MGI %>% filter(gene %in% cancer_list$oncoKB_gene, pathway %in% enrich_Edge_ORA$pathway) %>% pull(gene))
      n_core_gene = length(core_gene_Edge_ORA)
      n_oncoKB = length(gene)
      n_oncoKB_core = length(intersect(gene, core_gene_Edge_ORA))
      prop_oncoKB_inter_core = n_oncoKB_core / n_core_gene
      prop_inter_oncoKB = n_oncoKB_core / n_oncoKB

      gene <- unique(pathway_MGI %>% filter(gene %in% cancer_list$COSMIC_gene, pathway %in% enrich_Edge_ORA$pathway) %>% pull(gene))
      n_COSMIC = length(gene)
      n_COSMIC_core = length(intersect(gene, core_gene_Edge_ORA))
      prop_COSMIC_inter_core = n_COSMIC_core / n_core_gene
      prop_inter_COSMIC = n_COSMIC_core / n_COSMIC

      Edge_ORA_row <- data.frame(
        pathway = path_name, cancer = cancer, method = "Edge_ORA", n_enrich_pw = nrow(enrich_Edge_ORA),
        n_core_miR = n_core_miR, n_CMC = n_CMC, n_CMC_core = n_CMC_core,
        prop_inter_core = prop_inter_core, prop_inter_CMC = prop_inter_CMC,
        n_core_gene = n_core_gene, n_oncoKB = n_oncoKB, n_oncoKB_core = n_oncoKB_core,
        prop_oncoKB_inter_core = prop_oncoKB_inter_core, prop_inter_oncoKB = prop_inter_oncoKB,
        n_COSMIC = n_COSMIC, n_COSMIC_core = n_COSMIC_core,
        prop_COSMIC_inter_core = prop_COSMIC_inter_core, prop_inter_COSMIC = prop_inter_COSMIC, stringsAsFactors = FALSE
      )
    } else {
      Edge_ORA_row <- data.frame(
        pathway = path_name, cancer = cancer, method = "Edge_ORA", n_enrich_pw = nrow(enrich_Edge_ORA),
        n_core_miR = NA, n_CMC = NA, n_CMC_core = NA,
        prop_inter_core = NA, prop_inter_CMC = NA,
        n_core_gene = NA, n_oncoKB = NA, n_oncoKB_core = NA,
        prop_oncoKB_inter_core = NA, prop_inter_oncoKB = NA,
        n_COSMIC = NA, n_COSMIC_core = NA,
        prop_COSMIC_inter_core = NA, prop_inter_COSMIC = NA, stringsAsFactors = FALSE
      )
    }


    # Edge_Score
    if (nrow(enrich_Edge_Score) > 0){
      split_list <- strsplit(core_Edge_Score, ":")
      all_elements <- unlist(split_list)
      core_miR_Edge_Score  <-unique(all_elements[grepl("^hsa-", all_elements)])
      core_gene_Edge_Score <- unique(all_elements[!grepl("^hsa-", all_elements)])

      CMC_miR <- unique(pathway_MGI %>% filter(miRNA %in% cancer_list$CMC_miR, pathway %in% enrich_Edge_Score$pathway) %>% pull(miRNA))
      n_core_miR = length(core_miR_Edge_Score)
      n_CMC = length(CMC_miR)
      n_CMC_core = length(intersect(CMC_miR, core_miR_Edge_Score))
      prop_inter_core = length(intersect(CMC_miR, core_miR_Edge_Score))/length(core_miR_Edge_Score)
      prop_inter_CMC = length(intersect(CMC_miR, core_miR_Edge_Score))/length(CMC_miR)

      gene <- unique(pathway_MGI %>% filter(gene %in% cancer_list$oncoKB_gene, pathway %in% enrich_Edge_Score$pathway) %>% pull(gene))
      n_core_gene = length(core_gene_Edge_Score)
      n_oncoKB = length(gene)
      n_oncoKB_core = length(intersect(gene, core_gene_Edge_Score))
      prop_oncoKB_inter_core = n_oncoKB_core / n_core_gene
      prop_inter_oncoKB = n_oncoKB_core / n_oncoKB

      gene <- unique(pathway_MGI %>% filter(gene %in% cancer_list$COSMIC_gene, pathway %in% enrich_Edge_Score$pathway) %>% pull(gene))
      n_COSMIC = length(gene)
      n_COSMIC_core = length(intersect(gene, core_gene_Edge_Score))
      prop_COSMIC_inter_core = n_COSMIC_core / n_core_gene
      prop_inter_COSMIC = n_COSMIC_core / n_COSMIC

      Edge_Score_row <- data.frame(
        pathway = path_name, cancer = cancer, method = "Edge_Score", n_enrich_pw = nrow(enrich_Edge_Score),
        n_core_miR = n_core_miR, n_CMC = n_CMC, n_CMC_core = n_CMC_core,
        prop_inter_core = prop_inter_core, prop_inter_CMC = prop_inter_CMC,
        n_core_gene = n_core_gene, n_oncoKB = n_oncoKB, n_oncoKB_core = n_oncoKB_core,
        prop_oncoKB_inter_core = prop_oncoKB_inter_core, prop_inter_oncoKB = prop_inter_oncoKB,
        n_COSMIC = n_COSMIC, n_COSMIC_core = n_COSMIC_core,
        prop_COSMIC_inter_core = prop_COSMIC_inter_core, prop_inter_COSMIC = prop_inter_COSMIC, stringsAsFactors = FALSE
      )
    } else {
      Edge_Score_row <- data.frame(
        pathway = path_name, cancer = cancer, method = "Edge_Score", n_enrich_pw = nrow(enrich_Edge_Score),
        n_core_miR = NA, n_CMC = NA, n_CMC_core = NA,
        prop_inter_core = NA, prop_inter_CMC = NA,
        n_core_gene = NA, n_oncoKB = NA, n_oncoKB_core = NA,
        prop_oncoKB_inter_core = NA, prop_inter_oncoKB = NA,
        n_COSMIC = NA, n_COSMIC_core = NA,
        prop_COSMIC_inter_core = NA, prop_inter_COSMIC = NA, stringsAsFactors = FALSE
      )
    }


    # Edge_Network
    if (nrow(enrich_Edge_Network) > 0){
      split_list <- strsplit(core_Edge_Network, ":")
      all_elements <- unlist(split_list)
      core_miR_Edge_Network  <- unique(all_elements[grepl("^hsa-", all_elements)])
      core_gene_Edge_Network <- unique(all_elements[!grepl("^hsa-", all_elements)])

      CMC_miR <- unique(pathway_MGI %>% filter(miRNA %in% cancer_list$CMC_miR, pathway %in% enrich_Edge_Network$pathway) %>% pull(miRNA))
      n_core_miR = length(core_miR_Edge_Network)
      n_CMC = length(CMC_miR)
      n_CMC_core = length(intersect(CMC_miR, core_miR_Edge_Network))
      prop_inter_core = length(intersect(CMC_miR, core_miR_Edge_Network))/length(core_miR_Edge_Network)
      prop_inter_CMC = length(intersect(CMC_miR, core_miR_Edge_Network))/length(CMC_miR)

      gene <- unique(pathway_MGI %>% filter(gene %in% cancer_list$oncoKB_gene, pathway %in% enrich_Edge_Network$pathway) %>% pull(gene))
      n_core_gene = length(core_gene_Edge_Network)
      n_oncoKB = length(gene)
      n_oncoKB_core = length(intersect(gene, core_gene_Edge_Network))
      prop_oncoKB_inter_core = n_oncoKB_core / n_core_gene
      prop_inter_oncoKB = n_oncoKB_core / n_oncoKB

      gene <- unique(pathway_MGI %>% filter(gene %in% cancer_list$COSMIC_gene, pathway %in% enrich_Edge_Network$pathway) %>% pull(gene))
      n_COSMIC = length(gene)
      n_COSMIC_core = length(intersect(gene, core_gene_Edge_Network))
      prop_COSMIC_inter_core = n_COSMIC_core / n_core_gene
      prop_inter_COSMIC = n_COSMIC_core / n_COSMIC

      Edge_Network_row <- data.frame(
        pathway = path_name, cancer = cancer, method = "Edge_Network", n_enrich_pw = nrow(enrich_Edge_Network),
        n_core_miR = n_core_miR, n_CMC = n_CMC, n_CMC_core = n_CMC_core,
        prop_inter_core = prop_inter_core, prop_inter_CMC = prop_inter_CMC,
        n_core_gene = n_core_gene, n_oncoKB = n_oncoKB, n_oncoKB_core = n_oncoKB_core,
        prop_oncoKB_inter_core = prop_oncoKB_inter_core, prop_inter_oncoKB = prop_inter_oncoKB,
        n_COSMIC = n_COSMIC, n_COSMIC_core = n_COSMIC_core,
        prop_COSMIC_inter_core = prop_COSMIC_inter_core, prop_inter_COSMIC = prop_inter_COSMIC, stringsAsFactors = FALSE
      )
    } else {
      Edge_Network_row <- data.frame(
        pathway = path_name, cancer = cancer, method = "Edge_Network", n_enrich_pw = nrow(enrich_Edge_Network),
        n_core_miR = NA, n_CMC = NA, n_CMC_core = NA,
        prop_inter_core = NA, prop_inter_CMC = NA,
        n_core_gene = NA, n_oncoKB = NA, n_oncoKB_core = NA,
        prop_oncoKB_inter_core = NA, prop_inter_oncoKB = NA,
        n_COSMIC = NA, n_COSMIC_core = NA,
        prop_COSMIC_inter_core = NA, prop_inter_COSMIC = NA, stringsAsFactors = FALSE
      )
    }

    df <- bind_rows(TG_Score_row, Edge_ORA_row, Edge_Score_row, Edge_Network_row)
    summary <- rbind(summary, df)
  }
}

summary$sqrt_CMC = sqrt(summary$prop_inter_core * summary$prop_inter_CMC)
summary$sqrt_oncoKB = sqrt(summary$prop_oncoKB_inter_core * summary$prop_inter_oncoKB)
summary$sqrt_COSMIC = sqrt(summary$prop_COSMIC_inter_core * summary$prop_inter_COSMIC)
summary$F1_CMC = 2 * summary$prop_inter_core * summary$prop_inter_CMC / (summary$prop_inter_core + summary$prop_inter_CMC)
summary$F1_oncoKB = 2 * summary$prop_oncoKB_inter_core * summary$prop_inter_oncoKB / (summary$prop_oncoKB_inter_core + summary$prop_inter_oncoKB)
summary$F1_COSMIC = 2 * summary$prop_COSMIC_inter_core * summary$prop_inter_COSMIC / (summary$prop_COSMIC_inter_core + summary$prop_inter_COSMIC)

write.csv(summary, file = paste0(result_dir, "result.csv"), row.names = FALSE)

# 2. draw ----
summary <- read.csv(paste0(result_dir, "result.csv"))
# sqrt (geometric-median)
df_long <- summary %>%
  pivot_longer(
    cols = c(sqrt_CMC, sqrt_oncoKB, sqrt_COSMIC),
    names_to = "Metric",
    values_to = "Value"
  )
df_long$method <- factor(df_long$method, levels = c("TG_Score", "Edge_ORA", "Edge_Score", "Edge_Network"))
df_long$Metric <- factor(df_long$Metric,
                         levels = c("sqrt_CMC", "sqrt_oncoKB", "sqrt_COSMIC"),
                         labels = c("CMC (miRNAs)", "OncoKB (genes)", "COSMIC (genes)"))

g_median_plot <- ggplot(na.omit(df_long), aes(x = method, y = Value, fill = method)) +
  # geom_violin(alpha = 0.6, width = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 1, color = "black", size = 0.75, width = 0.6) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.75, color = "black") +
  scale_fill_manual(values = fill_col) +
  geom_point(size = 2) + # position = position_jitter(width = 0.2),
  facet_wrap(~ Metric, scales = "fixed") + # "free_y"
  scale_y_continuous(
    limits = c(0.03, 0.57),
    breaks = seq(0, 0.55, by = 0.1),
    expand = c(-0.05,0)
  ) +
  theme_bw() +
  labs(title = "Geometric Median for Cancer-related miRNAs/genes Identification", x = NULL, y = "Geometric Median") +
  theme(
    panel.grid = element_blank(),
    # axis.line = element_line(color = "black"), #, linewidth = 0.5
    # axis.ticks = element_line(color = "black"), # , linewidth = 0.5)
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(color = "black", linewidth = 0.75),
    strip.text = element_text(size = 12, face = "bold"),
    axis.ticks.length = unit(0.2, "cm"),
    # strip.text = element_text(size = 8),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.border = element_rect(color = "black", linewidth = 0.75),
    legend.position = "none"
  )

ggsave(paste0(result_dir, "summary.pdf"), g_median_plot, width = 8, height = 8)

