# This is the script to draw Figure 4 for BLCA case study:
# 1. (A-C) summary plot (plot/summary.pdf)
# 2. (D) heatmap (plot/heatmap_Edge_Network.pdf)
# 3. (E) sankey plot for the experimental evidence (plot/sankey.pdf).

setwd("/scratch/project_2011179/code/miREA/") # change your own directory here
library(dplyr)
library(ggalluvial)
library(ComplexHeatmap)
library(circlize) # colorRamp2()
library(ggnewscale) # new_scale_fill()
library(patchwork)
library(scales)

library(readxl)
library(ggplot2)
library(tidyr)

result_dir <- "analysis/4_case_study/"
plot_path <- "analysis/4_case_study/plot/"
cancer = "BLCA"
path_name = "hallmark"
penrichCutoff = 0.05
methods <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network")

### 1. summary and heatmap ----
# load result
enrich_dir <- paste0("results/", path_name, "/", cancer, "/")
load(paste0(enrich_dir, "RData/result.RData"))
invalid_methods <- setdiff(names(result$result), methods)
if (length(invalid_methods) > 0){
  result$result[[invalid_methods]] <- NULL
}

# load input data
input_dir <- paste0("data/input_data/", path_name, "/")
load(paste0(input_dir, cancer, "_", path_name, "_input_data.RData"))

# load gene and miRNA annotation
load("data/cancer_list/cancer_list.RData")

# load plot function
source("R/plot_heatmap_sankey.R")
source("R/plot_summary.R")
source("R/plot_heatmap.R")


gene_annot <- list("COSMIC" = cancer_list$COSMIC_gene,
                   "OncoKB" = cancer_list$oncoKB_gene)
mir_annot <- list("CMC" = cancer_list$CMC_miR) #

ht_plot <- plot_heatmap(method = "Edge_Network", result = result, input_data = input_data, n_mir_heatmap = 10, n_pathway = 10, plot_path = plot_path, penrichCutoff = penrichCutoff,
                        height = 8, width = 16, gene_annot = gene_annot, mir_annot = mir_annot, annot_color = "palegreen3")

ht_sankey_plot <- plot_heatmap_sankey(method = "Edge_Network", result = result, input_data = input_data, n_mir_heatmap = 10, n_mir_sankey = 5,  n_pathway = 10, plot_path = plot_path, penrichCutoff = penrichCutoff,
                                height = 12, width = 18, sankey_prop = 1.5, gene_annot = gene_annot, mir_annot = mir_annot, annot_color = "palegreen3") # #A6DDBA"


summary_plot <- plot_summary(result = result, penrichCutoff = penrichCutoff, plot_path = plot_path,
                             fill_col = c("TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
                                          "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_2Ddist" = "#626FB3",
                                          "Edge_Topology" = "#EAA5C2","Edge_Network" = "#FA8072"))

data <- input_data$data$data
data <- data %>% mutate(DEMGI = ifelse(cor < 0 & corPadj <= 0.05 & abs(mirLogFC) >= 1 & mirPadj <= 0.05 & abs(geneLogFC) >= 1 & genePadj <= 0.05, TRUE, FALSE))
write.csv(data, file = paste0(result_dir,"data.csv"), row.names = FALSE)

### sankey plot for references ----
table <- readxl::read_excel(paste0(result_dir, "sankey_table.xlsx"))
if ("miREA" %in% unique(table$Type)){
  table$Type[table$Type == "miREA"] = "Edge-Network"
}

expr_col_map <- c("upregulated" = "#ffb6c1", "downregulated" = "#add8e6")
cancer_col_map <- c("BLCA" = "#66c2a5", "LUAD" = "#fc8d62", "CRC" = "#e78ac3")
table2 <- table %>%
  mutate(miR_fill = expr_col_map[miR_type],
         gene_fill = expr_col_map[gene_type],
         cancer_fill = cancer_col_map[Cancer],
         pw_fill = "#fd8c8a",
         flow_value = 1)

table2$Cancer <- factor(table2$Cancer, levels = unique(table2$Cancer))
table2$MiRNA <- factor(table2$MiRNA, levels = unique(table2$MiRNA))
table2$Gene <- factor(table2$Gene, levels = unique(table2$Gene))
table2$Type <- factor(table2$Type, levels = c("Edge-Network", "reference"))

cancer_mir_gene_table <- table2 %>% select(Cancer, MiRNA, Gene, Type, miR_fill, gene_fill, cancer_fill, flow_value) %>% distinct()
gene_pathway_table <- table2 %>% select(Gene, Pathway, Type, gene_fill, pw_fill, flow_value) %>% distinct()

df_mir_gene_long <- ggalluvial::to_lodes_form(cancer_mir_gene_table, key = "Demographic", axes = 1:3) %>% distinct()
df_gene_path_long <- ggalluvial::to_lodes_form(gene_pathway_table, key = "Demographic", axes = 1:2) %>% distinct()

# eusure the order of the stratum is cancer-miRNA-gene-pathway
df_mir_gene_long$Demographic <- factor(df_mir_gene_long$Demographic,
                                       levels =c("Cancer", "MiRNA", "Gene", "Pathway"))
df_mir_gene_long <- df_mir_gene_long %>%
  mutate(fill = case_when(
    Demographic == "Cancer" ~ cancer_fill,
    Demographic == "MiRNA"  ~ miR_fill,
    Demographic == "Gene"   ~ gene_fill
  ))
df_gene_path_long$Demographic <- factor(df_gene_path_long$Demographic,
                                        levels =c("Cancer", "MiRNA", "Gene", "Pathway"))
df_gene_path_long <- df_gene_path_long %>%
  mutate(fill = case_when(
    Demographic == "Gene" ~ gene_fill,
    Demographic == "Pathway" ~ pw_fill
  ))

# To ensure the cancer-miR-gene panel and gene-pathway panel have the same height
gene_flow_from_miR <- df_mir_gene_long %>%
  filter(Demographic == "Gene") %>%
  group_by(stratum) %>%  # stratum 是 gene
  summarise(total_flow = sum(flow_value), .groups = "drop") %>%
  rename(Gene = stratum)
gene_path_counts <- df_gene_path_long %>%
  filter(Demographic == "Gene") %>%
  group_by(stratum) %>%
  summarise(n_pathway = n(), .groups = "drop") %>%
  rename(Gene = stratum)
gene_info <- gene_flow_from_miR %>%
  left_join(gene_path_counts, by = "Gene") %>%
  mutate(flow_value = total_flow / n_pathway)
df_gene_path_long <- df_gene_path_long %>%
  left_join(gene_info %>% select(Gene, gene_flow_value = flow_value),
            by = c("stratum" = "Gene"))

# pathway info
pathway_info <- df_gene_path_long %>%
  filter(Demographic == "Pathway") %>%
  select(alluvium, pathway = stratum)

# gene node and flow_value
gene_info_per_edge <- df_gene_path_long %>%
  filter(Demographic == "Gene") %>%
  select(alluvium, gene = stratum, gene_flow_value)

# edges to connect gene and pathway
pathway_flow <- pathway_info %>%
  left_join(gene_info_per_edge, by = "alluvium") %>%
  group_by(pathway) %>%
  summarise(flow_value = sum(gene_flow_value, na.rm = TRUE)/n(), .groups = "drop")

df_gene_path_long <- df_gene_path_long %>%
  mutate(flow_value = case_when(
    Demographic == "Gene" ~ gene_flow_value,
    TRUE ~ NA_real_
  )) %>%
  left_join(pathway_flow, by = c("stratum" = "pathway")) %>%
  mutate(flow_value = if_else(Demographic == "Pathway", flow_value.y, flow_value.x)) %>%
  select(-gene_flow_value, -flow_value.y, - flow_value.x)

sankey <- ggplot() +
  geom_flow(data = df_mir_gene_long,
            aes(x = Demographic,
                stratum = stratum,
                alluvium = alluvium,
                y = flow_value,
                fill = Type),
            stat = "alluvium",
            color = "black",
            alpha = 0.7,
            aes.flow = "forward") +
  scale_fill_manual(name = "Source", values = c("reference" = "#C4E6C3", "Edge-Network" = "#FDBF6F")) +
  new_scale_fill() +
  geom_flow(data = df_gene_path_long,
            aes(x = Demographic,
                stratum = stratum,
                alluvium = alluvium,
                y = flow_value,
                fill = Type
            ),
            color = "black",
            stat = "alluvium",
            alpha = 0.7,
            aes.flow = "forward") +
  scale_fill_manual(name = "Source", values = c("reference" = "#C4E6C3", "Edge-Network" = "#FDBF6F")) +
  new_scale_fill() +
  geom_stratum(data = df_mir_gene_long,
               aes(x = Demographic, stratum = stratum, y = flow_value, fill = fill)) +
  geom_stratum(data = df_gene_path_long,
               aes(x = Demographic, stratum = stratum, y = flow_value, fill = fill)) +
  scale_x_discrete(limits = c("Cancer", "MiRNA", "Gene", "Pathway"), expand = c(0.02, 0), position = "bottom") +
  geom_text(data = df_mir_gene_long,
            aes(x = Demographic, stratum = stratum, label = after_stat(stratum), y = flow_value), #
            stat = "stratum", size = 4, color = "black") +
  geom_text(data = subset(df_gene_path_long, Demographic == "Pathway"),
            aes(x = Demographic, stratum = stratum, y = flow_value, label = stringr::str_wrap(after_stat(stratum), width = 20)),
            stat = "stratum", size = 4, color = "black") +
  #scale_fill_identity() +
  scale_fill_identity(
    name = "Expression",
    breaks = c("#ffb6c1", "#add8e6"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend"
  ) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 16, face = "bold", color = "black") ,
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
    legend.position = "left",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10)
  ) +
  scale_y_continuous(expand = c(0.02,0.02))

ggsave(paste0(work_dir, "/plot/sankey.pdf"), sankey, width = 17, height = 8)


