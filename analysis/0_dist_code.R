# This is the script to summarize the distributions, including:
# 1. boxplot_genesize.pdf: Figure S2.
# 2. boxplot_MGI_score.pdf: Figure S3.
# 3. density_centrality_hallmark.pdf: Figure S4.
# 4. boxplot_pvalue_Topology.pdf: FIgure S5.

setwd("/scratch/project_2011179/code/miREA/")
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)

pw_dir <- "data/raw_data/pathway/"
result_dir <- "results/"
paths <- c("Reactome", "hallmark", "cancer-specific")
input_data_dir <- "data/input_data/hallmark/"

# 1. gene set size distributions for three pathways----

pw_hallmark <- read.csv(paste0(pw_dir, "hallmark/hallmark_gene.csv"))
pw_reactome <- read.csv(paste0(pw_dir, "Reactome/Reactome_gene.csv"))
pw_TP <- read.csv(paste0(pw_dir, "TP/TP_pathway_all_cancers.csv"))

size_hallmark <- pw_hallmark %>%
  group_by(pathway) %>%
  summarise(
    geneSize = n()
  ) %>%
  ungroup() %>%
  mutate(path_name = "hallmark")

size_reactome <- pw_reactome %>%
  group_by(pathway) %>%
  summarise(
    geneSize = n()
  ) %>%
  ungroup() %>%
  mutate(path_name = "Reactome")

size_TP <- pw_TP %>%
  group_by(pathway) %>%
  summarise(
    geneSize = n()
  ) %>%
  ungroup() %>%
  mutate(path_name = "cancer-specific")

size_df <- rbind(size_hallmark, size_reactome, size_TP)
size_df$path_name <- factor(size_df$path_name, levels = paths)

p_boxplot <- ggplot(size_df, aes(x = path_name, y = geneSize, fill = path_name)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  labs(
    title = "Gene set size distribution for pathway databases",
    x = "Pathway database",
    y = "Gene set size"
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(color = "black", linewidth = 0.75),
    strip.text = element_text(size = 12, face = "bold"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", linewidth = 0.75)
  )
ggsave("analysis/boxplot_genesize.pdf", p_boxplot, width = 6, height = 6)

# 2. MGI scores distribution for all cancers ----
files <- list.files(
  path = input_data_dir,
  pattern = "_hallmark_input_data\\.RData$",
  full.names = TRUE
)

MGI_list <- vector("list", length(files))

for (i in seq_along(files)) {

  file <- files[i]

  cancer <- sub("_hallmark_input_data\\.RData$", "",
                basename(file))
  cat(cancer,"\n")

  env <- new.env()
  load(file, envir = env)

  tmp_data <- env$input_data$data$data

  tmp_data$cancer <- cancer

  MGI_list[[i]] <- tmp_data
}

MGI_data <- do.call(rbind, MGI_list)
MGI_data$cancer <- factor(MGI_data$cancer,
                          levels = sort(unique(MGI_data$cancer)))

# draw
p <- ggplot(MGI_data, aes(x = cancer, y = strength)) +
  geom_violin(aes(fill = cancer), color = "black")+
  # scale_fill_brewer(palette = "Set3")+
  scale_fill_hue()+
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "MGI Scores Distribution Across Cancer Types",
    x = NULL,
    y = "MGI score"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(color = "black", linewidth = 0.75),
    strip.text = element_text(size = 12, face = "bold"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", linewidth = 0.75)
  )

ggsave("analysis/boxplot_MGI_score.pdf", p,
       width = 9,
       height = 6)
# 3. Edge betweenness distribution for 10 cancer hallmarks ----
load(paste0(input_data_dir, "BLCA_hallmark_input_data.RData"))

pathway_list <- split(input_data$pathway$Edge, input_data$pathway$Edge$pathway)
pathway_GGI_list <- split(input_data$pathway_GGI, input_data$pathway_GGI$pathway)

pathway_name <- unique(input_data$pathway$Edge$pathway)
n_pathway <- length(pathway_name)


min_max <- function(x) {
  diff <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (diff != 0){
    scaled <- (x - min(x, na.rm = TRUE)) / diff
  } else {
    scaled <- 0
  }
}

compute_centrality <- function(i){
  name <- pathway_name[i]
  cat("      ", i, "/", n_pathway, ",", name, ":")

  PW <- pathway_list[[name]]
  PW_GGI <- pathway_GGI_list[[name]]

  PW <- data.table::as.data.table(PW)
  PW$MGI = paste0(PW$miRNA, ":", PW$gene)
  data.table::setkey(PW, miRNA, gene, MGI)

  PW_name <- unique(PW$pathway)
  mirnas <- unique(PW$miRNA)

  PW_edge <- rbind(
    data.table(from = PW$miRNA, to = PW$gene),
    data.table(from = PW_GGI$from, to = PW_GGI$to)
  )
  PW_edge <- unique(PW_edge)

  g <- igraph::graph_from_data_frame(d = PW_edge, directed = TRUE)
  edge_eb <- igraph::edge_betweenness(g)
  edge_dt <- data.table(igraph::as_data_frame(g, what = "edges"))[, EB := edge_eb]
  edge_dt[, EB := fifelse(is.na(EB), 0, EB)]
  edge_dt[, centrality := min_max(EB) + 0.01]
  edge_dt[, pathway := PW_name]
  return(edge_dt)
}

centrality_result <- dplyr::bind_rows(lapply(seq_len(n_pathway), compute_centrality))


label_df <- centrality_result %>%
  group_by(pathway) %>%
  summarise(
    x = max(centrality, na.rm = TRUE),
    y = Inf
  )

p_dens <- ggplot(centrality_result, aes(x = centrality)) + # , fill = pathway

  # density
  geom_density(color = "black", alpha = 0.7, fill = "grey") +

  # points
  geom_point(
    aes(y = 0),
    shape = 21,
    size = 1.2,
    stroke = 0.3,
    color = "black",
    alpha = 0.8,
    position = position_jitter(height = 1)
  ) +

  # pathway name
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = pathway),
    inherit.aes = FALSE,
    hjust = 0.95,
    vjust = 1.5,
    size = 3
  ) +

  # scale_fill_brewer(palette = "Set3") +
  labs(
    y = "Density",
    x = "Centrality",
    title = "Normalized Edge Betweenness Centrality For 10 Cancer Hallmarks"
  ) +

  facet_wrap(~ pathway, ncol = 2, scales = "free_y") +

  theme_bw() +
  theme(
    panel.grid = element_blank(),

    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.text = element_text(color = "black"),

    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(5,10,5,5),

    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),

    panel.border = element_rect(color = "black", linewidth = 0.75)
  )

ggsave("analysis/density_centrality_hallmark.pdf", p_dens, width = 7, height = 10)

# 4. p-value distributions for Edge-Topology ----
df <- data.frame(
  path_name = character(),
  pathway = character(),
  cancer = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (path_name in paths){
  cat(path_name,"...\n")
  if (path_name == "cancer-specific"){
    path_name <- "TP"
  }

  base_dir <- file.path(result_dir, path_name)
  cancer_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

  df_list <- lapply(cancer_dirs, function(cancer_dir) {

    cancer <- basename(cancer_dir)
    cat(cancer,"\n")
    csv_file <- file.path(cancer_dir, "result",
                          paste0(cancer, "_Edge_Topology.csv"))

    if (!file.exists(csv_file)) return(NULL)

    tmp <- read.csv(csv_file, stringsAsFactors = FALSE)

    if (path_name == "TP"){
      path_name <- "cancer-specific"
    }
    data.frame(
      path_name = path_name,
      pathway   = tmp$pathway,
      cancer    = cancer,
      p_value   = tmp$p_value,
      stringsAsFactors = FALSE
    )
  })

  df <- rbind(df, do.call(rbind, df_list))
}

df <- df %>% left_join(size_df, by = c("path_name", "pathway"))
df$path_name <- factor(df$path_name, levels = paths)




p_p_boxplot <- ggplot(df, aes(x = path_name, y = p_value, fill = path_name)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  scale_fill_brewer(palette="Set3")+
  theme_bw() +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  labs(
    title = "P value distribution for enrichment across cancer types",
    x = "Pathway database",
    y = "P value"
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(color = "black", linewidth = 0.75),
    strip.text = element_text(size = 12, face = "bold"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", linewidth = 0.75)
  )
ggsave("analysis/boxplot_pvalue_Topology.pdf", p_p_boxplot, width = 6, height = 6)
