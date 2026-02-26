# This is the script to generate true positive rates heatmap across 16 cancers
# Generates TPR_ht.pdf: Figure 2A.

setwd("/scratch/project_2011179/code/miREA/") # change your own directory here
library(dplyr)
library(readr)
library(tidyr)

library(ComplexHeatmap)
library(ggplot2)

result_dir <- "analysis/2.1_positive_benchmark/"
if (!dir.exists(result_dir)){
  dir.create(result_dir)
}

plot_path <- paste0(result_dir, "TPR_ht.pdf")

path_name = "TP"
method_order <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score",
                  "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network")

fill_col = c("TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
             "Edge_ORA" = "#EDC194", "Edge_Score" = "#F09137", "Edge_2Ddist" = "#9AA4D6",
             "Edge_Topology" = "#EAA5C2","Edge_Network" = "#FA8072")

# 1. get summary table----
# (pathways with padj < 0.05 are regarded as enriched)

for (path_name in path_name){
  cat("Start processing", path_name, "...\n")
  base_dir <- paste0("results/", path_name, "/")
  cancer_dirs <- list.dirs(base_dir, recursive = FALSE)
  df_list <- list()

  for (cancer_dir in cancer_dirs) {
    cancer <- basename(cancer_dir)
    file_path <- file.path(cancer_dir, paste0(cancer, "_", path_name, "_result.csv"))

    if (file.exists(file_path)) {
      df <- read_csv(file_path, col_types = cols())
      df <- df %>% mutate(cancer = cancer)
      df_list[[cancer]] <- df
    } else {
      warning(paste("File not found:", file_path))
    }
  }

  summary_df <- bind_rows(df_list) %>%
    select(cancer, pathway, method, n_enrich, n_pathway, time, positive_rate)
  # summary_df <- summary_df %>% filter(method != "Edge_manova")
  method_order <- unique(summary_df$method)

  # cancer × method positive_rate dataframe
  summary <- summary_df %>%
    group_by(cancer, method) %>%
    summarize(positive_rate = mean(positive_rate, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = method, values_from = positive_rate) %>%
    select(cancer, all_of(method_order))

  write.csv(summary, file = paste0(result_dir, path_name, "_summary.csv"), row.names = FALSE)
  write.csv(summary_df, file = paste0(result_dir, path_name, "_summary_df.csv"), row.names = FALSE)
}

# 2. get TP cancer data annotation ----
cancers <- unique(summary$cancer)
TP_stat <- data.frame(cancer = character(), n_MGI = integer(), n_DEMGI = integer(), DEprop = numeric(), stringsAsFactors = FALSE)
for(cancer in cancers){
  load(paste0("data/input_data/", path_name, "/", cancer, "_", path_name, "_input_data.RData"))
  pathway_MGI <- input_data$pathway$Edge
  pathway_MGIs <- unique(paste0(pathway_MGI$miRNA, ":", pathway_MGI$gene))
  data <- input_data$data$data
  data <- data %>% mutate(MGI = paste0(miRNA, ":", gene)) %>% filter(MGI %in% pathway_MGIs)
  n_MGI <- length(unique(data$MGI))
  data_DEMGI <- data %>% filter(abs(geneLogFC) > 1, genePadj < 0.05, abs(mirLogFC) > 1, mirPadj < 0.05, cor < 0, corPadj < 0.05)
  n_DEMGI = length(unique(data_DEMGI$MGI))
  TP_row <- data.frame(cancer = cancer, n_MGI = n_MGI, n_DEMGI = n_DEMGI, DEprop = n_DEMGI/n_MGI, stringsAsFactors = FALSE)
  TP_stat <- rbind(TP_stat, TP_row)
  rm(input_data)
}

TP_MGI <- read.csv("data/raw_data/pathway/TP/TP_MGI_all_cancers.csv")
TP_size <- TP_MGI %>% dplyr::select(pathway, cancer) %>% dplyr::distinct() %>%
  dplyr::group_by(cancer) %>% dplyr::summarise(n_path = n()) %>% dplyr::ungroup()

TP_stat <- merge(TP_stat, TP_size, by = "cancer") %>%
  arrange(desc(n_path))
cancer_order <- TP_stat$cancer
TP_stat$cancer <- factor(TP_stat$cancer, levels = cancer_order)

write.csv(TP_stat, file = paste0(result_dir, path_name, "_data_stat.csv"), row.names = FALSE)

# 3. get TP summary across all cancer types----
TP_summary_df <- summary_df %>% group_by(method) %>%
  summarise(TP = sum(n_enrich), n_pos = sum(n_pathway)) %>%
  ungroup() %>%
  mutate(TPR = TP / n_pos)
write.csv(TP_summary_df, file = paste0(result_dir, path_name, "_overall_performance.csv"), row.names = FALSE)

# 4. draw heatmap plot----
## load data above ----
TP_input <- read.csv(paste0(result_dir, path_name, "_data_stat.csv"))
cancer_order <- TP_input %>% arrange(desc(n_path)) %>% pull(cancer)

TP_summary <- read.csv(paste0(result_dir, path_name, "_summary.csv"), header = TRUE)
rownames(TP_summary) <- TP_summary$cancer
TP_summary <- TP_summary[, -1]
TP_matrix <- as.matrix(t(TP_summary))
TP_matrix <- TP_matrix[,cancer_order]

TP_overall <- read.csv(paste0(result_dir, path_name, "_overall_performance.csv"))

## right annotation: TP overall performance----
df_long <- TP_overall %>%
  pivot_longer(cols = TPR, names_to = "Metric", values_to = "Value")

node_methods <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score")
edge_methods <- setdiff(unique(TP_overall$method), node_methods)

df_long$method <- factor(df_long$method, levels = rev(method_order))

split_x <- length(edge_methods) + 0.5

fmt_label <- function(x){
  ifelse(x == 0, "0",
         ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

p <- ggplot(df_long, aes(x = method, y = Value, fill = method)) + # fill = Metric
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7,
           color = "black",
           linewidth = 0.5) +
  geom_text(aes(label = fmt_label(Value), color = Metric),
            position = position_dodge(width = 0.75),
            vjust = 0,
            hjust = -0.2,
            size = 2.5,
            angle = 0,
            show.legend = FALSE) +
  geom_vline(xintercept = split_x, color = "grey", linewidth = 0.5) +
  scale_fill_manual(values = fill_col) +
  scale_color_manual(values = c(TPR = "black")) +
  # scale_fill_manual(values = c(TPR = "#F4A7B9")) +
  # scale_color_manual(values = c(TPR = "red3")) +
  scale_y_continuous(
    limits = c(0, max(df_long$Value, na.rm = TRUE) * 2),
    breaks = seq(0, max(df_long$Value, na.rm = TRUE) * 2, by = 0.05),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(title = NULL,
       x = NULL,
       y = NULL,
       fill = NULL) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length.x = unit(-0.1, "cm"),
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),  # Add border
    panel.grid = element_blank()
  )+
  coord_flip()


grob_pr <- ggplotGrob(p)
anno_ggplot <- AnnotationFunction(
  fun = function(index, k, n){
    pushViewport(viewport()) # width = unit(4, "cm"), height = unit(6, "cm")
    grid.draw(grob_pr)
    popViewport()
  },
  width = unit(4, "cm"),
  subsettable = FALSE,
  which = "row"
)
label_text <- "True positive rate"

row_anno <- rowAnnotation(
  PR = anno_ggplot,
  annotation_label = label_text,
  annotation_name_side = "bottom",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_offset = unit(0.3, "cm")
  # gap = unit(-1.5, "pt")
)

## column annotation for n_MGI, DEprop, and n_pathway ----
if (!is.null(TP_input) && all(c("n_MGI", "DEprop") %in% colnames(TP_input))){
  data_annot <- TP_input
  data_annot <- data_annot[data_annot[[1]] %in% colnames(TP_matrix), ]
  n_MGI_range <- c(min(data_annot$n_MGI, na.rm = TRUE), quantile(data_annot$n_MGI, 0.25, na.rm = TRUE), median(data_annot$n_MGI, na.rm = TRUE), quantile(data_annot$n_MGI, 0.75, na.rm = TRUE), max(data_annot$n_MGI, na.rm = TRUE))
  n_MGI_fun = circlize::colorRamp2(n_MGI_range, c("#F9F7E6", "#F7EBD2", "#F4E1AB", "#EDD594", "#E1C379") )
  DEprop_range <- c(min(data_annot$DEprop, na.rm = TRUE), quantile(data_annot$DEprop, 0.25, na.rm = TRUE), median(data_annot$DEprop, na.rm = TRUE), quantile(data_annot$DEprop, 0.75, na.rm = TRUE), max(data_annot$DEprop, na.rm = TRUE))
  DEprop_fun = circlize::colorRamp2(DEprop_range, c("#DFEDD3", "#CEE5B9", "#BFDDA3", "#B1D691", "#92C467") )

  if ("n_path" %in% colnames(data_annot)){
    ha_top <- HeatmapAnnotation(
      "Number of cancer-related pathways" = anno_barplot(
        data_annot$n_path,
        add_numbers = TRUE,
        height = unit(1.5, "cm"),
        gp = gpar(fill = "#C95F7A"),
        axis = FALSE,
        axis_param = list(at = NULL, labels = NULL, labels_rot = NULL)
      ),
      show_annotation_name = TRUE,
      show_legend = FALSE,
      annotation_name_gp = gpar(fontsize = 8)
    )
  }


  ha_bottom <- HeatmapAnnotation(
    "Number of detected TP MGIs" = data_annot$n_MGI,
    "Proportion of detected TP DEMGIs" = data_annot$DEprop,
    col = list("Number of detected TP MGIs" = n_MGI_fun, "Proportion of detected TP DEMGIs" = DEprop_fun),
    show_annotation_name = TRUE,
    show_legend = FALSE,
    annotation_name_gp = gpar(fontsize = 8),
    annotation_name_side = "left",
    annotation_legend_param = list(
      "Number of detected TP MGIs" = list(
        legend_height = unit(1.9, "cm"),
        # at = c(min(data_annot$n_MGI), as.integer(median(data_annot$n_MGI)), max(data_annot$n_MGI)),
        at = n_MGI_range, # pretty(range(data_annot$n_MGI, na.rm = TRUE), n = 3),
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8)
      ),
      "Proportion of detected TP DEMGIs" = list(
        at = pretty(range(data_annot$DEprop, na.rm = TRUE), n = 3),
        legend_height = unit(1.9, "cm"),
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8)
      )
    )
  )
  rm(data_annot)
}

## left annotation: method type ----
left_anno <- function(matrix){
  type_col <- c("Node" = "#589D58" , "MGI" = "#DB8041" ) #"darkgreen","#d95f02"
  data_col <- c("DE" = "cornflowerblue", "Score" = "#FF99C0")
  # left method annotation ----
  annot_df <- data.frame(
    method = c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA", "Edge_Score", "Edge_2Ddist", "Edge_Topology", "Edge_Network"),
    Type = c("Node", "Node","Node","Node","MGI","MGI","MGI","MGI","MGI"),
    Data = c("DE", "Score", "DE", "Score", "DE", "Score", "Score", "Score", "DE")
  )
  annot_df <- annot_df[match(rownames(matrix), annot_df$method), ]
  # method annotation bar
  ha_left <- rowAnnotation(
    Type = annot_df$Type,
    Data = annot_df$Data,
    col = list(
      Type = type_col,
      Data = data_col
    ),
    annotation_label = c("Method Type", "Input Data Type"),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    annotation_legend_param = list(
      Type = list(
        title_gp = gpar(fontsize = 9, fontface = "bold"),
        labels_gp = gpar(fontsize = 8)
      ),
      Data = list(
        title_gp = gpar(fontsize = 9, fontface = "bold"),
        labels_gp = gpar(fontsize = 8)
      )
    )
  )
  return(ha_left)
}
TP_left_anno <- left_anno(TP_matrix)

## main heatmap----
TP_col_fun = circlize::colorRamp2(seq(0,1, length.out = 5), c("#FBEBF0", "#F7D7DE", "#EDB2BF", "#E38B9F", "#C95F7A"))
TP_max_pos <- apply(round(TP_matrix,2), 2, function(col) {
  if (all(is.na(col)) || all(col == 0, na.rm = TRUE)) {
    return(integer(0))
  } else {
    return(which(col == max(col, na.rm = TRUE)))
  }
})

TP_args_list <- list(
  matrix = TP_matrix, name = "TPR", col = TP_col_fun,
  na_col = "grey90",
  column_title = paste0("True Positive Rates Across ", ncol(TP_matrix)," Cancer Types"),
  column_title_gp = gpar(fontface = "bold", fontsize = 12),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  row_names_side = "left",
  # row_split = annot_df$Type,
  left_annotation = TP_left_anno,
  right_annotation = row_anno,
  # top_annotation = ha_top,
  # bottom_annotation = ha_bottom,
  row_dend_side = "left",
  border_gp = gpar(col = "black", lwd = 1.5),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  width = ncol(TP_matrix) * unit(1, "cm"),
  height = nrow(TP_matrix) * unit(1, "cm"),
  show_heatmap_legend = FALSE,
  cell_fun =  function(j, i, x, y, width, height, fill) {
    # grid.rect(x,y, width = width,height = height,
    #            gp = gpar(lwd = 3,col = "#ff6555",fill = fill))
    if (!is.na(TP_matrix[i, j])){
      if(TP_matrix[i, j] == 0){
        grid.text(sprintf("%.0f", TP_matrix[i, j]), x, y, gp = gpar(fontsize = 7.5))
      } else if (sprintf("%.2f", TP_matrix[i, j]) < 0.01){
        grid.text("< 0.01", x, y, gp = gpar(fontsize = 7.5))
      } else{
        grid.text(sprintf("%.2f", TP_matrix[i, j]), x, y, gp = gpar(fontsize = 7.5))
      }
    }
    if (sum(i == as.integer(TP_max_pos[[j]])) > 0) {
      grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 2))
    }
  },
  column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 9),
  #column_names_rot = 45,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)
if (exists("ha_top") && !is.null(ha_top)) {
  TP_args_list$top_annotation <- ha_top
}
if (exists("ha_bottom") && !is.null(ha_bottom)) {
  TP_args_list$bottom_annotation <- ha_bottom
}
TP_ht <- do.call(ComplexHeatmap::Heatmap, TP_args_list)

## legend setting ----
type_col <- c("Node" = "#589D58" , "MGI" = "#DB8041" ) #"darkgreen","#d95f02"
data_col <- c("DE" = "cornflowerblue", "Score" = "#FF99C0")

lgd_method_type <- Legend(title = "Method Type", legend_gp = gpar(fill = c("#589D58", "#DB8041")), labels = c("Node", "MGI"),
                          labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 8, fontface = "bold"))
lgd_input_data <- Legend(title = "Input Data Type", legend_gp = gpar(fill = c("cornflowerblue", "#FF99C0")), labels = c("DE", "Score"),
                         labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 8, fontface = "bold"))
if (exists("ha_bottom") && !is.null(ha_bottom)) {
  lgd_n_MGI <- Legend(title = "num_MGIs", col_fun = n_MGI_fun,
                      labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 8, fontface = "bold"),
                      at = as.integer(n_MGI_range))
  lgd_DEprop <- Legend(title = "prop_DEMGIs", col_fun = DEprop_fun,
                       labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 8, fontface = "bold")
  )
}


legend_list <- c(list(lgd_method_type, lgd_input_data),
                 if (exists("lgd_n_MGI")) list(lgd_n_MGI) else list(),
                 if (exists("lgd_DEprop")) list(lgd_DEprop) else list()#,
                 ) # list(lgd_method)

## draw plot ----
pdf(plot_path, height = 6, width = 11)
draw(TP_ht, merge_legend = TRUE,
     annotation_legend_list = legend_list
     )
if (exists("ha_top")){
  decorate_annotation("Number of cancer-related pathways", {
    grid.rect(gp = gpar(col = "black", lwd = 1.5, fill = NA))
  })
}
dev.off()


