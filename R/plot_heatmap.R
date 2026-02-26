# only support one method each time
#' Description: codes for drawing heatmap plot.
#' @param method a character that specifies only one method which you have the corresponding enrichment result (only five edge-based methods are supported).
#' @param result the output of miREA() function. Just use the unprocessed result list as the input here.
#' @param input_data the output of get_all_input_data() function. Just use the unprocessed input_data list as the input here.
#' @param n_mir_heatmap The number of miRNAs included in the heatmap. Default is 10.
#' @param n_pathway The number of top enriched pathways included in the heatmap. Default it 10. If the number of enriched pathways is less than n_pathway, all enriched pathways will be included.
#' @param plot_path The directory to save the plot. If not null, it will generate a pdf file named heatmap_[method].pdf at the plot path. If null, the plot will be drawn at the current window.
#' @param penrichCutoff  The threshold for deciding enrichment of pathways. Default is 0.05.
#' @param height The height of the plot. If not, it will be set to default values.
#' @param width The width of the plot. If not, it will be set to default values.
#' @param gene_annot A list contains the annotations of genes, where each sub-element denotes a character vector containing a list of genes. It will generate a annotation bar on the top of the heatmap to show whether the genes hit in the annotation list. Default is NULL, which means no annotation bar will be generated.
#' @param mir_annot Same format as gene_annot. If not NULL, it will generate an annotation bar for miRNAs at the left side of the heatmap, to show whether miRNAs hit in the given list.
#' @param annot_color The color of the annotation bars for genes and miRNAs

#' @return A heatmap plot either in the current window, or saved as heatmap_[method].pdf in the plot_path.


plot_heatmap <- function(method, result, input_data, n_mir_heatmap = 10, n_pathway = 10,
                                plot_path = NULL, penrichCutoff = 0.05, height = NULL, width = NULL,
                                gene_annot = NULL, mir_annot = NULL, annot_color = "#A6DDBA"){
  if (length(method) != 1){
    stop("Please specify only ONE [method] each time!")
  }
  if (! method %in% names(result$result)){
    stop("No valid result for method:", method, ". Please check your result file or method parameter!")
  }
  result <- list(
    result = result$result[method],
    time = result$time,
    parameter = result$parameter
  )

  # methods <- names(result$result)
  if (sum(grepl("Edge", method)) == 0){
    stop("Only useful when the methods include Edge-based method.")
  }

  if (!is.null(gene_annot)) {

    if (!is.list(gene_annot)) {
      stop("gene_annot must be a list when not NULL.")
    }

    if (is.null(names(gene_annot)) || any(names(gene_annot) == "")) {
      stop("gene_annot must be a *named* list.")
    }

    if (!all(sapply(gene_annot, is.character))) {
      stop("Each element of gene_annot must be a character vector of gene names.")
    }
  }

  if (!is.null(mir_annot)) {

    if (!is.list(mir_annot)) {
      stop("mir_annot must be a list when not NULL.")
    }

    if (is.null(names(mir_annot)) || any(names(mir_annot) == "")) {
      stop("mir_annot must be a *named* list.")
    }

    if (!all(sapply(mir_annot, is.character))) {
      stop("Each element of mir_annot must be a character vector of miRNA names.")
    }
  }

  n_mir <- n_mir_heatmap
  pathway <- input_data$pathway$Edge
  all_paths <- unique(pathway$pathway)


  enriched_path <- result$result[[method]] %>% arrange(padj) %>% filter(pathway %in% all_paths, padj < penrichCutoff)
  if (nrow(enriched_path) > n_pathway){
    enriched_path <- enriched_path[1:n_pathway, ]
  } else if (nrow(enriched_path) == 0 || is.null(enriched_path)){
    warning("No enriched items with padj value smaller than", penrichCutoff, ". Plot the top ", n_pathway, " with smallest padj values instead.\n")
    enriched_path <- result$result[[method]] %>% arrange(padj) %>% slice_head(n = n_pathway)
  }
  all_paths <- unique(enriched_path$pathway)
  pathway <- pathway %>% filter(pathway %in% all_paths)


  # pathway <- input_data$pathway$Edge %>% filter(pathway %in% all_paths)
  data <- input_data$data$data
  # result <- result$result

  path_MGI <- table(factor(pathway$pathway, levels = all_paths))

  pathway_miRs <- unique(pathway$miRNA)
  pathway_TGs <- unique(pathway$gene)

  DEMGI <- data %>% filter(cor < 0, corPadj < 0.05, abs(mirLogFC) > 1, mirPadj < 0.05, abs(geneLogFC) > 1, genePadj < 0.05) %>%
    mutate(MGI = paste0(miRNA, ":", gene))
  DEG <- data %>% filter(abs(geneLogFC) > 1, genePadj < 0.05) %>% pull(gene)
  MGIList <- data %>% arrange(strength) %>% select(miRNA, gene, strength) %>% filter(miRNA %in% pathway_miRs, gene %in% pathway_TGs)

  DEMGI_path <- pathway %>% inner_join(DEMGI, by = c("miRNA", "gene")) %>% filter(pathway %in% all_paths)
  path_DEMGI <- table(factor(DEMGI_path$pathway, levels = all_paths))

  ## 1. Annotation -----
  #### 1.1 pathway-level annotation: path_order, path_MGI, path_DEMGI -----
  # filter MGI dataframe
  n_total <- length(unique(MGIList$miRNA))
  n_keep <- min(n_mir, n_total)
  DEmiRs <- data %>% filter(abs(mirLogFC) > 1, mirPadj < 0.05) %>% pull(miRNA) %>% unique()

  top_miRNA <- MGIList %>% mutate(MGI = paste0(miRNA, ":", gene)) %>%
    mutate(is_DE = ifelse(MGI %in% unique(DEMGI$MGI), 1, 0)) %>%
    # select(miRNA, is_DE) %>% distinct() %>%
    group_by(miRNA) %>% summarise(n_DE = sum(is_DE)) %>%
    ungroup() %>% arrange(desc(n_DE)) %>%
    slice_head(n = n_keep) %>%
    pull(miRNA)
  #包含在当前通路中的gene与miRNA组成的MGI中，差异表达MGI的个数最多的前n_keep个。
  # 即，gene必须包含在当前all_paths中，miRNA不一定要差异表达，但是只有当差异表达miRNA不足n_keep个时才会被画入图中（意味着没有DEMGI，因为DEMGI的前提时miRNA差异表达）

  MGIList <- MGIList %>% filter(miRNA %in% top_miRNA)


  MGIList <- MGIList %>% group_by(miRNA) %>%
    arrange(strength) %>%
    mutate(gene_count = n_distinct(gene)) %>%
    mutate(row_number = row_number()) %>%
    filter(gene_count <= 10 | row_number <= 10) %>%
    ungroup()
  # 如果一个miRNA调控超过10个基因，只保留strength最小的前10个
  unique_gene <- unique(MGIList$gene)

  #### 1.4 gene bottom annotation: gene_anno dataframe ----
  gene_anno <- data %>% select(gene, geneLogFC, genePadj) %>% distinct() %>%
    filter(gene %in% unique_gene) %>%
    arrange(geneLogFC)
  gene_order <- gene_anno$gene

  # MGIList <- get_top_mir(MGIList, min(n_mir, length(unique(MGIList$miRNA)))) # dataframe: miRNA, gene, strength, only n_mir unique miRNAs involved.
  # A. build miRNA-gene matrix: miR_mat, miR_DElabel_mat ----
  # unique_mir <- unique(MGIList$miRNA)
  # unique_gene <- unique(MGIList$gene)
  miR_mat <- matrix(NA, nrow = length(top_miRNA), ncol = length(gene_order),
                    dimnames = list(top_miRNA, gene_order))
  for (i in seq_len(nrow(MGIList))){
    miR_mat[MGIList$miRNA[i], MGIList$gene[i]] <- MGIList$strength[i]
  }

  miR_DElabel_mat <- matrix(NA, nrow = length(top_miRNA), ncol = length(gene_order),
                            dimnames = list(top_miRNA, gene_order))
  DEMGIList <- MGIList %>% mutate(MGI = paste0(miRNA, ":", gene)) %>% filter(MGI %in% DEMGI$MGI)
  for (i in seq_len(nrow(DEMGIList))){
    miR_DElabel_mat[DEMGIList$miRNA[i], DEMGIList$gene[i]] <- 1
  }

  round_up_to <- function(x, step = 0.1) {
    ceiling(x / step) * step
  }

  miR_vals <- miR_mat[!is.na(miR_mat)]
  miR_min <- min(miR_vals)
  miR_max <- max(miR_vals)

  # v <- ceiling(max(abs(miR_vals)))
  v <- round_up_to(max(abs(miR_vals), na.rm = TRUE))

  # unit color map
  if (sign(miR_min) != sign(miR_max)){
    mir_col <- colorRamp2(
      breaks = c(-v, 0, v),
      colors = c("#ff6700", "grey80", "#238B45") # lightgreen
    )
  } else if (miR_max < 0){
    mir_col <- colorRamp2(
      breaks = c(miR_min, 0),
      colors = c("#ff6700",  "#FFD7B5")#"#FFB38A"
    )
  } else{
    mir_col <- colorRamp2(
      breaks = c(0, miR_max),
      colors = c("#C7E9C0", "#238B45")
    )
  }


  #### 1.2 miRNA-level annotation: mir_target, mir_DEtarget ----
  bg <- input_data$background_MGI
  mir_target <- table(factor(bg$miRNA, levels = top_miRNA))
  # mir_target <- table(bg$miRNA)[top_miRNA]

  bg_DEG <- bg %>% filter(gene %in% DEG) %>% distinct()
  mir_DEtarget <- table(factor(bg_DEG$miRNA, levels = top_miRNA))

  bg_DEMGI <- bg %>% mutate(MGI = paste0(miRNA, ":", gene)) %>% filter(MGI %in% unique(DEMGI$MGI))
  mir_DEMGI <- table(factor(bg_DEMGI$miRNA, levels = top_miRNA))

  DEMGI_path <- DEMGI_path %>% select(miRNA, gene) %>% unique()
  mir_DEMGI_path <- table(factor(DEMGI_path$miRNA, levels = top_miRNA))

  #### 1.3 miRNA right annotation: mir_anno dataframe ----
  mir_anno <- data %>% select(miRNA, mirLogFC, mirPadj) %>% distinct() %>%
    filter(miRNA %in% top_miRNA)
  mir_anno$miRNA <- factor(mir_anno$miRNA, levels = top_miRNA)
  mir_anno <- mir_anno[order(mir_anno$miRNA), ]
  mir_max <- round_up_to(max(abs(mir_anno$mirLogFC)))

  #### 1.4 gene bottom annotation: gene_anno dataframe ----
  gene_max <- round_up_to(max(abs(gene_anno$geneLogFC)))

  #### 1.5 pathway method padj right annotation: method_anno list ----
  padj_df <- result$result[[method]][, c("pathway", "padj")] %>% filter(pathway %in% all_paths)
  padj_df$pathway <- factor(padj_df$pathway, levels = all_paths)
  method_anno <- na.omit(padj_df[order(padj_df$pathway),])

  #### 1.6 miRNA left cancer-lelated annotaion ----
  if (!is.null(mir_annot)){
    mir_anno_df <- lapply(mir_annot, function(gset) {
      top_miRNA %in% gset
    })
    mir_anno_df <- as.data.frame(mir_anno_df, check.names = FALSE)
    rownames(mir_anno_df) <- top_miRNA

    anno_col <- setNames(
      replicate(ncol(mir_anno_df), c("TRUE" = annot_color, "FALSE" = "white"), simplify = FALSE),
      colnames(mir_anno_df)
    )

    left_mir_annot <- rowAnnotation(
      mir_tar = anno_text(
        paste('(', as.integer(mir_DEMGI_path), "/", as.integer(mir_DEMGI), "/",
              as.integer(mir_DEtarget), '/', as.integer(mir_target), ')', sep = ''),
        gp = gpar(fontsize = 12)
      ),
      mir_annot_flag = anno_simple(
        as.data.frame(lapply(mir_anno_df, as.character)),
        col = if (ncol(mir_anno_df) == 1) anno_col[[1]] else anno_col,
        border = FALSE,
        width = unit(rep(4, ncol(mir_anno_df)), "mm")
      ),
      show_annotation_name = TRUE,
      annotation_name_side = "top",
      annotation_name_rot = 90,
      annotation_name_offset = unit(-3, "mm"),
      annotation_label = list(
        mir_annot_flag = colnames(mir_anno_df)          # 使用 mir_annot 的列名
      ),
      show_legend = FALSE
    )
  } else {
    left_mir_annot = rowAnnotation(
      mir_tar = anno_text(
        paste('(', as.integer(mir_DEMGI_path), "/", as.integer(mir_DEMGI), "/",
              as.integer(mir_DEtarget), '/', as.integer(mir_target), ')', sep = ''),
        gp = gpar(fontsize = 12)
      )
    )
  }

  #
  # B. build pathway-gene matrix: path_mat, path_DElabel_mat ----
  path_mat <- matrix(0, nrow = length(all_paths), ncol = length(gene_order),
                     dimnames = list(all_paths, gene_order))
  pathway_filter <- pathway %>% filter(pathway %in% all_paths, gene %in% gene_order)
  for (i in seq_len(nrow(pathway_filter))) {
    pw <- pathway_filter$pathway[i]
    gene <- pathway_filter$gene[i]
    if (pw %in% all_paths && gene %in% gene_order) {
      path_mat[pw, gene] <- 1
    }
  }

  #### 1.7 gene upper annotation:gene_annot (whether includes in cancer-related genes)----
  if (!is.null(gene_annot)){
    gene_anno_df <- lapply(gene_annot, function(gset) {
      gene_order %in% gset
    })
    gene_anno_df <- as.data.frame(gene_anno_df, check.names = FALSE)
    rownames(gene_anno_df) <- gene_order
    anno_col <- setNames(
      replicate(
        ncol(gene_anno_df),
        c("TRUE" = annot_color, "FALSE" = "white"),
        simplify = FALSE
      ),
      colnames(gene_anno_df)
    )

    top_gene_annot <- HeatmapAnnotation(
      df = gene_anno_df,
      col = anno_col,
      annotation_height = unit(rep(4, ncol(gene_anno_df)), "mm"),
      annotation_name_side = "left",
      show_legend = FALSE
    )

  } else {
    top_gene_annot = NULL
  }


  #### 1.8 Methods annotation (right)----


  anno_list <- anno_simple(
    method_anno$padj,
    col = colorRamp2(c(0, penrichCutoff, 1), c("#ff6666","mistyrose", "#f7f7f7")),
    which = "row"#, width = unit(1, "cm")
  )

  row_anno <- rowAnnotation(
    padj = anno_list,
    show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_rot = 90,
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 12)
  )

  ## 2. Legends ----
  ## miRNA log2FoldChange
  col_mir_fc <- colorRamp2(c(-mir_max,0,mir_max),c("#add8e6",  "white", "#ffb6c1")) # c("#0571b0", "#f7f7f7", "#ca0020") ,c("#add8e6",  "white", "#ffb6c1")
  range_mir_fc <- seq(-mir_max, mir_max, length = 5)
  lgd_mir_fc <- Legend(title = "miRNA Log2FC", col_fun = col_mir_fc, at = range_mir_fc, labels = paste(range_mir_fc))
  ## miRNA padj
  col_mir_padj <- colorRamp2(c(0, 0.05, 1),  c("#ff6666","mistyrose", "#f7f7f7")) # "#f03b20", "#feb24c", "#ffeda0"
  lgd_mir_padj <- Legend(title = "miRNA Padj", col_fun = col_mir_padj, at = c(0, 0.05, 1))
  ## method padj
  # lgd_method_list <- lapply(methods, function(method) {
  #   Legend(title = paste0(method, " Padj"), col_fun = colorRamp2(c(0, 0.05, 1), c(scales::alpha(col_fill[[method]], 0.7), "grey90", "grey90")), at = c(0, 0.05, 1))
  # })
  lgd_method <- Legend(title = paste0(method, " Padj"), col_fun = colorRamp2(c(0, penrichCutoff, 1), c("#ff6666", "mistyrose", "grey90")), at = c(0, penrichCutoff, 1))
  # names(lgd_method_list) <- methods
  ## gene log2FoldChange
  lgd_gene_fc_cat <- Legend(title = "Gene Expression \n(* shows DEG)", legend_gp = gpar(fill = c("blue", "red")), labels = c("downregulated", "upregulated"))
  col_gene_fc <- colorRamp2(c(-gene_max,0,gene_max),c("#add8e6",  "white", "#ffb6c1"))
  range_gene_fc <- seq(-gene_max, gene_max, length = 5)
  lgd_gene_fc <- Legend(title = "Gene Log2FC", col_fun = col_gene_fc, at = range_gene_fc, labels = paste(range_gene_fc))

  ## DE rectangle
  # lgd_DEG <- Legend(labels = "DEG", type = "points", pch = 22, legend_gp = gpar(fill = "grey50", col = "red", lwd = 2))
  # lgd_DEMGI <- Legend(labels = "DEMGI", type = "points", pch = 22, legend_gp = gpar(fill = "white", col = "#B266CC", lwd = 2))
  # lgd_DE <- Legend(title = "DE", labels = c("DEG", "DEMGI"), type = c("points", "points"), pch = c(22,22), legend_gp = gpar(fill = c("grey50", "white"), col = c("red", "#b266cc"), lwd = c(2,2)))
  lgd_DE <- Legend(title = "MGI Expression", labels = c("DEMGI"),  type = c("points"), pch = c(22),
                   legend_gp = gpar(fill = c("white"), col = c("#b266cc"), lwd = c(2)), size = unit(5, "mm"))

  # sankey plot annotation
  # lgd_expression <- Legend(title = "Expression", legend_gp = gpar(fill = c("#FFB3B3", "#99D2E7")), labels = c("upregulated", "downregulated"))
  # lgd_strength <- Legend(title = "Strength", legend_gp = gpar(lty = c("dashed", "solid"), col = c("black", "black")),type = "lines", labels = c("positive", "negative"))
  # lgd_enrich <- Legend(title = "Enriched", legend_gp = gpar(fill = c("#ff6666", "grey")), labels = c("TRUE", "FALSE"))

  # cancer-related gene
  if (!is.null(gene_annot)){
    gene_top_lgd <- Legend(title = "Cancer-related Genes", legend_gp = gpar(fill = c(annot_color, "white")), labels = c("TRUE", "FALSE"))
  } else {
    gene_top_lgd <- NULL
  }

  if (!is.null(mir_annot)){
    mir_left_lgd <- Legend(title = "Cancer-related miRNAs", legend_gp = gpar(fill = c(annot_color, "white")), labels = c("TRUE", "FALSE"))
  } else {
    mir_left_lgd <- NULL
  }

  ht_list <- c(gene_top_lgd, lgd_method, lgd_gene_fc_cat, mir_left_lgd, lgd_DE, lgd_gene_fc, lgd_mir_fc, lgd_mir_padj) #
  #sankey_list <- list(lgd_expression, lgd_strength) # , lgd_enrich
  #legend_list <- c(ht_list, sankey_list)
  legend_list <- ht_list


  # packed_legends <- packLegend(
  #   list = legend_list,
  #   direction = "vertical",
  #   gap = unit(5, "mm")
  # )

  ## 3. Heatmaps ----
  ht_unit <- min(c(0.0264*(800-50)/dim(path_mat)[1],0.5))
  wd_unit <- min(c(0.0264*(1700-100)/dim(path_mat)[2],0.5))

  #### heatmap 1 (gene-pathway)----
  h1_ht <- ceiling(ht_unit*dim(path_mat)[1])
  h1_wd <- ceiling(wd_unit*dim(path_mat)[2])

  ht1 <- ComplexHeatmap::Heatmap(
    path_mat,
    name = "ht1",
    column_title = paste0("Heatmap plot for ", method, " method"),
    column_title_gp = gpar(fontsize = 20, fontface = "bold", hjust = 0.5),
    col  = c("1" = "grey50", "0" = "white"),
    column_names_rot = 90,
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_heatmap_legend = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 12),
    width = unit(h1_wd, "cm"),
    height = unit(h1_ht, "cm"),
    column_gap = unit(c(2), "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,

    left_annotation = rowAnnotation(totgene_text = anno_text(
      paste('(', as.integer(path_DEMGI),'/',as.integer(path_MGI),')', sep = ''),
      gp = gpar(fontsize = 12)
    )),

    right_annotation = row_anno,
    top_annotation = top_gene_annot,
    bottom_annotation = columnAnnotation(
      geneLog2FC = anno_barplot(
        gene_anno$geneLogFC,
        border = FALSE,
        baseline = 0,
        ylim = c(-gene_max, gene_max),
        gp = gpar(
          fontsize = 12,
          col = ifelse(gene_anno$geneLogFC < 0, "blue", "red"),
          fill = ifelse(gene_anno$geneLogFC < 0, "blue", "red")
        ),
        height = unit(1.5, "cm")
      ),
      star = anno_text(
        ifelse(abs(gene_anno$geneLogFC) > 1, "*", ""),
        gp = gpar(fontsize = 18, col = ifelse(gene_anno$geneLogFC < 0, "blue", "red")),
        # location = unit(1, "npc"),
        # location = unit(gene_anno$geneLogFC + 0.2 * sign(gene_anno$geneLogFC), "native"),
        location = 0.5,
        just = "center",
        rot = 0
      ),
      show_legend = FALSE,
      show_annotation_name = TRUE,
      annotation_name_rot = 0,
      annotation_name_side = "left"
    )
  )

  #### heatmap 2 (gene-miRNA) ----

  h2_ht <- ceiling(ht_unit*dim(miR_mat)[1])
  h2_wd <- ceiling(wd_unit*dim(miR_mat)[2])

  ht2 <- ComplexHeatmap::Heatmap(
    miR_mat,
    name = "ht2",
    na_col = "white",
    col = mir_col,
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = 12),
    show_row_names = TRUE,
    row_names_side = 'left',
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_gp = gpar(fontsize = 12),
    # bottom_annotation = sankey_anno,
    width = unit(h2_wd, "cm"),
    height = unit(h2_ht, "cm"),
    column_gap = unit(c(2), "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      # grid.rect(x, y, width, height, gp = gpar(col = "grey70", lwd = 0.5, fill = NA))
      if (!is.na(miR_DElabel_mat[i, j]))
      {
        grid.rect(x,y,width = width,height = height,
                  gp = gpar(lwd = 3, col = "#B266CC",fill = fill))# #B19CD9
      }
    },
    heatmap_legend_param = list(
      direction = "vertical",
      title = "MGI Strength",
      col = mir_col,
      na_col = "white"
    ),
    left_annotation = left_mir_annot,
    # left_annotation = rowAnnotation(mir_tar = anno_text(
    #   paste('(', as.integer(mir_DEMGI_path), "/", as.integer(mir_DEMGI), "/", as.integer(mir_DEtarget), '/', as.integer(mir_target), ')', sep = ''),
    #   gp = gpar(fontsize = 12)
    # )),

    right_annotation = rowAnnotation(
      mir_fc = mir_anno$mirLogFC,
      mir_padj = mir_anno$mirPadj,
      col = list(
        mir_fc = col_mir_fc,
        mir_padj = col_mir_padj
      ),
      show_legend = FALSE,
      show_annotation_name = TRUE,
      annotation_name_rot = 90,
      annotation_name_side = "top",
      annotation_width = unit(c(1, 1), "cm")  # 每列宽度设置为1cm
    )
  )

  ## 4. Combine ----
  ht_list <- ht1 %v% ht2

  if (!is.null(plot_path)){
    if (is.null(height)){height = 2* h1_ht + 2*h2_ht} # unit((h1_ht + h2_ht + n_mir_sankey) * 0.8, "cm")
    if (is.null(width)){width = h1_wd*1.2} # unit(h1_wd * 0.9, "cm")
    pdf(paste0(plot_path, "/heatmap_", method,".pdf"), height =  height, width = width)
    draw(ht_list, merge_legend = TRUE,
         annotation_legend_list = legend_list,
         heatmap_legend_side = "bottom",
         annotation_legend_side = "bottom",
         align_heatmap_legend = "heatmap_center")
    decorate_heatmap_body("ht1", {
      grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
    })

    decorate_heatmap_body("ht2", {
      grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
    })
    dev.off()
    cat("    Please check the heatmap plot for", method, "at:", plot_path, "\n")
    return(paste0("Please check the heatmap plot for", method, "at: ", plot_path))
  } else {
    draw(ht_list, merge_legend = TRUE,
         annotation_legend_list = legend_list,
         column_title = paste0("Heatmap for ", method, " method"),
         column_title_gp = gpar(fontsize = 20, fontface = "bold"),
         heatmap_legend_side = "right",
         annotation_legend_side = "right")
    decorate_heatmap_body("ht1", {
      grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
    })

    decorate_heatmap_body("ht2", {
      grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
    })
    cat("  Visualization finished!\n***If you want to save them to pdf files, please specify the [plot_path] parameter!***")
    return(plot = list(ht_list, legend_list))
  }

}


