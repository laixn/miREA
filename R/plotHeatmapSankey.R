#' Visualize an enrichment result as a combined heatmap-Sankey plot
#'
#' Draws a heatmap of the miRNA-gene-pathway enrichment result for one
#' edge-based method, combined with a Sankey diagram showing the
#' miRNA-gene-pathway flow. Only supports one method each time.
#'
#' @param method a character that specifies only one method which you have the corresponding enrichment result (only five edge-based methods are supported).
#' @param result the output of miREA() function, used here as the unprocessed result list.
#' @param input_data the output of getAllInputData() function, used here as the unprocessed input_data list.
#' @param n_mir_heatmap The number of miRNAs included in the heatmap. Default is 20.
#' @param n_mir_sankey The number of miRNAs included in the Sankey diagram. Default is 10.
#' @param n_pathway The number of top enriched pathways included in the heatmap. Default it 10. If the number of enriched pathways is less than n_pathway, all enriched pathways will be included.
#' @param plot_path The directory to save the plot. If not null, it will generate a pdf file named heatmap_sankey_[method].pdf at the plot path. If null, the plot will be drawn at the current window.
#' @param penrichCutoff  The threshold for deciding enrichment of pathways. Default is 0.05.
#' @param height The height of the plot. If not, it will be set to default values.
#' @param width The width of the plot. If not, it will be set to default values.
#' @param sankey_prop The ratio of the height of Sankey diagram and heatmap. Default is 1.5.
#' @param gene_annot A list contains the annotations of genes, where each sub-element denotes a character vector containing a list of genes. It will generate a annotation bar on the top of the heatmap to show whether the genes hit in the annotation list. Default is NULL, which means no annotation bar will be generated.
#' @param mir_annot Same format as gene_annot. If not NULL, it will generate an annotation bar for miRNAs at the left side of the heatmap, to show whether miRNAs hit in the given list.
#' @param annot_color The color of the annotation bars for genes and miRNAs
#' @param DElogFCCutoff optional, the absolute log2FC threshold used to call a miRNA/gene differentially
#'        expressed for the heatmap/Sankey annotation. If NULL (the default), falls back to
#'        `input_data$threshold$DElogFCCutoff` (set by [getAllInputData()]), or 1 if that is absent.
#' @param DEpadjCutoff optional, the adjusted p-value threshold used to call a miRNA/gene differentially
#'        expressed for the heatmap/Sankey annotation. If NULL (the default), falls back to
#'        `input_data$threshold$DEpadjCutoff`, or 0.05 if that is absent.
#' @param corCutoff optional, the miRNA-gene expression correlation threshold for the heatmap/Sankey
#'        annotation. If NULL (the default), falls back to `input_data$threshold$corCutoff`, or 0 if
#'        that is absent.
#' @param corPadjCutoff optional, the adjusted p-value threshold for the miRNA-gene expression
#'        correlation, for the heatmap/Sankey annotation. If NULL (the default), falls back to
#'        `input_data$threshold$corPadjCutoff`, or 0.05 if that is absent.

#' @return A heatmap-sankey plot either in the current window, or saved as heatmap_sankey_[method].pdf in the plot_path.
#' @examples
#' data(toy_pathway); data(toy_gene_DEdata); data(toy_mir_DEdata)
#' data(toy_gene_mat); data(toy_mir_mat); data(background_MGI)
#' # Edge_Score2D is included so input_data carries the MGI "strength" scores
#' # plotHeatmapSankey() needs for its heatmap cells, without running it in miREA().
#' input_data <- getAllInputData(
#'   methods = c("Edge_ORA", "Edge_Score2D"), pathway = toy_pathway,
#'   mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
#'   background_MGI = background_MGI,
#'   gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
#' )
#' background <- list(
#'   Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))
#' )
#' result <- miREA(
#'   methods = "Edge_ORA", input_data = input_data, background = background,
#'   minSize = NULL, maxSize = NULL, pvalueCutoff = 1
#' )
#' plotHeatmapSankey(method = "Edge_ORA", result = result, input_data = input_data, penrichCutoff = 1)
#' @export
plotHeatmapSankey <- function(method, result, input_data, n_mir_heatmap = 20, n_mir_sankey = 10,  n_pathway = 10,
                                plot_path = NULL, penrichCutoff = 0.05, height = NULL, width = NULL, sankey_prop = 1.5,
                                gene_annot = NULL, mir_annot = NULL, annot_color = "#A6DDBA",
                                DElogFCCutoff = NULL, DEpadjCutoff = NULL,
                                corCutoff = NULL, corPadjCutoff = NULL){
  if (length(method) != 1){
    stop("Please specify only ONE [method] each time!")
  }

  thr <- input_data[["threshold"]]
  DElogFCCutoff <- DElogFCCutoff %||% thr[["DElogFCCutoff"]] %||% 1
  DEpadjCutoff  <- DEpadjCutoff  %||% thr[["DEpadjCutoff"]]  %||% 0.05
  corCutoff     <- corCutoff     %||% thr[["corCutoff"]]     %||% 0
  corPadjCutoff <- corPadjCutoff %||% thr[["corPadjCutoff"]] %||% 0.05
  if (! method %in% names(result$result)){
    stop("No valid result for method:", method, ". Please check your result file or method parameter!")
  }
  result <- list(
    result = result$result[method],
    time = result$time,
    parameter = result$parameter
  )

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

    if (!all(vapply(gene_annot, is.character, logical(1)))) {
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

    if (!all(vapply(mir_annot, is.character, logical(1)))) {
      stop("Each element of mir_annot must be a character vector of miRNA names.")
    }
  }

  n_mir <- n_mir_heatmap
  pathway <- input_data$pathway$Edge
  all_paths <- unique(pathway$pathway)


  enriched_path <- result$result[[method]] %>% arrange(padj) %>% filter(pathway %in% all_paths, padj < penrichCutoff)
  if (nrow(enriched_path) > n_pathway){
    enriched_path <- enriched_path[seq_len(n_pathway), ]
  } else if (nrow(enriched_path) == 0 || is.null(enriched_path)){
    warning("No enriched items with padj value smaller than", penrichCutoff, ". Plot the top ", n_pathway, " with smallest padj values instead.\n")
    enriched_path <- result$result[[method]] %>% arrange(padj) %>% slice_head(n = n_pathway)
  }
  all_paths <- unique(enriched_path$pathway)
  pathway <- pathway %>% filter(pathway %in% all_paths)

  data <- input_data$data$data

  path_MGI <- table(factor(pathway$pathway, levels = all_paths))

  pathway_miRs <- unique(pathway$miRNA)
  pathway_TGs <- unique(pathway$gene)

  DEMGI <- data %>% filter(cor < corCutoff, corPadj < corPadjCutoff, abs(mirLogFC) > DElogFCCutoff, mirPadj < DEpadjCutoff, abs(geneLogFC) > DElogFCCutoff, genePadj < DEpadjCutoff) %>%
    mutate(MGI = paste0(miRNA, ":", gene))
  DEG <- data %>% filter(abs(geneLogFC) > DElogFCCutoff, genePadj < DEpadjCutoff) %>% pull(gene)
  MGIList <- data %>% arrange(strength) %>% select(miRNA, gene, strength) %>% filter(miRNA %in% pathway_miRs, gene %in% pathway_TGs)

  DEMGI_path <- pathway %>% inner_join(DEMGI, by = c("miRNA", "gene")) %>% filter(pathway %in% all_paths)
  path_DEMGI <- table(factor(DEMGI_path$pathway, levels = all_paths))

  ## 1. Annotation -----
  #### 1.1 pathway-level annotation: path_order, path_MGI, path_DEMGI -----
  n_total <- length(unique(MGIList$miRNA))
  n_keep <- min(n_mir, n_total)
  DEmiRs <- data %>% filter(abs(mirLogFC) > DElogFCCutoff, mirPadj < DEpadjCutoff) %>% pull(miRNA) %>% unique()

  top_miRNA <- MGIList %>% mutate(MGI = paste0(miRNA, ":", gene)) %>%
    mutate(is_DE = ifelse(MGI %in% unique(DEMGI$MGI), 1, 0)) %>%
    group_by(miRNA) %>% summarise(n_DE = sum(is_DE)) %>%
    ungroup() %>% arrange(desc(n_DE)) %>%
    slice_head(n = n_keep) %>%
    pull(miRNA)
  # Among the MGIs formed by genes/miRNAs within the current pathway, take the top n_keep
  # by number of DE MGIs. I.e. the gene must be in the current all_paths; the miRNA need not
  # be DE, but is only drawn when there are fewer than n_keep DE miRNAs (meaning no DEMGI
  # exists for it, since a DEMGI requires a DE miRNA).

  MGIList <- MGIList %>% filter(miRNA %in% top_miRNA)


  MGIList <- MGIList %>% group_by(miRNA) %>%
    arrange(strength) %>%
    mutate(gene_count = n_distinct(gene)) %>%
    mutate(row_number = row_number()) %>%
    filter(gene_count <= 10 | row_number <= 10) %>%
    ungroup()
  # if a miRNA regulates more than 10 genes, keep only the 10 with the smallest strength
  unique_gene <- unique(MGIList$gene)

  #### 1.4 gene bottom annotation: gene_anno dataframe ----
  gene_anno <- data %>% select(gene, geneLogFC, genePadj) %>% distinct() %>%
    filter(gene %in% unique_gene) %>%
    arrange(geneLogFC)
  gene_order <- gene_anno$gene

  # A. build miRNA-gene matrix: miR_mat, miR_DElabel_mat ----
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
      colors = c("#ff6700",  "#FFD7B5")
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
      replicate(ncol(mir_anno_df),
                c("TRUE" = annot_color, "FALSE" = "white"),
                simplify = FALSE),
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
        mir_annot_flag = colnames(mir_anno_df)
      ),
      show_legend = FALSE
    )

  } else {
    left_mir_annot <- rowAnnotation(
      mir_tar = anno_text(
        paste('(', as.integer(mir_DEMGI_path), "/", as.integer(mir_DEMGI), "/",
              as.integer(mir_DEtarget), '/', as.integer(mir_target), ')', sep = ''),
        gp = gpar(fontsize = 12)
      )
    )
  }

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
    top_gene_annot <- NULL
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
  col_mir_fc <- colorRamp2(c(-mir_max,0,mir_max),c("#add8e6",  "white", "#ffb6c1"))
  range_mir_fc <- seq(-mir_max, mir_max, length = 5)
  lgd_mir_fc <- Legend(title = "miRNA Log2FC", col_fun = col_mir_fc, at = range_mir_fc, labels = paste(range_mir_fc))
  ## miRNA padj
  col_mir_padj <- colorRamp2(c(0, 0.05, 1),  c("#ff6666","mistyrose", "#f7f7f7"))
  lgd_mir_padj <- Legend(title = "miRNA Padj", col_fun = col_mir_padj, at = c(0, 0.05, 1))
  ## method padj
  lgd_method <- Legend(title = paste0(method, " Padj"), col_fun = colorRamp2(c(0, penrichCutoff, 1), c("#ff6666", "mistyrose", "grey90")), at = c(0, penrichCutoff, 1))
  ## gene log2FoldChange
  lgd_gene_fc_cat <- Legend(title = "Gene Expression \n(* shows DEG)", legend_gp = gpar(fill = c("blue", "red")), labels = c("downregulated", "upregulated"))
  col_gene_fc <- colorRamp2(c(-gene_max,0,gene_max),c("#add8e6",  "white", "#ffb6c1"))
  range_gene_fc <- seq(-gene_max, gene_max, length = 5)
  lgd_gene_fc <- Legend(title = "Gene Log2FC", col_fun = col_gene_fc, at = range_gene_fc, labels = paste(range_gene_fc))

  ## DE rectangle
  lgd_DE <- Legend(title = "MGI Expression", labels = c("DEMGI"),  type = c("points"), pch = c(22),
                   legend_gp = gpar(fill = c("white"), col = c("#b266cc"), lwd = c(2)), size = unit(5, "mm"))

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

  ht_list <- c(gene_top_lgd, lgd_method, lgd_gene_fc_cat, mir_left_lgd, lgd_DE, lgd_gene_fc, lgd_mir_fc, lgd_mir_padj)
  legend_list <- ht_list

  ## 3. Heatmaps ----
  ht_unit <- min(c(0.0264*(800-50)/dim(path_mat)[1],0.5))
  wd_unit <- min(c(0.0264*(1700-100)/dim(path_mat)[2],0.5))

  #### heatmap 1 (gene-pathway)----
  h1_ht <- ceiling(ht_unit*dim(path_mat)[1])
  h1_wd <- ceiling(wd_unit*dim(path_mat)[2])

  ht1 <- ComplexHeatmap::Heatmap(
    path_mat,
    name = "ht1",
    column_title = paste0("Heatmap and Sankey plot for ", method, " method"),
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
        ifelse(abs(gene_anno$geneLogFC) > DElogFCCutoff, "*", ""),
        gp = gpar(fontsize = 18, col = ifelse(gene_anno$geneLogFC < 0, "blue", "red")),
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

  ##### sankey plot ----
  sankey <- plot_sankey(dict = data, pathway = pathway, enriched_path = enriched_path, all_paths = all_paths, top_miRNA = top_miRNA, mir_max = mir_max, gene_max = gene_max, strength_max = v , penrichCutoff = penrichCutoff, n_mir = n_mir_sankey, n_pathway = n_pathway,
                        DElogFCCutoff = DElogFCCutoff, DEpadjCutoff = DEpadjCutoff, corCutoff = corCutoff, corPadjCutoff = corPadjCutoff)
  sankey <- ggplotGrob(sankey)
  anno_ggplot <- AnnotationFunction(
    fun = function(index, k, n){
      pushViewport(viewport()) # width = unit(4, "cm"), height = unit(6, "cm")
      grid.draw(sankey)
      popViewport()
    },
    height = unit((h2_ht+h1_ht)*sankey_prop,"cm"),# unit(n_mir_sankey * 1.5, "cm"),#  unit((n_mir_heatmap/n_mir_sankey) * 8, "cm"),
    subsettable = FALSE,
    which = "column"
  )

  sankey_anno <- HeatmapAnnotation(
    PR = anno_ggplot,
    show_annotation_name = FALSE
  )

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
    bottom_annotation = sankey_anno,
    width = unit(h2_wd, "cm"),
    height = unit(h2_ht, "cm"),
    column_gap = unit(c(2), "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(miR_DElabel_mat[i, j]))
      {
        grid.rect(x,y,width = width,height = height,
                  gp = gpar(lwd = 3, col = "#B266CC",fill = fill))
      }
    },
    heatmap_legend_param = list(
      direction = "vertical",
      title = "MGI Strength",
      col = mir_col,
      na_col = "white"
    ),
    left_annotation = left_mir_annot,
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
      annotation_width = unit(c(1, 1), "cm")
    )
  )

  ## 4. Combine ----
  ht_list <- ht1 %v% ht2

  if (!is.null(plot_path)){
    if (is.null(height)){height <- 2* h1_ht + 2*h2_ht} # unit((h1_ht + h2_ht + n_mir_sankey) * 0.8, "cm")
    if (is.null(width)){width <- h1_wd*1.2} # unit(h1_wd * 0.9, "cm")
    pdf(paste0(plot_path, "/heatmap_sankey_", method,".pdf"), height =  height, width = width)
    draw(ht_list, merge_legend = TRUE,
         annotation_legend_list = legend_list,
         heatmap_legend_side = "right",
         annotation_legend_side = "right")
    decorate_heatmap_body("ht1", {
      grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
    })

    decorate_heatmap_body("ht2", {
      grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
    })
    dev.off()
    message("    Please check the heatmap & sankey plot for ", method, " at: ", plot_path)
    return(paste0("Please check the heatmap & sankey plot for", method, "at: ", plot_path))
  } else {
    draw(ht_list, merge_legend = TRUE,
         annotation_legend_list = legend_list,
         column_title = paste0("Heatmap and Sankey plot for ", method, " method"),
         column_title_gp = gpar(fontsize = 20, fontface = "bold"),
         heatmap_legend_side = "right",
         annotation_legend_side = "right")
    decorate_heatmap_body("ht1", {
      grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
    })

    decorate_heatmap_body("ht2", {
      grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
    })
    message("  Visualization finished!\n***If you want to save them to pdf files, please specify the [plot_path] parameter!***")
    return(plot = list(ht_list, legend_list))
  }

}

# data (dict), pathway
# all_paths, top_miRNA，
# col_mir_fc()
# col_gene_fc()
# max_strength #(for strength)

# stringr::str_wrap() only breaks on whitespace, so underscore-joined names
# (e.g. MSigDB pathway IDs like "LINDGREN_BLADDER_CANCER_CLUSTER_2B") are
# treated as a single unbreakable word and never wrap. Break on underscores
# as well as spaces, inserting a newline in place of the delimiter without
# otherwise altering the original characters (so real spaces in a label are
# left as spaces, not turned into underscores).
wrap_underscored_label <- function(x, width) {
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    tokens <- regmatches(s, gregexpr("[^_ ]*[_ ]|[^_ ]+$", s))[[1]]
    lines <- character(0)
    line <- ""
    for (tok in tokens) {
      candidate <- paste0(line, tok)
      if (nchar(line) > 0 && nchar(sub("[_ ]$", "", candidate)) > width) {
        lines <- c(lines, sub("[_ ]$", "", line))
        line <- tok
      } else {
        line <- candidate
      }
    }
    lines <- c(lines, sub("[_ ]$", "", line))
    paste(lines, collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

plot_sankey <- function(dict, pathway, enriched_path, all_paths, top_miRNA, mir_max, gene_max, strength_max , penrichCutoff = 0.05, n_mir = 25, n_pathway = 10,
                         DElogFCCutoff = 1, DEpadjCutoff = 0.05, corCutoff = 0, corPadjCutoff = 0.05){

  pathway <- pathway %>% filter(pathway %in% all_paths)
  dict_DE <- dict %>%
    filter(cor < corCutoff, corPadj < corPadjCutoff, genePadj < DEpadjCutoff, abs(geneLogFC) >= DElogFCCutoff, mirPadj < DEpadjCutoff, abs(mirLogFC) >= DElogFCCutoff) %>%
    arrange(desc(abs(strength)))
  path_DEMGI <- pathway %>% inner_join(dict_DE, by = c("miRNA", "gene")) %>% filter(pathway %in% all_paths)

  pw_MGI <- unique(paste0(pathway$miRNA, ":", pathway$gene))
  dict_DEMGI <- dict_DE %>%
    filter(paste0(miRNA,":", gene) %in% pw_MGI)
  n_total <- length(top_miRNA)
  n_keep <- min(n_mir, n_total)
  unique_mir <- top_miRNA[seq_len(n_keep)]
  MGI <- dict_DEMGI %>% filter(miRNA %in% unique_mir)
  # if a miRNA has more than 5 target genes, keep only the 5 with the smallest strength
  MGI <- MGI %>% group_by(miRNA) %>%
    arrange(strength) %>%
    mutate(gene_count = n_distinct(gene)) %>%
    mutate(row_number = row_number()) %>%
    filter(gene_count <= 5 | row_number <= 5) %>%
    ungroup()
  unique_gene <- unique(MGI$gene)

  MGI <- MGI %>% mutate(sign = ifelse(strength > 0, "pos", "neg"), flow_value = 1)

  ## sankey plot ----
  # mutate MGI expression profile, pathway member profile, and pathway enrichment profile.
  df_sankey <- pathway %>%
    inner_join(MGI, by = c("miRNA", "gene"))
  df_sankey <- df_sankey %>%
    inner_join(enriched_path, by = "pathway") %>%
    arrange(desc(abs(strength)))

  df_mir_gene <- df_sankey %>%
    select(miRNA, gene, flow_value,  sign,mirLogFC, geneLogFC, strength) %>%
    distinct()
  df_gene_path <- df_sankey %>%
    select(gene, pathway, geneLogFC, padj) %>%
    distinct()


  # order by log2FC
  gene_order <- df_mir_gene %>% arrange(desc(geneLogFC)) %>% pull(gene) %>% unique()
  mir_order <- df_mir_gene %>% arrange(desc(mirLogFC)) %>% pull(miRNA) %>% unique()
  pw_order <- df_gene_path %>% arrange(padj) %>% pull(pathway) %>% unique()
  df_mir_gene$gene <- factor(df_mir_gene$gene, levels = gene_order)
  df_mir_gene$miRNA <- factor(df_mir_gene$miRNA, levels = mir_order)
  df_gene_path$gene <- factor(df_gene_path$gene, levels = gene_order)
  df_gene_path$pathway <- factor(df_gene_path$pathway, levels = pw_order)

  df_mir_gene_long <- ggalluvial::to_lodes_form(df_mir_gene, key = "Demographic", axes = seq_len(2)) %>% distinct()
  df_gene_path_long <- ggalluvial::to_lodes_form(df_gene_path, key = "Demographic", axes = seq_len(2)) %>% distinct()


  gene_flow_from_miR <- df_mir_gene_long %>%
    filter(Demographic == "gene") %>%
    group_by(stratum) %>%
    summarise(total_flow = sum(flow_value), .groups = "drop") %>%
    rename(gene = stratum)
  gene_path_counts <- df_gene_path_long %>%
    filter(Demographic == "gene") %>%
    group_by(stratum) %>%
    summarise(n_pathway = n(), .groups = "drop") %>%
    rename(gene = stratum)
  gene_info <- gene_flow_from_miR %>%
    left_join(gene_path_counts, by = "gene") %>%
    mutate(flow_value = total_flow / n_pathway)
  df_gene_path_long <- df_gene_path_long %>%
    left_join(gene_info %>% select(gene, gene_flow_value = flow_value),
              by = c("stratum" = "gene"))
  pathway_info <- df_gene_path_long %>%
    filter(Demographic == "pathway") %>%
    select(alluvium, pathway = stratum)

  gene_info_per_edge <- df_gene_path_long %>%
    filter(Demographic == "gene") %>%
    select(alluvium, gene = stratum, gene_flow_value)

  pathway_flow <- pathway_info %>%
    left_join(gene_info_per_edge, by = "alluvium") %>%
    group_by(pathway) %>%
    summarise(flow_value = sum(gene_flow_value, na.rm = TRUE)/n(), .groups = "drop")

  df_gene_path_long <- df_gene_path_long %>%
    mutate(flow_value = case_when(
      Demographic == "gene" ~ gene_flow_value,
      TRUE ~ NA_real_
    )) %>%
    left_join(pathway_flow, by = c("stratum" = "pathway")) %>%
    mutate(flow_value = if_else(Demographic == "pathway", flow_value.y, flow_value.x)) %>%
    select(-gene_flow_value, -flow_value.y, - flow_value.x)

  miRNA_fill_df <- df_mir_gene_long %>%
    filter(Demographic == "miRNA") %>%
    distinct(stratum, mirLogFC)
  gene_fill_df <- df_mir_gene_long %>%
    filter(Demographic == "gene") %>%
    distinct(stratum, geneLogFC)
  pathway_fill_df <- df_gene_path_long %>%
    filter(Demographic == "pathway") %>%
    distinct(stratum, padj)

  col_mir_fc <- colorRamp2(c(-mir_max,0,mir_max),c("#add8e6",  "white", "#ffb6c1"))
  miRNA_fill_df <- miRNA_fill_df %>%
    mutate(fill_color = col_mir_fc(mirLogFC))

  # gene: geneLogFC mapped light-blue (down) -> white -> light-pink (up)
  col_gene_fc <- colorRamp2(c(-gene_max,0,gene_max),c("#add8e6",  "white", "#ffb6c1"))
  gene_fill_df <- gene_fill_df %>%
    mutate(fill_color = col_gene_fc(geneLogFC))

  color_func <- scales::col_numeric(
    palette = c("#ff6666", "mistyrose"),
    domain = c(0, penrichCutoff)
  )

  pathway_fill_df <- pathway_fill_df %>%
    mutate(
      fill_color = NA_character_,
      fill_color = replace(fill_color, padj <= penrichCutoff, color_func(padj[padj <= penrichCutoff])),
      fill_color = replace(fill_color, padj > penrichCutoff, "#f7f7f7")
    )

  stratum_fill_df <- bind_rows(
    miRNA_fill_df %>% mutate(Demographic = "miRNA"),
    gene_fill_df %>% mutate(Demographic = "gene"),
    pathway_fill_df %>% mutate(Demographic = "pathway")
  ) %>% select(stratum, Demographic, fill_color)


  df_mir_gene_long <- df_mir_gene_long %>%
    left_join(stratum_fill_df, by = c("stratum", "Demographic"))

  df_gene_path_long <- df_gene_path_long %>%
    left_join(stratum_fill_df, by = c("stratum", "Demographic"))

  # eusure the order of the stratum is miRNA-gene-pathway
  df_mir_gene_long$Demographic <- factor(df_mir_gene_long$Demographic,
                                         levels = c("miRNA", "gene", "pathway"))
  df_gene_path_long$Demographic <- factor(df_gene_path_long$Demographic,
                                          levels = c("miRNA", "gene", "pathway"))


  gene_colors <- setNames(rep("#B3D49D90", length(unique(df_mir_gene$gene))), unique(df_mir_gene$gene))
  sankey <- ggplot() +
    geom_flow(data = df_mir_gene_long,
              aes(x = Demographic,
                  stratum = stratum,
                  alluvium = alluvium,
                  y = flow_value,
                  fill = strength),
              stat = ggalluvial::StatAlluvium,
              color = "purple",
              alpha = 0.7,
              aes.flow = "forward") +
    scale_fill_gradientn(
      colors = c("#ff6700", "grey80", "#238B45"), limits = c(-strength_max, strength_max),
      values = scales::rescale(c(-strength_max, 0, strength_max)), breaks = c(-strength_max, -strength_max/2, 0, strength_max/2, strength_max), name = "Strength", guide = "none"
    ) +
    new_scale_fill() +
    geom_flow(data = df_gene_path_long,
              aes(x = Demographic,
                  stratum = stratum,
                  alluvium = alluvium,
                  fill = stratum,
                  y = flow_value
              ),
              color = "black",
              stat = ggalluvial::StatAlluvium,
              alpha = 0.7,
              aes.flow = "forward") +
    scale_fill_manual(values = gene_colors, guide = "none") +
    new_scale_fill() +
    geom_stratum(data = df_mir_gene_long,
                 aes(x = Demographic, stratum = stratum, y = flow_value, fill = fill_color)) + #
    geom_stratum(data = df_gene_path_long,
                 aes(x = Demographic, stratum = stratum, y = flow_value, fill = fill_color)) + #
    scale_x_discrete(limits = c("miRNA", "gene", "pathway"), expand = c(0.02, 0.15), position = "bottom") +
    geom_text(data = df_mir_gene_long,
              aes(x = Demographic, stratum = stratum, label = after_stat(stratum), y = flow_value), #
              stat = ggalluvial::StatStratum, size = 3, color = "black") +
    geom_text(data = subset(df_gene_path_long, Demographic == "pathway"),
              aes(x = Demographic, stratum = stratum, y = flow_value, label = wrap_underscored_label(after_stat(stratum), width = 22)),
              stat = ggalluvial::StatStratum, size = 3, color = "black") +
    scale_fill_identity() +

    theme_void() +
    theme(
      axis.text.x = element_text(size = 16, face = "bold", color = "black") ,
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    ) +
    scale_y_continuous(expand = c(0.02,0))


  return(sankey)

}

