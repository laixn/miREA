#' Build the shared matrices/annotations/legends for the heatmap plots
#'
#' Shared by [plotHeatmap()] and [plotHeatmapSankey()], which otherwise differ
#' only in the Sankey diagram, the heatmap `column_title`, and the `draw()`
#' legend placement.
#' @noRd
build_heatmap_components <- function(method, result, input_data, n_mir_heatmap, n_pathway,
                                     penrichCutoff, gene_annot, mir_annot, annot_color,
                                     DElogFCCutoff, DEpadjCutoff, corCutoff, corPadjCutoff) {
    if (length(method) != 1) {
        stop("Please specify only ONE [method] each time!")
    }

    thr <- input_data[["threshold"]]
    DElogFCCutoff <- DElogFCCutoff %||% thr[["DElogFCCutoff"]] %||% 1
    DEpadjCutoff <- DEpadjCutoff %||% thr[["DEpadjCutoff"]] %||% 0.05
    corCutoff <- corCutoff %||% thr[["corCutoff"]] %||% 0
    corPadjCutoff <- corPadjCutoff %||% thr[["corPadjCutoff"]] %||% 0.05
    if (!method %in% names(result$result)) {
        stop("No valid result for method:", method, ". Please check your result file or method parameter!")
    }
    result <- list(
        result = result$result[method],
        time = result$time,
        parameter = result$parameter
    )

    if (sum(grepl("Edge", method)) == 0) {
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

    enriched_path <- result$result[[method]] %>%
        arrange(padj) %>%
        filter(pathway %in% all_paths, padj < penrichCutoff)
    if (nrow(enriched_path) > n_pathway) {
        enriched_path <- enriched_path[seq_len(n_pathway), ]
    } else if (nrow(enriched_path) == 0 || is.null(enriched_path)) {
        warning("No enriched items with padj value smaller than", penrichCutoff, ". Plot the top ", n_pathway, " with smallest padj values instead.\n")
        enriched_path <- result$result[[method]] %>%
            arrange(padj) %>%
            slice_head(n = n_pathway)
    }
    all_paths <- unique(enriched_path$pathway)
    pathway <- pathway %>% filter(pathway %in% all_paths)

    data <- input_data$data$data

    path_MGI <- table(factor(pathway$pathway, levels = all_paths))

    pathway_miRs <- unique(pathway$miRNA)
    pathway_TGs <- unique(pathway$gene)

    DEMGI <- data %>%
        filter(cor < corCutoff, corPadj < corPadjCutoff, abs(mirLogFC) > DElogFCCutoff, mirPadj < DEpadjCutoff, abs(geneLogFC) > DElogFCCutoff, genePadj < DEpadjCutoff) %>%
        mutate(MGI = paste0(miRNA, ":", gene))
    DEG <- data %>%
        filter(abs(geneLogFC) > DElogFCCutoff, genePadj < DEpadjCutoff) %>%
        pull(gene)
    MGIList <- data %>%
        arrange(strength) %>%
        select(miRNA, gene, strength) %>%
        filter(miRNA %in% pathway_miRs, gene %in% pathway_TGs)

    DEMGI_path <- pathway %>%
        inner_join(DEMGI, by = c("miRNA", "gene")) %>%
        filter(pathway %in% all_paths)
    path_DEMGI <- table(factor(DEMGI_path$pathway, levels = all_paths))

    ## 1. Annotation -----
    #### 1.1 pathway-level annotation: path_order, path_MGI, path_DEMGI -----
    n_total <- length(unique(MGIList$miRNA))
    n_keep <- min(n_mir, n_total)
    DEmiRs <- data %>%
        filter(abs(mirLogFC) > DElogFCCutoff, mirPadj < DEpadjCutoff) %>%
        pull(miRNA) %>%
        unique()

    top_miRNA <- MGIList %>%
        mutate(MGI = paste0(miRNA, ":", gene)) %>%
        mutate(is_DE = ifelse(MGI %in% unique(DEMGI$MGI), 1, 0)) %>%
        group_by(miRNA) %>%
        summarise(n_DE = sum(is_DE)) %>%
        ungroup() %>%
        arrange(desc(n_DE)) %>%
        slice_head(n = n_keep) %>%
        pull(miRNA)
    # Among the MGIs formed by genes/miRNAs within the current pathway, take the top n_keep
    # by number of DE MGIs. I.e. the gene must be in the current all_paths; the miRNA need not
    # be DE, but is only drawn when there are fewer than n_keep DE miRNAs (meaning no DEMGI
    # exists for it, since a DEMGI requires a DE miRNA).

    MGIList <- MGIList %>% filter(miRNA %in% top_miRNA)

    MGIList <- MGIList %>%
        group_by(miRNA) %>%
        arrange(strength) %>%
        mutate(gene_count = n_distinct(gene)) %>%
        mutate(row_number = row_number()) %>%
        filter(gene_count <= 10 | row_number <= 10) %>%
        ungroup()
    # if a miRNA regulates more than 10 genes, keep only the 10 with the smallest strength
    unique_gene <- unique(MGIList$gene)

    #### 1.4 gene bottom annotation: gene_anno dataframe ----
    gene_anno <- data %>%
        select(gene, geneLogFC, genePadj) %>%
        distinct() %>%
        filter(gene %in% unique_gene) %>%
        arrange(geneLogFC)
    gene_order <- gene_anno$gene

    # A. build miRNA-gene matrix: miR_mat, miR_DElabel_mat ----
    miR_mat <- matrix(NA,
        nrow = length(top_miRNA), ncol = length(gene_order),
        dimnames = list(top_miRNA, gene_order)
    )
    for (i in seq_len(nrow(MGIList))) {
        miR_mat[MGIList$miRNA[i], MGIList$gene[i]] <- MGIList$strength[i]
    }

    miR_DElabel_mat <- matrix(NA,
        nrow = length(top_miRNA), ncol = length(gene_order),
        dimnames = list(top_miRNA, gene_order)
    )
    DEMGIList <- MGIList %>%
        mutate(MGI = paste0(miRNA, ":", gene)) %>%
        filter(MGI %in% DEMGI$MGI)
    for (i in seq_len(nrow(DEMGIList))) {
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
    if (sign(miR_min) != sign(miR_max)) {
        mir_col <- colorRamp2(
            breaks = c(-v, 0, v),
            colors = c("#ff6700", "grey80", "#238B45") # lightgreen
        )
    } else if (miR_max < 0) {
        mir_col <- colorRamp2(
            breaks = c(miR_min, 0),
            colors = c("#ff6700", "#FFD7B5")
        )
    } else {
        mir_col <- colorRamp2(
            breaks = c(0, miR_max),
            colors = c("#C7E9C0", "#238B45")
        )
    }

    #### 1.2 miRNA-level annotation: mir_target, mir_DEtarget ----
    bg <- input_data$background_MGI
    mir_target <- table(factor(bg$miRNA, levels = top_miRNA))

    bg_DEG <- bg %>%
        filter(gene %in% DEG) %>%
        distinct()
    mir_DEtarget <- table(factor(bg_DEG$miRNA, levels = top_miRNA))

    bg_DEMGI <- bg %>%
        mutate(MGI = paste0(miRNA, ":", gene)) %>%
        filter(MGI %in% unique(DEMGI$MGI))
    mir_DEMGI <- table(factor(bg_DEMGI$miRNA, levels = top_miRNA))

    DEMGI_path <- DEMGI_path %>%
        select(miRNA, gene) %>%
        unique()
    mir_DEMGI_path <- table(factor(DEMGI_path$miRNA, levels = top_miRNA))

    #### 1.3 miRNA right annotation: mir_anno dataframe ----
    mir_anno <- data %>%
        select(miRNA, mirLogFC, mirPadj) %>%
        distinct() %>%
        filter(miRNA %in% top_miRNA)
    mir_anno$miRNA <- factor(mir_anno$miRNA, levels = top_miRNA)
    mir_anno <- mir_anno[order(mir_anno$miRNA), ]
    mir_max <- round_up_to(max(abs(mir_anno$mirLogFC)))

    #### 1.4 gene bottom annotation: gene_max for color scale ----
    gene_max <- round_up_to(max(abs(gene_anno$geneLogFC)))

    #### 1.5 pathway method padj right annotation: method_anno list ----
    padj_df <- result$result[[method]][, c("pathway", "padj")] %>% filter(pathway %in% all_paths)
    padj_df$pathway <- factor(padj_df$pathway, levels = all_paths)
    method_anno <- na.omit(padj_df[order(padj_df$pathway), ])

    #### 1.6 miRNA left cancer-lelated annotaion ----
    if (!is.null(mir_annot)) {
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
                paste("(", as.integer(mir_DEMGI_path), "/", as.integer(mir_DEMGI), "/",
                    as.integer(mir_DEtarget), "/", as.integer(mir_target), ")",
                    sep = ""
                ),
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
                paste("(", as.integer(mir_DEMGI_path), "/", as.integer(mir_DEMGI), "/",
                    as.integer(mir_DEtarget), "/", as.integer(mir_target), ")",
                    sep = ""
                ),
                gp = gpar(fontsize = 12)
            )
        )
    }

    # B. build pathway-gene matrix: path_mat, path_DElabel_mat ----
    path_mat <- matrix(0,
        nrow = length(all_paths), ncol = length(gene_order),
        dimnames = list(all_paths, gene_order)
    )
    pathway_filter <- pathway %>% filter(pathway %in% all_paths, gene %in% gene_order)
    for (i in seq_len(nrow(pathway_filter))) {
        pw <- pathway_filter$pathway[i]
        gene <- pathway_filter$gene[i]
        if (pw %in% all_paths && gene %in% gene_order) {
            path_mat[pw, gene] <- 1
        }
    }

    #### 1.7 gene upper annotation:gene_annot (whether includes in cancer-related genes)----
    if (!is.null(gene_annot)) {
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
        col = colorRamp2(c(0, penrichCutoff, 1), c("#ff6666", "mistyrose", "#f7f7f7")),
        which = "row" # , width = unit(1, "cm")
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
    col_mir_fc <- colorRamp2(c(-mir_max, 0, mir_max), c("#add8e6", "white", "#ffb6c1"))
    range_mir_fc <- seq(-mir_max, mir_max, length = 5)
    lgd_mir_fc <- Legend(title = "miRNA Log2FC", col_fun = col_mir_fc, at = range_mir_fc, labels = paste(range_mir_fc))
    ## miRNA padj
    col_mir_padj <- colorRamp2(c(0, 0.05, 1), c("#ff6666", "mistyrose", "#f7f7f7"))
    lgd_mir_padj <- Legend(title = "miRNA Padj", col_fun = col_mir_padj, at = c(0, 0.05, 1))
    ## method padj
    lgd_method <- Legend(title = paste0(method, " Padj"), col_fun = colorRamp2(c(0, penrichCutoff, 1), c("#ff6666", "mistyrose", "grey90")), at = c(0, penrichCutoff, 1))
    ## gene log2FoldChange
    lgd_gene_fc_cat <- Legend(title = "Gene Expression \n(* shows DEG)", legend_gp = gpar(fill = c("blue", "red")), labels = c("downregulated", "upregulated"))
    col_gene_fc <- colorRamp2(c(-gene_max, 0, gene_max), c("#add8e6", "white", "#ffb6c1"))
    range_gene_fc <- seq(-gene_max, gene_max, length = 5)
    lgd_gene_fc <- Legend(title = "Gene Log2FC", col_fun = col_gene_fc, at = range_gene_fc, labels = paste(range_gene_fc))

    ## DE rectangle
    lgd_DE <- Legend(
        title = "MGI Expression", labels = c("DEMGI"), type = c("points"), pch = c(22),
        legend_gp = gpar(fill = c("white"), col = c("#b266cc"), lwd = c(2)), size = unit(5, "mm")
    )

    # cancer-related gene; legend title derived from the annotation list's own names
    # rather than a hard-coded "Cancer-related ..." string, so the label reflects
    # whatever gene_annot/mir_annot the caller actually supplied.
    gene_top_lgd <- if (!is.null(gene_annot)) {
        Legend(
            title = paste(names(gene_annot), collapse = " / "),
            legend_gp = gpar(fill = c(annot_color, "white")), labels = c("TRUE", "FALSE")
        )
    } else {
        NULL
    }

    mir_left_lgd <- if (!is.null(mir_annot)) {
        Legend(
            title = paste(names(mir_annot), collapse = " / "),
            legend_gp = gpar(fill = c(annot_color, "white")), labels = c("TRUE", "FALSE")
        )
    } else {
        NULL
    }

    legend_list <- c(gene_top_lgd, lgd_method, lgd_gene_fc_cat, mir_left_lgd, lgd_DE, lgd_gene_fc, lgd_mir_fc, lgd_mir_padj)

    ## 3. Heatmap dimensions ----
    ht_unit <- min(c(0.0264 * (800 - 50) / dim(path_mat)[1], 0.5))
    wd_unit <- min(c(0.0264 * (1700 - 100) / dim(path_mat)[2], 0.5))

    h1_ht <- ceiling(ht_unit * dim(path_mat)[1])
    h1_wd <- ceiling(wd_unit * dim(path_mat)[2])
    h2_ht <- ceiling(ht_unit * dim(miR_mat)[1])
    h2_wd <- ceiling(wd_unit * dim(miR_mat)[2])

    list(
        result = result, data = data, pathway = pathway, all_paths = all_paths, enriched_path = enriched_path,
        top_miRNA = top_miRNA, gene_order = gene_order,
        path_mat = path_mat, path_MGI = path_MGI, path_DEMGI = path_DEMGI,
        miR_mat = miR_mat, miR_DElabel_mat = miR_DElabel_mat, mir_col = mir_col, v = v,
        mir_anno = mir_anno, mir_max = mir_max, gene_anno = gene_anno, gene_max = gene_max,
        left_mir_annot = left_mir_annot, top_gene_annot = top_gene_annot, row_anno = row_anno,
        col_mir_fc = col_mir_fc, col_mir_padj = col_mir_padj,
        legend_list = legend_list,
        h1_ht = h1_ht, h1_wd = h1_wd, h2_ht = h2_ht, h2_wd = h2_wd,
        DElogFCCutoff = DElogFCCutoff, DEpadjCutoff = DEpadjCutoff,
        corCutoff = corCutoff, corPadjCutoff = corPadjCutoff
    )
}

#' Visualize an enrichment result as a heatmap
#'
#' Draws a heatmap of the miRNA-gene-pathway enrichment result for one
#' edge-based method, showing expression fold-changes and significance.
#' Only supports one method each time.
#'
#' @param method a character that specifies only one method which you have the
#'        corresponding enrichment result (only five edge-based methods are
#'        supported).
#' @param result the output of miREA() function, used here as the unprocessed
#'        result list.
#' @param input_data the output of getAllInputData() function, used here as the
#'        unprocessed input_data list.
#' @param n_mir_heatmap The number of miRNAs included in the heatmap. Default is
#'        10.
#' @param n_pathway The number of top enriched pathways included in the heatmap.
#'        Default it 10. If the number of enriched pathways is less than
#'        n_pathway, all enriched pathways will be included.
#' @param plot_path The directory to save the plot. If not null, it will
#'        generate a pdf file named heatmap_[method].pdf at the plot path. If
#'        null, the plot will be drawn at the current window.
#' @param penrichCutoff  The threshold for deciding enrichment of pathways.
#'        Default is 0.05.
#' @param height The height of the plot, in cm. If not, it will be set to
#'        default values.
#' @param width The width of the plot, in cm. If not, it will be set to default
#'        values.
#' @param gene_annot A list contains the annotations of genes, where each
#'        sub-element denotes a character vector containing a list of genes. It
#'        will generate a annotation bar on the top of the heatmap to show
#'        whether the genes hit in the annotation list. Default is NULL, which
#'        means no annotation bar will be generated.
#' @param mir_annot Same format as gene_annot. If not NULL, it will generate an
#'        annotation bar for miRNAs at the left side of the heatmap, to show
#'        whether miRNAs hit in the given list.
#' @param annot_color The color of the annotation bars for genes and miRNAs
#' @param DElogFCCutoff optional, the absolute log2FC threshold used to call a
#'        miRNA/gene differentially
#'        expressed for the heatmap annotation. If NULL (the default), falls
#'        back to
#'        `input_data$threshold$DElogFCCutoff` (set by [getAllInputData()]), or
#'        1 if that is absent.
#' @param DEpadjCutoff optional, the adjusted p-value threshold used to call a
#'        miRNA/gene differentially
#'        expressed for the heatmap annotation. If NULL (the default), falls
#'        back to
#'        `input_data$threshold$DEpadjCutoff`, or 0.05 if that is absent.
#' @param corCutoff optional, the miRNA-gene expression correlation threshold
#'        for the heatmap annotation.
#'        If NULL (the default), falls back to `input_data$threshold$corCutoff`,
#'        or 0 if that is absent.
#' @param corPadjCutoff optional, the adjusted p-value threshold for the
#'        miRNA-gene expression correlation,
#'        for the heatmap annotation. If NULL (the default), falls back to
#'        `input_data$threshold$corPadjCutoff`, or 0.05 if that is absent.

#' @return Invisibly returns `list(ht_list, legend_list)` (the drawn
#'        [ComplexHeatmap::HeatmapList-class] and its annotation legends),
#'        whether or not
#'        `plot_path` is supplied. If `plot_path` is supplied, the heatmap is
#'        additionally
#'        saved as `heatmap_[method].pdf` at that path and a message reports the
#'        file location.
#' @examples
#' data(toy_pathway)
#' data(toy_gene_DEdata)
#' data(toy_mir_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
#' data(background_MGI)
#' # Edge_Score2D is included so input_data carries the MGI "strength" scores
#' # plotHeatmap() needs for its heatmap cells, without running it in miREA().
#' input_data <- getAllInputData(
#'     methods = c("Edge_ORA", "Edge_Score2D"), pathway = toy_pathway,
#'     mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
#'     background_MGI = background_MGI,
#'     gene_mat = toy_gene_mat, mir_mat = toy_mir_mat, scoreFun = "rank"
#' )
#' background <- list(
#'     Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))
#' )
#' result <- miREA(
#'     methods = "Edge_ORA", input_data = input_data, background = background,
#'     minSize = NULL, maxSize = NULL, pvalueCutoff = 1
#' )
#' plotHeatmap(method = "Edge_ORA", result = result, input_data = input_data,
#' penrichCutoff = 1)
#' @export
plotHeatmap <- function(method, result, input_data, n_mir_heatmap = 10, n_pathway = 10,
                        plot_path = NULL, penrichCutoff = 0.05, height = NULL, width = NULL,
                        gene_annot = NULL, mir_annot = NULL, annot_color = "#A6DDBA",
                        DElogFCCutoff = NULL, DEpadjCutoff = NULL,
                        corCutoff = NULL, corPadjCutoff = NULL) {
    comp <- build_heatmap_components(
        method, result, input_data, n_mir_heatmap, n_pathway,
        penrichCutoff, gene_annot, mir_annot, annot_color,
        DElogFCCutoff, DEpadjCutoff, corCutoff, corPadjCutoff
    )

    ht1 <- ComplexHeatmap::Heatmap(
        comp$path_mat,
        name = "ht1",
        col = c("1" = "grey50", "0" = "white"),
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
        width = unit(comp$h1_wd, "cm"),
        height = unit(comp$h1_ht, "cm"),
        column_gap = unit(c(2), "mm"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        left_annotation = rowAnnotation(totgene_text = anno_text(
            paste("(", as.integer(comp$path_DEMGI), "/", as.integer(comp$path_MGI), ")", sep = ""),
            gp = gpar(fontsize = 12)
        )),
        right_annotation = comp$row_anno,
        top_annotation = comp$top_gene_annot,
        bottom_annotation = columnAnnotation(
            geneLog2FC = anno_barplot(
                comp$gene_anno$geneLogFC,
                border = FALSE,
                baseline = 0,
                ylim = c(-comp$gene_max, comp$gene_max),
                gp = gpar(
                    fontsize = 12,
                    col = ifelse(comp$gene_anno$geneLogFC < 0, "blue", "red"),
                    fill = ifelse(comp$gene_anno$geneLogFC < 0, "blue", "red")
                ),
                height = unit(1.5, "cm")
            ),
            star = anno_text(
                ifelse(abs(comp$gene_anno$geneLogFC) > comp$DElogFCCutoff, "*", ""),
                gp = gpar(fontsize = 18, col = ifelse(comp$gene_anno$geneLogFC < 0, "blue", "red")),
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

    ht2 <- ComplexHeatmap::Heatmap(
        comp$miR_mat,
        name = "ht2",
        na_col = "white",
        col = comp$mir_col,
        column_names_rot = 90,
        column_names_gp = gpar(fontsize = 12),
        show_row_names = TRUE,
        row_names_side = "left",
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_gp = gpar(fontsize = 12),
        width = unit(comp$h2_wd, "cm"),
        height = unit(comp$h2_ht, "cm"),
        column_gap = unit(c(2), "mm"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(comp$miR_DElabel_mat[i, j])) {
                grid.rect(x, y,
                    width = width, height = height,
                    gp = gpar(lwd = 3, col = "#B266CC", fill = fill)
                )
            }
        },
        heatmap_legend_param = list(
            direction = "vertical",
            title = "MGI Strength",
            col = comp$mir_col,
            na_col = "white"
        ),
        left_annotation = comp$left_mir_annot,
        right_annotation = rowAnnotation(
            mir_fc = comp$mir_anno$mirLogFC,
            mir_padj = comp$mir_anno$mirPadj,
            col = list(
                mir_fc = comp$col_mir_fc,
                mir_padj = comp$col_mir_padj
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
    legend_list <- comp$legend_list

    if (!is.null(plot_path)) {
        if (is.null(height)) {
            height <- 2 * comp$h1_ht + 2 * comp$h2_ht
        }
        if (is.null(width)) {
            width <- comp$h1_wd * 1.2
        }
        pdf_path <- file.path(plot_path, paste0("heatmap_", method, ".pdf"))
        # height/width above are in cm; pdf() expects inches.
        pdf(pdf_path, height = height / 2.54, width = width / 2.54)
        on.exit(dev.off(), add = TRUE)
        draw(ht_list,
            merge_legend = TRUE,
            annotation_legend_list = legend_list,
            column_title = paste0("Heatmap for ", method, " method"),
            column_title_gp = gpar(fontsize = 20, fontface = "bold"),
            heatmap_legend_side = "bottom",
            annotation_legend_side = "bottom",
            align_heatmap_legend = "heatmap_center"
        )
        decorate_heatmap_body("ht1", {
            grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
        })

        decorate_heatmap_body("ht2", {
            grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
        })
        message("Heatmap for ", method, " saved at: ", pdf_path)
    } else {
        draw(ht_list,
            merge_legend = TRUE,
            annotation_legend_list = legend_list,
            column_title = paste0("Heatmap for ", method, " method"),
            column_title_gp = gpar(fontsize = 20, fontface = "bold"),
            heatmap_legend_side = "right",
            annotation_legend_side = "right"
        )
        decorate_heatmap_body("ht1", {
            grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
        })

        decorate_heatmap_body("ht2", {
            grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
        })
        message("  Visualization finished!\n***If you want to save them to pdf files, please specify the [plot_path] parameter!***")
    }

    invisible(list(ht_list = ht_list, legend_list = legend_list))
}
