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
#' @param height The height of the plot, in cm. If not, it will be set to default values.
#' @param width The width of the plot, in cm. If not, it will be set to default values.
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

#' @return Invisibly returns `list(ht_list, legend_list)` (the drawn
#'        [ComplexHeatmap::HeatmapList-class] and its annotation legends), whether or not
#'        `plot_path` is supplied. If `plot_path` is supplied, the plot is additionally saved as
#'        `heatmap_sankey_[method].pdf` at that path and a message reports the file location.
#' @examples
#' data(toy_pathway)
#' data(toy_gene_DEdata)
#' data(toy_mir_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
#' data(background_MGI)
#' # Edge_Score2D is included so input_data carries the MGI "strength" scores
#' # plotHeatmapSankey() needs for its heatmap cells, without running it in miREA().
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
#' plotHeatmapSankey(method = "Edge_ORA", result = result, input_data = input_data, penrichCutoff = 1)
#' @export
plotHeatmapSankey <- function(method, result, input_data, n_mir_heatmap = 20, n_mir_sankey = 10, n_pathway = 10,
                              plot_path = NULL, penrichCutoff = 0.05, height = NULL, width = NULL, sankey_prop = 1.5,
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
        column_title = paste0("Heatmap and Sankey plot for ", method, " method"),
        column_title_gp = gpar(fontsize = 20, fontface = "bold", hjust = 0.5),
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

    ##### sankey plot ----
    sankey <- plot_sankey(
        dict = comp$data, pathway = comp$pathway, enriched_path = comp$enriched_path,
        all_paths = comp$all_paths, top_miRNA = comp$top_miRNA, mir_max = comp$mir_max,
        gene_max = comp$gene_max, strength_max = comp$v, penrichCutoff = penrichCutoff,
        n_mir = n_mir_sankey, n_pathway = n_pathway,
        DElogFCCutoff = comp$DElogFCCutoff, DEpadjCutoff = comp$DEpadjCutoff,
        corCutoff = comp$corCutoff, corPadjCutoff = comp$corPadjCutoff
    )
    sankey <- ggplotGrob(sankey)
    anno_ggplot <- AnnotationFunction(
        fun = function(index, k, n) {
            pushViewport(viewport())
            grid.draw(sankey)
            popViewport()
        },
        height = unit((comp$h2_ht + comp$h1_ht) * sankey_prop, "cm"),
        subsettable = FALSE,
        which = "column"
    )

    sankey_anno <- HeatmapAnnotation(
        PR = anno_ggplot,
        show_annotation_name = FALSE
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
        bottom_annotation = sankey_anno,
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
        pdf_path <- file.path(plot_path, paste0("heatmap_sankey_", method, ".pdf"))
        # height/width above are in cm; pdf() expects inches.
        pdf(pdf_path, height = height / 2.54, width = width / 2.54)
        on.exit(dev.off(), add = TRUE)
        draw(ht_list,
            merge_legend = TRUE,
            annotation_legend_list = legend_list,
            heatmap_legend_side = "right",
            annotation_legend_side = "right"
        )
        decorate_heatmap_body("ht1", {
            grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
        })

        decorate_heatmap_body("ht2", {
            grid.rect(gp = gpar(col = "black", lwd = 3, fill = NA))
        })
        message("Heatmap & Sankey plot for ", method, " saved at: ", pdf_path)
    } else {
        draw(ht_list,
            merge_legend = TRUE,
            annotation_legend_list = legend_list,
            column_title = paste0("Heatmap and Sankey plot for ", method, " method"),
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
        if (is.na(s)) {
            return(NA_character_)
        }
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

plot_sankey <- function(dict, pathway, enriched_path, all_paths, top_miRNA, mir_max, gene_max, strength_max, penrichCutoff = 0.05, n_mir = 25, n_pathway = 10,
                        DElogFCCutoff = 1, DEpadjCutoff = 0.05, corCutoff = 0, corPadjCutoff = 0.05) {
    pathway <- pathway %>% filter(pathway %in% all_paths)
    dict_DE <- dict %>%
        filter(cor < corCutoff, corPadj < corPadjCutoff, genePadj < DEpadjCutoff, abs(geneLogFC) >= DElogFCCutoff, mirPadj < DEpadjCutoff, abs(mirLogFC) >= DElogFCCutoff) %>%
        arrange(desc(abs(strength)))
    path_DEMGI <- pathway %>%
        inner_join(dict_DE, by = c("miRNA", "gene")) %>%
        filter(pathway %in% all_paths)

    pw_MGI <- unique(paste0(pathway$miRNA, ":", pathway$gene))
    dict_DEMGI <- dict_DE %>%
        filter(paste0(miRNA, ":", gene) %in% pw_MGI)
    n_total <- length(top_miRNA)
    n_keep <- min(n_mir, n_total)
    unique_mir <- top_miRNA[seq_len(n_keep)]
    MGI <- dict_DEMGI %>% filter(miRNA %in% unique_mir)
    # if a miRNA has more than 5 target genes, keep only the 5 with the smallest strength
    MGI <- MGI %>%
        group_by(miRNA) %>%
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
        select(miRNA, gene, flow_value, sign, mirLogFC, geneLogFC, strength) %>%
        distinct()
    df_gene_path <- df_sankey %>%
        select(gene, pathway, geneLogFC, padj) %>%
        distinct()


    # order by log2FC
    gene_order <- df_mir_gene %>%
        arrange(desc(geneLogFC)) %>%
        pull(gene) %>%
        unique()
    mir_order <- df_mir_gene %>%
        arrange(desc(mirLogFC)) %>%
        pull(miRNA) %>%
        unique()
    pw_order <- df_gene_path %>%
        arrange(padj) %>%
        pull(pathway) %>%
        unique()
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
            by = c("stratum" = "gene")
        )
    pathway_info <- df_gene_path_long %>%
        filter(Demographic == "pathway") %>%
        select(alluvium, pathway = stratum)

    gene_info_per_edge <- df_gene_path_long %>%
        filter(Demographic == "gene") %>%
        select(alluvium, gene = stratum, gene_flow_value)

    pathway_flow <- pathway_info %>%
        left_join(gene_info_per_edge, by = "alluvium") %>%
        group_by(pathway) %>%
        summarise(flow_value = sum(gene_flow_value, na.rm = TRUE) / n(), .groups = "drop")

    df_gene_path_long <- df_gene_path_long %>%
        mutate(flow_value = case_when(
            Demographic == "gene" ~ gene_flow_value,
            TRUE ~ NA_real_
        )) %>%
        left_join(pathway_flow, by = c("stratum" = "pathway")) %>%
        mutate(flow_value = if_else(Demographic == "pathway", flow_value.y, flow_value.x)) %>%
        select(-gene_flow_value, -flow_value.y, -flow_value.x)

    miRNA_fill_df <- df_mir_gene_long %>%
        filter(Demographic == "miRNA") %>%
        distinct(stratum, mirLogFC)
    gene_fill_df <- df_mir_gene_long %>%
        filter(Demographic == "gene") %>%
        distinct(stratum, geneLogFC)
    pathway_fill_df <- df_gene_path_long %>%
        filter(Demographic == "pathway") %>%
        distinct(stratum, padj)

    col_mir_fc <- colorRamp2(c(-mir_max, 0, mir_max), c("#add8e6", "white", "#ffb6c1"))
    miRNA_fill_df <- miRNA_fill_df %>%
        mutate(fill_color = col_mir_fc(mirLogFC))

    # gene: geneLogFC mapped light-blue (down) -> white -> light-pink (up)
    col_gene_fc <- colorRamp2(c(-gene_max, 0, gene_max), c("#add8e6", "white", "#ffb6c1"))
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
        levels = c("miRNA", "gene", "pathway")
    )
    df_gene_path_long$Demographic <- factor(df_gene_path_long$Demographic,
        levels = c("miRNA", "gene", "pathway")
    )


    gene_colors <- setNames(rep("#B3D49D90", length(unique(df_mir_gene$gene))), unique(df_mir_gene$gene))
    sankey <- ggplot() +
        geom_flow(
            data = df_mir_gene_long,
            aes(
                x = Demographic,
                stratum = stratum,
                alluvium = alluvium,
                y = flow_value,
                fill = strength
            ),
            stat = ggalluvial::StatAlluvium,
            color = "purple",
            alpha = 0.7,
            aes.flow = "forward"
        ) +
        scale_fill_gradientn(
            colors = c("#ff6700", "grey80", "#238B45"), limits = c(-strength_max, strength_max),
            values = scales::rescale(c(-strength_max, 0, strength_max)), breaks = c(-strength_max, -strength_max / 2, 0, strength_max / 2, strength_max), name = "Strength", guide = "none"
        ) +
        new_scale_fill() +
        geom_flow(
            data = df_gene_path_long,
            aes(
                x = Demographic,
                stratum = stratum,
                alluvium = alluvium,
                fill = stratum,
                y = flow_value
            ),
            color = "black",
            stat = ggalluvial::StatAlluvium,
            alpha = 0.7,
            aes.flow = "forward"
        ) +
        scale_fill_manual(values = gene_colors, guide = "none") +
        new_scale_fill() +
        geom_stratum(
            data = df_mir_gene_long,
            aes(x = Demographic, stratum = stratum, y = flow_value, fill = fill_color)
        ) + #
        geom_stratum(
            data = df_gene_path_long,
            aes(x = Demographic, stratum = stratum, y = flow_value, fill = fill_color)
        ) + #
        scale_x_discrete(limits = c("miRNA", "gene", "pathway"), expand = c(0.02, 0.15), position = "bottom") +
        geom_text(
            data = df_mir_gene_long,
            aes(x = Demographic, stratum = stratum, label = after_stat(stratum), y = flow_value), #
            stat = ggalluvial::StatStratum, size = 3, color = "black"
        ) +
        geom_text(
            data = subset(df_gene_path_long, Demographic == "pathway"),
            aes(x = Demographic, stratum = stratum, y = flow_value, label = wrap_underscored_label(after_stat(stratum), width = 22)),
            stat = ggalluvial::StatStratum, size = 3, color = "black"
        ) +
        scale_fill_identity() +
        theme_void() +
        theme(
            axis.text.x = element_text(size = 16, face = "bold", color = "black"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()
        ) +
        scale_y_continuous(expand = c(0.02, 0))


    return(sankey)
}
