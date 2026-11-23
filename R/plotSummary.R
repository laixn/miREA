# result: the output of miREA
# plot_path: file directory
# pvalueCutoff: the threshold to distinguish enriched or not. If null, don't filter padj.
# output: 2 pdf files, summary.pdf and pvalue_density.pdf

#' Visualize and compare enrichment results across methods
#'
#' Draws summary plots comparing enrichment results across the different
#' miREA methods: a positive-rate/PCA/pairwise-overlap panel, and a p-value
#' density plot.
#'
#' @param result the output of miREA() function, used here as the unprocessed result list.
#' @param plot_path The directory to save the plot. If not null, it will generate two pdf file named summary.pdf and pvalue_density.pdf at the plot path. If null, the plot will be drawn at the current window.
#' @param penrichCutoff  The threshold for deciding enrichment of pathways. Default is 0.05.
#' @param fill_col a named vector to specify fill colors for different methods.
#' @return Two plots, either in the current window if plot_path is NULL, or saved as summary.pdf and pvalue_density.pdf at the plot_path.
#' @examples
#' data(toy_pathway)
#' data(toy_gene_DEdata)
#' data(toy_mir_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
#' data(background_MGI)
#' methods <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score", "Edge_ORA")
#' input_data <- getAllInputData(
#'   methods = methods, pathway = toy_pathway,
#'   mir_DEdata = toy_mir_DEdata, gene_DEdata = toy_gene_DEdata,
#'   background_MGI = background_MGI,
#'   gene_mat = toy_gene_mat, mir_mat = toy_mir_mat
#' )
#' background <- list(
#'   TG = unique(toy_pathway$gene),
#'   MiR = unique(input_data$pathway$MiR$miRNA),
#'   Edge = unique(paste0(background_MGI$miRNA, ":", background_MGI$gene))
#' )
#' result <- miREA(
#'   methods = methods, input_data = input_data, background = background,
#'   minSize = NULL, maxSize = NULL,
#'   pvalueType = list(TG = "std", MiR = "std", Edge = "neg"),
#'   pvalueCutoff = 1
#' )
#' plotSummary(result, penrichCutoff = 1)
#' @export
plotSummary <- function(result, penrichCutoff = NULL, plot_path = NULL,
                        fill_col = c(
                          "TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3", "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
                          "Edge_ORA" = "#EDC194", "Edge_Score1D" = "#F09137", "Edge_Score2D" = "#626FB3",
                          "Edge_Topology" = "#EAA5C2", "Edge_Network" = "#FA8072"
                        )) {
  enrich <- list()
  all_pathways <- unique(unlist(lapply(result$result, function(df) df$pathway)))
  n_pathway <- length(all_pathways)
  summary <- data.frame(method = character(), n_pathway = integer(), positive_rate = numeric(), stringsAsFactors = FALSE)
  if (is.null(penrichCutoff)) {
    penrichCutoff <- 1.01
  }
  for (method in names(result$result)) {
    enriched_path <- result$result[[method]] %>%
      filter(padj < penrichCutoff) %>%
      select(pathway) %>%
      pull()
    enrich[[method]] <- enriched_path
    summary <- rbind(summary, data.frame(method = method, n_pathway = n_pathway, positive_rate = length(enriched_path) / n_pathway))
  }
  result <- result$result

  # summary.pdf
  rate_plot <- plot_positive_rate(summary, fill_col)
  pca_plot <- plot_pca(result, fill_col)
  overlap_plot <- plot_overlap_ht(enrich)

  combined <- (rate_plot + overlap_plot + pca_plot) +
    plot_layout(widths = c(1, 1, 1), heights = 1) & # , guides = "collect"
    theme(plot.margin = margin(2, 2, 2, 2))


  # pvalue_density.pdf
  pvalue_plot <- plot_pvalue_dens(result, fill_col) # patchwork ggplot combine

  if (!is.null(plot_path)) {
    ggsave(paste0(plot_path, "/summary.pdf"), combined, width = 16, height = 5)
    ggsave(paste0(plot_path, "/pvalue_density.pdf"), pvalue_plot, width = 10, height = 10)
    message("    Please check the summary plot and pvalue density plot at: ", plot_path)
    return(paste0("Please check the summary plot and pvalue density plot at: ", plot_path, "\n"))
  } else {
    plot(combined)
    plot(pvalue_plot)
    message("  Visualization finished!\n***If you want to save them to pdf files directly, please specify the [plot_path] parameter!***")
    return(list(summary = combined, pvalue = pvalue_plot))
  }
}


plot_positive_rate <- function(summary, fill_col) { # plot_path
  colnames(summary)[seq_len(3)] <- c("method", "n_pathway", "positive_rate")
  methods <- summary$method
  n_path <- summary$n_pathway[1]
  x_range <- range(as.numeric(factor(summary$method)))
  y_range <- range(summary$positive_rate, na.rm = TRUE)
  # diff(y_range) is 0 whenever every method has the same positive rate (e.g.
  # all-enriched or all-unenriched), which would otherwise make ratio Inf/NaN
  # and coord_fixed() error out.
  ratio <- if (diff(y_range) == 0) 1 else diff(x_range) / diff(y_range)

  summary$method <- factor(summary$method, levels = unique(summary$method))
  fill_col <- fill_col[names(fill_col) %in% methods]

  p <- ggplot(summary, aes(x = method, y = positive_rate, fill = method)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = fill_col, guide = "none") +
    labs(
      x = NULL, y = paste0("Positive Rate (", n_path, " pathways)"),
      title = "Positive Rates for Different Methods"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      axis.title.x = element_text(color = "black", size = 10),
      axis.title.y = element_text(color = "black", size = 10),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks = element_line(color = "black"),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank()
    ) +
    coord_fixed(ratio = ratio, clip = "off") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_range[2] * 1.05))

  return(p)
}


plot_pca <- function(result, fill_col) {
  all_pathways <- unique(unlist(lapply(result, function(df) df$pathway)))
  methods <- names(result)
  pval_matrix <- matrix(1,
    nrow = length(all_pathways), ncol = length(result),
    dimnames = list(all_pathways, methods)
  )
  for (method in methods) {
    sub_matrix <- result[[method]]
    pval_matrix[sub_matrix$pathway, method] <- sub_matrix$p_value
  }

  matrix <- t(pval_matrix)
  matrix <- matrix[, colSums(is.na(matrix)) == 0] # drop columns with NA
  matrix <- matrix[, which(apply(matrix, 2, var) != 0)] # remove zero variance columns
  fill_col <- fill_col[names(fill_col) %in% methods]

  pca_results <- prcomp(matrix, center = TRUE, scale. = TRUE)


  pca_scores <- as_tibble(pca_results$x) # Keep patient IDs

  pca_summary <- summary(pca_results)

  variance_explained <- pca_summary$importance[2, ] * 100 # Importance[2,] contains the proportion of variance
  pc1_variance <- round(variance_explained[1], 1)
  pc2_variance <- round(variance_explained[2], 1)

  pca_scores$methods <- factor(methods) # , levels = methods
  pca_scores$shape_group <- ifelse(grepl("Edge_", pca_scores$methods), "Edge", "Other")
  shape_values <- ifelse(grepl("Edge_", methods), 16, 17)
  names(shape_values) <- methods

  pca_plot <- ggplot(
    pca_scores,
    aes(
      x = PC1, y = PC2,
      color = factor(methods, levels = methods),
      shape = factor(methods, levels = methods)
    )
  ) +
    geom_point(size = 3, alpha = 1, stroke = 1.5) + # aes(shape = grade)
    scale_color_manual(values = fill_col) +
    scale_shape_manual(values = shape_values) +
    labs(
      title = "PCA Based on P-values for different methods",
      x = paste0("PC1 (", pc1_variance, "% variance)"),
      y = paste0("PC2 (", pc2_variance, "% variance)"),
      colour = "Methods",
      shape = "Methods"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.caption = element_text(hjust = 0),
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_text(color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      axis.title.x = element_text(color = "black", size = 10),
      axis.title.y = element_text(color = "black", size = 10),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks = element_line(color = "black"),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7)
    ) +
    coord_fixed()

  return(pca_plot)
}

plot_overlap_ht <- function(enrich) {
  enrich <- Filter(function(x) length(x) > 0, enrich)
  if (length(enrich) < 2) {
    warning("The number of methods with at least one enriched pathways are below 2. Skip the overlap plot.\n")
    return(NULL)
  }
  n <- length(enrich)
  matrix <- matrix(0, n, n, dimnames = list(names(enrich), names(enrich)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      intersection <- length(intersect(enrich[[i]], enrich[[j]]))
      if (i >= j) { # lower plot: overlap
        matrix[i, j] <- intersection
      } else if (i < j) {
        # upper plot: jaccard
        union <- length(base::union(enrich[[i]], enrich[[j]]))
        if (union == 0) {
          matrix[i, j] <- 0
        } else {
          matrix[i, j] <- intersection / union
        }
      }
    }
  }
  df <- reshape2::melt(matrix, varnames = c("Method1", "Method2"), value.name = "Value")
  df$Value <- as.numeric(df$Value)

  df$Type <- ifelse(as.numeric(df$Method1) < as.numeric(df$Method2), "Upper",
    ifelse(df$Method1 == df$Method2, "Diag", "Lower")
  )
  df <- df %>% mutate(text = ifelse(Type %in% c("Diag", "Lower"), as.integer(Value), NA))

  df$Method2 <- factor(df$Method2, levels = rev(unique(df$Method2)))

  overlap <- ggplot(df, aes(x = Method1, y = Method2, label = text)) +
    # upper plot: jaccard index
    geom_point(
      data = subset(df, Type %in% "Upper"),
      aes(size = Value),
      shape = 21, fill = "grey60", color = "black"
    ) +
    # lower part: overlap count
    geom_point(
      data = subset(df, Type %in% c("Lower", "Diag")),
      aes(color = Value),
      shape = 15, size = 8
    ) +
    geom_tile(color = "black", linewidth = 0.5, linetype = 1, fill = NA) +
    {
      # When every method's overlap count is identical (e.g. all pathways
      # enriched under a lenient pvalueCutoff), the gradient breakpoints
      # below collapse to a single repeated value and scales::rescale()
      # can't build an interpolatable color ramp from them. Fall back to a
      # plain two-color gradient in that case.
      lower_diag_vals <- unique(subset(df, Type %in% c("Lower", "Diag"))$Value)
      if (length(lower_diag_vals) > 1) {
        scale_color_gradientn(
          colours = c("#F7FCF0", "#CCE8C5", "#A8DDB5", "#4EB3D3", "#084081"),
          values = scales::rescale(c(
            min(lower_diag_vals), quantile(lower_diag_vals, 0.25),
            median(lower_diag_vals), quantile(lower_diag_vals, 0.75),
            max(lower_diag_vals)
          ))
        )
      } else {
        scale_color_gradientn(colours = c("#F7FCF0", "#084081"))
      }
    } +
    geom_text(
      data = subset(df, !is.na(text)),
      aes(label = text),
      size = 0.5 * length(enrich)
    ) +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Enrichment Overlap for Different Methods", x = NULL, y = NULL, color = "Overlap Count", size = "Jaccard Index") +
    theme(
      axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      panel.grid = element_blank(),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size  = guide_legend(order = 2)
    )
  return(overlap)
}


plot_pvalue_dens <- function(result, fill_col) {
  plist <- list()
  for (i in names(result)) {
    plist[[i]] <- result[[i]]$p_value
  }
  methods <- names(result)
  fill_col <- fill_col[names(fill_col) %in% methods]

  plot_list <- lapply(names(plist), function(method) {
    df <- data.frame(pvalue = plist[[method]])
    ggplot(df, aes(x = pvalue)) +
      geom_density(aes(y = after_stat(density / max(density))),
        fill = fill_col[method], alpha = 0.6
      ) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = label_number(accuracy = 0.01)) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        labels = scales::percent_format(accuracy = 1) # accuracy: only save integer part
      ) +
      ggtitle(method) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
      ) +
      labs(x = "pvalue", y = "density") +
      guides(fill = "none")
  })
  names(plot_list) <- names(plist)
  p_all <- wrap_plots(plot_list, nrow = 3, guides = "collect") +
    plot_annotation(
      title = "Pvalue Density Distribution for Different Methods",
      subtitle = NULL,
      caption = NULL,
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)
      )
    ) &
    labs(x = "pvalue", y = "density")
  return(p_all)
}
