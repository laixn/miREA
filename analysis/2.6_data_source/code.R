# Figure S10: paired TPR bars for each enrichment method.
# Altered-input variants use separate '+' and '-' symbol patterns.
# Original input data remain unchanged; outputs are written to work_dir.

suppressWarnings(suppressMessages({
  library(dplyr)
}))

work_dir <- "C:/Users/rtd648/Desktop/project/miREA/submission/s3_ploscb/revise/FigureS10-2.6_data_source/"
plot_path  <- file.path(work_dir, "FigureS10.pdf")

combined <- read.csv(file.path(work_dir, "TPR_data.csv"), stringsAsFactors = FALSE)

method_order <- c("TG_ORA", "TG_Score", "MiR_ORA", "MiR_Score",
                   "Edge_ORA", "Edge_Score1D", "Edge_Score2D", "Edge_Topology", "Edge_Network")
group_levels <- c("only miRNA data", "miRNA and gene data")

# same base fill color per method family, regardless of input-type group
fill_col_base <- c(
  "TG_ORA" = "#97D7F2", "TG_Score" = "#07AEE3",
  "MiR_ORA" = "#B3D49D", "MiR_Score" = "#35B257",
  "Edge_ORA" = "#EDC194", "Edge_Score1D" = "#F09137",
  "Edge_Score2D" = "#9AA4D6", "Edge_Topology" = "#EAA5C2", "Edge_Network" = "#FA8072"
)

# Two classes of altered-input method variants.
plus_methods <- c(
  "TG_ORA (DEmiR-neg DETG)",
  "MiR_ORA (pathway neg miR set)",
  "MiR_Score (pathway neg miR set)"
)
minus_methods <- c(
  "Edge_ORA (DEmiR-TG)",
  "Edge_Network (DEmiR-TG)"
)

# 1. build the 18-row (9 method x 2 group) plotting table, filling missing combos with 0 ----
combined <- combined %>% mutate(base_method = sub(" \\(.*\\)$", "", method))

full_grid <- expand.grid(base_method = method_order, group = group_levels, stringsAsFactors = FALSE)

plot_df <- full_grid %>%
  left_join(combined %>% select(base_method, group, method, TPR), by = c("base_method", "group")) %>%
  mutate(
    TPR = ifelse(is.na(TPR), 0, TPR),
    pattern = case_when(
      method %in% plus_methods ~ "+",
      method %in% minus_methods ~ "-",
      TRUE ~ "none"
    ),
    base_method = factor(base_method, levels = method_order),
    group = factor(group, levels = group_levels)
  ) %>%
  arrange(base_method, group)

write.csv(plot_df, file = file.path(work_dir, "FigureS10_data.csv"), row.names = FALSE)

# 2. bar geometry: 2 bars per method, small inner gap, one tick per method at the midpoint ----
pos <- seq_along(method_order)
half_width <- 0.4
gap_inner  <- 0.05

bar_pos <- match(plot_df$base_method, method_order)
plot_df$xl <- ifelse(plot_df$group == "only miRNA data", bar_pos - half_width, bar_pos + gap_inner / 2)
plot_df$xr <- ifelse(plot_df$group == "only miRNA data", bar_pos - gap_inner / 2, bar_pos + half_width)

fmt_label <- function(x){
  if (x == 0) return("0")
  if (x < 0.001) return("<0.001")
  sprintf("%.3f", x)
}

# fixed at 0.3 per reviewer-figure convention (data max is 0.267, so it still clears the
# tallest bar with headroom for its value label)
y_max <- 0.3

# Input-group transparency is independent of the symbol pattern.
alpha_miR_only <- 0.55

# Draw isolated '+' / '-' glyphs on a regular grid inside a bar.
# Spacing and margins are defined in physical inches so that the symbols
# do not stretch when the x/y axis scales differ.
draw_symbol_pattern <- function(xl, xr, height, symbol,
                                ink = "grey40", spacing_x = 0.115,
                                spacing_y = 0.115, cex = 0.6) {
  if (!symbol %in% c("+", "-") || !is.finite(height) || height <= 0) {
    return(invisible(NULL))
  }

  # Convert bar dimensions to device inches.
  left_in   <- grconvertX(xl, from = "user", to = "inches")
  right_in  <- grconvertX(xr, from = "user", to = "inches")
  bottom_in <- grconvertY(0, from = "user", to = "inches")
  top_in    <- grconvertY(height, from = "user", to = "inches")

  # Leave sufficient margins to prevent glyphs crossing bar borders.
  margin_x <- 0.045
  margin_y <- 0.050
  if (right_in - left_in <= 2 * margin_x ||
      top_in - bottom_in <= 2 * margin_y) {
    return(invisible(NULL))
  }

  # nx <- max(1L, floor((right_in - left_in - 2 * margin_x) / spacing_x) + 1L)
  # ny <- max(1L, floor((top_in - bottom_in - 2 * margin_y) / spacing_y) + 1L)
  # xs <- (left_in + right_in) / 2 + (seq_len(nx) - (nx + 1) / 2) * spacing_x
  # ys <- (bottom_in + top_in) / 2 + (seq_len(ny) - (ny + 1) / 2) * spacing_y
  
  nx <- 3L  # exactly three symbol columns per bar
  ny <- max(1L, floor((top_in - bottom_in - 2 * margin_y) / spacing_y) + 1L)
  xs <- seq(left_in + margin_x, right_in - margin_x, length.out = nx)
  ys <- (bottom_in + top_in) / 2 + (seq_len(ny) - (ny + 1) / 2) * spacing_y

  for (xx in xs) {
    for (yy in ys) {
      text(grconvertX(xx, from = "inches", to = "user"),
           grconvertY(yy, from = "inches", to = "user"),
           labels = symbol, col = ink, cex = cex, font = 2)
    }
  }
  invisible(NULL)
}

# 3. draw ----
# family = "Helvetica": true Arial is not installed on this system (checked via
# systemfonts::system_fonts() / fc-list) -- Helvetica is the metric-compatible design basis
# for Arial, ships as one of the PDF base-14 fonts (always renders correctly, stays true
# vector text, no font embedding needed), and is the documented fallback for "Arial" already
# used elsewhere in this project (see revise/miRNet/R/DataUtils.R's CairoFonts() call).
pdf(plot_path, width = 8.5, height = 6.5, family = "Helvetica")
par(mar = c(5, 5, 1.5, 1.5), family = "Helvetica") # c(6, 5, 4, 1.5)

plot(NA, xlim = c(1 - half_width - 0.3, length(method_order) + half_width + 0.3), ylim = c(0, y_max),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", yaxs = "i")

for (i in seq_len(nrow(plot_df))) {
  row <- plot_df[i, ]
  col <- fill_col_base[[as.character(row$base_method)]]
  if (row$group == "only miRNA data") col <- adjustcolor(col, alpha.f = alpha_miR_only)
  rect(row$xl, 0, row$xr, row$TPR, col = col, border = "black", lwd = 1)
  if (row$pattern != "none" && row$TPR > 0) {
    draw_symbol_pattern(row$xl, row$xr, row$TPR, symbol = row$pattern)
  }
  if (row$TPR > 0) {
    text((row$xl + row$xr) / 2, row$TPR + y_max * 0.018, fmt_label(row$TPR), cex = 0.55, xpd = TRUE)
  }
}

box(bty = "o")
axis(2, las = 1, cex.axis = 0.9)
mtext("True Positive Rate", side = 2, line = 3, cex = 1)

axis(1, at = pos, labels = FALSE)
usr <- par("usr")
text(x = pos, y = usr[3] - 0.035 * (usr[4] - usr[3]), labels = method_order,
     srt = 30, adj = 1, xpd = TRUE, cex = 0.85)

# mtext("Overall True Positive Rates for Methods", side = 3, line = 2.3, cex = 1.15, font = 2)
# mtext("Split by Input Data Type", side = 3, line = 1.1, cex = 0.9, font = 1)

# Separate mini legends: bar opacity encodes input type;
# plus/minus glyphs encode the two altered-input classes.

legend_xl <- 0.73 - 0.6
legend_xr <- 0.97 - 0.6
legend_key_height <- 0.011
legend_y_top <- c(0.286 + 0.005, 0.268 + 0.005, 0.250 + 0.005, 0.232 + 0.005)
legend_labels <- c("Left: only miRNA data", "Right: miRNA and gene data",
                   "Add information", "Reduce information")
legend_fills <- c(adjustcolor("grey40", alpha.f = alpha_miR_only),
                  "grey40", "white", "white")
legend_symbols <- c("", "", "+", "-")

for (j in seq_along(legend_labels)) {
  yt <- legend_y_top[j]
  yb <- yt - legend_key_height
  rect(legend_xl, yb, legend_xr, yt,
       col = legend_fills[j], border = "black",
       lwd = 1)  #if (j >= 3L) 1.6 else 0.8)
  if (j >= 3L) {
    text((legend_xl + legend_xr) / 2, (yb + yt) / 2,
         labels = legend_symbols[j], col = "grey30", cex = 0.75, font = 2)
  }
  text(legend_xr + 0.10, (yb + yt) / 2,
       labels = legend_labels[j], adj = c(0, 0.5), cex = 0.72)
}


dev.off()

cat("Saved:\n  ", plot_path, "\n  ", file.path(work_dir, "FigureS10_data.csv"), "\n")
