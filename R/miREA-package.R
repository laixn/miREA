#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr mutate select filter arrange rename distinct group_by
#' @importFrom dplyr summarise left_join inner_join bind_rows case_when pull
#' @importFrom dplyr ungroup if_else %>% semi_join n n_distinct desc row_number
#' @importFrom dplyr slice_head
#' @importFrom ggplot2 ggplot aes geom_point geom_tile geom_text geom_bar
#' @importFrom ggplot2 geom_density scale_fill_manual scale_color_manual
#' @importFrom ggplot2 scale_fill_gradientn scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_x_discrete labs theme theme_minimal theme_void
#' @importFrom ggplot2 element_text element_blank element_line element_rect
#' @importFrom ggplot2 ggtitle guides guide_legend margin ggsave after_stat
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 scale_shape_manual scale_size scale_linetype_manual
#' @importFrom ggplot2 scale_color_gradientn guide_colorbar ggplotGrob
#' @importFrom ggplot2 scale_fill_identity
#' @importFrom grid gpar unit grid.draw grid.rect pushViewport popViewport
#' @importFrom grid viewport unit.c
#' @importFrom circlize colorRamp2
#' @importFrom patchwork wrap_plots plot_layout plot_annotation
#' @importFrom tibble as_tibble tibble
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggalluvial geom_stratum geom_flow
#' @importFrom ComplexHeatmap Heatmap Legend HeatmapAnnotation rowAnnotation
#' @importFrom ComplexHeatmap columnAnnotation draw anno_text anno_simple
#' @importFrom ComplexHeatmap packLegend
#' @importFrom ComplexHeatmap anno_barplot AnnotationFunction
#' @importFrom ComplexHeatmap decorate_heatmap_body "%v%"
#' @importFrom igraph V E "V<-" "E<-" betweenness keeping_degseq ecount vcount
#' @importFrom igraph rewire
#' @importFrom igraph graph_from_data_frame edge_betweenness as_data_frame
#' @importFrom igraph union make_line_graph as_adjacency_matrix
#' @importFrom Matrix sparseMatrix Matrix
#' @importFrom stats setNames na.omit p.adjust sd cor quantile median pt
#' @importFrom stats complete.cases prcomp var pnorm
#' @importFrom scales label_number
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' @importFrom methods is as
#' @importFrom utils object.size
#' @importFrom data.table rbindlist setkey setnames setDT fifelse data.table
#' @importFrom data.table as.data.table
#' @importFrom BiocParallel bplapply bpnworkers SerialParam SnowParam
## usethis namespace: end
NULL

## Required for data.table's `[.data.table` NSE (e.g. dict[, ratio := ...])
## to work correctly when data.table is only used via `::` (never attached).
## Without this, data.table's cedta() check fails and `[` silently falls back
## to `[.data.frame` semantics, breaking unquoted column references.
## See data.table's "datatable-importing" vignette.
.datatable.aware <- TRUE

## These are data.table/dplyr/ggplot2 non-standard-evaluation column names and
## computed aesthetics (e.g. dict[, ratio := ...], aes(x = Method1)), not
## undefined globals -- see the data.table "importing" vignette for why this
## declaration is needed to silence R CMD check's static-analysis NOTE.
utils::globalVariables(c(
    ":=", ".", "alluvium", "centrality", "corPadj", "Demographic", "density",
    "EB", "edge", "fill_color", "flow_value", "flow_value.x", "flow_value.y",
    "from", "gene", "gene_count", "gene_flow_value", "geneLogFC", "genePadj",
    "geneStat", "i.strength", "ID", "is_DE", "method", "Method1", "Method2",
    "MGI", "mirLogFC", "miRNA", "mirPadj", "mirStat", "n_DE", "n_pathway", "normratio",
    "p_value", "padj", "pathway", "PC1", "PC2", "positive_rate", "pvalue",
    "rand_score_dist", "ratio", "score", "stratum",
    "strength", "text", "to", "total_flow", "Type", "Value", "x", "y"
))
