options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "dplyr", "stringr", "conflicted",
  "igraph", "V8", "data.table", "Matrix", "reticulate",
  "circlize", "patchwork", "scales",
  "tibble", "ggalluvial", "ggplot2", "ggnewscale", "reshape2"
)
bioc_pkgs <- c("clusterProfiler", "ComplexHeatmap", "miRBaseConverter")

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

library(dplyr)
library(stringr) # deal with strings
library(conflicted) # deal with conflicted functions
library(miRBaseConverter) # transform miRNA name version
# library(HGNChelper) # transform HGNC symbol to the newest one.

library(clusterProfiler) # enrichment analysis

library(igraph)#network analysis

library(V8) # for javascript

library(data.table)
library(Matrix) # for Edge-Network
library(parallel)# parallel computing
library(reticulate) # load python

## visualization
library(ComplexHeatmap)
library(grid)
library(circlize)
library(patchwork)
library(scales)
library(tibble)

library(ggalluvial)
library(ggplot2)
library(ggnewscale)

library(reshape2)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::coalesce)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(base::`%*%`)
