library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr) # deal with strings
library(conflicted) # deal with conflicted functions
# library(miRBaseConverter) # transform miRNA name version
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
library(RColorBrewer)
library(grid)
library(circlize)
library(gridExtra)
library(patchwork)
library(scales)
library(tibble)

library(ggalluvial)
library(ggplot2)
library(ggnewscale)
library(colorspace)

library(reshape2)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::coalesce)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(base::`%*%`)
