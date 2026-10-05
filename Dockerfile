# miREA analysis environment
# Base image ships R 4.4.x + Bioconductor 3.20 + RStudio Server, matching the
# package versions pinned in README.md (ComplexHeatmap 2.22.0 / R 4.4.2 era).
FROM bioconductor/bioconductor_docker:3.20

# --- system deps -------------------------------------------------------
# python3 + numpy/scipy are required by reticulate::py_require() calls in
# R/get_input_data.R (calc_MGI_corr uses scipy.stats.spearmanr).
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir --break-system-packages numpy scipy

ENV RETICULATE_PYTHON=/usr/bin/python3

# --- R deps --------------------------------------------------------------
# Bioconductor packages first (binary repo configured by the base image).
RUN R -e "BiocManager::install(c('clusterProfiler', 'ComplexHeatmap'), update = FALSE, ask = FALSE)"

# CRAN packages used by R/lib.R.
RUN R -e "install.packages(c( \
    'dplyr', 'tidyr', 'tidyverse', 'stringr', 'conflicted', \
    'igraph', 'V8', 'data.table', 'Matrix', \
    'RColorBrewer', 'circlize', 'gridExtra', 'patchwork', 'scales', 'tibble', \
    'ggalluvial', 'ggplot2', 'ggnewscale', 'colorspace', 'reshape2', 'reticulate' \
    ), repos = BiocManager::repositories())"

# --- project code ---------------------------------------------------------
# Only the library code + tutorials are baked in. data/ (~530MB) and
# analysis/ (~90MB, mostly precomputed benchmark outputs for the paper) are
# intentionally NOT copied into the image — mount them at runtime instead,
# e.g.:
#   docker run -v /path/to/miREA:/home/rstudio/miREA ...
WORKDIR /home/rstudio/miREA

COPY R/ R/
COPY tutorial/ tutorial/
COPY function.RData .
COPY NAMESPACE .
COPY README.md .

RUN chown -R rstudio:rstudio /home/rstudio/miREA

EXPOSE 8787
