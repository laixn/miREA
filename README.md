# miREA: a network-based tool for edge-based miRNA-oriented enrichment analysis

We present miRNA-oriented Enrichment Analysis (**miREA**) framework, which integrates miRNA–gene interaction (**MGI**) networks with miRNA and gene transcriptomic data to improve functional interpretation of miRNAs.

miREA is designed to address the intrinsic many-to-many regulatory interactions of MGIs and reduce biases introduced by conventional node-centric enrichment aprroaches, by explicitly modeling **miRNA–gene regulatory networks**, **pathway topology**, and **expression-informed MGI edge scores**.

---

## Key Features

* **MGI edge-based** miRNA enrichment analysis
* **Edge scores** to measure and differentiate regulatory strengths of numerous MGIs
* Five novel edge-based enrichment algorithms that include **over-representation**, **scoring-based**, **topology-aware**, and **network-propagation** apporaches.
* **Comprehansive benchmarking framework** that evaluates sensitivity, specificity, discriminative ability, robustness, biological relevance, and computational efficiency across multiple cancer types.
* **Advanced visualization modules** to enhance the interpretability for intuitive exploration of miRNA-gene-pathway regulatory mechanisms.

---

## Usage
miREA are implemented based on R. Required R packages include (but are not limited to): 
* dplyr (1.14), tidyr (1.3.1), tidyverse (2.0.0), stringr (1.5.1), conflicted (1.2.0), igraph (2.1.4), V8 (6.0.1), data.table (1.16.4), Matrix (1.7.2), parallel (4.4.2), reticulate (1.40.0), ComplexHeatmap (2.22.0), RColorBrewer (1.1.3), grid (4.4.2), circlize (0.4.16), gridExtra (2.3), patchwork (1.3.0), scales (1.3.0), tibble (3.2.1), ggalluvial (0.12.5), ggplot2 (3.5.1), ggnewscale (0.5.0), colorspace (2.1.1), reshape2 (1.4.4)

Please follow the following steps to get enrichment results (see details at analysis/3_case_study/example.R):
1. load all functions in R/ and library packages needed.
   ```r
   setwd("[directory of miREA]") # set working directory
   load(function.RData) # it contains all functions
   source(R/lib.R)
   ```
3. Get raw data ready (see details at ```get_all_input_data()``` function).
2. Run ```get_all_input_data()``` to get prepared for data required for miREA, which will be saved at `data/input_data/[pathway]/`.
3. Run ```miREA()``` to conduct enrichment analysis under prefered methods, which will be save at `results/[pathway]/[cancer]/`.
4. Run ```plot_summary()``` and ```plot_heatmap_sankey()``` to get visualization plots, which will be save at `results/[pathway]/[cancer]/plot/`.

---

## Output

miREA returns results containing:

* a list for the enrichment result, including parameters, valid methods, time for different methods, and results containing a dataframe for each approach.
* visualization plots for methods comparison, and heatmap-sankey plot for each approach to discovery miRNA-gene-pathway regulatory mechanisms.

---

## Repository Structure

```
miREA/
├── function.RData     # all functions for miREA, integrating functions in R/
├── R/                 # all R functions
├── Javascript/        # Javascript codes for estimating p-values for permutation test

├── data/              # processed data
  ├── raw_data/
    ├── background/    # global miRNA-gene interaction and gene-gene interaction background
    ├── cancer_data/   # data needed for enrichment
      ├── DEG/         # gene differential expression analysis result
      ├── DEmiR/       # miRNA differential expression analysis result
      └── paired/      # gene and miRNA expression profiles and paired information
    └── pathway/       # contains extracted pathway gene sets
  ├── input_data/      # input data generated from get_all_input_data() function
  └── cancer_list/     # contains cancer-related genes and miRNAs from well-established research

├── results/           # enrichment results

├── analysis/          # codes, results, and plots generated for analysis
  ├── 2.1_positive_benchmark  # corresponds to Figure 2A
  ├── 2.2_negative_benchmark  # corresponds to Figure 2B, Figure 2C, and Supplementary Fig. S2
  ├── 2.3_distinguish_ability # corresponds to Figure 3A
  ├── 2.4_cancer_identify_ability # corresponds to Figure 3B
  ├── 2.5_robustness  # corresponds to Supplementary Fig. S3 and Fig. S4
  ├── 2.6_data_source  # corresponds to Supplementary Fig. S5
  ├── 2.7_time_test  # corresponds to Figure 3C and Supplementary Fig. S6
  ├── 3_case_study  # corresponds to Figure 4

└── README.md
```

---

## Citation

---

## Contact

