# miREA: a network-based tool for edge-based miRNA-oriented enrichment analysis

We present miRNA-oriented Enrichment Analysis (**miREA**) framework, which integrates miRNA–gene interaction (**MGI**) networks with transcriptomic data to improve functional interpretation of miRNAs.

miREA is designed to address the intrinsic many-to-many regulatory interactions of MGIs and reduce biases introduced by conventional node-centric enrichment aprroaches, by explicitly modeling **miRNA–gene regulatory networks**, **pathway topology**, and **expression-informed MGI edge scores**.

---

## Key Features

* **MGI edge-based** miRNA enrichment analysis
* **Edge scores** to measure and differentiate regulatory strengths of numerous MGIs
* Five novel edge-based enrichment algorithms that include **over-representation**, **scoring-based**, **topology-aware**, and **network-propagation** apporaches.
* **Comprehansive benchmarking framework** that evaluates sensitivity, specificity, discriminative ability, robustness, biological relevance, and computational efficiency across multiple cancer types.
* **Advanced visualization modules** to enhance the interpretability for intuitive exploration of miRNA-gene-pathway regulatory mechanisms.

---

## Installation

Currently, MiREA is under active development.

```bash
git clone https://github.com/your-username/MiREA.git
```

Required R packages include (but are not limited to):

* data.table
* igraph
* parallel
* fgsea
* ComplexHeatmap

---

## Input Data

MiREA requires the following inputs:

1. **Gene expression matrix**

   * Rows: genes
   * Columns: samples

2. **miRNA expression matrix**

   * Rows: miRNAs
   * Columns: samples

3. **miRNA–gene interaction network**

   * Directed edges from miRNA to gene
   * Optionally weighted by correlation or binding confidence

4. **Pathway definitions**

   * Gene sets or pathway graphs (e.g., KEGG, Reactome)

---

## Usage Example

```r
# Load MiREA functions
source("R/MiREA.R")

# Run MiREA
result <- MiREA(
  expr_gene   = gene_expr,
  expr_mirna  = mirna_expr,
  mgi_network = mgi_edges,
  pathways    = pathway_list,
  nperm       = 1000
)
```

---

## Output

MiREA returns a data frame containing:

* Enrichment score
* Normalized enrichment score (NES)
* Empirical p-value
* Adjusted p-value (BH)
* Pathway size and edge statistics

---

## Reproducibility

All analyses are fully reproducible given the same input data and random seed. Permutation procedures support parallel execution on HPC environments.

---

## Project Structure

```
MiREA/
├── R/                 # Core R functions
├── data/              # Example or processed data
├── example/           # Minimal working examples
├── scripts/           # Analysis and pipeline scripts
├── README.md
└── LICENSE
```

---

## Citation

If you use MiREA in your research, please cite:

> Author(s). *MiREA: a network-based framework for miRNA enrichment analysis*. (Manuscript in preparation)

---

## License

This project is released under the MIT License.

---

## Contact

For questions, suggestions, or collaborations, please contact:

**[Your Name]**
Doctoral Researcher, Tampere University
Email: [your.email@tuni.fi](mailto:your.email@tuni.fi)

