## Regenerates the toy_gene_mat/toy_gene_DEdata/toy_mir_mat/toy_mir_DEdata example
## datasets shipped in data/. Not run during R CMD build/check; re-run manually if
## the toy data needs to be rebuilt. Requires data.table and access to the
## manuscript's deposited raw data (see `zenodo_dir` below).
##
## This revision replaces an earlier version of this script (recovered from git
## history at commit 5bd0f45; deleted at c60cd6e) that re-ran DESeq2 on a fresh
## 80-tumor-sample subsample of the raw TCGA counts. Instead, it directly subsets
## the manuscript's own already-computed, already-deposited DE results and raw
## expression matrices for BLCA (`raw_data/cancer_data/{DEG,DEmiR,paired}/`, the
## same Zenodo deposit `toy_pathway`'s 10 Reactome pathways were themselves
## selected from -- see R/data.R and the enrichment run under
## `data/input_data/{hallmark,Reactome}/BLCA_*_input_data.RData`, which share
## this same underlying DE/expression data regardless of which pathway database
## was tested downstream -- verified by exact value match on sampled
## genes/miRNAs/samples against `BLCA_hallmark_input_data.RData$raw_data`).
## This is strictly better than re-deriving our own DESeq2 run on a hand-picked
## sample subset:
##  - no new sample-selection decision to make or justify;
##  - the DE results are numerically identical to the ones `toy_pathway` was
##    chosen from, so there is no risk of the toy pathway/DE data disagreeing;
##  - the raw expression matrices already contain all 403 tumor-only samples
##    (TCGA sample-type code "01A"), giving far more statistical power for the
##    correlation/permutation-based methods than the old 80-sample subset,
##    without having to pick a new sample count ourselves;
##  - no DESeq2 re-run needed at all (much faster, and not dependent on access
##    to the raw per-base-pair TCGA count files any more).
##
## toy_pathway.rda/background_MGI.rda/background_GGI.rda are NOT regenerated here
## (unchanged; toy_pathway is an independently-curated, fixed set of 10 Reactome
## pathways, and background_MGI/GGI are shipped in full, undownsampled form).

suppressMessages(library(data.table))

zenodo_dir <- "/scratch/project_2011179/code/miREA_zenodo/data/raw_data/cancer_data"

## DEG.csv/DEmiR.csv: gene/miRNA, baseMean, log2FoldChange, lfcSE, stat, pvalue,
## padj, then per-sample normalized counts (TP+NT) -- only the first 7 (summary)
## columns are needed here. Rows with NA padj (DESeq2 independent filtering) are
## already excluded, matching `results(dds)` with NAs dropped.
full_gene_DEdata <- fread(file.path(zenodo_dir, "DEG/BLCA_DEG.csv"),
  select = c("gene", "log2FoldChange", "stat", "padj")
)
full_mir_DEdata <- fread(file.path(zenodo_dir, "DEmiR/BLCA_DEmiR.csv"),
  select = c("miRNA", "log2FoldChange", "stat", "padj")
)

## paired/BLCA_{gene,mir}_TP.csv: raw counts, tumor samples only (unnamed first
## column holds the gene/miRNA symbol).
full_gene_mat <- fread(file.path(zenodo_dir, "paired/BLCA_gene_TP.csv"))
full_gene_mat <- as.data.frame(full_gene_mat, check.names = FALSE)
rownames(full_gene_mat) <- full_gene_mat[[1]]
full_gene_mat[[1]] <- NULL

full_mir_mat <- fread(file.path(zenodo_dir, "paired/BLCA_mir_TP.csv"))
full_mir_mat <- as.data.frame(full_mir_mat, check.names = FALSE)
rownames(full_mir_mat) <- full_mir_mat[[1]]
full_mir_mat[[1]] <- NULL

load(file.path("data", "toy_pathway.rda")) # provides `toy_pathway` (unchanged)
load(file.path("data", "background_MGI.rda")) # provides `background_MGI` (unchanged, full)

pathway_genes <- unique(toy_pathway$gene)

## Pathway-relevant miRNAs: anything in the (full, undownsampled) background_MGI
## that targets a pathway gene, intersected with what's actually DE-testable in
## the manuscript's own miRNA universe.
pathway_mirs <- unique(background_MGI$miRNA[background_MGI$gene %in% pathway_genes])
pathway_mirs <- intersect(pathway_mirs, full_mir_DEdata$miRNA)

## A handful of pathway genes can be absent from gene_DEdata (DESeq2's
## independent filtering assigns NA padj to very-low-count genes, dropping them
## from `results()`, even though they're still present in the raw gene_mat) --
## expected and harmless: those genes just won't carry DE info downstream.
missing_from_DEdata <- setdiff(pathway_genes, full_gene_DEdata$gene)
if (length(missing_from_DEdata) > 0) {
  message(
    length(missing_from_DEdata), " pathway gene(s) have no DESeq2 result (independent filtering): ",
    paste(missing_from_DEdata, collapse = ", ")
  )
}

toy_gene_DEdata <- as.data.frame(full_gene_DEdata[full_gene_DEdata$gene %in% pathway_genes, ])
colnames(toy_gene_DEdata) <- c("gene", "log2FoldChange", "stat", "padj")

toy_mir_DEdata <- as.data.frame(full_mir_DEdata[full_mir_DEdata$miRNA %in% pathway_mirs, ])
colnames(toy_mir_DEdata) <- c("miRNA", "log2FoldChange", "stat", "padj")

toy_gene_mat <- full_gene_mat[rownames(full_gene_mat) %in% pathway_genes, , drop = FALSE]
toy_mir_mat <- full_mir_mat[rownames(full_mir_mat) %in% pathway_mirs, , drop = FALSE]

message(nrow(toy_gene_DEdata), " genes / ", ncol(toy_gene_mat), " tumor samples in toy_gene_mat/toy_gene_DEdata.")
message(nrow(toy_mir_DEdata), " miRNAs / ", ncol(toy_mir_mat), " tumor samples in toy_mir_mat/toy_mir_DEdata.")

## ---- Save ------------------------------------------------------------------
usethis::use_data(toy_gene_DEdata, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_mir_DEdata, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_gene_mat, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_mir_mat, overwrite = TRUE, compress = "xz")
