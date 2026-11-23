#' Toy BLCA example data for miREA
#'
#' A small, self-contained example dataset derived from the real BLCA
#' (bladder urothelial carcinoma) data used in the miREA manuscript. Gene and
#' miRNA differential expression are computed (DESeq2, tumor vs. normal) from
#' 80 tumor + 19 normal TCGA samples. The resulting DE-testable gene/miRNA
#' universe is used to subset the 10 Reactome pathways most enriched in the
#' full manuscript BLCA run and the manuscript's background_MGI/
#' background_GGI. 80 tumor samples is the
#' smallest size found to give every method, including the edge-based
#' methods' stricter correlation-significance filters, real non-empty
#' results, so every exported function in the package can run on real (not
#' synthetic) data. `toy_gene_mat`/`toy_mir_mat` use only the tumor samples
#' (correlation is computed within the tumor population, as in the
#' manuscript). For the full datasets used in the manuscript, see the
#' Zenodo deposit at \url{https://zenodo.org/records/18803209}.
#'
#' @format \describe{
#'   \item{toy_pathway}{A data.frame with columns `pathway`, `gene`: pathway-gene
#'     membership for the 10 Reactome pathways most enriched in the full
#'     manuscript BLCA run.}
#'   \item{toy_gene_DEdata}{A data.frame with columns `gene`, `log2FoldChange`,
#'     `stat`, `padj`: differential expression results for the gene universe.}
#'   \item{toy_mir_DEdata}{A data.frame with columns `miRNA`, `log2FoldChange`,
#'     `stat`, `padj`: differential expression results for the miRNA universe.}
#'   \item{toy_gene_mat}{A gene-by-sample expression matrix (data.frame, genes
#'     as rownames), the 80 tumor TCGA samples.}
#'   \item{toy_mir_mat}{A miRNA-by-sample expression matrix (data.frame,
#'     miRNAs as rownames), the same 80 tumor TCGA samples as `toy_gene_mat`.}
#' }
#' @source \describe{
#'   \item{toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat}{TCGA
#'     BLCA RNA-seq/miRNA-seq read counts (same source data as the
#'     manuscript, see Methods 2.3, \doi{10.64898/2026.02.27.708509}),
#'     subset to 80 tumor + 19 normal samples: every BLCA patient with
#'     paired tumor+normal RNA-seq and miRNA-seq is used for the normal
#'     group (19 patients); 80 tumor samples is the smallest tumor-group
#'     size found empirically to give every method, including the
#'     edge-based methods' stricter correlation-significance filters, a
#'     non-empty result. Genes/miRNAs with read counts > 10 in at least
#'     five samples were tested with DESeq2 (design `~ type`, tumor vs.
#'     normal), as in the manuscript's own DE pipeline, run independently
#'     on this toy sample subset (not a subset of the manuscript's own DE
#'     results). `toy_gene_mat`/`toy_mir_mat` are the raw counts for the
#'     tumor samples only (matching the manuscript's correlation step,
#'     computed within the tumor population).}
#'   \item{toy_pathway}{The 10 Reactome pathways with the lowest (padj,
#'     p_value) in the Edge_Score1D output of the full manuscript BLCA run
#'     (Reactome pathway set as in Methods 2.2, 2,766 pathways, v94;
#'     \doi{10.64898/2026.02.27.708509}; full enrichment results at
#'     \url{https://zenodo.org/records/18803209}), restricted to genes
#'     present in `toy_gene_DEdata`. Reactome's data is
#'     CC0 (public domain), so this carries no usage restriction, unlike
#'     the KEGG-derived MSigDB sets used here previously.}
#' }
#' @note These toy objects are paired with the full, real background_MGI
#'   and background_GGI_Reactome datasets (see [miREA_background_data]) for
#'   examples/tests that need a background_MGI/background_GGI argument —
#'   there is no separate "toy" subset of those, since the full datasets are
#'   already small enough to use directly.
#' @name miREA_example_data
#' @aliases toy_pathway toy_gene_DEdata toy_mir_DEdata toy_gene_mat toy_mir_mat
#' @docType data
#' @keywords datasets
#' @examples
#' data(toy_pathway)
#' data(toy_gene_DEdata)
#' data(toy_mir_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
NULL

#' Full background MGI/GGI data for miREA
#'
#' The complete (non-toy) miRNA-gene interaction (MGI) background, and the
#' Reactome gene-gene interaction (GGI) background, used in the miREA
#' manuscript. These are the datasets [getAllInputData()] falls back to when
#' `background_MGI` is not supplied, or `background_GGI` is not supplied and
#' `GGI_source = "Reactome"`. There is no bundled Omnipath GGI background:
#' `GGI_source = "Omnipath"` instead fetches Omnipath's full interaction set
#' live, with no license/source filtering, since Omnipath aggregates many
#' underlying databases with their own (sometimes restrictive) licenses,
#' unlike Reactome's data, which is CC0. Users who need a license-restricted
#' subset should filter it themselves and supply it via `background_GGI`
#' (see the "Third-party data licensing" section of the package README).
#'
#' @format \describe{
#'   \item{background_MGI}{A data.frame with columns `miRNA`, `gene`:
#'     306,150 miRNA-gene interactions covering 2,562 miRNAs and 5,135
#'     protein-coding genes.}
#'   \item{background_GGI_Reactome}{A data.frame with columns
#'     `source_gene`, `target_gene`, `direction`: a frozen snapshot of
#'     464,630 Reactome Functional Interaction network edges, covering
#'     13,640 genes.}
#' }
#' @source \describe{
#'   \item{background_MGI}{The most comprehensive MGI dataset available in
#'     the literature: 251,888,701 putative MGIs derived from miRNA target
#'     prediction algorithms, filtered using Argonaute (AGO) CLIP-seq data
#'     to retain experimentally supported interactions, with miRNA and gene
#'     names unified via miRBase and HGNC gene symbols v115 (09/2025). See
#'     Hua, X., Li, Y., Pentaparthi, S.R., et al. (2023). Landscape of
#'     MicroRNA Regulatory Network Architecture and Functional Rerouting in
#'     Cancer. Cancer Res, 83, 59-73.}
#'   \item{background_GGI_Reactome}{Reconstructed from the Reactome
#'     Functional Interaction network's predefined gene-gene interactions,
#'     with undirected/mixed-sign edges expanded into directed pairs, then
#'     harmonized to current HGNC gene symbols. This is a frozen snapshot
#'     of a specific dated Reactome FI release; fetching Reactome live via
#'     [getAllInputData()]'s `GGI_source = "Reactome"` option reconstructs
#'     the network from a pinned (but not HGNC-harmonized) release instead.
#'     Reactome's data is licensed CC0 (public domain dedication), so this
#'     snapshot carries no usage restrictions.}
#' }
#' @name miREA_background_data
#' @aliases background_MGI background_GGI_Reactome
#' @docType data
#' @keywords datasets
#' @examples
#' data(background_MGI)
#' data(background_GGI_Reactome)
NULL
