#' Toy BLCA example data for miREA
#'
#' A small, self-contained example dataset derived from the real BLCA
#' (bladder urothelial carcinoma) data used in the miREA manuscript. Gene and
#' miRNA differential expression are computed (DESeq2, tumor vs. normal) from
#' 80 tumor + 19 normal TCGA samples. The resulting DE-testable gene/miRNA
#' universe is used to subset 16 bladder-cancer-specific MSigDB gene sets
#' (e.g. KEGG_BLADDER_CANCER, LINDGREN_BLADDER_CANCER_*) and the
#' manuscript's background_MGI/background_GGI. 80 tumor samples is the
#' smallest size found to give every method, including the edge-based
#' methods' stricter correlation-significance filters, real non-empty
#' results, so every exported function in the package can run on real (not
#' synthetic) data. `toy_gene_mat`/`toy_mir_mat` use only the tumor samples
#' (correlation is computed within the tumor population, as in the
#' manuscript). See `data-raw/generate_toy_data.R` in the package's GitHub
#' repository (\url{https://github.com/laixn/miREA}; not included in the
#' installed package, since it depends on local TCGA files that are not
#' redistributed) for the generation script. For the full datasets used in
#' the manuscript, see the Zenodo deposit at
#' \url{https://zenodo.org/records/18803209}.
#'
#' @format \describe{
#'   \item{toy_pathway}{A data.frame with columns `pathway`, `gene`: pathway-gene
#'     membership for 16 bladder-cancer-specific MSigDB gene sets.}
#'   \item{toy_gene_DEdata}{A data.frame with columns `gene`, `log2FoldChange`,
#'     `stat`, `padj`: differential expression results for the gene universe.}
#'   \item{toy_mir_DEdata}{A data.frame with columns `miRNA`, `log2FoldChange`,
#'     `stat`, `padj`: differential expression results for the miRNA universe.}
#'   \item{toy_gene_mat}{A gene-by-sample expression matrix (data.frame, genes
#'     as rownames), the 80 tumor TCGA samples.}
#'   \item{toy_mir_mat}{A miRNA-by-sample expression matrix (data.frame,
#'     miRNAs as rownames), the same 80 tumor TCGA samples as `toy_gene_mat`.}
#'   \item{toy_background_MGI}{A data.frame with columns `miRNA`, `gene`:
#'     curated miRNA-gene interactions restricted to the toy universe.}
#'   \item{toy_background_GGI}{A data.frame with columns `source_genesymbol`,
#'     `target_genesymbol`, `interact`: Omnipath gene-gene interactions
#'     restricted to the toy gene universe.}
#' }
#' @source \describe{
#'   \item{toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat}{Derived
#'     from TCGA BLCA RNA-seq/miRNA-seq data, as processed for the miREA
#'     manuscript (\doi{10.64898/2026.02.27.708509}).}
#'   \item{toy_pathway}{Subset of MSigDB bladder-cancer-specific gene sets
#'     (16 sets, including `KEGG_BLADDER_CANCER`). KEGG-derived gene sets
#'     carry their own redistribution restrictions (see
#'     \url{https://www.genome.jp/kegg/legal.html}), separate from MSigDB's
#'     own license terms; this has not been independently verified for the
#'     bundled set.}
#'   \item{toy_background_MGI}{Subset of a miRNA-gene target-site prediction
#'     table from Hua, X., Li, Y., Pentaparthi, S.R., et al. (2023).
#'     Landscape of MicroRNA Regulatory Network Architecture and Functional
#'     Rerouting in Cancer. Cancer Res, 83, 59-73, deduplicated to
#'     (miRNA, gene) pairs and harmonized to current HGNC gene symbols.}
#'   \item{toy_background_GGI}{Subset of Omnipath interactions
#'     (\code{OmnipathR::import_all_interactions()}), filtered to
#'     interactions with experimental support (\code{n_references != 0}),
#'     with \code{interact} derived from \code{is_stimulation -
#'     is_inhibition}, harmonized to current HGNC gene symbols. Omnipath
#'     aggregates many underlying interaction databases, some of which
#'     (e.g. SPIKE, NetPath) impose non-commercial-use restrictions; which
#'     resources contributed to this particular subset has not been
#'     independently verified.}
#' }
#' @name miREA_example_data
#' @aliases toy_pathway toy_gene_DEdata toy_mir_DEdata toy_gene_mat toy_mir_mat toy_background_MGI toy_background_GGI
#' @docType data
#' @keywords datasets
#' @examples
#' data(toy_pathway)
#' data(toy_gene_DEdata)
#' data(toy_mir_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
#' data(toy_background_MGI)
#' data(toy_background_GGI)
NULL

#' Full background MGI/GGI data for miREA
#'
#' The complete (non-toy) miRNA-gene interaction (MGI) and gene-gene
#' interaction (GGI) backgrounds used in the miREA manuscript. Unlike the
#' `toy_background_MGI`/`toy_background_GGI` subsets above, these are not
#' restricted to any particular gene/miRNA universe, and are the datasets
#' [getAllInputData()] falls back to when `background_MGI`/`background_GGI`
#' is not supplied.
#'
#' @format \describe{
#'   \item{background_MGI}{A data.frame with columns `miRNA`, `gene`:
#'     306,150 miRNA-gene interactions covering 2,562 miRNAs and 5,135
#'     protein-coding genes.}
#'   \item{background_GGI_Omnipath}{A data.frame with columns
#'     `source_genesymbol`, `target_genesymbol`, `interact`: a frozen
#'     snapshot of 103,316 Omnipath gene-gene interactions.}
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
#'   \item{background_GGI_Omnipath}{282,504 experimentally validated,
#'     Homo-sapiens-only Omnipath interactions (\code{n_references != 0}),
#'     with undirected interactions duplicated as bidirectional edges, then
#'     harmonized to current HGNC gene symbols. This is a frozen snapshot:
#'     fetching Omnipath live via [getAllInputData()]'s `GGI_source =
#'     "Omnipath"` option reflects Omnipath's current (larger, still
#'     growing) database contents instead. Omnipath aggregates many
#'     underlying interaction databases, some of which (e.g. SPIKE,
#'     NetPath) impose non-commercial-use restrictions; which resources
#'     contributed to this particular snapshot has not been independently
#'     verified.}
#'   \item{background_GGI_Reactome}{Reconstructed from the Reactome
#'     Functional Interaction network's predefined gene-gene interactions,
#'     with undirected/mixed-sign edges expanded into directed pairs, then
#'     harmonized to current HGNC gene symbols. This is a frozen snapshot
#'     of a specific dated Reactome FI release; fetching Reactome live via
#'     [getAllInputData()]'s `GGI_source = "Reactome"` option reconstructs
#'     the network from a pinned (but not HGNC-harmonized) release instead.}
#' }
#' @name miREA_background_data
#' @aliases background_MGI background_GGI_Omnipath background_GGI_Reactome
#' @docType data
#' @keywords datasets
#' @examples
#' data(background_MGI)
#' data(background_GGI_Omnipath)
#' data(background_GGI_Reactome)
NULL
