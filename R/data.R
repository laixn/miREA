#' Toy BLCA example data for miREA
#'
#' A small, self-contained example dataset derived from the real BLCA
#' (bladder urothelial carcinoma) data used in the miREA manuscript.
#' `toy_pathway` is a fixed set of 10 Reactome pathways selected from the full
#' manuscript BLCA run (see `@source` below). `toy_gene_DEdata`,
#' `toy_mir_DEdata`,
#' `toy_gene_mat`, and `toy_mir_mat` are the manuscript's own differential
#' expression results and raw expression matrices (403 tumor and 19 normal TCGA
#' samples). The normal samples are only used to compute differential
#' expression, not correlation (see below). These objects are restricted to
#' the genes/miRNAs needed to exercise `toy_pathway`, so every exported
#' function in the package can run on real, not synthetic, data without
#' shipping the full ~35,000-gene DE-testable universe. `toy_gene_mat` and
#' `toy_mir_mat` use only the tumor samples, since correlation is computed
#' within the tumor population, as in the manuscript. For the full datasets
#' used in the manuscript, see the Zenodo deposit at
#' \url{https://zenodo.org/records/18803209}.
#'
#' @format \describe{
#'   \item{toy_pathway}{A data.frame with columns `pathway` and `gene`, giving
#'     pathway-gene membership for the 10 Reactome pathways most enriched in
#'     the full manuscript BLCA run.}
#'   \item{toy_gene_DEdata}{A data.frame with columns `gene`, `log2FoldChange`,
#'     `stat`, and `padj`, giving differential expression results for the
#'     genes in `toy_pathway`. 1,142 of 1,145 pathway genes have a result.
#'     The other 3 were dropped by DESeq2's independent filtering and simply
#'     carry no DE info downstream.}
#'   \item{toy_mir_DEdata}{A data.frame with columns `miRNA`, `log2FoldChange`,
#'     `stat`, and `padj`, giving differential expression results for the
#'     miRNAs that target a `toy_pathway` gene (via `background_MGI`) and are
#'     DE-testable.}
#'   \item{toy_gene_mat}{A gene-by-sample expression matrix (data.frame, genes
#'     as rownames), restricted to 1,145 genes in `toy_pathway` and 403 tumor
#'     TCGA samples.}
#'   \item{toy_mir_mat}{A miRNA-by-sample expression matrix (data.frame,
#'     miRNAs as rownames), restricted to the miRNAs in `toy_mir_DEdata`, the
#'     same 403 tumor TCGA samples as `toy_gene_mat`.}
#' }
#' @source \describe{
#'   \item{toy_gene_DEdata, toy_mir_DEdata, toy_gene_mat, toy_mir_mat}{Subset,
#'     by gene and miRNA, of the manuscript's own full TCGA BLCA differential
#'     expression results and raw expression matrices (403 tumor and 19 normal
#'     samples, DESeq2, design `~ type`, tumor vs. normal, Methods 2.3,
#'     \doi{10.64898/2026.02.27.708509}). These are the exact DE results
#'     `toy_pathway` below was itself selected from, not an independently
#'     re-derived DE run on a different sample subset. Restricted to the genes
#'     in `toy_pathway`, and the miRNAs that target one of those genes via
#'     `background_MGI` and are DE-testable. `toy_gene_mat` and `toy_mir_mat`
#'     are the raw counts for the tumor samples only (matching the
#'     manuscript's correlation step, computed within the tumor population).}
#'   \item{toy_pathway}{The 10 Reactome pathways with the lowest (padj,
#'     p_value) in the Edge_Score1D output of the full manuscript BLCA run
#'     (Reactome pathway set as in Methods 2.2, 2,766 pathways, v94,
#'     \doi{10.64898/2026.02.27.708509}). Full enrichment results are
#'     available at \url{https://zenodo.org/records/18803209}. Reactome's
#'     data is CC0 (public domain), so this carries no usage restriction,
#'     unlike the KEGG-derived MSigDB sets used here previously.}
#' }
#' @note These toy objects are paired with the full, real background_MGI
#'   and background_GGI_Reactome datasets (see [miREA_background_data]) for
#'   examples/tests that need a background_MGI/background_GGI argument.
#'   There is no separate "toy" subset of those, since the full datasets are
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
#' `GGI_source = "Reactome"`. There is no bundled Omnipath GGI background.
#' `GGI_source = "Omnipath"` instead fetches Omnipath's full interaction set
#' live, with no license/source filtering, since Omnipath aggregates many
#' underlying databases with their own (sometimes restrictive) licenses,
#' unlike Reactome's data, which is CC0. Users who need a license-restricted
#' subset should filter it themselves and supply it via `background_GGI`
#' (see the "Third-party data licensing" section of the package README).
#'
#' @format \describe{
#'   \item{background_MGI}{A data.frame with columns `miRNA` and `gene`,
#'     giving 306,150 miRNA-gene interactions covering 2,562 miRNAs and 5,135
#'     protein-coding genes.}
#'   \item{background_GGI_Reactome}{A data.frame with columns
#'     `source_gene`, `target_gene`, and `direction`, giving a frozen
#'     snapshot of 464,630 Reactome Functional Interaction network edges,
#'     covering 13,640 genes.}
#' }
#' @source \describe{
#'   \item{background_MGI}{The most comprehensive MGI dataset available in
#'     the literature, built from 251,888,701 putative MGIs derived from
#'     miRNA target prediction algorithms, filtered using Argonaute (AGO)
#'     CLIP-seq data
#'     to retain experimentally supported interactions, with miRNA and gene
#'     names unified via miRBase and HGNC gene symbols v115 (09/2025). See
#'     Hua, X., Li, Y., Pentaparthi, S.R., et al. (2023). Landscape of
#'     MicroRNA Regulatory Network Architecture and Functional Rerouting in
#'     Cancer. Cancer Res, 83, 59-73.}
#'   \item{background_GGI_Reactome}{Reconstructed from the Reactome
#'     Functional Interaction network's predefined gene-gene interactions,
#'     with undirected/mixed-sign edges expanded into directed pairs, then
#'     harmonized to current HGNC gene symbols. This is a frozen snapshot
#'     of a specific dated Reactome FI release. Fetching Reactome live via
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
