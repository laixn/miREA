## Reconstructs how background_GGI_Omnipath/background_GGI_Reactome (bundled
## in data/, see data-raw/generate_background_data.R) were originally built
## from Omnipath and Reactome. Not run during R CMD build/check; this is
## provenance documentation, not part of the installed package. The package
## itself no longer fetches these live at runtime (see data(background_MGI),
## data(background_GGI_Omnipath), data(background_GGI_Reactome)) -- re-run
## the functions below only if you want a fresh snapshot instead of the
## bundled one (Omnipath's database grows over time, so a fresh pull will
## not match the bundled snapshot row-for-row; Reactome's dated release URLs
## are pinned and immutable, so re-running with the same date reproduces the
## bundled snapshot exactly).
##
## After the steps below, the manuscript's background_GGI_Omnipath.csv and
## background_GGI_Reactome.csv additionally went through a gene-symbol
## harmonization pass (check_gene_names(), which queries biomaRt for current
## HGNC symbols/aliases and is not reproduced here, since it requires a live
## Ensembl/biomaRt connection).

## ---- Omnipath ---------------------------------------------------------
## Requires the OmnipathR package (not a runtime dependency of miREA itself).
build_background_GGI_omnipath <- function() {
  interactions <- OmnipathR::import_all_interactions()
  interactions <- interactions[, c("source_genesymbol", "target_genesymbol", "is_directed",
                                    "is_stimulation", "is_inhibition", "n_references")]

  ## keep only interactions with experimental support
  interactions <- interactions[interactions$n_references != 0, ]

  ## drop rows describing protein complexes (comma-separated) or other
  ## compound identifiers (underscore-separated)
  interactions <- interactions[!grepl(",", interactions$source_genesymbol) &
                                  !grepl("_", interactions$source_genesymbol) &
                                  !grepl("_", interactions$target_genesymbol), ]

  ## undirected edges (is_directed == 0) are duplicated with source/target
  ## swapped, so they are represented in both directions
  undirected <- interactions[interactions$is_directed == 0, c(2, 1, 3:6)]
  names(undirected) <- names(interactions)
  interactions <- rbind(undirected, interactions)

  ## encode the interaction sign: +1 stimulation, -1 inhibition, 0 unclear
  interactions$interact <- interactions$is_stimulation - interactions$is_inhibition
  interactions <- interactions[, c("source_genesymbol", "target_genesymbol", "interact")]

  ## drop residual miRNA nodes (named like "hsa-...")
  interactions[!grepl("^hsa-", interactions$source_genesymbol) &
                 !grepl("^hsa-", interactions$target_genesymbol), ]
}

## ---- Reactome ----------------------------------------------------------
## Input: the Reactome Functional Interaction (FI) network flat file,
## distributed as "FIsInGene_<date>_with_annotations.txt" from Reactome's
## download page (reactome.org/download/tools/ReactomeFIs/), pinned here to
## the exact dated release used to build background_GGI_Reactome. Reactome
## hosts past dated releases permanently at stable URLs, so this will not
## drift over time (unlike Omnipath's continuously growing database).
## Column `Direction` encodes the relationship type: "-" (undirected),
## "->"/"<-" (single direction), "|-|"/"|->"/"<-|" (both-ways with
## activating/inhibiting flavors).
build_background_GGI_reactome <- function(
    url = "https://reactome.org/download/tools/ReactomeFIs/FIsInGene_04142025_with_annotations.txt.zip") {
  tmp_zip <- tempfile(fileext = ".zip")
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(c(tmp_zip, tmp_dir), recursive = TRUE), add = TRUE)

  utils::download.file(url, destfile = tmp_zip, quiet = TRUE, mode = "wb")
  utils::unzip(tmp_zip, exdir = tmp_dir)
  txt_file <- list.files(tmp_dir, pattern = "_with_annotations\\.txt$", full.names = TRUE)[1]
  fi <- utils::read.table(txt_file, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)

  expand_direction(fi)
}

expand_direction <- function(fi) {
  ## each direction code maps to how source/target are swapped on the
  ## second (reverse) row; single-direction rows ("->"/"<-") are kept as-is
  split_one <- function(gene1, gene2, dir) {
    switch(dir,
      "-"   = list(c(gene1, gene2, "-"), c(gene2, gene1, "-")),
      "<->" = list(c(gene1, gene2, "<-"), c(gene2, gene1, "->")),
      "|-|" = list(c(gene1, gene2, "|-"), c(gene2, gene1, "-|")),
      "|->" = list(c(gene1, gene2, "|-"), c(gene2, gene1, "->")),
      "<-|" = list(c(gene1, gene2, "<-"), c(gene2, gene1, "-|")),
      list(c(gene1, gene2, dir))
    )
  }
  expanded <- do.call(rbind, unlist(
    Map(split_one, fi$Gene1, fi$Gene2, fi$Direction),
    recursive = FALSE
  ))
  result <- as.data.frame(expanded, stringsAsFactors = FALSE, row.names = seq_len(nrow(expanded)))
  colnames(result) <- c("source_gene", "target_gene", "direction")
  result
}
