## Regenerates the full (non-toy) background_MGI/background_GGI_Omnipath/
## background_GGI_Reactome datasets shipped in data/. These are the complete
## manuscript backgrounds (not the toy-universe subsets in
## generate_toy_data.R). Not run during R CMD build/check; re-run manually
## if the bundled backgrounds need to be refreshed. Requires access to the
## raw CSVs under miREA_zenodo/data/raw_data/background/ (project root, not
## part of this package).

raw_dir <- "/scratch/project_2011179/code/miREA_zenodo/data/raw_data/background"

background_MGI <- read.csv(file.path(raw_dir, "background_MGI.csv"), stringsAsFactors = FALSE)
colnames(background_MGI) <- c("miRNA", "gene")

background_GGI_Omnipath <- read.csv(file.path(raw_dir, "background_GGI_Omnipath.csv"), stringsAsFactors = FALSE)
colnames(background_GGI_Omnipath) <- c("source_genesymbol", "target_genesymbol", "interact")

background_GGI_Reactome <- read.csv(file.path(raw_dir, "background_GGI_Reactome.csv"), stringsAsFactors = FALSE)
colnames(background_GGI_Reactome) <- c("source_gene", "target_gene", "direction")

usethis::use_data(background_MGI, overwrite = TRUE, compress = "xz")
usethis::use_data(background_GGI_Omnipath, overwrite = TRUE, compress = "xz")
usethis::use_data(background_GGI_Reactome, overwrite = TRUE, compress = "xz")
