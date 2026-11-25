# miREA 0.99.5

* Fixed a broken `@examples` block in `plotHeatmapSankey()`: a multi-line
  comment was missing a leading `#` on its last line, which R CMD check
  parsed as code and failed on with a syntax error (also affecting the R
  code problem check).
* `getAllInputData()` now always computes the per-MGI `strength` and `MGI`
  columns for any `Edge_` method, not only `Edge_Score1D`/`Edge_Score2D`/
  `Edge_Topology`. Previously, building `input_data` from `Edge_ORA` or
  `Edge_Network` alone produced a `data` table missing `strength`, which
  made `plotHeatmap()`/`plotHeatmapSankey()` error when called on such
  input_data. `scoreFun`'s default/validation now applies whenever any
  `Edge_` method is requested, matching the above.
* Standardized the differential-expression threshold comparison to be
  boundary-inclusive (`>=`, not `>`) everywhere it is used: the node-based
  and edge-based DE/DEMGI definitions in `getAllInputData()`, the
  DEMGI/DEG/DEmiRs filters and heatmap `*` DE-star annotations in
  `plotHeatmap()`/`plotHeatmapSankey()`. This fixes a prior inconsistency
  where the internal Sankey-building helper already used `>=` while every
  other DE check used strict `>`, which could make a boundary-value MGI
  count as DE in one plot panel but not the other within the same
  `plotHeatmapSankey()` call. Note this also changes the DE-calling rule
  used by `TG_ORA`/`TG_Score`/`MiR_ORA`/`Edge_ORA`/`Edge_Network`'s
  enrichment analyses themselves (not only plotting), so results for MGIs
  exactly at a cutoff boundary may differ from previous releases.
* Removed the unused, dead-code-only internal `fetch_background_GGI_reactome()`
  helper: `getAllInputData()`'s `GGI_source = "Reactome"` has always loaded
  the bundled `background_GGI_Reactome` snapshot directly and never called
  this function. Corrected the accompanying `background_GGI_Reactome`
  documentation, which incorrectly described a live Reactome re-fetch path.
* Removed dead default argument values from the internal `plot_sankey()`
  helper (its single call site always passed every argument explicitly);
  also removed its `n_pathway` parameter, which the function never used
  (the top-N-enriched-pathway restriction is already applied once, upstream,
  in `build_heatmap_components()`).
* Added regression tests for the above: `Edge_ORA`/`Edge_Network`-only
  `input_data` now carries `strength` and plots successfully; the DEMGI
  classification in `plotHeatmap()`/`plotHeatmapSankey()` is independently
  re-derived and checked against the documented inclusive (`>=`) rule.

# miREA 0.99.4

* Restored consistent 4-space indentation across the R source, fixing an
  accidental partial reformat to 2-space indentation introduced while
  preparing the previous release.
* Wrapped long roxygen documentation lines, vignette prose, and the
  `DESCRIPTION` `biocViews` field to at most 80 characters, addressing the
  BiocCheck line-length and indentation NOTEs. 

# miREA 0.99.3

* Addressed comments from the first round of Bioconductor review
  (bioconductor/BiocContributions#244).
* Fixed `getAllInputData()`'s `background_GGI` fallback: a broken `load()`
  call to a path that did not exist in the installed package was replaced with a live `OmnipathR` fetch or a bundled Reactome snapshot.
* Fixed RNG handling in `Edge_Score2D()`, `Edge_Topology()`, and
  `Edge_Network()`: removed a global-environment side effect from seed
  snapshotting, draw permutations from a single seeded stream instead of
  reseeding per iteration (via `BiocParallel`'s `RNGseed` in parallel
  runs), and give `Edge_Topology()` independent per-pathway seeds.
* Fixed `Edge_Network()` silently dropping zero-valued permutation scores
  from each pathway's null distribution, which had biased the mean, SD,
  and z-score, and could drop pathways from the output entirely.
* Fixed an `NA`/zero-SD branch condition in `Edge_Score2D()`/
  `Edge_Topology()` that produced `NA` p-values.
* Fixed `Edge_Score2D()`/`Edge_Topology()` mutating caller-supplied
  data.frames in place via `data.table::setDT()`.
* Added `p_empirical`, the empirical permutation p-value, alongside the
  existing normal-approximation `p_value`, and documented
  `Edge_Topology()`'s adaptive early stopping (the actual number of
  permutations run for a pathway can be lower than `iter`).
* Clarified data provenance and licensing: `toy_pathway` now uses Reactome
  pathways (CC0) instead of the previous MSigDB/KEGG-derived set, and
  `?miREA_background_data` documents the source database, version, and
  license of `background_MGI` and `background_GGI_Reactome`.
* Added a vignette section on supplying custom `background_MGI`/
  `background_GGI`, and fixed a documented `background = NULL` fallback
  that did not match `miREA()`'s actual behavior (it always errors on a
  missing background).
* Minor documentation fixes: removed a stray reference to a nonexistent
  `Edge_TopoScore` method, corrected the `pvalueCutoff` documentation
  (pathways above, not below, the cutoff are filtered out), removed an
  invalid `decrease` argument from an internal `arrange()` call, and fixed
  several typos.
* Switched parallel computing in `Edge_Score2D()`, `Edge_Topology()`, and
  `Edge_Network()` from `parallel::makeCluster()`/`parLapply()` to
  `BiocParallel::bplapply()` with a `BPPARAM` argument, removing the
  `RhpcBLASctl` dependency.
* Added `DElogFCCutoff`, `DEpadjCutoff`, `corCutoff`, and `corPadjCutoff`
  arguments to `getAllInputData()` for thresholds that were previously
  hard-coded.
* `getAllInputData()` now accepts `DESeq2::results()`/`limma::topTable()`
  output directly as differential expression input, matched by column
  name instead of a fixed 4-column position, and its internal
  `coerce_expr_mat()` gained an `assay` argument for `SummarizedExperiment`
  input.
* Replaced several `warning()`/`stop(paste(...))` calls with `message()`/
  `stop(...)`.
* `miREAResult` gained a `setValidity()` method and slot-name validation
  in `$`/`[[`, and now records methods that failed, with their error
  messages, instead of only emitting a warning.
* Fixed `plotSummary()` returning a message string instead of its plots
  invisibly when saving to `plot_path`.
* Removed unused internal functions `ranknorm_pos()` and `sigmoid_pos()`.
* Reduced the size of the bundled toy datasets to the 1,142 genes and 790
  miRNAs actually used by `toy_pathway`.
* Added tests for seed reproducibility (the same seed gives identical
  results, `seed = NULL` gives different results across runs), input
  immutability, and the `background_GGI = NULL` fallback path. 
* Unit test coverage is now 92.4% (`covr::package_coverage()`).
* Noted the recommended `iter` (>= 1000) and the expected runtime of
  `Edge_Topology()`/`Edge_Network()` on real-sized data in the vignette.
* Added nine additional vignettes, one per enrichment method, alongside
  the existing Quick Start vignette.

# miREA 0.99.0

* Initial Bioconductor submission.
* Four node-based and five edge-based miRNA-oriented enrichment methods:
  `TG_ORA`, `TG_Score`, `MiR_ORA`, `MiR_Score`, `Edge_ORA`, `Edge_Score1D`,
  `Edge_Score2D`, `Edge_Topology`, `Edge_Network`.
* `getAllInputData()` to prepare method-specific input data from raw
  expression and differential-expression tables.
* `plotSummary()`, `plotHeatmap()`, and `plotHeatmapSankey()` for
  visualizing and comparing enrichment results.
