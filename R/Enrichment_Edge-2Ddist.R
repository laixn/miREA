#' Title: Edge-2Ddist
#' Description:
#' @param MGI2D A five-column dataframe, including MGI miRNA gene cor(x) ranknormratio(y).
#' @param pathway A three-column dataframe contains pathway miRNA gene.
#'@param background A character vector contains background MGIs in the format of [miRNA]:[gene].
#' @param pAdjMethod The method used for multiple testing correction to adjust p-values and control the false discovery rate.
#'        Use default.pAdjMethod to find all possible values.
#'        If not specified, Benjamini & Hochberg (BH) method will be automatically used.
#'        If you don't want any adjustment, please set pAdjMethod = "none".
#' @param pvalueCutoff The threshold for statistical significance, filtering out pathways with padj values below the specified pvalueCutoff. Default pvalueCutoff is 0.05.
#' @param minSize The minimum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs less than the minSize threshold.
#' @param maxSize The maximum MGI set sizes allowed for analysis, filtering out pathways with number of member MGIs more than the maxSize threshold.
#' @param iter Number of iterations for permutation test
#' @param ncores Number of cores for calculating randomized r during the permutation test.
#'
#' @return a dataframe contains enrichment result for Edge-2Ddist
#' @examples
#'

Edge_2Ddist <- function(MGI2D, pathway, background = NULL, pAdjMethod = "BH", pvalueCutoff = 0.05, minSize = NULL, maxSize = NULL, iter = 1000, ncores = NULL){
  cat("\n  Start Edge_2Ddist analysis ...\n")

  data.table::setDT(MGI2D)
  MGI2D <- na.omit(MGI2D)
  if (is.null(MGI2D)||!is.data.frame(MGI2D)){
    stop("Please make sure the input MGI2D is a dataframe!")
  }
  if (nrow(MGI2D) == 0){
    stop("No valid MGI pairs in MGI2D!")
  }

  if (is.null(background)){
    cat("    No background data is given, use the global background list instead. \n")
    background <- read.csv("data/raw_data/background/background_MGI.csv", header = TRUE)
    background <- paste0(background$miRNA, ":", background$gene)
  }

  # filtering MGI in the background
  if (all(MGI2D$MGI %in% background) == FALSE) {
    MGI2D <- MGI2D[MGI %in% background]
    if (nrow(MGI2D) == 0) {
      stop("No valid MGI input from background list!")
    }
  }
  colnames(MGI2D) <- c("MGI", "miRNA", "gene", "x", "y")

  r_all <- sqrt(((1 - MGI2D$x)/2)^2 + (MGI2D$y)^2)


  # filtering pathway in the background
  data.table::setDT(pathway)
  pathway[, MGI := paste0(miRNA, ":", gene)]
  pathway <- pathway[MGI %in% MGI2D$MGI]
  pathway <- pathway[MGI2D, on = .(MGI), nomatch = 0L]


  ### filtering pathway size
  pathway_counts <- table(pathway$pathway)
  if (!is.null(minSize)){
    pathways_to_remove_min <- names(pathway_counts[pathway_counts < minSize])
    pathway <- pathway[!pathway$pathway %in% pathways_to_remove_min, ]
  }
  if (!is.null(maxSize)){
    pathways_to_remove_max <- names(pathway_counts[pathway_counts > maxSize])
    pathway <- pathway[!pathway$pathway %in% pathways_to_remove_max, ]
  }

  n_pathway <- length(unique(pathway$pathway))
  if (n_pathway == 0) {
    warning("No pathways meet the size criteria.")
    #return(NULL)
    return(data.frame(
      pathway = character(), MGISize = integer(),
      r = numeric(), p_value = numeric(), padj = numeric(),
      nMoreExtreme = integer(), top10_MGI = character(),
      r_random = character(), stringsAsFactors = FALSE
    ))}

  if (is.null(ncores) || ncores == 1){
    ncores = 1
    # warning("You are use unparallel computing. It could cause huge amount of time!!!
    #         \nFor parallel computing, please specify [ncores] parameter to be larger than 1 (We highly recommend that)!!!")
  }

  # generate random pathways
  cat("    Generate random MGI set radius on", ncores, "cores...\n")
  setDT(pathway)
  setDT(MGI2D)

  # 按 pathway 大小分类
  pathway_sizes <- table(pathway$pathway)
  unique_sizes <- sort(unique(pathway_sizes))

  set.seed(12345)
  seeds_null = sample(100000, iter, replace = FALSE)

  compute_random_r <- function(r_all, size, i){
    set.seed(seeds_null[i])
    r <- sum(sample(r_all, size))
  }
  if (ncores == 1){
    r_random_list <- list()
    for (size in unique_sizes) {
      r_random_list[[as.character(size)]] <- unlist(lapply(1:iter, compute_random_r, size = size, r_all = r_all))
    }} else {
      random_core <- min(ncores, length(unique_sizes))
      cl <- parallel::makeCluster(random_core, type = "SOCK")
      parallel::clusterExport(cl, c("compute_random_r", "r_all", "iter", "seeds_null"), envir = environment())
      size_groups <- split(unique_sizes,
                           rep(1:random_core, length.out = length(unique_sizes)))
      seeds_split <- split(seeds_null, rep(1:random_core, length.out = length(seeds_null)))

      r_random_list <-  parallel::parLapply(
        cl,
        X = split(unique_sizes, rep(1:random_core, length.out = length(unique_sizes))),
        fun = function(size_group) {
          sub_result <- list()
          for (size in size_group) {
            r_values <- sapply(1:iter, function(i) compute_random_r(r_all = r_all, size = size, i))
            sub_result[[as.character(size)]] <- r_values
          }
          return(sub_result)
        }
      )
      parallel::stopCluster(cl)
      r_random_list <- unlist(r_random_list, recursive = FALSE)
      names(r_random_list) <- sub("^\\d+\\.", "", names(r_random_list))
    }
  r_random_list <- lapply(r_random_list, sort)

  cat("    Calculate pathway radius and siginificant level...\n")
  result <- Edge_2Ddist_simple(MGI2D = MGI2D, r_random_list = r_random_list, pathway = pathway, iter = iter)
  result <- result %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = pAdjMethod)) %>%
    dplyr::relocate(padj, .after = p_value) %>%
    dplyr::filter(
      p_value <= pvalueCutoff,
      padj <= pvalueCutoff) %>%
    arrange(padj, decrease = FALSE)

  cat("  Edge-2Ddist analysis have finished!\n")
  return(result)
}


Edge_2Ddist_simple <- function(MGI2D, r_random_list, pathway, iter){
  pathway_names <- unique(pathway$pathway)
  pathway_list <- split(pathway, pathway$pathway)
  ctx <- V8::v8()
  decimal_code <- paste(readLines("Javascript/decimal.js"), collapse = "\n")
  getP_code     <- paste(readLines("Javascript/GetpValueFromZ.js"), collapse = "\n")
  ctx$eval(decimal_code)
  ctx$eval(getP_code)
  compute_pathway <- function(i){
    pathway_name <- pathway_names[i]
    pw_single <- pathway_list[[pathway_name]]
    pw_size <- nrow(pw_single)
    r_random <- r_random_list[[as.character(pw_size)]]

    cat("      ", i, "/", length(pathway_names), ",", pathway_name, ":")

    result <- Edge_2Ddist_single(MGI2D = MGI2D, PW = pw_single, r_random = r_random, ctx = ctx)

    cat(round(result$p_value, 4), "\n")

    return(as.data.frame(result))
  }

  pathway_result <- bind_rows(lapply(seq_len(length(pathway_names)), compute_pathway))
  return(pathway_result)
}


Edge_2Ddist_single <- function(MGI2D, PW, r_random, ctx){
  MGI2D <- na.omit(MGI2D)
  if (is.null(MGI2D)||!is.data.frame(MGI2D)){
    stop("Please make sure the input MGI2D is a dataframe!")
  }
  if (nrow(MGI2D) == 0){
    stop("No valid MGI pairs in MGI2D!")
  }
  if (length(unique(PW$pathway)) != 1){
    warning("The pathway is not calculated one-by-one, please check the PW. Calculate the first one instead.")
    name <- unique(PW$pathway)[1]
    PW <- PW %>% dplyr::filter(pathway == name)
  }

  points_data <- semi_join(MGI2D, PW, by = c("miRNA", "gene"))

  if (nrow(points_data) == 0) {return(NULL)}


  points_data <- points_data %>%
    mutate(r = sqrt(((1-x)/2)^2 + y^2))

  r <- sum(points_data$r)
  count_leq <- findInterval(r, r_random) # count <= r
  nMoreExtreme <- length(r_random) - count_leq
  # p_value <- (nMoreExtreme + 1) / (length(r_random) + 1)
  # z <- (r - mean(r_random)) / sd(r_random)
  # p_value <- ctx$call("GetpValueFromZ", z, "right")

  if (is.na(sd(r_random)) || sd(r_random) != 0){
    z <- (r - mean(r_random)) / sd(r_random)
    p_value <- ctx$call("GetpValueFromZ", z, "right")
  } else {
    p_value <- 1.0000
  }


  threshold <- quantile(points_data$r, 0.9)
  top_pairs <- points_data %>%
    filter(r >= threshold) %>%
    mutate(MGI = paste0(miRNA,":", gene)) %>%
    pull(MGI)

  ## permutation test
  PW_name <- unique(PW$pathway)
  MGI_list <- PW$MGI
  size <- nrow(PW)

  random_ds <- round(quantile(r_random, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 4)
  random_mean <- round(mean(r_random, na.rm = TRUE), 4)
  random_ds <- c(
    min = random_ds[1],
    q25 = random_ds[2],
    median = random_ds[3],
    mean = random_mean,
    q75 = random_ds[4],
    max = random_ds[5]
  )

  result <- data.frame(
    "pathway" = PW_name,
    "MGISize" = size,
    "score" = r,
    "rand_score_dist" = paste(random_ds, collapse = ","),
    "p_value" = p_value,
    "nMoreExtreme" = nMoreExtreme,
    "TopMGI" = paste(top_pairs, collapse = "/"),
    # "rand_score_dist" = paste(round(r_random, 4), collapse = ","),
    stringsAsFactors = FALSE
  )
  return(result)

}

