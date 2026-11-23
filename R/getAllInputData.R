# if TG_method, initially use the differential expressed miRNAs (interested miRNAs), and link to genes (as interested genes, or geneList (need gene_data)).
# if MiR_method, use mir_DEdata directly
# if Edge_method, need both mir_DEdata and gene_DEdata

#' Generate input data for miREA
#'
#' Functions to generate input data for miREA, based on the original expression and differential expression results for miRNAs and genes.
#'
#' @param methods compulsory, see default_method_list for all supportive methods. You could either input one method, or a vector contains several methods you want to try.
#' @param pathway compulsory. It should be a two-column pathway, containing pathway and gene, or a GSEABase GeneSetCollection object.
#' @param mir_DEdata compulsory. A dataframe of differential expression results, with columns matched
#'        by name (not position), so extra columns are fine. Recognized logFC columns:
#'        `log2FoldChange`/`logFC`; stat columns: `stat`/`t`; padj columns: `padj`/`adj.P.Val`
#'        (i.e. `DESeq2::results()` or `limma::topTable()` output can be passed directly). The miRNA
#'        ID is read from a `miRNA`/`ID`/`SYMBOL`/`Row.names` column if present, otherwise from
#'        `rownames(mir_DEdata)` (the default for `DESeq2::results()`/`limma::topTable()` output).
#' @param gene_DEdata optional, except for TG-ORA, MiR-ORA, and MiR-Score methods. Accepts the same
#'        shapes as `mir_DEdata` (with `gene` instead of `miRNA` as the ID).
#' @param background_MGI optional. It should be either a character vector in the format of [miRNA]:[gene], or a two-column dataframe, containing miRNA gene. If NULL, miREA falls back to its bundled `background_MGI` dataset (the full manuscript background).
#' @param background_GGI optional, it should be a two-column data frame that contains source gene and target gene. If supplied, GGI_source is ignored (a message notes this). If NULL and GGI_source = "Omnipath" (or GGI_source = NULL, which is treated the same as "Omnipath"), miREA fetches Omnipath's full interaction set live (requires the OmnipathR package and an internet connection), with no license/source filtering (see the "Third-party data licensing" section of the package README for why, and how to supply a license-restricted subset yourself). If NULL and GGI_source = "Reactome", miREA falls back to its bundled `background_GGI_Reactome` dataset (a frozen CC0-licensed snapshot, no commercial-use restriction).
#' @param GGI_source optional, select from "Omnipath" or "Reactome", or NULL. Selects how background_GGI is obtained when not provided directly; NULL behaves like "Omnipath". The default is "Omnipath".
#' @param gene_mat optional, needed for Edge_ method, used for calculating correlations between miRNAs and genes. It should be a raw gene expression matrix (each row a gene, each column a sample), either as a dataframe/matrix or as a SummarizedExperiment object.
#' @param mir_mat optional, needed for Edge_ method, used for calculating correlations between miRNAs and genes. It should be a raw miRNA expression matrix (each row a miRNA, each column a sample), either as a dataframe/matrix or as a SummarizedExperiment object.
#' @param gene_assay optional, only used when `gene_mat` is a SummarizedExperiment. Selects which assay to
#'        extract, either by integer index or assay name (passed as `i` to [SummarizedExperiment::assay()]).
#'        Default is 1 (the first assay).
#' @param mir_assay optional, only used when `mir_mat` is a SummarizedExperiment. Selects which assay to
#'        extract, either by integer index or assay name. Default is 1 (the first assay).
#' @param scoreFun optional, needed for Edge_Score1D, Edge_Score2D, Edge_Topology. See default_scoreFun for all supportive methods.
#' @param DElogFCCutoff optional, the absolute log2FC threshold used to call a miRNA/gene differentially
#'        expressed when building the DE sets for the ORA methods and the plotting functions. Default is 1.
#' @param DEpadjCutoff optional, the adjusted p-value threshold used to call a miRNA/gene differentially
#'        expressed when building the DE sets for the ORA methods and the plotting functions. Default is 0.05.
#' @param corCutoff optional, the miRNA-gene expression correlation threshold (MGIs with `cor < corCutoff`
#'        are considered) used when building the DE sets for the ORA methods and the plotting functions.
#'        Default is 0.
#' @param corPadjCutoff optional, the adjusted p-value threshold for the miRNA-gene expression correlation,
#'        used when building the DE sets for the ORA methods and the plotting functions. Default is 0.05.
#' @note gene_mat and mir_mat should have same sample name format. If the number of shared samples below 50 % of the minimum of gene samples and miRNA samples, a warning would happen to remind you recheck them.
#' @return A list that contains all the data formats needed for miREA, including a `threshold` element
#'        recording the `DElogFCCutoff`/`DEpadjCutoff`/`corCutoff`/`corPadjCutoff` values actually used.
#' @references Zhang Z, Lai X (2026). miREA: a network-based tool for
#' microRNA-oriented enrichment analysis. bioRxiv. \doi{10.64898/2026.02.27.708509}
#' @examples
#' data(toy_pathway)
#' data(toy_gene_DEdata)
#' data(toy_mir_DEdata)
#' data(toy_gene_mat)
#' data(toy_mir_mat)
#' data(background_MGI)
#' input_data <- getAllInputData(
#'     methods = default_method_list,
#'     pathway = toy_pathway,
#'     mir_DEdata = toy_mir_DEdata,
#'     gene_DEdata = toy_gene_DEdata,
#'     background_MGI = background_MGI,
#'     GGI_source = "Reactome",
#'     gene_mat = toy_gene_mat,
#'     mir_mat = toy_mir_mat,
#'     scoreFun = "rank"
#' )
#' names(input_data)
#' @export
getAllInputData <- function(methods, pathway,
                            mir_DEdata, gene_DEdata = NULL,
                            background_MGI = NULL, background_GGI = NULL,
                            GGI_source = "Omnipath",
                            gene_mat = NULL, mir_mat = NULL, gene_assay = 1, mir_assay = 1,
                            scoreFun = NULL,
                            DElogFCCutoff = 1, DEpadjCutoff = 0.05,
                            corCutoff = 0, corPadjCutoff = 0.05) {
    message("Generate all input data...")
    pathway <- coerce_pathway(pathway)
    gene_mat <- coerce_expr_mat(gene_mat, assay = gene_assay)
    mir_mat <- coerce_expr_mat(mir_mat, assay = mir_assay)
    sum_TG <- sum(grepl("TG_", methods))
    sum_MiR <- sum(grepl("MiR_", methods))
    sum_Edge <- sum(grepl("Edge_", methods))
    target_method <- c("Edge_Topology", "Edge_Network")
    score_method <- c("Edge_Score1D", "Edge_Score2D", "Edge_Topology")
    colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
    # check scoreFun
    if (length(intersect(score_method, methods)) > 0 && is.null(scoreFun)) {
        message("No [scoreFun] parameter has been specified. Use the default [rank] instead!\n")
        scoreFun <- "rank"
    }
    if (length(intersect(score_method, methods)) > 0 && (!scoreFun %in% default_scoreFun)) {
        stop("The [scoreFun] you input is not valid. Please select one from the default_scoreFun!")
    }

    result <- list()

    # check background_MGI
    # NOTE: nearly every method (all but MiR_Score) derives its input via
    # background_MGI somewhere below (e.g. TG_ORA maps significant miRNAs to
    # target genes through it), so this is required unconditionally rather
    # than only for Edge_/MiR_ methods.
    if (is.null(background_MGI)) {
        message("No background_MGI is given. Using the bundled default background_MGI instead.\n")
        utils::data("background_MGI", package = "miREA", envir = environment())
    }
    colnames(background_MGI)[seq_len(2)] <- c("miRNA", "gene")

    # check background_GGI
    if (length(intersect(target_method, methods)) != 0) {
        if (!is.null(background_GGI)) {
            message("background_GGI is supplied directly; GGI_source is ignored.\n")
        } else {
            # GGI_source = NULL means "use the default source", i.e. Omnipath.
            if (is.null(GGI_source)) GGI_source <- "Omnipath"
            if (!GGI_source %in% default_GGI_source) {
                stop("You tend to extract default background_GGI for ", GGI_source, ". But [GGI_source] you have specified are not involved. Please use default_GGI_source to see all supporting options!")
            }
            if (GGI_source == "Omnipath") {
                message("No background_GGI is given. Fetching Omnipath's full interaction set live (no license/source filtering).\n")
                background_GGI <- fetch_background_GGI_omnipath()
            } else {
                bundled_env <- new.env()
                message("No background_GGI is given. Using the bundled default background_GGI_Reactome instead.\n")
                utils::data("background_GGI_Reactome", package = "miREA", envir = bundled_env)
                background_GGI <- bundled_env[["background_GGI_Reactome"]]
            }
        }
        colnames(background_GGI)[seq_len(2)] <- c("from", "to")
    }

    # check expression data
    if (!is.null(gene_DEdata)) {
        gene_DEdata <- standardize_DEdata(gene_DEdata, "gene")
        gene_DEdata <- gene_DEdata %>% distinct(gene, .keep_all = TRUE)
        gene_DEdata <- data.table::as.data.table(gene_DEdata)
    }

    mir_DEdata <- standardize_DEdata(mir_DEdata, "miRNA")
    mir_DEdata <- mir_DEdata %>% distinct(miRNA, .keep_all = TRUE)
    mir_DEdata <- data.table::as.data.table(mir_DEdata)

    # ensure that mir_mat and gene_mat have their names as rownames
    if (!is.null(mir_mat) && !is.null(gene_mat)) {
        aligned <- align_expr_mat_samples(mir_mat, gene_mat, rename_id_cols = TRUE, dedup = TRUE, check_samples = FALSE)
        mir_mat <- aligned$mir_mat
        gene_mat <- aligned$gene_mat
    }


    raw_data <- list(
        methods = methods,
        mir_DEdata = mir_DEdata, gene_DEdata = gene_DEdata,
        gene_mat = gene_mat, mir_mat = mir_mat,
        scoreFun = scoreFun
        # check_mir_names = check_mir_names, mir_version = mir_version,
        # check_gene_names = check_gene_names
    )
    result[["raw_data"]] <- raw_data
    result[["background_MGI"]] <- background_MGI
    result[["background_GGI"]] <- background_GGI


    message("  Start generate expression data...")
    data <- get_data(
        methods = methods, mir_DEdata = mir_DEdata, gene_DEdata = gene_DEdata,
        background_MGI = background_MGI, gene_mat = gene_mat, mir_mat = mir_mat, scoreFun = scoreFun,
        DElogFCCutoff = DElogFCCutoff, DEpadjCutoff = DEpadjCutoff,
        corCutoff = corCutoff, corPadjCutoff = corPadjCutoff
    )
    result[["data"]] <- data
    result[["threshold"]] <- list(
        DElogFCCutoff = DElogFCCutoff, DEpadjCutoff = DEpadjCutoff,
        corCutoff = corCutoff, corPadjCutoff = corPadjCutoff
    )

    message("  Start generate pathway data...")
    pathway_data <- get_pathway_data(methods = methods, pathway = pathway, background_MGI = background_MGI)
    result[["pathway"]] <- pathway_data

    if (length(intersect(target_method, methods)) != 0) {
        message("  Start generate pathway_GGI data...")
        pathway_GGI <- get_pathway_GGI_data(methods = methods, pathway = pathway, background_GGI = background_GGI)
        result[["pathway_GGI"]] <- pathway_GGI

        if ("Edge_Network" %in% methods) {
            message("  Start generate GGI data for whole network...")
            GGI <- get_GGI_data(methods = methods, DEMGI = result[["data"]][["Edge_Network"]], pathway = pathway, background_GGI = background_GGI)
            result[["GGI"]] <- GGI
            message("  Start generate MGI data for whole network...")
            MGI <- get_MGI_data(
                methods = methods, DEMGI = result[["data"]][["Edge_Network"]], pathway = pathway,
                pathway_MGI = result[["pathway"]][["Edge"]], background_MGI = background_MGI
            )
            result[["MGI"]] <- MGI
        }

        if ("Edge_Topology" %in% methods) {
            message("  Start generate dict_MGI data for permutation...")
            dict_MGI <- get_dict_MGI(
                methods = methods, mir_mat = mir_mat, gene_mat = gene_mat, data = data$data,
                pathway_MGI = pathway_data$Edge, scoreFun = scoreFun
            )
            result[["dict_MGI"]] <- dict_MGI
        }
    }

    message("All the input data is ready!")
    return(result)
}

## get all input data, including:
### data, pathway, pathway_GGI, dict_MGI, dict_GGI
## all the original input include:
### method, gene_DEdata, mir_DEdata, pathway, gene_mat, mir_mat, background_MGI, background_GGI, scoreFun
### Notice: gene_mat, mir_mat should only conclude tumor samples


# gene_mat: Edge_ORA
# mir_mat: Edge_ORA
# gene_DEdata: TG_ORA, TG_Score
# mir_DEdata: MiR_ORA, MiR_Score
# background_MGI: could be either two-column dataframe, or a character vector with format [miRNA]:[gene]
#                 Edge_ORA, Edge_Score1D, Edge_Topology, Edge_Network
# background_GGI
# pathway
# scoreFun: Edge_Score1D

# gene_DEdata: except for TG_ORA, MiR_ORA and MiR_Score, should be four-column dataframe, containing gene, log2FC, stat, padj
# mir_DEdata: compulsory, should be four-column dataframe, containing miRNA, log2FC, stat, padj
# gene_mat: only needed for Edge_ method
# mir_mat: only needed for Edge_ method
# scoreFun: only needed for Edge_Score1D, Edge_Score2D, Edge_Topology, choose one from rank, sigmoid, normalize
# background_MGI: only needed for Edge_ method
# will return a list, contains the input of each method
get_data <- function(methods, mir_DEdata, gene_DEdata = NULL,
                     background_MGI, gene_mat = NULL, mir_mat = NULL, scoreFun,
                     DElogFCCutoff = 1, DEpadjCutoff = 0.05,
                     corCutoff = 0, corPadjCutoff = 0.05) {
    scoreFun <- get_data_validate_inputs(methods, gene_DEdata, gene_mat, mir_mat, background_MGI, scoreFun)

    coerced <- get_data_coerce_expr(gene_DEdata, mir_DEdata)
    gene_DEdata <- coerced$gene_DEdata
    mir_DEdata <- coerced$mir_DEdata
    background_MGI <- get_data_coerce_background_mgi(background_MGI)

    result <- get_data_node_based(methods, mir_DEdata, gene_DEdata, background_MGI, DElogFCCutoff, DEpadjCutoff)

    if (sum(grepl("Edge_", methods)) > 0) {
        result <- c(result, get_data_edge_based(
            methods, mir_DEdata, gene_DEdata, background_MGI, mir_mat, gene_mat, scoreFun,
            DElogFCCutoff, DEpadjCutoff, corCutoff, corPadjCutoff
        ))
    }
    result
}

#' @noRd
get_data_validate_inputs <- function(methods, gene_DEdata, gene_mat, mir_mat, background_MGI, scoreFun) {
    invalid_method <- base::setdiff(methods, default_method_list)
    if (length(invalid_method) > 0) {
        stop("The method(s):", invalid_method, "you are input is(are) not one of the supportive methods. Please revise the method(s)!")
    }
    sum_Edge <- sum(grepl("Edge_", methods))
    no_gene_DEdata <- c("TG_ORA", "MiR_ORA", "MiR_Score")
    score_method <- c("Edge_Score1D", "Edge_Score2D", "Edge_Topology")
    if (length(setdiff(methods, no_gene_DEdata)) > 0 && is.null(gene_DEdata)) {
        stop("Please make sure you have input gene_DEdata!")
    }
    if (sum_Edge != 0 && (is.null(gene_mat) || is.null(mir_mat))) {
        stop("Please make sure you have input gene_mat and mir_mat")
    }
    if (any(score_method %in% methods) && is.null(scoreFun)) {
        message("You didn't specify scoreFun for transform expression score, set 'rank' automatically.\n")
        scoreFun <- "rank"
    }
    # NOTE: unreachable in practice — getAllInputData() resolves background_MGI
    # (falling back to the bundled default) before calling get_data(), which is
    # the only caller of this function. Kept as a defensive safety net.
    if (is.null(background_MGI)) {
        stop(
            "background_MGI must be supplied: a two-column data.frame of miRNA, gene ",
            "pairs, or use the bundled default background_MGI via getAllInputData()."
        )
    }
    scoreFun
}

#' Standardize a differential-expression table to miREA's internal column layout
#'
#' Columns are always matched by name, never by position. Recognizes miREA's own canonical
#' column names (so this is idempotent on an already-standardized table), the column names
#' produced by `DESeq2::results()` (`log2FoldChange`, `stat`, `padj`), and those produced by
#' `limma::topTable()` (`logFC`, `t`, `adj.P.Val`). `df` may have any number of extra columns
#' (e.g. `baseMean`/`lfcSE` from DESeq2, `AveExpr`/`B` from limma) -- they are ignored. If no ID
#' column is recognized, `rownames(df)` is used instead (this is how `as.data.frame(DESeq2::
#' results(dds))` and `limma::topTable()` return their output by default). Errors with a message
#' listing the recognized names if logFC/stat/padj can't be matched.
#' @param df A data.frame of DE results for `entity`.
#' @param entity Either `"gene"` or `"miRNA"`.
#' @return A four-column data.frame named `c(entity, paste0(entity, c("LogFC","Stat","Padj")))`.
#' @noRd
standardize_DEdata <- function(df, entity) {
    if (!is.data.frame(df)) {
        stop("Please make sure the input ", entity, "_DEdata is a dataframe.")
    }

    # miREA's internal convention prefixes the logFC/stat/padj columns with "mir" for
    # the miRNA entity (not "miRNA"), e.g. mirLogFC/mirStat/mirPadj; "gene" uses itself.
    prefix <- if (entity == "miRNA") "mir" else entity
    id_col <- entity
    logFC_col <- paste0(prefix, "LogFC")
    stat_col <- paste0(prefix, "Stat")
    padj_col <- paste0(prefix, "Padj")

    # each candidate list also includes miREA's own canonical output name, so that
    # standardize_DEdata() is idempotent when applied to an already-standardized table
    # (it is called a second time inside get_data_coerce_expr()).
    id_candidates <- c(id_col, "ID", "SYMBOL", "Row.names")
    logFC_candidates <- c(logFC_col, "log2FoldChange", "logFC") # DESeq2, limma
    stat_candidates <- c(stat_col, "stat", "t") # DESeq2, limma
    padj_candidates <- c(padj_col, "padj", "adj.P.Val") # DESeq2, limma

    matched_logFC <- intersect(logFC_candidates, colnames(df))
    matched_stat <- intersect(stat_candidates, colnames(df))
    matched_padj <- intersect(padj_candidates, colnames(df))

    if (length(matched_logFC) == 0 || length(matched_stat) == 0 || length(matched_padj) == 0) {
        stop(
            "Please make sure the input ", entity, "_DEdata has recognizable logFC/stat/padj columns ",
            "(tried logFC: ", paste(logFC_candidates, collapse = ", "), "; ",
            "stat: ", paste(stat_candidates, collapse = ", "), "; ",
            "padj: ", paste(padj_candidates, collapse = ", "), "). ",
            "Columns found: ", paste(colnames(df), collapse = ", "), "."
        )
    }

    matched_id <- intersect(id_candidates, colnames(df))
    if (length(matched_id) >= 1) {
        id_values <- df[[matched_id[1]]]
    } else if (!is.null(rownames(df)) && !identical(rownames(df), as.character(seq_len(nrow(df))))) {
        id_values <- rownames(df)
    } else {
        stop(
            "Please make sure the input ", entity, "_DEdata has an ID column (tried: ",
            paste(id_candidates, collapse = ", "), ") or informative rownames."
        )
    }

    out <- data.frame(id_values, df[[matched_logFC[1]]], df[[matched_stat[1]]], df[[matched_padj[1]]],
        stringsAsFactors = FALSE
    )
    colnames(out) <- c(id_col, logFC_col, stat_col, padj_col)
    out
}

#' @noRd
get_data_coerce_expr <- function(gene_DEdata, mir_DEdata) {
    if (!is.null(gene_DEdata)) {
        gene_DEdata <- standardize_DEdata(gene_DEdata, "gene")
        gene_DEdata <- data.table::as.data.table(gene_DEdata)
    }

    mir_DEdata <- standardize_DEdata(mir_DEdata, "miRNA")
    mir_DEdata <- data.table::as.data.table(mir_DEdata)

    list(gene_DEdata = gene_DEdata, mir_DEdata = mir_DEdata)
}

#' @noRd
get_data_coerce_background_mgi <- function(background_MGI) {
    if (is.data.frame(background_MGI)) {
        if (ncol(background_MGI) != 2) {
            stop("Please make sure the input background_MGI is a two-column dataframe contains miRNA and gene")
        }
        colnames(background_MGI) <- c("miRNA", "gene")
    } else if (is.character(background_MGI) && is.vector(background_MGI)) {
        background_MGI <- data.frame(
            miRNA = sub(":.*", "", background_MGI),
            gene = sub(".*:", "", background_MGI),
            stringsAsFactors = FALSE
        )
    } else {
        stop("Please make sure the input background_MGI is either a two-column dataframe contains miRNA and gene,
         or a character vector with each element in the format of [miRNA]:[gene]!")
    }
    colnames(background_MGI) <- c("miRNA", "gene")
    data.table::as.data.table(background_MGI)
}

#' @noRd
get_data_node_based <- function(methods, mir_DEdata, gene_DEdata, background_MGI,
                                DElogFCCutoff = 1, DEpadjCutoff = 0.05) {
    result <- list()
    if (sum(methods %in% c("TG_ORA", "TG_Score", "MiR_ORA")) > 0) {
        interested_miRNAs <- mir_DEdata[abs(mirLogFC) > DElogFCCutoff & mirPadj < DEpadjCutoff, unique(miRNA)]
    }

    if ("TG_ORA" %in% methods) {
        result[["TG_ORA"]] <- background_MGI[miRNA %in% interested_miRNAs, unique(gene)]
    }
    if ("TG_Score" %in% methods) {
        genes <- background_MGI[miRNA %in% interested_miRNAs, unique(gene)]
        data <- gene_DEdata[gene %in% genes, setNames(geneStat, gene)]
        result[["TG_Score"]] <- sort(data, decreasing = TRUE)
    }
    if ("MiR_ORA" %in% methods) {
        result[["MiR_ORA"]] <- interested_miRNAs
    }
    if ("MiR_Score" %in% methods) {
        data <- mir_DEdata[, setNames(mirStat, miRNA)]
        result[["MiR_Score"]] <- sort(data, decreasing = TRUE)
    }
    result
}

#' @noRd
get_data_edge_based <- function(methods, mir_DEdata, gene_DEdata, background_MGI, mir_mat, gene_mat, scoreFun,
                                DElogFCCutoff = 1, DEpadjCutoff = 0.05,
                                corCutoff = 0, corPadjCutoff = 0.05) {
    DE_method <- c("Edge_ORA", "Edge_Network")
    score_method <- c("Edge_Score1D", "Edge_Score2D", "Edge_Topology")
    result <- list()

    dict <- merge(background_MGI, mir_DEdata[, .(miRNA, mirLogFC, mirPadj)], by = "miRNA", all = FALSE)
    dict <- merge(dict, gene_DEdata[, .(gene, geneLogFC, genePadj)], by = "gene", all = FALSE)
    cor <- data.table::as.data.table(calc_MGI_corr(dict = dict, mir_mat = mir_mat, gene_mat = gene_mat))
    dict <- merge(dict, cor, by = c("miRNA", "gene"), all = FALSE)
    dict <- na.omit(dict)
    data.table::setcolorder(dict, c("miRNA", "gene", "cor", "corPadj", "mirLogFC", "mirPadj", "geneLogFC", "genePadj"))
    result[["data"]] <- dict

    if (any(DE_method %in% methods)) {
        sel <- dict[
            cor < corCutoff & corPadj < corPadjCutoff &
                abs(mirLogFC) > DElogFCCutoff & mirPadj < DEpadjCutoff &
                abs(geneLogFC) > DElogFCCutoff & genePadj < DEpadjCutoff
        ]
        data <- if (nrow(sel) > 0) paste0(sel$miRNA, ":", sel$gene) else character(0)
        for (sub_method in intersect(DE_method, methods)) {
            result[[sub_method]] <- data
        }
    }
    if (any(score_method %in% methods)) {
        dict[, ratio := (geneLogFC / mirLogFC) * ((1 - genePadj) / (1 + mirPadj))]
        dict[, normratio := apply_scoreFun(ratio, scoreFun)]
        dict[, `:=`( # `:=`() form modifies multiple columns; infix `:=` modifies a single column
            strength = cor * normratio,
            MGI = paste0(miRNA, ":", gene)
        )]
        data.table::setorder(dict, strength)
        dict <- unique(dict, by = "MGI")

        if ("Edge_Score1D" %in% methods) {
            result[["Edge_Score1D"]] <- setNames(as.numeric(dict$strength), dict$MGI)
        }
        if ("Edge_Topology" %in% methods) {
            result[["Edge_Topology"]] <- as.data.frame(dict[, .(MGI, miRNA, gene, strength)])
        }
        if ("Edge_Score2D" %in% methods) {
            result[["Edge_Score2D"]] <- as.data.frame(dict[, .(MGI, miRNA, gene, cor, normratio)])
        }
    }
    result
}


get_pathway_data <- function(methods, pathway, background_MGI) {
    colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
    sum_TG <- sum(grepl("TG_|^Edge_Network$", methods)) # Edge-Network needs pathway gene set to map GGI.
    sum_MiR <- sum(grepl("MiR_", methods))
    sum_Edge <- sum(grepl("Edge_", methods))

    result <- list()

    if (sum_TG > 0) {
        result[["TG"]] <- pathway
    }
    if (sum_MiR + sum_Edge > 0) {
        colnames(background_MGI)[seq_len(2)] <- c("miRNA", "gene")
        background_MGI <- data.table::as.data.table(background_MGI)
        pathway <- data.table::as.data.table(pathway)

        pathway_genes <- pathway[, .(genes = list(gene)), by = pathway]
        result_list <- lapply(seq_len(nrow(pathway_genes)), function(i) {
            pw <- pathway_genes$pathway[i]
            genes <- pathway_genes$genes[[i]]
            background_MGI[gene %in% genes, .(pathway = pw, miRNA, gene)]
        })
        pathway_MGI <- unique(rbindlist(result_list))
        if (sum_MiR > 0) {
            result[["MiR"]] <- unique(pathway_MGI[, .(pathway, miRNA)])
        }
        if (sum_Edge > 0) {
            result[["Edge"]] <- pathway_MGI
        }
    }
    return(result)
}


# only useful for Edge_Topology, Edge_Network
get_pathway_GGI_data <- function(methods, pathway, background_GGI) {
    target_method <- c("Edge_Network", "Edge_Topology")
    if (length(intersect(target_method, methods)) == 0) {
        message("The methods you have selected don't need gene-gene interactions data. Skip the pathway_GGI process.\n")
        return(NULL)
    } else {
        colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
        colnames(background_GGI)[seq_len(2)] <- c("from", "to")
        background_GGI <- data.table::as.data.table(background_GGI)
        pathway <- data.table::as.data.table(pathway)

        pathway_genes <- pathway[, .(genes = list(gene)), by = pathway]
        result_list <- vector("list", nrow(pathway_genes))

        result_list <- lapply(seq_len(nrow(pathway_genes)), function(i) {
            pw <- pathway_genes$pathway[i]
            genes <- pathway_genes$genes[[i]]

            background_GGI[
                from %in% genes & to %in% genes,
                .(pathway = pw, from, to)
            ] # , type = "GGI"
        })

        pathway_GGI <- unique(rbindlist(result_list))
        return(pathway_GGI)
    }
}

# DEMGI: result[["data"]][["Edge_Network"]]
get_GGI_data <- function(methods, DEMGI, pathway, background_GGI) { # gene_mat, gene_DEdata,
    if (!"Edge_Network" %in% methods) {
        message("The methods you have selected don't need global GGI data. Skip the get_GGI_data process.\n")
        return(NULL)
    } else {
        DEMGI <- data.frame(
            miRNA = sub(":.*", "", DEMGI),
            gene = sub(".*:", "", DEMGI),
            stringsAsFactors = FALSE
        )

        colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
        colnames(background_GGI)[seq_len(2)] <- c("from", "to")
        background_GGI <- data.table::as.data.table(background_GGI)
        pathway <- data.table::as.data.table(pathway)
        DEMGI <- data.table::as.data.table(DEMGI)

        pathway_genes <- unique(pathway$gene)
        DEMGI_genes <- unique(DEMGI$gene)
        global_genes <- unique(c(pathway_genes, DEMGI_genes))
        global_GGI <- background_GGI[
            from %in% global_genes & to %in% global_genes,
            .(from = from, to = to)
        ]
        global_GGI <- as.data.frame(global_GGI)
        return(global_GGI)
    }
}

get_MGI_data <- function(methods, DEMGI, pathway, pathway_MGI, background_MGI) {
    if (!"Edge_Network" %in% methods) {
        message("The methods you have selected don't need global MGI data. Skip the get_MGI_data process.\n")
        return(NULL)
    } else {
        DEMGI <- data.frame(
            miRNA = sub(":.*", "", DEMGI),
            gene = sub(".*:", "", DEMGI),
            stringsAsFactors = FALSE
        )

        colnames(pathway)[seq_len(2)] <- c("pathway", "gene")
        colnames(background_MGI)[seq_len(2)] <- c("miRNA", "gene")
        colnames(pathway_MGI)[seq_len(3)] <- c("pathway", "miRNA", "gene")
        background_MGI <- data.table::as.data.table(background_MGI)
        pathway <- data.table::as.data.table(pathway)
        pathway_MGI <- data.table::as.data.table(pathway_MGI)
        DEMGI <- data.table::as.data.table(DEMGI)

        pathway_genes <- unique(pathway$gene)
        pathway_mirs <- unique(pathway_MGI$miRNA)
        DEMGI_genes <- unique(DEMGI$gene)
        DEMGI_mirs <- unique(DEMGI$miRNA)
        global_genes <- unique(c(pathway_genes, DEMGI_genes))
        global_mirs <- unique(c(pathway_mirs, DEMGI_mirs))

        global_MGI <- background_MGI[
            gene %in% global_genes & miRNA %in% global_mirs,
            .(miRNA = miRNA, gene = gene)
        ]
        global_MGI <- as.data.frame(global_MGI)
        return(global_MGI)
    }
}


# useful in Edge_Topology and Edge_Network
# data is data$data get in the get_data() function
get_dict_MGI <- function(methods, mir_mat, gene_mat, data, pathway_MGI, scoreFun) {
    if (!"Edge_Topology" %in% methods) {
        message("The methods you have selected don't need dict_MGI data. Skip the get_dict_MGI process.\n")
        return(NULL)
    } else {
        aligned <- align_expr_mat_samples(mir_mat, gene_mat)
        mir_mat <- aligned$mir_mat
        gene_mat <- aligned$gene_mat
        common_samples <- aligned$common_samples

        colnames(pathway_MGI)[seq_len(3)] <- c("pathway", "miRNA", "gene")
        pathway_MGI <- data.table::as.data.table(pathway_MGI)

        data <- data.table::as.data.table(data)

        # build mat that have both miRNA and gene with only shared samples.
        mir_sub <- mir_mat[, common_samples, drop = FALSE]
        gene_sub <- gene_mat[, common_samples, drop = FALSE]

        # build dict dataframe: all possible combinations between miRNAs and genes,
        # who are existed in pathway and meanwhile have valid expresison profiles.
        unique_miRNAs <- intersect(unique(pathway_MGI$miRNA), unique(data$miRNA))
        unique_genes <- intersect(unique(pathway_MGI$gene), unique(data$gene))
        dict <- data.table::CJ(miRNA = unique_miRNAs, gene = unique_genes)

        # merge with gene and miRNA differential expressed data
        mir_info <- unique(data[, .(miRNA, mirLogFC, mirPadj)])
        gene_info <- unique(data[, .(gene, geneLogFC, genePadj)])

        dict <- merge(dict, mir_info, by = "miRNA", all.x = TRUE)
        dict <- merge(dict, gene_info, by = "gene", all.x = TRUE)

        results_mat <- calc_MGI_corr(dict = dict, mir_mat = mir_sub, gene_mat = gene_sub)
        results_mat <- data.table::as.data.table(results_mat)

        dict <- merge(dict, results_mat[, .(miRNA, gene, cor, corPadj)], by = c("miRNA", "gene"), all.x = TRUE)
        dict <- dict[complete.cases(dict)]

        # dict[, corPadj := p.adjust(pvalue, method = "BH")]

        dict[, ratio := (geneLogFC / mirLogFC) * ((1 - genePadj) / (1 + mirPadj))]
        dict[, normratio := apply_scoreFun(ratio, scoreFun)]

        dict[, strength := cor * normratio]
        dict[, MGI := paste0(miRNA, ":", gene)]
        data.table::setorder(dict, -strength)
        dict <- unique(dict)
        dict <- dict[, .(
            MGI, miRNA, gene, cor, corPadj, mirLogFC, mirPadj, geneLogFC, genePadj,
            ratio, normratio, strength
        )]
        dict <- as.data.frame(dict)
        return(dict)
    }
}


#' @noRd
align_expr_mat_rownames <- function(mat, id_col_name = NULL, dedup = FALSE) {
    if (!is.numeric(mat[[1]]) && !is.integer(mat[[1]])) {
        # we suppose the first column is the miRNA/gene name in that case.
        if (!is.null(id_col_name)) colnames(mat)[1] <- id_col_name
        if (dedup) mat <- mat %>% distinct(.[1], .keep_all = TRUE) # only keep the first appearance if duplicated values.
        rownames(mat) <- mat[[1]]
        mat <- mat[, -1, drop = FALSE]
    }
    mat
}

#' @noRd
align_expr_mat_samples <- function(mir_mat, gene_mat, rename_id_cols = FALSE, dedup = FALSE, check_samples = TRUE) {
    mir_mat <- align_expr_mat_rownames(mir_mat, id_col_name = if (rename_id_cols) "miRNA" else NULL, dedup = dedup)
    gene_mat <- align_expr_mat_rownames(gene_mat, id_col_name = if (rename_id_cols) "gene" else NULL, dedup = dedup)

    common_samples <- NULL
    if (check_samples) {
        common_samples <- intersect(colnames(mir_mat), colnames(gene_mat))
        if (length(common_samples) == 0) {
            stop("No matched samples between miRNA expression and gene expression profiles,
         please ensure mir_mat and gene_mat have shared colnames as sample names!")
        }
        if (length(common_samples) < min(ncol(mir_mat), ncol(gene_mat)) * 0.5) {
            warning(
                "The number of matched samples between mir_mat and gene_mat is:", length(common_samples),
                ", which is below 50% of the total number of samples.
            We suggest to double check the colnames of miRNA and gene expression profiles have shared colnames as sample names!"
            )
        }
    }

    list(mir_mat = mir_mat, gene_mat = gene_mat, common_samples = common_samples)
}

#' @noRd
apply_scoreFun <- function(ratio, scoreFun) {
    if (scoreFun == "rank") {
        ranknorm(ratio)
    } else if (scoreFun == "sigmoid") {
        sigmoid(ratio)
    } else if (scoreFun == "normalize") {
        max_min(ratio)
    } else {
        NA_real_
    }
}

#' @param dict data.frame/data.table with at least `miRNA` and `gene` columns naming the pairs to correlate.
#' @param mir_mat miRNA-by-sample expression matrix (or data.frame with miRNA names in the first column).
#' @param gene_mat gene-by-sample expression matrix (or data.frame with gene names in the first column).
#' @return A data.frame with `miRNA`, `gene`, `cor` (Spearman rho), and `corPadj` (BH-adjusted p-value).
#' @details Computes Spearman correlation for every (miRNA, gene) pair in `dict` by rank-transforming
#' each expression matrix once, then evaluating all pairs with matched-index matrix algebra
#' (no per-pair loop). Pairs are processed in chunks of 2e5 to bound peak memory on very large
#' `dict` inputs (e.g. whole-pathway background correlation for `Edge_Topology`).
#' @noRd
calc_MGI_corr <- function(dict, mir_mat, gene_mat) {
    aligned <- align_expr_mat_samples(mir_mat, gene_mat)
    mir_mat <- aligned$mir_mat
    gene_mat <- aligned$gene_mat
    common_samples <- aligned$common_samples

    # build mat that have both miRNA and gene with only shared samples.
    mir_sub <- mir_mat[, common_samples, drop = FALSE]
    gene_sub <- gene_mat[, common_samples, drop = FALSE]

    # rank-transform each row once; Spearman correlation is Pearson correlation on ranks.
    rank_mir <- t(apply(as.matrix(mir_sub), 1, rank))
    rank_gene <- t(apply(as.matrix(gene_sub), 1, rank))

    idx_mir_all <- match(dict$miRNA, rownames(mir_sub))
    idx_gene_all <- match(dict$gene, rownames(gene_sub))
    n_samples <- ncol(rank_mir)
    n_pairs <- nrow(dict)

    chunk_size <- 2e5
    rho_all <- numeric(n_pairs)
    pval_all <- numeric(n_pairs)
    for (start in seq(1, n_pairs, by = chunk_size)) {
        idx <- start:min(start + chunk_size - 1, n_pairs)
        Rm <- rank_mir[idx_mir_all[idx], , drop = FALSE]
        Rg <- rank_gene[idx_gene_all[idx], , drop = FALSE]
        mean_m <- rowMeans(Rm)
        mean_g <- rowMeans(Rg)
        Rm_c <- Rm - mean_m
        Rg_c <- Rg - mean_g
        cov_mg <- rowSums(Rm_c * Rg_c)
        sd_m <- sqrt(rowSums(Rm_c^2))
        sd_g <- sqrt(rowSums(Rg_c^2))
        rho <- cov_mg / (sd_m * sd_g)
        tstat <- rho * sqrt((n_samples - 2) / (1 - rho^2))
        rho_all[idx] <- rho
        pval_all[idx] <- 2 * stats::pt(-abs(tstat), df = n_samples - 2)
    }

    data.frame(
        miRNA = dict$miRNA,
        gene = dict$gene,
        cor = rho_all,
        corPadj = stats::p.adjust(pval_all, method = "BH")
    )
}


min_max <- function(x) {
    diff <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (diff != 0) {
        scaled <- (x - min(x, na.rm = TRUE)) / diff
    } else {
        scaled <- 0
    }
}

ranknorm <- function(x) {
    ranks <- rank(x)
    result <- (max(ranks) + 1 - ranks) / length(ranks)
    return(result)
}

sigmoid <- function(x) {
    1 / (1 + exp(x))
}

max_min <- function(x) {
    diff <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (diff != 0) {
        scaled <- 1 - (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        # scaled * 0.98 + 0.01
    } else {
        scaled <- 0
    }
}
