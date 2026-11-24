#' S4 class for miREA() results
#'
#' A lightweight container bundling the per-method enrichment result tables
#' together with the timing/validity/parameter metadata returned by
#' [miREA()]. Supports `$` and `[[` access for backward compatibility with
#' the function's previous plain-list return value.
#'
#' @slot result Named list, one enrichment result data.frame per requested method.
#' @slot time A data.frame recording the elapsed time (in seconds) for each method.
#' @slot valid_methods Character vector of methods that completed without error.
#' @slot failed_methods Named list mapping each method that errored to its error message.
#' @slot parameter A list of the parameters [miREA()] was called with.
#'
#' @examples
#' res <- miREAResult(
#'   result = list(TG_ORA = data.frame(pathway = "PW1", padj = 0.01)),
#'   time = data.frame(method = "TG_ORA", time = 0.1),
#'   valid_methods = "TG_ORA",
#'   parameter = list(n_pathway = 1)
#' )
#' res
#' res$valid_methods
#' res[["time"]]
#' as.data.frame(res)
#'
#' @name miREAResult-class
#' @rdname miREAResult-class
#' @exportClass miREAResult
methods::setClass("miREAResult", slots = c(
  result = "list",
  time = "data.frame",
  valid_methods = "character",
  failed_methods = "list",
  parameter = "list"
), prototype = methods::prototype(failed_methods = list()))

methods::setValidity("miREAResult", function(object) {
  errors <- character()
  if (!all(object@valid_methods %in% names(object@result))) {
    errors <- c(errors, "valid_methods must be a subset of names(result)")
  }
  if (!all(object@time$method %in% object@valid_methods)) {
    errors <- c(errors, "time$method must be a subset of valid_methods")
  }
  if (!is.list(object@failed_methods) ||
    (length(object@failed_methods) > 0 && is.null(names(object@failed_methods)))) {
    errors <- c(errors, "failed_methods must be a named list")
  }
  if (length(errors) == 0) TRUE else errors
})

#' @param result Named list, one enrichment result data.frame per requested method.
#' @param time A data.frame recording the elapsed time (in seconds) for each method.
#' @param valid_methods Character vector of methods that completed without error.
#' @param failed_methods Named list mapping each method that errored to its error message.
#'        Defaults to an empty list.
#' @param parameter A list of the parameters [miREA()] was called with.
#' @return A [miREAResult-class] object.
#' @rdname miREAResult-class
#' @export
miREAResult <- function(result, time, valid_methods, failed_methods = list(), parameter) {
  methods::new("miREAResult",
    result = result, time = time,
    valid_methods = valid_methods, failed_methods = failed_methods,
    parameter = parameter
  )
}

#' @param x A [miREAResult-class] object.
#' @param name The slot to access, one of `"result"`, `"time"`,
#'        `"valid_methods"`, `"failed_methods"`, `"parameter"`.
#' @rdname miREAResult-class
#' @export
methods::setMethod("$", "miREAResult", function(x, name) {
  if (!name %in% methods::slotNames(x)) {
    stop(
      "'", name, "' is not a valid miREAResult slot. Available: ",
      paste(methods::slotNames(x), collapse = ", ")
    )
  }
  methods::slot(x, name)
})

#' @param i The slot to access, one of `"result"`, `"time"`,
#'        `"valid_methods"`, `"failed_methods"`, `"parameter"`.
#' @param j Not used.
#' @param ... Not used.
#' @rdname miREAResult-class
#' @export
methods::setMethod("[[", "miREAResult", function(x, i, j, ...) {
  if (!i %in% methods::slotNames(x)) {
    stop(
      "'", i, "' is not a valid miREAResult slot. Available: ",
      paste(methods::slotNames(x), collapse = ", ")
    )
  }
  methods::slot(x, i)
})

#' @param object A [miREAResult-class] object.
#' @rdname miREAResult-class
#' @importFrom methods show
#' @export
methods::setMethod("show", "miREAResult", function(object) {
  cat("miREAResult object\n")
  cat("  Methods completed:", paste(object@valid_methods, collapse = ", "), "\n")
  cat("  Pathways tested:", object@parameter$n_pathway, "\n")
  if (length(object@failed_methods) > 0) {
    cat("  Failed methods:", paste(names(object@failed_methods), collapse = ", "), "\n")
  }
  invisible(object)
})

#' @param row.names Not used.
#' @param optional Not used.
#' @rdname miREAResult-class
#' @importFrom BiocGenerics as.data.frame
#' @export
methods::setMethod("as.data.frame", "miREAResult", function(x, row.names = NULL, optional = FALSE, ...) {
  dplyr::bind_rows(x@result, .id = "method")
})
