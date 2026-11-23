# Memoizes calls that are repeated verbatim across test files, keyed on the
# literal argument expressions. test-getAllInputData.R and
# test-input-immutability.R test getAllInputData() itself and must keep
# calling it directly, not through this cache.
make_cached_call <- function(fn) {
    force(fn)
    cache <- new.env(parent = emptyenv())
    function(...) {
        key <- paste(deparse(match.call()), collapse = " ")
        if (!exists(key, envir = cache, inherits = FALSE)) {
            assign(key, fn(...), envir = cache)
        }
        get(key, envir = cache, inherits = FALSE)
    }
}

cached_getAllInputData <- make_cached_call(getAllInputData)
cached_miREA <- make_cached_call(miREA)
