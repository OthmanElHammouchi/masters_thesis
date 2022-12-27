library(data.table)

genConfig <- function(...) {

    args <- list(...)
    nargs <- length(args)

    indices.args <- list(rep(0, nargs))

    for (idx in seq_len(nargs)) {
        if (is.vector(args[[idx]])) {
            indices.args[[idx]] <- seq_len(length(args[[idx]]))
        } else if (is.array(args[[idx]])) {
            indices.args[[idx]] <- seq_len(nrow(args[[idx]]))
        }
    }

    indices <- do.call(expand.grid, indices.args)

    config.args <- list(rep(0, nargs))

    for (idx in seq_len(nargs)) {
        if (is.vector(args[[idx]])) {
            config.args[[idx]] <- args[[idx]][indices[, idx]]
        } else if (is.array(args[[idx]])) {
            config.args[[idx]] <- args[[idx]][indices[, idx], ]
        }
    }

    config <- do.call(cbind.data.frame, config.args)

    return(as.data.table(config))
}