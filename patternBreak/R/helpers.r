#' Generate configuration for Mack bootstrap simulations.
#'
#' FUNCTION_DESCRIPTION
#'
#' @param ndev DESCRIPTION.
#' @param factors DESCRIPTION.
#' @param resids.type DESCRIPTION.
#' @param boot.type DESCRIPTION.
#' @param dist DESCRIPTION.
#' @param type DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
mackConfig <- function(ndev,
    factors,
    resids.type = c("raw", "scaled", "parametric"),
    boot.type = c("conditional", "unconditional"),
    dist = c("normal", "gamma"),
    type = "single") {

    if (type == "single") {

        npts <- (ndev ** 2 - ndev) / 2  # First column can't be outlier.
        outlier.points <- matrix(rep(0, 2 * npts), ncol = 2)
        k <- 1
        for (i in seq_len(ndev - 1)) {
            for (j in seq(2, ndev + 1 - i)) {
                outlier.points[k, ] <- c(i, j)
                k <- k + 1
            }
        }

        excl.points <- outlier.points

        indices <- do.call(expand.grid,
            list(
                seq_len(nrow(outlier.points)),
                seq_along(factors),
                seq_len(nrow(excl.points)),
                seq_along(resids.type),
                seq_along(boot.type),
                seq_along(dist)
            ))


        config <- do.call(cbind.data.frame, 
            list(
                outlier.points[indices[, 1], ],
                factors[indices[, 2]],
                excl.points[indices[, 3], ],
                resids.type[indices[, 4]],
                boot.type[indices[, 5]],
                dist[indices[, 6]]
            ))

        config <- data.table::as.data.table(config)

        names(config) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx", "resids.type", "boot.type", "dist")

        return(config)

    } else if (type == "calendar" || type == "origin") {

        excl.diags <- outlier.diags <- seq_len(ndev - 1)

        indices <- do.call(expand.grid,
            list(
                seq_along(outlier.diags),
                seq_along(excl.diags),
                seq_along(factors),
                seq_along(resids.type),
                seq_along(boot.type),
                seq_along(dist)
            ))


        config <- do.call(cbind.data.frame, 
            list(
                outlier.diags[indices[, 1]],
                factors[indices[, 2]],
                excl.diags[indices[, 3]],
                resids.type[indices[, 4]],
                boot.type[indices[, 5]],
                dist[indices[, 6]]
            ))

        config <- data.table::as.data.table(config)

        names(config) <- c("outlier.diagidx", "factor", "excl.diagidx", "resids.type", "boot.type", "dist")

        return(config)
    }
}

#' Generate configuration for GLM bootstrap simulations.
#'
#' FUNCTION_DESCRIPTION
#'
#' @param ndev DESCRIPTION.
#' @param factors DESCRIPTION.
#' @param resids.type DESCRIPTION.
#' @param boot.type DESCRIPTION.
#' @param dist DESCRIPTION.
#' @param type DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
glmConfig <- function(ndev,
    factors,
    type = "single") {

    if (type == "single") {

        npts <- (ndev ** 2 + ndev) / 2
        outlier.points <- matrix(rep(0, 2 * npts), ncol = 2)
        k <- 1
        for (i in seq_len(ndev)) {
            for (j in seq_len(ndev + 1 - i)) {
                outlier.points[k, ] <- c(i, j)
                k <- k + 1
            }
        }

        excl.points <- outlier.points

        indices <- do.call(expand.grid,
            list(
                seq_len(nrow(outlier.points)),
                seq_along(factors),
                seq_len(nrow(excl.points))
            ))

        config <- do.call(cbind.data.frame, 
            list(
                outlier.points[indices[, 1], ],
                factors[indices[, 2]],
                excl.points[indices[, 3], ]
            ))

        config <- data.table::as.data.table(config)

        names(config) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx")

        return(config)

    } else if (type == "calendar" || type == "origin") {

        excl.diags <- outlier.diags <- seq_len(ndev - 1)

        indices <- do.call(expand.grid,
            list(
                seq_along(outlier.diags),
                seq_along(excl.diags),
                seq_along(factors)
            ))

        config <- do.call(cbind.data.frame, 
            list(
                outlier.diags[indices[, 1]],
                factors[indices[, 2]],
                excl.diags[indices[, 3]]
            ))

        config <- data.table::as.data.table(config)

        names(config) <- c("outlier.diagidx", "factor", "excl.diagidx")

        return(config)
    }
}

#' Convert cumulative triangle to incremental one.
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
cum2incr <- function(triangle) {

    ndev <- ncol(triangle)

    for (j in seq(2, ndev)) {
        for (i in seq(1, ndev + 1 - j)) {
            triangle[i, j] <- triangle[i, j] - sum(triangle[i, 1:(j - 1)])
        }
    }

    return(triangle)
}


#' Convert incremental triangle to cumulative one.
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
incr2cum <- function(triangle) {

    ndev <- ncol(triangle)

    for (j in seq(ndev, 2)) {
        for (i in seq(1, ndev + 1 - j)) {
            triangle[i, j] <- triangle[i, j] + sum(triangle[i, 1:(j - 1)])
        }
    }

    return(triangle)
}
