reserveBoot <- function(triangle, nboot, ...,
    resids.type = "raw",
    boot.type = "conditional",
    dist = "normal",
    excl.resids = NULL) {

    if (dist == "gamma" && resids.type == "parametric") {
        stop("Parametric residuals can only be used with normal distribution")
    }

    triangle <- unname(triangle)
    ndev <- ncol(triangle)

    devfacs <- rep(0, ndev - 1)
    sigmas <- rep(0, ndev - 1)
    indiv.devfacs <- matrix(nrow = ndev - 1, ncol = ndev - 1)
    resids <- matrix(nrow = ndev - 1, ncol = ndev - 1)

    for (colidx in 2:ndev) {

        nrows <- ndev + 1 - colidx

        model <- lm(
            y ~ x + 0,
            weights = 1 / triangle[1:nrows, colidx - 1],
            data = data.frame(
                x = triangle[1:nrows, colidx - 1],
                y = triangle[1:nrows, colidx]
            )
        )

        devfacs[colidx - 1] <- unname(model$coefficients)

        if (colidx < ndev) {
            sigmas[colidx - 1] <- summary(model)$sigma
        } else {
            sigmas[colidx - 1] <- sqrt(min(
                sigmas[(colidx - 2)]**2,
                sigmas[colidx - 3]**2,
                sigmas[colidx - 2]**4 / sigmas[colidx - 3]**2
            ))
        }

        indiv.devfacs[1:nrows, colidx - 1] <- triangle[1:nrows, colidx] / triangle[1:nrows, colidx - 1]

        if (resids.type == "raw") {

            resids[1:nrows, colidx - 1] <-
                (indiv.devfacs[1:nrows, colidx - 1] - devfacs[colidx - 1]) * sqrt(triangle[1:nrows, colidx - 1]) /
                    sigmas[colidx - 1]

        } else if (resids.type == "scaled") {

            if (colidx < ndev) {
                scale.factors <- sqrt(1 - triangle[1:nrows, colidx - 1] / sum(triangle[, colidx - 1]))

                resids[1:nrows, colidx - 1] <- (indiv.devfacs[1:nrows, colidx - 1] - devfacs[colidx - 1]) *
                    sqrt(triangle[1:nrows, colidx - 1]) /
                        (sigmas[colidx - 1] * scale.factors)

            } else {
                resids[1:nrows, colidx - 1] <- (indiv.devfacs[1:nrows, colidx - 1] - devfacs[colidx - 1]) *
                        sqrt(triangle[1:nrows, colidx - 1]) / sigmas[colidx - 1]
                }
        }

    }

    # remove residuals for the excluded points if there are any
    if (!is.null(excl.resids)) {
        for (rowidx in seq_len(nrow(excl.resids))) {
            resids[excl.resids[rowidx, 1], excl.resids[rowidx, 2] - 1] <- NA
        }
    }

    # nresids <- (ndev ** 2 + ndev) / 2 - nrow(excl.resids)
    resids[1, ndev - 1] <- NA
    flat.resids <- as.vector(resids)
    flat.resids <- flat.resids[!is.na(flat.resids)]


    reserve <- rep(0, nboot)
    nNaN <- 0

    for (iboot in 1:nboot) {

        if (!(resids.type == "parametric")) {

            resids.boot <- matrix(nrow = ndev - 1, ncol = ndev - 1)
            for (colidx in 2:ndev) {
                nrows <- ndev + 1 - colidx
                resids.boot[1:nrows, colidx - 1] <- sample(flat.resids, replace = TRUE, size = nrows)
            }

        } else {

            resids.boot <- matrix(nrow = ndev - 1, ncol = ndev - 1)
            for (colidx in 2:ndev) {
                nrows <- ndev + 1 - colidx
                resids.boot[1:nrows, colidx - 1] <- rnorm(nrows)
            }

        }

        # compute bootstrapped quantities
        indiv.devfacs.boot <- matrix(nrow = ndev - 1, ncol = ndev  - 1)
        devfacs.boot <- rep(0, ndev - 1)
        sigmas.boot <- rep(0, ndev - 1)

        if (boot.type == "conditional") {

            for (colidx in 2:ndev) {

                nrows <- ndev + 1 - colidx

                indiv.devfacs.boot[1:nrows, colidx - 1] <- devfacs[colidx - 1] +
                    (resids.boot[1:nrows, colidx - 1] * sigmas[colidx - 1]) / sqrt(triangle[1:nrows, colidx  - 1])

                devfacs.boot[colidx - 1] <-
                sum(indiv.devfacs.boot[1:nrows, colidx - 1] * triangle[1:nrows, colidx - 1]) /
                    sum(triangle[1:nrows, colidx - 1])

                if (colidx < ndev) {

                    sigmas.boot[colidx - 1] <- sqrt(sum(triangle[1:nrows, colidx - 1] *
                    (indiv.devfacs.boot[1:nrows, colidx - 1] - devfacs.boot[colidx - 1])**2) / (nrows - 1))

                } else {

                    sigmas.boot[colidx - 1] <- sqrt(min(
                        sigmas.boot[colidx - 2]**2,
                        sigmas.boot[colidx - 3]**2,
                        sigmas.boot[colidx - 2]**4 / sigmas.boot[colidx - 3]**2
                    ))
                }
            }

        if (any(sigmas.boot == 0)) browser()

        } else if (boot.type == "unconditional") {

            resampled.triangle <- matrix(nrow = ndev, ncol = ndev)
            resampled.triangle[, 1] <- triangle[, 1]

            for (colidx in 2:ndev) {

                nrows <- ndev + 1 - colidx

                resampled.triangle[1:nrows, colidx] <- devfacs[colidx - 1] * resampled.triangle[1:nrows, colidx - 1] +
                 sigmas[colidx - 1] * sqrt(resampled.triangle[1:nrows, colidx - 1]) * resids.boot[1:nrows, colidx - 1]


                indiv.devfacs.boot[1:nrows, colidx - 1] <-
                    resampled.triangle[1:nrows, colidx] / resampled.triangle[1:nrows, colidx - 1]

                devfacs.boot[colidx - 1] <-
                    sum(resampled.triangle[1:nrows, colidx]) / sum(resampled.triangle[1:nrows, colidx - 1])

                if (colidx < ndev) {
                    sigmas.boot[colidx - 1] <- sqrt(sum(resampled.triangle[1:nrows, colidx - 1] *
                    (indiv.devfacs.boot[1:nrows, colidx - 1] - devfacs.boot[colidx - 1])**2) / (nrows - 1))
                } else {
                    sigmas.boot[colidx - 1] <- sqrt(min(
                        sigmas.boot[colidx - 2]**2,
                        sigmas.boot[colidx - 3]**2,
                        sigmas.boot[colidx - 2]**4 / sigmas.boot[colidx - 3]**2
                    ))
                }
            }
        }

        if (dist == "normal") {

            # complete the triangle
            triangle.boot <- triangle

            for (diagidx in 1:(ndev - 1)) {
                for (rowidx in (diagidx + 1):ndev) {

                    colidx <- ndev + diagidx + 1 - rowidx
                    # off-diagonal elements satisfy i + j = I + 1 + (diagonal number)
                    triangle.boot[rowidx, colidx] <-
                        rnorm(
                            1,
                            triangle.boot[rowidx, colidx - 1] * devfacs.boot[colidx - 1],
                            sqrt(triangle.boot[rowidx, colidx - 1]) * sigmas.boot[colidx - 1]
                        )

                }
            }
        } else if (dist == "gamma") {

            triangle.boot <- triangle

            for (diagidx in 1:(ndev - 1)) {
                for (rowidx in (diagidx + 1):ndev) {
                    colidx <- ndev + diagidx + 1 - rowidx
                    alpha <- devfacs.boot[colidx - 1]**2 * triangle.boot[rowidx, colidx - 1] /
                    sigmas.boot[colidx - 1]**2
                    beta <- devfacs.boot[colidx - 1] / (sigmas.boot[colidx - 1]**2)
                    triangle.boot[rowidx, colidx] <-
                        rgamma(1, shape = alpha, rate = beta)
                    if (triangle.boot[rowidx, colidx] <= 0) browser()
                }
            }
        }

        # if the triangle contains undefined values, the simulated reserve became negative at some point,
        # and we throw away the triangle
        if (!(NaN %in% triangle.boot)) {
            latest <- triangle.boot[col(triangle.boot) + row(triangle.boot) == ndev + 1]
            reserve[iboot] <- sum(triangle.boot[, ndev] - latest)
        } else {
            nNaN <- nNaN + 1
        }
    }
    attr(reserve, "nNaN") <- nNaN
    return(reserve)
}

reserveBootFortran <- function(triangle, nboot,
    resids.type = "raw",
    boot.type = "conditional",
    dist = "normal",
    excl.resids = NULL) {

    dyn.load("fortran/build/reserve_boot.so")

    triangle[is.na(triangle)] <- 0
    ndev <- ncol(triangle)
    reserve <- rep(0, nboot)

    dist.options <- list(normal = 1L, gamma = 2L)
    resids.type.options <- list(raw = 1L, scaled = 2L, parametric = 3L)
    boot.type.options <- list(conditional = 1L, unconditional = 2L)

    if (is.null(excl.resids)) {
        excl.resids <- matrix(c(0L, 0L), nrow = 2)
    }

    res <- .Fortran(
        "reserve_boot",
        nboot = as.integer(nboot),
        ndev = ndev,
        triangle = triangle,
        reserve = reserve,
        dist = dist.options[[dist]],
        resids.type = resids.type.options[[resids.type]],
        boot.type = boot.type.options[[boot.type]],
        excl.resids = excl.resids,
        excl.resids.ncols = ncol(excl.resids)
    )

    return(res$reserve)
}

singleOutlier <- function(outlier.rowidx, outlier.colidx, factor, ...,
    initcol,
    devfacs,
    sigmas,
    dist = "normal") {

    ndev <- length(initcol)

    if (outlier.colidx == 1) {
        stop("Outlier column index must be greater than 1.")
    }

    if (dist == "normal") {
        triangle <- matrix(ncol = ndev, nrow = ndev)
        triangle[, 1] <- initcol

        for (colidx in 2:ndev) {
            for (rowidx in setdiff(1:(ndev + 1 - colidx), outlier.rowidx)) {
                prevc <- triangle[rowidx, colidx - 1]
                triangle[rowidx, colidx] <-
                    rnorm(1, devfacs[colidx - 1] * prevc, sigmas[colidx - 1] * sqrt(prevc))
            }
        }

        if (outlier.colidx > 2) {
            for (colidx in 2:(outlier.colidx - 1)) {
                prevc <- triangle[outlier.rowidx, colidx - 1]
                triangle[outlier.rowidx, colidx] <-
                    rnorm(1, devfacs[colidx - 1] * prevc, sigmas[colidx - 1] * sqrt(prevc))
            }
        }

        prevc <- triangle[outlier.rowidx, outlier.colidx - 1]
        triangle[outlier.rowidx, outlier.colidx] <-
            rnorm(
                1,
                factor * devfacs[outlier.colidx - 1] * prevc,
                sigmas[outlier.colidx - 1] * sqrt(prevc)
            )

        if (outlier.colidx < ndev) {
            for (colidx in (outlier.colidx + 1):(ndev + 1 - outlier.rowidx)) {
                prevc <- triangle[outlier.rowidx, colidx - 1]
                triangle[outlier.rowidx, colidx] <-
                    rnorm(1, devfacs[colidx - 1] * prevc, sigmas[colidx - 1] * sqrt(prevc))
            }
        }

        return(triangle)
    } else {
        triangle <- matrix(ncol = ndev, nrow = ndev)
        triangle[, 1] <- initcol

        for (colidx in 2:ndev) {
            for (rowidx in setdiff(1:(ndev + 1 - colidx), outlier.rowidx)) {
                prevc <- triangle[rowidx, colidx - 1]
                alpha <- devfacs[colidx - 1]**2 * prevc / sigmas[colidx - 1]**2
                beta <- devfacs[colidx - 1] / (sigmas[colidx - 1]**2)

                triangle[rowidx, colidx] <-
                    rgamma(1, shape = alpha, rate = beta)
            }
        }

        if (outlier.colidx > 2) {
            for (colidx in 2:(outlier.colidx - 1)) {
                prevc <- triangle[outlier.rowidx, colidx - 1]
                alpha <- devfacs[colidx - 1]**2 * prevc / sigmas[colidx - 1]**2
                beta <- devfacs[colidx - 1] / (sigmas[colidx - 1]**2)

                triangle[outlier.rowidx, colidx] <-
                    rgamma(1, shape = alpha, rate = beta)
            }
        }

        prevc <- triangle[outlier.rowidx, outlier.colidx - 1]
        alpha <- factor * devfacs[outlier.colidx - 1]**2 * prevc / sigmas[outlier.colidx - 1]**2
        beta <- devfacs[outlier.colidx - 1] / (sigmas[outlier.colidx - 1]**2)

        triangle[outlier.rowidx, outlier.colidx] <-
            rgamma(1, shape = alpha, rate = beta)

        if (outlier.colidx < ndev) {
            for (colidx in (outlier.colidx + 1):(ndev + 1 - outlier.rowidx)) {
                prevc <- triangle[outlier.rowidx, colidx - 1]
                alpha <- devfacs[colidx - 1]**2 * prevc / sigmas[colidx - 1]**2
                beta <- devfacs[colidx - 1] / (sigmas[colidx - 1]**2)

                triangle[outlier.rowidx, colidx] <-
                    rgamma(1, shape = alpha, rate = beta)
            }
        }

        return(triangle)
    }
}

calendarOutlier <- function(outlier.diagidx, factor, ..., initcol, devfacs, sigmas) {
    ndev <- length(initcol)

    triangle <- matrix(ncol = ndev, nrow = ndev)
    triangle[, 1] <- initcol

    for (colidx in 2:ndev) {
        # points belonging to the diagidx'th antidiagonal satisfy nDev + 1 - diagidx == colidx + rowidx
        diag.rowidx <- ndev + 1 - outlier.diagidx - colidx
        for (rowIdx in setdiff(1:(ndev + 1 - colidx), diag.rowidx)) {
            prevc <- triangle[rowIdx, colidx - 1]
            triangle[rowIdx, colidx] <-
                rnorm(1, devfacs[colidx - 1] * prevc, sigmas[colidx - 1] * sqrt(prevc))
        }
        if (diag.rowidx > 0) {
            prevc <- triangle[diag.rowidx, colidx - 1]
            triangle[diag.rowidx, colidx] <-
                rnorm(1, factor * devfacs[colidx - 1] * prevc, sigmas[colidx - 1] * sqrt(prevc))
        }
    }
    return(triangle)
}

accidentOutlier <- function(outlier.rowidx, factor, ..., initcol, devfacs, sigmas, distribution = "normal") {
    nDev <- length(initcol)

    claimsTriangle <- matrix(ncol = nDev, nrow = nDev)
    claimsTriangle[, 1] <- initcol

    for (colidx in (2:nDev)) {
        for (rowidx in setdiff(1:(nDev + 1 - colidx), outlier.rowidx)) {
            prevc <- claimsTriangle[rowidx, colidx - 1]

            if (distribution == "gamma") {
                alpha <- devfacs[colidx - 1]**2 * prevc / sigmas[colidx - 1]**2
                beta <- devfacs[colidx - 1] / (sigmas[colidx - 1]**2)

                claimsTriangle[rowidx, colidx] <-
                    rgamma(1, shape = alpha, rate = beta)
            } else {
                claimsTriangle[rowidx, colidx] <-
                    rnorm(1, devfacs[colidx - 1] * prevc, sigmas[colidx - 1] * sqrt(prevc))
            }
        }

        prevc <- claimsTriangle[outlier.rowidx, colidx - 1]
        claimsTriangle[outlier.rowidx, colidx] <-
            rnorm(1, factor * devfacs[colidx - 1] * prevc, sigmas[colidx - 1] * sqrt(prevc))
    }
    return(claimsTriangle)
}

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

    return(config)
}