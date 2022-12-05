reserveBoot <- function(triangle, nboot, ...,
    resids.type = "raw",
    bootstrap.type = "conditional",
    distribution = "normal",
    exclude.residuals = NULL) {

    triangle <- unname(triangle)
    ndev <- ncol(triangle)

    devfacs <- c()
    sigmas <- c()
    indiv.devfacs <- list()
    prev.cs <- list()

    purrr::map(1:(ndev - 1), function(colidx) {
        model <- lm(
            y ~ x + 0,
            weights = 1 / triangle[1:(ndev - colidx), colidx],
            data = data.frame(
                x = triangle[1:(ndev - colidx), colidx],
                y = triangle[1:(ndev - colidx), colidx + 1]
            )
        )
        devfacs[colidx] <<- unname(model$coefficients)

        if (colidx != (ndev - 1)) {
            sigmas[colidx] <<- summary(model)$sigma
        } else {
            sigmas[colidx] <<-
            sqrt(min(
                sigmas[(length(sigmas) - 1)]**2,
                sigmas[length(sigmas)]**2, 
                sigmas[length(sigmas)]**4 / sigmas[(length(sigmas) - 1)]**2
            ))
        }

        indiv.devfacs[[colidx]] <<- triangle[1:(ndev - colidx), colidx + 1] / triangle[1:(ndev - colidx), colidx]
        prev.cs[[colidx]] <<- triangle[1:(ndev - colidx), colidx]
    })

    if (resids.type == "raw") {
        resids <- list()

        purrr::map(1:(ndev - 1), function(colidx) {
            resids[[colidx]] <<- (indiv.devfacs[[colidx]] - devfacs[colidx]) * sqrt(prev.cs[[colidx]]) / sigmas[colidx]
        })
    } else if (resids.type == "scaled") {
        scaleFactors <- list()
        resids <- list()

        purrr::map(1:(ndev - 1), function(colIdx) {
            factors <- sqrt(1 - prev.cs[[colIdx]] / sum(prev.cs[[colIdx]]))

            scaleFactors[[colIdx]] <<- if (length(factors) > 1) factors else 1

            resids[[colIdx]] <<-
             (indiv.devfacs[[colIdx]] - devfacs[colIdx]) * sqrt(prev.cs[[colIdx]]) / (sigmas[colIdx] * scaleFactors[[colIdx]])
        })
    }


    # remove residuals for the excluded points if there are any
    if (!is.null(exclude.residuals)) {
        for (point in exclude.residuals) {
            resids[[point[2] - 1]] <- resids[[point[2] - 1]][-point[1]]
        }
    }

    reserve <- c()
    nNaN <- 0
    for (iBoot in 1:nboot) {
        if (!(resids.type == "parametric")) {
            residBoot <- list()
            resids <- unlist(resids)
            residSample <- sample(resids, replace = TRUE, size = length(unlist(prev.cs)))
            for (j in 1:(ndev - 1)) {
                residBoot[[j]] <- residSample[1:(ndev - j)]
                residSample <- residSample[-1:-(ndev - j)]
            }
        } else {
            residBoot <- list()
            residSample <- rnorm(length(unlist(prev.cs)))
            for (j in 1:(ndev - 1)) {
                residBoot[[j]] <- residSample[1:(ndev - j)]
                residSample <- residSample[-1:-(ndev - j)]
            }
        }

        # compute bootstrapped quantities
        indivDevFacBoot <- list()
        devFacBoot <- c()
        sigmaBoot <- c()

        if (bootstrap.type == "conditional") {
            purrr::map(1:(ndev - 1), function(colIdx) {
                if (resids.type == "scaled") {
                    indivDevFacBoot[[colIdx]] <<- devfacs[colIdx] + (scaleFactors[[colIdx]] * residBoot[[colIdx]] * sigmas[colIdx]) / sqrt(prev.cs[[colIdx]])
                } else {
                    indivDevFacBoot[[colIdx]] <<- devfacs[colIdx] + (residBoot[[colIdx]] * sigmas[colIdx]) / sqrt(prev.cs[[colIdx]])
                }

                devFacBoot[colIdx] <<- sum(indivDevFacBoot[[colIdx]] * prev.cs[[colIdx]]) / sum(prev.cs[[colIdx]])

                if (colIdx != ndev - 1) {
                    df <- ndev - colIdx - 1
                    sigmaBoot[colIdx] <<- sqrt(sum((prev.cs[[colIdx]] * (indivDevFacBoot[[colIdx]] - devfacs[colIdx])**2) / df))
                } else {
                    sigmaBoot[colIdx] <<- sqrt(min(
                        sigmaBoot[(length(sigmaBoot) - 1)]**2,
                        sigmaBoot[length(sigmaBoot)]**2,
                        sigmaBoot[length(sigmaBoot)]**4 / sigmaBoot[(length(sigmaBoot) - 1)]**2
                    ))
                }
            })
        } else if (bootstrap.type == "unconditional") {
            if (resids.type == "scaled") {
                indivDevFacBoot[[1]] <- devfacs[1] + (scaleFactors[[1]] * residBoot[[1]] * sigmas[1]) / sqrt(prev.cs[[1]])
            } else {
                indivDevFacBoot[[1]] <- devfacs[1] + (residBoot[[1]] * sigmas[1]) / sqrt(prev.cs[[1]])
            }

            devFacBoot[1] <- sum(indivDevFacBoot[[1]] * prev.cs[[1]]) / sum(prev.cs[[1]])

            df <- ndev - 2
            sigmaBoot[1] <- sqrt(sum((prev.cs[[1]] * (indivDevFacBoot[[1]] - devfacs[1])**2) / df))

            new.cs <- list(triangle[, 1])

            for (colidx in 2:(ndev - 1)) {
                new.cs[[colidx]] <- new.cs[[colidx - 1]][-length(new.cs[[colidx - 1]])] * devfacs[colidx - 1] +
                    sigmas[colidx - 1] * sqrt(new.cs[[colidx - 1]][-length(new.cs[[colidx - 1]])]) * residBoot[[colidx - 1]]

                if (resids.type == "scaled") {
                    indivDevFacBoot[[colidx]] <- devfacs[colidx] + (scaleFactors[[colidx]] * residBoot[[colidx]] * sigmas[colidx]) / sqrt(prev.cs[[colidx]])
                } else {
                    indivDevFacBoot[[colidx]] <- devfacs[colidx] + (residBoot[[colidx]] * sigmas[colidx]) / sqrt(prev.cs[[colidx]])
                }

                devFacBoot[colidx] <- sum(indivDevFacBoot[[colidx]] * prev.cs[[colidx]]) / sum(prev.cs[[colidx]])

                if (colidx != ndev - 1) {
                    df <- ndev - colidx - 1
                    sigmaBoot[colidx] <- sqrt(sum((prev.cs[[colidx]] * (indivDevFacBoot[[colidx]] - devfacs[colidx])**2) / df))
                } else {
                    sigmaBoot[colidx] <- sqrt(min(
                        sigmaBoot[(length(sigmaBoot) - 1)]**2,
                        sigmaBoot[length(sigmaBoot)]**2,
                        sigmaBoot[length(sigmaBoot)]**4 / sigmaBoot[(length(sigmaBoot) - 1)]**2
                    ))
                }
            }
        }


        if (distribution == "normal") {
            # complete the triangle
            triangleBoot <- triangle
            for (diagIdx in 1:(ndev - 1)) {
                for (rowIdx in (diagIdx + 1):ndev) {
                    colIdx <- ndev + diagIdx + 1 - rowIdx
                    # off-diagonal elements satisfy i + j = I + 1 + (diagonal number)
                    triangleBoot[rowIdx, colIdx] <-
                        rnorm(
                            1,
                            triangleBoot[rowIdx, colIdx - 1] * devFacBoot[colIdx - 1],
                            sqrt(triangleBoot[rowIdx, colIdx - 1]) * sigmaBoot[colIdx - 1]
                        )
                }
            }
        } else if (distribution == "gamma") {
            triangleBoot <- triangle
            for (diagIdx in 1:(ndev - 1)) {
                for (rowIdx in (diagIdx + 1):ndev) {
                    colIdx <- ndev + diagIdx + 1 - rowIdx
                    alpha <- devfacs[colIdx - 1]**2 * triangleBoot[rowIdx, colIdx - 1] / sigmas[colIdx - 1]**2
                    beta <- devfacs[colIdx - 1] / (sigmas[colIdx - 1]**2)
                    triangleBoot[rowIdx, colIdx] <-
                        rgamma(1, shape = alpha, rate = beta)
                }
            }
        }

        # if the triangle contains undefined values, the simulated reserve became negative at some point, and we throw away the triangle
        if (!(NaN %in% triangleBoot)) {
            latest <- triangleBoot[col(triangleBoot) + row(triangleBoot) == ndev + 1]
            reserve <- c(reserve, sum(triangleBoot[, ndev] - latest))
        } else {
            nNaN <- nNaN + 1
        }
    }
    attr(reserve, "nNaN") <- nNaN
    return(reserve)
}

reserveBootFortran <- function(triangle, nboot, distribution = "normal", resids.type = "raw", bootstrap.type = "conditional", excl.resids = NULL) {
    dyn.load("fortran/build/reserve_boot.so")

    triangle[is.na(triangle)] <- 0
    ndev <- ncol(triangle)
    reserve <- rep(0, nboot)

    distribution.options <- list(normal = 1L, gamma = 2L)
    resids_type.options <- list(raw = 1L, scaled = 2L, parametric = 3L)
    bootstrap.type.options <- list(conditional = 1L, unconditional = 2L)

    excl.resids.array <- matrix(c(0L, 0L), nrow = 2)

    if (!is.null(excl.resids)) {
        excl.resids.array <- matrix(nrow = 2, ncol = length(excl.resids))
        for (listidx in seq_along(excl.resids)) {
            excl.resids.array[, listidx] <- excl.resids[[listidx]]
        }
    }

    excl.resids <- excl.resids.array

    browser()

    res <- .Fortran(
        "reserve_boot",
        nboot = as.integer(nboot),
        ndev = ndev,
        triangle = triangle,
        reserve = reserve,
        distribution = distribution.options[[distribution]],
        resids.type = resids_type.options[[resids.type]],
        bootstrap.type = bootstrap.type.options[[bootstrap.type]],
        excl.resids = excl.resids,
        excl.resids.ncols = ncol(excl.resids)
    )

    return(res$reserve)
}

singleOutlier <- function(outlier.rowidx, outlier.colidx, factor, ..., initcol, devfacs, sigmas, distribution = "normal") {
    ndev <- length(initcol)

    if (outlier.colidx == 1) {
        stop("Outlier column index must be greater than 1.")
    }

    if (distribution == "normal") {
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
