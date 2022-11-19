reserveBoot <- function(inputTriangle, nBoot, ..., resids_type = "raw", bootstrap_type = "conditional", distribution = "normal", exclude_residuals = NULL) {

    inputTriangle <- unname(inputTriangle)
    nDev <- ncol(inputTriangle)

    devFacs <- c()
    sigma <- c()
    indivDevFacs <- list()
    prevCs <- list()

    purrr::map(1:(nDev - 1), function(colIdx) {

        model <- lm(
            y ~ x + 0,
            weights = 1 / inputTriangle[1:(nDev - colIdx), colIdx],
            data = data.frame(x = inputTriangle[1:(nDev - colIdx), colIdx],
                y = inputTriangle[1:(nDev - colIdx), colIdx + 1])
        )
        devFacs[colIdx] <<- unname(model$coefficients)

        if (colIdx != (nDev - 1)) sigma[colIdx] <<- summary(model)$sigma
        else sigma[colIdx] <<- min(sigma[(length(sigma) - 1)]**2, sigma[length(sigma)]**2, sigma[length(sigma)]**4 / sigma[(length(sigma) - 1)]**2)

        indivDevFacs[[colIdx]] <<- inputTriangle[1:(nDev - colIdx), colIdx + 1] / inputTriangle[1:(nDev - colIdx), colIdx]
        prevCs[[colIdx]] <<- inputTriangle[1:(nDev - colIdx), colIdx]
    })


    if (resids_type == "raw") {

        resids <- list()

        purrr::map(1:(nDev - 1), function(colIdx) {

            resids[[colIdx]] <<- (indivDevFacs[[colIdx]] - devFacs[colIdx]) * sqrt(prevCs[[colIdx]]) / sigma[colIdx]
        })

    } else if (resids_type == "scaled") {

        scaleFactors <- list()
        resids <- list()

        purrr::map(1:(nDev - 1), function(colIdx) {

            factors <- sqrt(1 - prevCs[[colIdx]] / sum(prevCs[[colIdx]]))

            scaleFactors[[colIdx]] <<- if (length(factors) > 1) factors else 1

            resids[[colIdx]] <<- (indivDevFacs[[colIdx]] - devFacs[colIdx]) * sqrt(prevCs[[colIdx]]) / (sigma[colIdx] * scaleFactors[[colIdx]])
        })

    }


    # remove residuals for the excluded points if there are any
    if (!is.null(exclude_residuals)) {
        for (point in exclude_residuals) {
            resids[[point[2]]] <- resids[[point[2]]][-point[1]]
        }
    }

    reserve <- c()
    nNaN <- 0
    for (iBoot in 1:nBoot) {

        if (!(resids_type == "parametric")) {
            residBoot <- list()
            resids <- unlist(resids)
            residSample <- sample(resids, replace = TRUE, size = length(unlist(prevCs)))
            for (j in 1:(nDev - 1)) {
                residBoot[[j]] <- residSample[1:(nDev - j)]
                residSample <- residSample[-1:-(nDev - j)]
            }
        } else {
            residBoot <- list()
            residSample <- rnorm(length(unlist(prevCs)))
            for (j in 1:(nDev - 1)) {
                residBoot[[j]] <- residSample[1:(nDev - j)]
                residSample <- residSample[-1:-(nDev - j)]
            }
        }

        # compute bootstrapped quantities
        indivDevFacBoot <- list()
        devFacBoot <- c()
        sigmaBoot <- c()

        if (bootstrap_type == "conditional") {
            purrr::map(1:(nDev - 1), function(colIdx) {
                if (resids_type == "scaled") {
                    indivDevFacBoot[[colIdx]] <<- devFacs[colIdx] + (scaleFactors[[colIdx]] * residBoot[[colIdx]] * sigma[colIdx]) / sqrt(prevCs[[colIdx]])
                } else {
                    indivDevFacBoot[[colIdx]] <<- devFacs[colIdx] + (residBoot[[colIdx]] * sigma[colIdx]) / sqrt(prevCs[[colIdx]])
                }

                devFacBoot[colIdx] <<- sum(indivDevFacBoot[[colIdx]] * prevCs[[colIdx]]) / sum(prevCs[[colIdx]])

                if (colIdx != nDev - 1) {
                    df <- nDev - colIdx - 1
                    sigmaBoot[colIdx] <<- sqrt(sum((prevCs[[colIdx]] * (indivDevFacBoot[[colIdx]] - devFacs[colIdx])**2) / df))
                } else {
                    sigmaBoot[colIdx] <<- min(
                        sigmaBoot[(length(sigmaBoot) - 1)]**2,
                        sigmaBoot[length(sigmaBoot)]**2,
                        sigmaBoot[length(sigmaBoot)]**4 / sigmaBoot[(length(sigmaBoot) - 1)]**2
                    )
                }
            })
        } else if (bootstrap_type == "unconditional") {

            if (resids_type == "scaled") {
                indivDevFacBoot[[1]] <- devFacs[1] + (scaleFactors[[1]] * residBoot[[1]] * sigma[1]) / sqrt(prevCs[[1]])
            } else {
            indivDevFacBoot[[1]] <- devFacs[1] + (residBoot[[1]] * sigma[1]) / sqrt(prevCs[[1]])
            }

            devFacBoot[1] <- sum(indivDevFacBoot[[1]] * prevCs[[1]]) / sum(prevCs[[1]])

            df <- nDev - 2
            sigmaBoot[1] <- sqrt(sum((prevCs[[1]] * (indivDevFacBoot[[1]] - devFacs[1])**2) / df))

            newCs <- list(inputTriangle[, 1])

            for (colIdx in 2:(nDev - 1)) {

                newCs[[colIdx]] <- newCs[[colIdx - 1]][-length(newCs[[colIdx - 1]])] * devFacs[colIdx] +
                    sigma[colIdx] * sqrt(newCs[[colIdx - 1]][-length(newCs[[colIdx - 1]])]) * residBoot[[colIdx - 1]]
                if (resids_type == "scaled") {
                    indivDevFacBoot[[colIdx]] <- devFacs[colIdx] + (scaleFactors[[colIdx]] * residBoot[[colIdx]] * sigma[colIdx]) / sqrt(prevCs[[colIdx]])
                } else {
                    indivDevFacBoot[[colIdx]] <- devFacs[colIdx] + (residBoot[[colIdx]] * sigma[colIdx]) / sqrt(prevCs[[colIdx]])
                }

                devFacBoot[colIdx] <- sum(indivDevFacBoot[[colIdx]] * prevCs[[colIdx]]) / sum(prevCs[[colIdx]])

                if (colIdx != nDev - 1) {
                    df <- nDev - colIdx - 1
                    sigmaBoot[colIdx] <- sqrt(sum((prevCs[[colIdx]] * (indivDevFacBoot[[colIdx]] - devFacs[colIdx])**2) / df))
                } else {
                    sigmaBoot[colIdx] <- min(
                        sigmaBoot[(length(sigmaBoot) - 1)]**2,
                        sigmaBoot[length(sigmaBoot)]**2,
                        sigmaBoot[length(sigmaBoot)]**4 / sigmaBoot[(length(sigmaBoot) - 1)]**2
                    )
                }
            }
        }

        if (distribution == "normal") {
            # complete the triangle
            triangleBoot <- inputTriangle
            for (diagIdx in 1:(nDev - 1)) {
                for (rowIdx in (diagIdx + 1):nDev) {
                    colIdx <- nDev + diagIdx + 1 - rowIdx
                    # off-diagonal elements satisfy i + j = I + 1 + (diagonal number)
                    triangleBoot[rowIdx, colIdx] <-
                        rnorm(
                            1,
                            triangleBoot[rowIdx, colIdx - 1] * devFacBoot[colIdx - 1],
                            sqrt(triangleBoot[rowIdx, colIdx - 1]) * sigmaBoot[colIdx - 1])
                }
            }
        } else if (distribution == "gamma") {
            triangleBoot <- inputTriangle
            for (diagIdx in 1:(nDev - 1)) {
                for (rowIdx in (diagIdx + 1):nDev) {
                    colIdx <- nDev + diagIdx + 1 - rowIdx
                    alpha <- devFacs[colIdx - 1]**2 * triangleBoot[rowIdx, colIdx - 1] / sigma[colIdx - 1]**2
                    beta <- devFacs[colIdx - 1] / (sigma[colIdx - 1]**2)
                    triangleBoot[rowIdx, colIdx] <-
                        rgamma(1, shape = alpha, rate = beta)
                }
            }
        }

        # if the triangle contains undefined values, the simulated reserve became negative at some point, and we throw away the triangle
        if (!(NaN %in% triangleBoot)) {
            latest <- triangleBoot[col(triangleBoot) + row(triangleBoot) == nDev + 1]
            reserve[iBoot] <- sum(triangleBoot[, nDev] - latest)
        } else {
            reserve[iBoot] <- NaN
            nNaN <- nNaN + 1
        }
    }
    attr(reserve, "nNaN") <- nNaN
    return(reserve)
}

singleOutlier <- function(outlierColIdx, outlierRowIdx, factor, ..., initCol, devFacs, sigma) {

    nDev <- length(initCol)

    if (outlierColIdx == 1) {
        stop("Outlier column index must be greater than 1.")
    }

    claimsTriangle <- matrix(ncol = nDev, nrow = nDev)
    claimsTriangle[, 1] <- initCol

    for (colIdx in 2:nDev) {
        for (rowIdx in setdiff(1:(nDev + 1 - colIdx), outlierRowIdx)) {
            prevC <- claimsTriangle[rowIdx, colIdx - 1]
            claimsTriangle[rowIdx, colIdx] <-
                rnorm(1, devFacs[colIdx - 1] * prevC, sigma[colIdx - 1] * sqrt(prevC))
        }
    }

    if (outlierColIdx > 2) {
        for (colIdx in 2:(outlierColIdx - 1)) {
            prevC <- claimsTriangle[outlierRowIdx, colIdx - 1]
            claimsTriangle[outlierRowIdx, colIdx] <-
                rnorm(1, devFacs[colIdx - 1] * prevC, sigma[colIdx - 1] * sqrt(prevC))
        }
    }

    prevC <- claimsTriangle[outlierRowIdx, outlierColIdx - 1]
    claimsTriangle[outlierRowIdx, outlierColIdx] <-
        rnorm(1,
            factor * devFacs[outlierColIdx - 1] * prevC,
            sigma[outlierColIdx - 1] * sqrt(prevC))

    if (outlierColIdx < nDev) {
        for (colIdx in (outlierColIdx + 1):(nDev + 1 - outlierRowIdx)) {
            prevC <- claimsTriangle[outlierRowIdx, colIdx - 1]
            claimsTriangle[outlierRowIdx, colIdx] <-
                rnorm(1, devFacs[colIdx - 1] * prevC, sigma[colIdx - 1] * sqrt(prevC))
        }
    }
    return(claimsTriangle)
}

# bootReserveGamma <- function(inputTriangle, nBoot) {

#     nDev <- ncol(inputTriangle)

#     mackResults <- suppressWarnings(MackChainLadder(inputTriangle))

#     devFacs <- mackResults$f[-nDev]
#     sigma <- mackResults$sigma

#     indivDevFacs <- list()
#     prevCs <- list()
#     resids <- list()

#     for (idx in 1:(nDev - 1)) {
#         indivDevFacs[[idx]] <- unname(inputTriangle[1:(nDev - idx), idx + 1] / inputTriangle[1:(nDev - idx), idx])
#         prevCs[[idx]] <- unname(inputTriangle[1:(nDev - idx), idx])

#         mean <- devFacs[idx - 1] ** 2 * prevCs[[idx]] / sigma[idx - 1] ** 2
#         std <- devFacs[idx - 1] * prevCs[[idx]] / sigma[idx - 1] ** 2

#         resids[[idx]] <- (indivDevFacs[[idx]] - mean) / std
#     }

#     #sample residuals and put them in proper list format
#     reserveBoot <- c()
#     for (iBoot in 1:nBoot) {
#         residBoot <- list()
#         residSample <- sample(unlist(resids), replace = TRUE)
#         for (j in 1:(nDev - 1)) {
#             residBoot[[j]] <- residSample[1:(nDev - j)]
#             residSample <- residSample[-1:-(nDev - j)]
#         }

#         #compute bootstrapped quantities
#         FBoot <- list()
#         fBoot <- c()
#         sigmaBoot <- c()
#         for (j in 1:(nDev - 1)) {
#             mean <- devFacs[idx - 1] ** 2 * prevCs[[idx]] / sigma[idx - 1] ** 2
#             std <- devFacs[idx - 1] * prevCs[[idx]] / sigma[idx - 1]**2

#             FBoot[[j]] <- residBoot[[j]] * std + mean
#             fBoot[j] <- sum(FBoot[[j]] * prevCs[[j]]) / sum(prevCs[[j]])
#             sigmaBoot[j] <- mean(prevCs[[j]] * (FBoot[[j]] - fBoot[j])**2)
#         }

#         #complete the triangle
#         bootTriangle <- inputTriangle
#         for (diagIdx in 1:(nDev-1)) {
#             for (rowIdx in (diagIdx + 1):nDev){
#                 colIdx <- nDev + diagIdx + 1 - rowIdx
#                 bootTriangle[rowIdx, colIdx] <-
#                     rgamma(1,
#                     devFacs[colIdx - 1]**2 * bootTriangle[rowIdx, colIdx - 1] / sigma[colIdx - 1]**2,
#                     devFacs[colIdx - 1] / (sigma[colIdx - 1]**2))
#                 }
#             }
#         #if the triangle contains undefined values, the simulated reserve became negative at some point,
#         # and we throw away the triangle
#         if (!(NA %in% bootTriangle) && !(NaN %in% bootTriangle)) {
#             latest <- bootTriangle[col(bootTriangle) + row(bootTriangle) == nDev + 1]
#             reserve <- sum(bootTriangle[, nDev] - latest)
#             reserveBoot <- c(reserveBoot, reserve)
#         }
#     }
#     return(reserveBoot)
# }

# singleOutlierGamma <- function(nDev, initMean, initStd, devFac, sigma, outlierColIdx, outlierRowIdx, pert = 1.1) {

#     initCol <- rgamma(nDev, initMean**2 / initStd**2, initMean / initStd**2)

#     claimsTriangle <- matrix(ncol = nDev, nrow = nDev)
#     claimsTriangle[, 1] <- initCol

#     for (colIdx in 2:nDev) {
#         for (rowIdx in 1:(nDev + 1 - colIdx)) {
#             prevC <- claimsTriangle[rowIdx, colIdx - 1]
#             claimsTriangle[rowIdx, colIdx] <-
#                 rgamma(1,
#                     devFac[colIdx - 1]**2 * claimsTriangle[rowIdx, colIdx - 1] / sigma[colIdx - 1]**2,
#                     devFac[colIdx - 1] / (sigma[colIdx - 1]**2))
#         }
#     }

#     claimsTriangle[outlierRowIdx, -1] <- NA

#     if (outlierColIdx > 2) {
#         for (colIdx in 2:(outlierColIdx - 1)) {
#             prevC <- claimsTriangle[outlierRowIdx, colIdx - 1]
#             claimsTriangle[outlierRowIdx, colIdx] <-
#                 rgamma(1,
#                     devFac[colIdx - 1]**2 * prevC / sigma[colIdx - 1]**2,
#                     devFac[colIdx - 1] / sigma[colIdx - 1]**2)
#         }
#     }

#     prevC <- claimsTriangle[outlierRowIdx, outlierColIdx - 1]
#     claimsTriangle[outlierRowIdx, outlierColIdx] <-
#         rgamma(1,
#             (pert * devFac[outlierColIdx - 1])**2 * prevC / sigma[outlierColIdx - 1]**2,
#             (pert * devFac[outlierColIdx - 1]) / sigma[outlierColIdx - 1]**2)

#     if (outlierColIdx < nDev) {
#         for (colIdx in (outlierColIdx + 1):(nDev + 1 - outlierRowIdx)) {
#             prevC <- claimsTria13560.28ngle[outlierRowIdx, colIdx - 1]
#             claimsTriangle[outlierRowIdx, colIdx] <-
#                 rgamma(1,
#                     devFac[colIdx - 1]**2 * prevC / sigma[colIdx - 1]**2,
#                     devFac[colIdx - 1] / sigma[colIdx - 1]**2)
#         }
#     }
#     return(claimsTriangle)
# }
