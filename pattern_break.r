reserveBoot <- function(inputTriangle, nBoot) {
    inputTriangle %<>% unname(inputTriangle)
    nDev <- ncol(inputTriangle)

    mackResults <- suppressWarnings(MackChainLadder(inputTriangle))

    devFacs <- mackResults$f[-nDev]
    sigma <- mackResults$sigma

    indivDevFacs <- list()
    prevCs <- list()
    resids <- list()

    for (idx in 1:(nDev - 1)) {
        indivDevFacs[[idx]] <- inputTriangle[1:(nDev - idx), idx + 1] / inputTriangle[1:(nDev - idx), idx]
        prevCs[[idx]] <- inputTriangle[1:(nDev - idx), idx]
        resids[[idx]] <- (indivDevFacs[[idx]] - devFacs[idx]) * sqrt(prevCs[[idx]]) / sigma[idx]
    }

    results <- tibble(
        indivDevFacBoot = list(),
        devFacBoot = list(),
        sigmaBoot = list(),
        triangleBoot = list(),
        reserve = numeric()
    )

    for (iBoot in 1:nBoot) {

        #sample residuals and put them in proper list format
        residBoot <- list()
        residSample <- sample(unlist(resids), replace = TRUE)
        for (j in 1:(nDev - 1)) {
            residBoot[[j]] <- residSample[1:(nDev - j)]
            residSample <- residSample[-1:-(nDev - j)]
        }

        #compute bootstrapped quantities
        indivDevFacBoot <- list()
        devFacBoot <- c()
        sigmaBoot <- c()
        for (j in 1:(nDev - 1)) {
            indivDevFacBoot[[j]] <- devFacs[j] + (residBoot[[j]] * sigma[j]) / sqrt(prevCs[[j]])
            devFacBoot[j] <- sum(indivDevFacBoot[[j]] * prevCs[[j]]) / sum(prevCs[[j]])
            sigmaBoot[j] <- mean(prevCs[[j]] * (indivDevFacBoot[[j]] - devFacBoot[j])**2)
        }

        #complete the triangle
        triangleBoot <- inputTriangle
        for (diagIdx in 1:(nDev - 1)) {
            for (rowIdx in (diagIdx + 1):nDev){
                colIdx <- nDev + diagIdx + 1 - rowIdx
                #off-diagonal elements satisfy i + j = I + 1 + (diagonal number)
                triangleBoot[rowIdx, colIdx] <-
                        rnorm(
                            1,
                            triangleBoot[rowIdx, colIdx - 1] * devFacBoot[colIdx - 1],
                            sqrt(triangleBoot[rowIdx, colIdx - 1]) * sigmaBoot[colIdx - 1])
                if (triangleBoot[rowIdx, colIdx] < 0) {
                    browser()
                }
                }
            }
        #if the triangle contains undefined values, the simulated reserve became negative at some point,
        # and we throw away the triangle
        if (!(NA %in% triangleBoot) && !(NaN %in% triangleBoot)) {
            latest <- triangleBoot[col(triangleBoot) + row(triangleBoot) == nDev + 1]
            reserve <- sum(triangleBoot[, nDev] - latest)
        } else {
            reserve <- NA
        }
    # browser()
    results %<>% add_row(
        indivDevFacBoot = list(indivDevFacBoot),
        devFacBoot = list(devFacBoot),
        sigmaBoot = list(sigmaBoot),
        triangleBoot = list(triangleBoot),
        reserve = reserve)
  
    }

    return(results)
}

singleOutlier <- function(nDev, initMean, initStd, devFac, sigma, outlierColIdx, outlierRowIdx, pert = 1.1) {
    initCol <- rnorm(nDev, initMean, initStd)

    claimsTriangle <- matrix(ncol = nDev, nrow = nDev)
    claimsTriangle[, 1] <- initCol

    for (colIdx in 2:nDev) {
        for (rowIdx in 1:(nDev + 1 - colIdx)) {
            prevC <- claimsTriangle[rowIdx, colIdx - 1]
            claimsTriangle[rowIdx, colIdx] <-
                rnorm(1, devFac[colIdx - 1] * prevC, sigma[colIdx - 1] * sqrt(prevC))
        }
    }
    claimsTriangle[outlierRowIdx, -1] <- NA

    if (outlierColIdx > 2) {
        for (colIdx in 2:(outlierColIdx - 1)) {
            prevC <- claimsTriangle[outlierRowIdx, colIdx - 1]
            claimsTriangle[outlierRowIdx, colIdx] <-
                rnorm(1, devFac[colIdx - 1] * prevC, sigma[colIdx - 1] * sqrt(prevC))
        }
    }
    prevC <- claimsTriangle[outlierRowIdx, outlierColIdx - 1]
    claimsTriangle[outlierRowIdx, outlierColIdx] <-
        rnorm(1,
            pert * devFac[outlierColIdx - 1] * prevC,
            sigma[outlierColIdx - 1] * sqrt(prevC))

    if (outlierColIdx < nDev) {
        for (colIdx in (outlierColIdx + 1):(nDev + 1 - outlierRowIdx)) {
            prevC <- claimsTriangle[outlierRowIdx, colIdx - 1]
            claimsTriangle[outlierRowIdx, colIdx] <-
                rnorm(1, devFac[colIdx - 1] * prevC, sigma[colIdx - 1] * sqrt(prevC))
        }
    }
    return(claimsTriangle)
}

bootReserveGamma <- function(inputTriangle, nBoot) {

    nDev <- ncol(inputTriangle)

    mackResults <- suppressWarnings(MackChainLadder(inputTriangle))

    devFacs <- mackResults$f[-nDev]
    sigma <- mackResults$sigma

    indivDevFacs <- list()
    prevCs <- list()
    resids <- list()

    for (idx in 1:(nDev - 1)) {
        indivDevFacs[[idx]] <- unname(inputTriangle[1:(nDev - idx), idx + 1] / inputTriangle[1:(nDev - idx), idx])
        prevCs[[idx]] <- unname(inputTriangle[1:(nDev - idx), idx])

        mean <- devFacs[idx - 1] ** 2 * prevCs[[idx]] / sigma[idx - 1] ** 2
        std <- devFacs[idx - 1] * prevCs[[idx]] / sigma[idx - 1] ** 2

        resids[[idx]] <- (indivDevFacs[[idx]] - mean) / std
    }

    #sample residuals and put them in proper list format
    reserveBoot <- c()
    for (iBoot in 1:nBoot) {
        residBoot <- list()
        residSample <- sample(unlist(resids), replace = TRUE)
        for (j in 1:(nDev - 1)) {
            residBoot[[j]] <- residSample[1:(nDev - j)]
            residSample <- residSample[-1:-(nDev - j)]
        }

        #compute bootstrapped quantities
        FBoot <- list()
        fBoot <- c()
        sigmaBoot <- c()
        for (j in 1:(nDev - 1)) {
            mean <- devFacs[idx - 1] ** 2 * prevCs[[idx]] / sigma[idx - 1] ** 2
            std <- devFacs[idx - 1] * prevCs[[idx]] / sigma[idx - 1]**2

            FBoot[[j]] <- residBoot[[j]] * std + mean
            fBoot[j] <- sum(FBoot[[j]] * prevCs[[j]]) / sum(prevCs[[j]])
            sigmaBoot[j] <- mean(prevCs[[j]] * (FBoot[[j]] - fBoot[j])**2)
        }

        #complete the triangle
        bootTriangle <- inputTriangle
        for (diagIdx in 1:(nDev-1)) {
            for (rowIdx in (diagIdx + 1):nDev){
                colIdx <- nDev + diagIdx + 1 - rowIdx
                bootTriangle[rowIdx, colIdx] <-
                    rgamma(1,
                    devFacs[colIdx - 1]**2 * bootTriangle[rowIdx, colIdx - 1] / sigma[colIdx - 1]**2,
                    devFacs[colIdx - 1] / (sigma[colIdx - 1]**2))
                }
            }
        #if the triangle contains undefined values, the simulated reserve became negative at some point,
        # and we throw away the triangle
        if (!(NA %in% bootTriangle) && !(NaN %in% bootTriangle)) {
            latest <- bootTriangle[col(bootTriangle) + row(bootTriangle) == nDev + 1]
            reserve <- sum(bootTriangle[, nDev] - latest)
            reserveBoot <- c(reserveBoot, reserve)
        }
    }
    return(reserveBoot)
}

singleOutlierGamma <- function(nDev, initMean, initStd, devFac, sigma, outlierColIdx, outlierRowIdx, pert = 1.1) {

    initCol <- rgamma(nDev, initMean**2 / initStd**2, initMean / initStd**2)

    claimsTriangle <- matrix(ncol = nDev, nrow = nDev)
    claimsTriangle[, 1] <- initCol

    for (colIdx in 2:nDev) {
        for (rowIdx in 1:(nDev + 1 - colIdx)) {
            prevC <- claimsTriangle[rowIdx, colIdx - 1]
            claimsTriangle[rowIdx, colIdx] <-
                rgamma(1,
                    devFac[colIdx - 1]**2 * claimsTriangle[rowIdx, colIdx - 1] / sigma[colIdx - 1]**2,
                    devFac[colIdx - 1] / (sigma[colIdx - 1]**2))
        }
    }

    claimsTriangle[outlierRowIdx, -1] <- NA

    if (outlierColIdx > 2) {
        for (colIdx in 2:(outlierColIdx - 1)) {
            prevC <- claimsTriangle[outlierRowIdx, colIdx - 1]
            claimsTriangle[outlierRowIdx, colIdx] <-
                rgamma(1,
                    devFac[colIdx - 1]**2 * prevC / sigma[colIdx - 1]**2,
                    devFac[colIdx - 1] / sigma[colIdx - 1]**2)
        }
    }

    prevC <- claimsTriangle[outlierRowIdx, outlierColIdx - 1]
    claimsTriangle[outlierRowIdx, outlierColIdx] <-
        rgamma(1,
            (pert * devFac[outlierColIdx - 1])**2 * prevC / sigma[outlierColIdx - 1]**2,
            (pert * devFac[outlierColIdx - 1]) / sigma[outlierColIdx - 1]**2)

    if (outlierColIdx < nDev) {
        for (colIdx in (outlierColIdx + 1):(nDev + 1 - outlierRowIdx)) {
            prevC <- claimsTriangle[outlierRowIdx, colIdx - 1]
            claimsTriangle[outlierRowIdx, colIdx] <-
                rgamma(1,
                    devFac[colIdx - 1]**2 * prevC / sigma[colIdx - 1]**2,
                    devFac[colIdx - 1] / sigma[colIdx - 1]**2)
        }
    }
    return(claimsTriangle)
}
