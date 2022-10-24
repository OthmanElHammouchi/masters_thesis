library(tidyverse)
library(ggplot2)
library(magrittr)
library(ggfortify)
library(ChainLadder)
library(patchwork)
library(MASS)

nDev <- 9

sigma <- rnorm(nDev - 1, 1, 0.1)
devFac <- c(1.4, 1.35, 1.25, 1.2, 1.1, 1.07, 1.05, 1.01)
initSigma <- 100

genMackData <- function(nDev, sigma, devFac, nCoh = nDev, initDist = rnorm, ...) {
    # generate synthetic claims triangle satisfying assumptions of Mack's method
    distParams <- list(...)
    initCol <- initDist(nCoh, distParams[[1]], distParams[[2]])

    claimsTriangle <- matrix(ncol = nDev, nrow = nCoh) %>%
        data.frame()
    claimsTriangle[, 1] <- initCol
    for (idx in 2:nDev) {
        prevClaims <- claimsTriangle[1:(nDev - idx + 1), idx - 1]
        newCol <- mvrnorm(
            1,
            devFac[idx - 1] * prevClaims,
            diag(prevClaims * sigma[idx - 1]**2, nrow = length(prevClaims)))
        claimsTriangle[1:(nDev - idx + 1), idx] <- newCol
    }
    return(claimsTriangle)
}

claimsTriangle <- genMackData(nDev, sigma, devFac, nCoh = nDev, initDist = rnorm, 2000, initSigma)

#we want to study what happens if our assumptions are grossly violated
#first, we need a function to generate a bootstrapped reserve distribution

mackEst <- function(inputTriangle) {
    mackResults <- suppressWarnings(MackChainLadder(inputTriangle))

    f <- mackResults$f[-nDev]
    sigma <- mackResults$sigma

    F_list <- list()
    prevCs <- list()
    resids <- list()

    for (idx in 1:(nDev - 1)) {
        F_list[[idx]] <- unname(claimsTriangle[1:(nDev - idx), idx + 1] / claimsTriangle[1:(nDev - idx), idx])
        prevCs[[idx]] <- unname(claimsTriangle[1:(nDev - idx), idx])
        resids[[idx]] <- (F_list[[idx]] - f[idx]) * sqrt(prevCs[[idx]]) / sigma[idx]
    }
    return(list(devFacs = f, stdevs = sigma, indivDevFacs = F_list, prevCs = prevCs, resids = resids))
}

bootReserve <- function(inputTriangle, nBoot) {

    res <- mackEst(inputTriangle)

    resids <- res$resids
    prevCs <- res$prevCs
    #indivDevFacs <- res$indivDevFacs
    devFacs <- res$devFacs

    #sample residuals and put them in proper list format
    reserveBoot <- c()
    for (i in 1:nBoot) {
        residBoot <- list()
        residSample <- sample(unlist(resids), replace = TRUE)
        for (j in 1:(nDev - 1)) {
            print(j)
            print(resids)
            residBoot[[j]] <- residSample[1:(nDev - j)]
            residSample <- residSample[-1:-(nDev - j)]
        }

        #compute bootstrapped quantities
        FBoot <- list()
        fBoot <- c()
        sigmaBoot <- c()
        for (j in 1:(nDev - 1)) {
            FBoot[[j]] <- devFacs[j] + (residBoot[[j]] * sigma[j]) / sqrt(prevCs[[j]])
            fBoot[j] <- sum(FBoot[[j]] * prevCs[[j]]) / sum(prevCs[[j]])
            sigmaBoot[j] <- mean(prevCs[[j]] * (FBoot[[j]] - fBoot[j])**2)
        }

        #complete triangle
        bootTriangle <- inputTriangle
        for (idx in 1:(nDev - 1)) {
            latest <- bootTriangle[row(bootTriangle) + col(bootTriangle) == nDev + idx]
            newDiag <- -1
                while (!prod(newDiag > 0)) {
                    newDiag <- mvrnorm(
                        1,
                        fBoot[idx:(nDev - 1)] * latest[1:(nDev - idx)],
                        diag(latest[1:(nDev - idx)] * sigmaBoot[idx:(nDev - 1)]**2, nrow = length(latest[1:(nDev - idx)])))
                }
            bootTriangle[row(bootTriangle) + col(bootTriangle) == nDev + idx + 1] <- newDiag
        }
        reserveBoot <- c(reserveBoot, sum(bootTriangle[, nDev]) - sum(inputTriangle[row(inputTriangle) + col(inputTriangle) == nDev + 1]))
    }
    return(reserveBoot)
}

res <- bootReserve(claimsTriangle, 1e3)

#what happens if the stationarity assumption is violated by a single cohort

claimsTriangle[5, 2] <- 2*devFac[1]*claimsTriangle[5, 1] + sigma[1]*sqrt(claimsTriangle[5, 1])*rnorm(1)

for (idx in 3:(nDev - 5 + 1)) {
    claimsTriangle[5, idx] <- devFac[idx - 1]*claimsTriangle[5, idx - 1] + sigma[idx - 1] * sqrt(claimsTriangle[5, idx - 1])*rnorm(1)
}

res <- bootReserve(claimsTriangle, 1e3)

# p <- ggplot(data = tibble(residuals = resid[[1]])) +
#      geom_path(mapping = aes(x = 1:length(residuals), y = residuals)) +
#      geom_point(aes(1:length(residuals), residuals))

# for (idx in 2:(length(resid) - 1)) {
#  p <- p + ggplot(data = tibble(residuals = resid[[idx]])) +
#      geom_path(mapping = aes(x = 1:length(residuals), y = residuals)) +
#      geom_point(aes(1:length(residuals), residuals))
# }

# p <- p + plot_layout(ncol = 4)

# p

# nBoot <- 1e3
