suppressPackageStartupMessages({
    library(MASS)
    library(ggplot2)
    library(ChainLadder)
    library(knitr)
    library(data.table)
    library(iterators)
    library(parallel)
    library(doParallel)
    library(foreach)
    library(moments)

    registerDoParallel()

    source("pattern_break.r")
})

suppressMessages({
    knitr::opts_chunk$set(echo = FALSE)
})

set.seed(5)

# extract simulation parameters from benchmark data
ndev <- ncol(UKMotor)

benchmark <- sapply(1:(ndev - 1), function(i) {

    model <- lm(y ~ x + 0,
        data = data.frame(x = UKMotor[1:(ndev - i), i], y = UKMotor[1:(ndev - i), i + 1]),
        weights = 1 / UKMotor[1:(ndev - i), i])

    return(c(model$coefficients, summary(model)$sigma))
})

devfacs <- benchmark[1, ]
sigmas <- benchmark[2, ]
sigmas[length(sigmas)] <- sqrt(min(sigmas[(length(sigmas) - 2)]**2, sigmas[length(sigmas) - 1]**2, sigmas[length(sigmas) - 1]**4 / sigmas[(length(sigmas) - 2)]**2))

# simulate triangle
initcol <- unname(UKMotor[, 1])

claims.triangle <- matrix(ncol = ndev, nrow = ndev)
claims.triangle[, 1] <- initcol

for (colidx in 2:ndev) {
    for (rowIdx in 1:(ndev + 1 - colidx)) {
        prevC <- claims.triangle[rowIdx, colidx - 1]
        claims.triangle[rowIdx, colidx] <-
            rnorm(1, devfacs[colidx - 1] * prevC, sigmas[colidx - 1] * sqrt(prevC))
    }
}

# create table of experiment configurations
points <- c()

for (colidx in 2:ndev) {
    nrows <- ndev + 1 - colidx
    points <- rbind(points, cbind(1:nrows, rep(colidx, nrows)))
}

pertfac <- seq(0.5, 1.5, by = 0.25)

residual.types <- c("parametric", "raw", "scaled")
bootstrap.types <- c("conditional", "unconditional")
distributions <- c("normal", "gamma")

indices <- expand.grid(
    seq_len(nrow(points)),
    seq_len(length(pertfac)),
    seq_len(nrow(points)),
    seq_len(length(residual.types)),
    seq_len(length(bootstrap.types)),
    seq_len(length(distributions))
)

config <- cbind.data.frame(
    points[indices[, 1], ],
    pertfac[indices[, 2]],
    points[indices[, 3], ],
    residual.types[indices[, 4]],
    bootstrap.types[indices[, 5]],
    distributions[indices[, 6]]
)

names(config) <- c(
    "outlier.rowidx",
    "outlier.colidx",
    "perfac",
    "excl.rowidx",
    "excl.colidx",
    "resids.type",
    "bootstrap.type",
    "distribution"
)

config <- config[!(config$distribution == "gamma" & config$resids.type == "parametric"), ]

# experiment loop
results <- foreach(config = iter(config, by = "row"), .combine = "rbind", .errorhandling = "stop") %do% {

    outlier.rowidx <- config[[1]]
    outlier.colidx <- config[[2]]
    pertfac <- config[[3]]
    excl.rowidx <- config[[4]]
    excl.colidx <- config[[5]]
    resids.type <- config[[6]]
    bootstrap.type <- config[[7]]
    distribution <- config[[8]]

    distortedTriangle <- singleOutlier(outlier.rowidx, outlier.colidx, pertfac,
        initcol = initcol,
        devfacs = devfacs,
        sigmas = sigmas,
    distribution = distribution)

    exclude.residuals <- list(c(excl.rowidx, excl.colidx))

    reserve <- reserveBootFortran(distortedTriangle, 1e3,
        excl.resids = exclude.residuals,
        distribution = distribution,
        resids.type = resids.type,
        bootstrap.type = bootstrap.type
    )

    row <- data.table(
        outlier.rowidx = outlier.rowidx,
        outlier.colidx = outlier.colidx,
        pertfac = pertfac,
        excl.rowidx = excl.rowidx,
        excl.colidx = excl.colidx,
        distribution = distribution,
        resids.type = resids.type,
        bootstrap.type = bootstrap.type,
        reserve = list(reserve))
}

##########################################

# triangle <- unname(UKMotor)
# initcol <- triangle[, 1]
# devfacs <- MackChainLadder(triangle)$f
# sigmas <- MackChainLadder(triangle)$sigmas

# singleOutlier(1, 2, 1.5,
#     initcol = initcol,
#     devfacs = devfacs,
#     sigmas = sigmas, 
#     distribution = "gamma")

# reserve <- reserveBoot(triangle, 1e3,
#     distribution = "gamma",
#     bootstrap.type = "unconditional",
#     resids.type = "raw")