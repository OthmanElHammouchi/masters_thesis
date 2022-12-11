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

suppressWarnings({
    triangle <- unname(UKMotor)
    initcol <- triangle[, 1]
    devfacs <- MackChainLadder(triangle)$f
    sigmas <- MackChainLadder(triangle)$sigma
})

points <- cbind(row(triangle)[!is.na(triangle)], col(triangle)[!is.na(triangle)])

factor <- seq(0.5, 1.5, by = 0.25)
resids.type <- c("parametric", "raw", "scaled")
boot.type <- c("conditional", "unconditional")
dist <- c("normal", "gamma")

config <- genConfig(points, factor, points, resids.type, boot.type, dist)

names(config) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx", "resids.type", "boot.type", "dist")

config <- config[!(config$dist == "gamma" &  config$resids.type == "parametric"), ]

config <- config[config$outlier.colidx != 1 & config $excl.colidx != 1, ]

nconfig <- nrow(config)
nboot <- 1e3

results <- config
results[, "reserve"] <- list(rep(0, nconfig))

progress.bar <- txtProgressBar(min = 0, max = nconfig, initial = 0, style = 3)

for (rowidx in seq_len(nconfig)) {

    setTxtProgressBar(progress.bar, rowidx)

    row <- config[rowidx, ]

    triangle <- singleOutlier(row$outlier.rowidx, row$outlier.colidx, row$factor,
        initcol = initcol,
        devfacs = devfacs,
        sigmas = sigmas,
        dist = row$dist)

    results$reserve[rowidx] <- list(reserveBootFortran(triangle, nboot,
        resids.type = row$resids.type,
        boot.type = row$boot.type,
        dist = row$dist,
        excl.resids = matrix(c(row$excl.rowidx, row$excl.colidx), nrow = 2)))
}

close(progress.bar)

saveRDS(results, "results/data_objects/single_outlier.RDS")

means <- sapply(results$reserve, mean, na.rm = TRUE)
