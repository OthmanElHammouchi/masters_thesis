library(data.table)
library(iterators)
library(parallel)
library(doParallel)
library(foreach)
library(moments)
suppressPackageStartupMessages({
library(ChainLadder)
})
registerDoParallel()

source("R/pattern_break.r")

triangle <- unname(UKMotor)
ndev <- ncol(triangle)

initcol <- triangle[, 1]

devfacs <- rep(0, ndev - 1)
sigmas <- rep(0, ndev - 1)

for (colidx in 2:ndev) {

    nrows <- ndev + 1 - colidx

    model <- lm(y ~ x + 0,
        data = data.frame(x = UKMotor[1:nrows, colidx - 1], y = UKMotor[1:nrows, colidx]),
        weights = 1 / UKMotor[1:nrows, colidx - 1])

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

}

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

    triangle <- with(row, singleOutlier(outlier.rowidx, outlier.colidx, factor,
        initcol = initcol,
        devfacs = devfacs,
        sigmas = sigmas,
        dist = dist))

    results$reserve[rowidx] <- with(row, list(reserveBootFortran(triangle, nboot,
        resids.type = resids.type,
        boot.type = boot.type,
        dist = dist,
        excl.resids = matrix(c(excl.rowidx, excl.colidx), nrow = 2))))
}

results <- as.data.table(results)[, .(reserve = unlist(reserve)), by = setdiff(names(results), "reserve")]

close(progress.bar)

saveRDS(results, "results/data_objects/single_outlier.RDS")