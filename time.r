library(ggplot2)
library(data.table)
suppressPackageStartupMessages(library(ChainLadder))
source("R/pattern_break.r")


triangle <- unname(UKMotor)
ndev <- ncol(triangle)
nboot <- 1e3

initcol <- triangle[, 1]
devfacs <- rep(0, ndev - 1)
sigmas <- rep(0, ndev - 1)
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
}

# comparing single run of R code vs Fortran

exclude.resids <- list(c(1, 1), c(1, 2))

res <- reserveBoot(gamma.triangle, 1e3, distribution = "gamma", bootstrap.type = "conditional", resids.type = "scaled")

times.r <- rep(0, 10)
times.fortran <- rep(0, 10)

for (i in seq_len(length(times.r))) {

start <- Sys.time()
reserve <- reserveBoot(triangle, 1e3, bootstrap.type = "unconditional")
end <- Sys.time()

times.r[i] <- end - start

start <- Sys.time()
reserveFortran <- reserveBootFortran(triangle, 1e3, distribution = "gamma")
end <- Sys.time()

times.fortran[i] <- end - start
}

result <- mean(times.r)
result.fortran <- mean(times.fortran)

#comparing simulation loop in R vs Fortran

########################## R ##############################

points <- cbind(row(triangle)[!is.na(triangle)], col(triangle)[!is.na(triangle)])

factor <- seq(0.5, 1.5, by = 0.25)
resids.type <- c(1, 2, 3)
boot.type <- c(1, 2)
dist <- c(1, 2)

config <- genConfig(points, factor, points, resids.type, boot.type, dist)

names(config) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx", "resids.type", "boot.type", "dist")

config <- config[!(config$dist == 2 & config$resids.type == 3), ]
config <- config[config$outlier.colidx != 1 & config$excl.colidx != 1, ]

nconfig <- nrow(config)

results <- config
results[, "reserve"] <- list(rep(0, nconfig))

start <- Sys.time()

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

end <- Sys.time()

time.r <- start - end

####################### Fortran ##########################

points <- cbind(row(triangle)[!is.na(triangle)], col(triangle)[!is.na(triangle)])

factor <- seq(0.5, 1.5, by = 0.25)
resids.type <- c(1, 2, 3)
boot.type <- c(1, 2)
dist <- c(1, 2)

config <- genConfig(points, factor, points, resids.type, boot.type, dist)

names(config) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx", "resids.type", "boot.type", "dist")

config <- config[!(config$dist == 2 & config$resids.type == 3), ]
config <- config[config$outlier.colidx != 1 & config$excl.colidx != 1, ]

start <- Sys.time()

results <- reserveSimFortran(triangle, nboot, config)


end <- Sys.time()

time.fortran <- start - end

names(results) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx", "resids.type", "boot.type", "dist", "reserve")
