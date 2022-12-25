library(data.table)
source("R/pattern_break.r")
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor

nboot <- 1e3

points <- cbind(row(triangle)[!is.na(triangle)], col(triangle)[!is.na(triangle)])

factor <- seq(0.5, 1.5, by = 0.25)
resids.type <- c(1, 2, 3)
boot.type <- c(1, 2)
dist <- c(1, 2)

config <- genConfig(points, factor, points, resids.type, boot.type, dist)

names(config) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx", "resids.type", "boot.type", "dist")

config <- config[!(config$dist == 2 & config$resids.type == 3), ]
config <- config[config$outlier.colidx != 1 & config$excl.colidx != 1, ]

results <- reserveSimFortran(triangle, 1e3, config)

saveRDS(results, "results/data_objects/single_outlier.RDS")
