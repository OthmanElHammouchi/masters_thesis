library(data.table)
source("R/helpers.r")
source("R/fortran_wrappers.r")
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
ndev <- nrow(triangle)
nboot <- 1e3

#single outlier

config <- patternBreak:::mackConfig(ndev, factors = seq(0.5, 1.5, by = 0.25))

results <- reserveSim(triangle, 1e3, config, "single")

# saveRDS(results, "results/data_objects/single_outlier.RDS")

# #calendar year outlier

config <- patternBreak:::mackConfig(ndev, factors = seq(0.5, 1.5, by = 0.25))

results <- reserveSim(triangle, 1e3, config, "calendar")

# saveRDS(results, "results/data_objects/calendar_outlier.RDS")

# #accident year outlier

rowidx <- 1:ndev

factor <- seq(0.5, 1.5, by = 0.25)
resids.type <- c("raw", "scaled", "parametric")
boot.type <- c("conditional", "unconditional")
dist <- c("normal", "gamma")

config <- genConfig(rowidx, factor, rowidx, resids.type, boot.type, dist)

names(config) <- c("outlier.rowidx", "factor", "excl.rowidx", "resids.type", "boot.type", "dist")

config <- config[!(dist == "gamma" & resids.type == "parametric")]

results <- reserveSim(triangle, 1e3, config, "origin")

# saveRDS(results, "results/data_objects/origin_outlier.RDS")
