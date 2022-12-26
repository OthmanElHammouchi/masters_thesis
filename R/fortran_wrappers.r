library(data.table)

dist.key <- list(normal = 1, gamma = 2)
resids.type.key <- list(raw = 1, scaled = 2, parametric = 3)
boot.type.key <- list(conditional = 1, unconditional = 2)
config.type.key <- list(single = 1, calendar = 2, origin = 3)

reserveBootFortran <- function(triangle, nboot,
                               resids.type = "raw",
                               boot.type = "conditional",
                               dist = "normal",
                               excl.resids = NULL) {
    nboot <- as.integer(nboot)
    ndev <- ncol(triangle)

    if (is.null(excl.resids)) {
        excl.resids <- matrix(c(0L, 0L), ncol = 1)
    }

    dyn.load("fortran/build/reserve_boot.so")

    triangle[is.na(triangle)] <- 0
    reserve <- rep(0, nboot)

    res <- .Fortran(
        "reserve_boot",
        triangle = triangle,
        nboot = nboot,
        ndev = ndev,
        reserve = reserve,
        resids.type = resids.type.key[[resids.type]],
        boot.type = boot.type.key[[boot.type]],
        dist = dist.key[[dist]],
        excl.resids = excl.resids,
        excl.resids.ncols = ncol(excl.resids)
    )

    results <- as.data.table(res$reserve)

    return(results)
}

reserveSim <- function(triangle, nboot, config, config.type) {
    dyn.load("fortran/build/reserve_sim.so")

    config[, c("resids.type", "boot.type", "dist") :=
        list(
            sapply(resids.type, function(resids.type) resids.type.key[[resids.type]]),
            sapply(boot.type, function(boot.type) boot.type.key[[boot.type]]),
            sapply(dist, function(dist) dist.key[[dist]])
        )]

    ndev <- nrow(triangle)
    config.spec <- numeric(3)
    config.spec[1:2] <- dim(config)
    config.spec[3] <- config.type.key[[config.type]]

    res <- matrix(rep(0, config.spec[1] * (config.spec[2] + 1) * nboot),
        nrow = config.spec[1] * nboot,
        ncol = config.spec[2] + 1
    )

    triangle[is.na(triangle)] <- 0

    results <- .Fortran("reserve_sim",
        unname(triangle),
        as.integer(nboot),
        as.integer(ndev),
        unname(as.matrix(config)),
        as.integer(config.spec),
        res = res
    )

    results <- as.data.table(results$res)

    if (config.type == "single") {

        names(results) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx", "resids.type", "boot.type", "dist", "reserve")

    } else if (config.type == "calendar") {

        names(results) <- c("outlier.diagidx", "factor", "excl.diagidx", "resids.type", "boot.type", "dist", "reserve")

    } else if (config.type == "origin") {
        
        names(results) <- c("outlier.rowidx", "factor", "excl.rowidx", "resids.type", "boot.type", "dist", "reserve")

    }

    results[, c("resids.type", "boot.type", "dist") :=
        list(
            sapply(resids.type, function(key) names(resids.type.key)[key]),
            sapply(boot.type, function(key) names(boot.type.key)[key]),
            sapply(dist, function(key) names(dist.key)[key])
        )]

    return(results)
}
