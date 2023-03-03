# # reserveSim <- function(triangle, nboot, config, config.type) {
# #     dyn.load("fortran/build/reserve_sim.so")

# #     #

# #     config[, c("resids.type", "boot.type", "dist") :=
# #         list(
# #             sapply(resids.type, function(resids.type) resids.type.key[[resids.type]]),
# #             sapply(boot.type, function(boot.type) boot.type.key[[boot.type]]),
# #             sapply(dist, function(dist) dist.key[[dist]])
# #         )]

# #     ndev <- nrow(triangle)
# #     config.spec <- numeric(3)
# #     config.spec[1:2] <- dim(config)
# #     config.spec[3] <- config.type.key[[config.type]]

# #     res <- matrix(rep(0, config.spec[1] * (config.spec[2] + 1) * nboot),
# #         nrow = config.spec[1] * nboot,
# #         ncol = config.spec[2] + 1
# #     )

# #     triangle[is.na(triangle)] <- 0

# #     results <- .Fortran("reserve_sim",
# #         unname(triangle),
# #         as.integer(nboot),
# #         as.integer(ndev),
# #         unname(as.matrix(config)),
# #         as.integer(config.spec),
# #         res = res
# #     )

# #     results <- as.data.table(results$res)

# #     if (config.type == "single") {

# #         names(results) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx", "resids.type", "boot.type", "dist", "reserve")

# #     } else if (config.type == "calendar") {

# #         names(results) <- c("outlier.diagidx", "factor", "excl.diagidx", "resids.type", "boot.type", "dist", "reserve")

# #     } else if (config.type == "origin") {

# #         names(results) <- c("outlier.rowidx", "factor", "excl.rowidx", "resids.type", "boot.type", "dist", "reserve")

# #     }

# #     results[, c("resids.type", "boot.type", "dist") :=
# #         list(
# #             sapply(resids.type, function(key) names(resids.type.key)[key]),
# #             sapply(boot.type, function(key) names(boot.type.key)[key]),
# #             sapply(dist, function(key) names(dist.key)[key])
# #         )]

# #     return(results)
# # }


# #' Compute bootstrap reserve based on GLM model.
# #'
# #' FUNCTION_DESCRIPTION
# #'
# #' @param triangle Run-off triangle.
# #' @param nboot DESCRIPTION.
# #'
# #' @return RETURN_DESCRIPTION
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @export
# glmBoot <- function(triangle, nboot) {

#     if (!is.double(triangle)) storage.mode(triangle) <- "double"
#     if (!is.integer(nboot)) storage.mode(nboot) <- "integer"

#     .Call("glm_boot_wrapper", triangle, nboot)
# }

#' Bootstrap simulation based on Mack's model.
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#' @param nboot DESCRIPTION.
#' @param config DESCRIPTION.
#' @param type DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
glmSim <- function(triangle, nboot, config, type) {

    key <- list(
        single = 1L,
        calendar =  2L,
        origin = 3L,
        normal = 1L,
        gamma = 2L,
        raw = 1L,
        scaled = 2L,
        parametric = 3L,
        conditional = 1L,
        unconditional = 2L
    )

    type <- key[[type]]

    config <- as.matrix(config)

    if (any(is.na(triangle))) {
        triangle[is.na(triangle)] <- 0
    }

    if (!is.double(triangle)) storage.mode(triangle) <- "double"
    if (!is.integer(nboot)) storage.mode(nboot) <- "integer"
    if (!is.double(config)) storage.mode(config) <- "double"

    result <- .glmSim(triangle, nboot, config, type)

    return(result)

}
