# #' Compute bootstrap reserve based on Mack's model.
# #'
# #' FUNCTION_DESCRIPTION
# #'
# #' @param triangle Run-off triangle.
# #' @param nboot Number of bootstrap samples.
# #' @param resids.type Type of residuals.
# #' @param boot.type Type of bootstrap.
# #' @param dist Distribution.
# #'
# #' @return A vector of bootstrap reserve estimates.
# #' @examples
# #' # ADD_EXAMPLES_HERE
# #' @export
# mackBoot <- function(triangle, nboot, resids.type, boot.type, dist) {

#     key <- list(
#         normal = 1L,
#         gamma = 2L,
#         raw = 1L,
#         scaled = 2L,
#         parametric = 3L,
#         conditional = 1L,
#         unconditional = 2L)

#     if (any(is.na(triangle))) {
#         triangle[is.na(triangle)] <- 0
#     }

#     if (!is.double(triangle)) storage.mode(triangle) <- "double"
#     if (!is.integer(nboot)) storage.mode(nboot) <- "integer"

#     resids.type <- key[[resids.type]]
#     boot.type <- key[[boot.type]]
#     dist <- key[[dist]]

#     .Call("mack_boot_cpp", triangle, nboot, resids.type, boot.type, dist)
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
mackSim <- function(triangle, nboot, config, type) {

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

    config[, c("resids.type", "boot.type", "dist") :=
        list(
            sapply(resids.type, function(resids.type) key[[resids.type]]),
            sapply(boot.type, function(boot.type) key[[boot.type]]),
            sapply(dist, function(dist) key[[dist]])
        )]

    config <- as.matrix(config)

    if (any(is.na(triangle))) {
        triangle[is.na(triangle)] <- 0
    }

    if (!is.double(triangle)) storage.mode(triangle) <- "double"
    if (!is.integer(nboot)) storage.mode(nboot) <- "integer"
    if (!is.double(config)) storage.mode(config) <- "double"
    
    result <- .mackSim(triangle, nboot, config, type)

    return(result)

}
