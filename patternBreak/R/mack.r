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

    type <- patternBreak:::.global$key[[type]]

    config[, c("resids.type", "boot.type", "dist") :=
        list(
            sapply(resids.type, function(resids.type) .global$key[[resids.type]]),
            sapply(boot.type, function(boot.type) .global$key[[boot.type]]),
            sapply(dist, function(dist) .global$key[[dist]])
        )]

    names.keep <- names(config)
    config <- as.matrix(config)

    if (any(is.na(triangle))) {
        triangle[is.na(triangle)] <- 0
    }

    if (!is.double(triangle)) storage.mode(triangle) <- "double"
    if (!is.integer(nboot)) storage.mode(nboot) <- "integer"
    if (!is.double(config)) storage.mode(config) <- "double"

    result <- .mackSim(triangle, nboot, config, type)

    result <- as.data.table(result)
    names(result) <- c(names.keep, "reserve")

    return(result)

}

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
mackBoot <- function(triangle, nboot, resids.type, boot.type, dist) {

    resids.type <- patternBreak:::.global$key[[resids.type]]
    boot.type <- patternBreak:::.global$key[[boot.type]]
    dist <- patternBreak:::.global$key[[dist]]

    if (!is.double(triangle)) storage.mode(triangle) <- "double"
    if (!is.integer(nboot)) storage.mode(nboot) <- "integer"

    result <- .mackBoot(triangle, nboot, resids.type, boot.type, dist)

}