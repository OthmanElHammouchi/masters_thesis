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

    type <- patternBreak:::.global$key[[type]]

    names.keep <- names(config)
    config <- as.matrix(config)

    if (any(is.na(triangle))) {
        triangle[is.na(triangle)] <- 0
    }

    if (!is.double(triangle)) storage.mode(triangle) <- "double"
    if (!is.integer(nboot)) storage.mode(nboot) <- "integer"
    if (!is.double(config)) storage.mode(config) <- "double"
    if (!is.integer(type)) storage.mode(type) <- "integer"

    result <- .glmSim(triangle, nboot, config, type)
    result <- data.table::as.data.table(result)
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
glmBoot <- function(triangle, nboot) {

    if (!is.double(triangle)) storage.mode(triangle) <- "double"
    if (!is.integer(nboot)) storage.mode(nboot) <- "integer"

    result <- .glmBoot(triangle, nboot)

}