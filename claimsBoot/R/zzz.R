.onLoad <- function(libname, pkgname) {
    if (is.na(Sys.getenv("COLUMNS", unset = NA))) {
        Sys.setenv(COLUMNS = getOption("width"))
    }
}

.global <- new.env(parent = emptyenv())
setPackageName("claimsBoot", .global)

.global$key <- list(
    single = 1L,
    calendar =  2L,
    origin = 3L,
    normal = 1L,
    gamma = 2L,
    poisson = 3L,
    parametric = 1L,
    residuals = 2L,
    pairs = 3L,
    standardised = 1L,
    studentised = 2L,
    lognormal = 3L,
    conditional = 1L,
    unconditional = 2L
)