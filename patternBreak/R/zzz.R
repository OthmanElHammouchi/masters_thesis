.onLoad <- function(libname, pkgname) {
    if (is.na(Sys.getenv("COLUMNS", unset = NA))) {
        Sys.setenv(COLUMNS = getOption("width"))
    }
}

.global <- new.env(parent = emptyenv())
setPackageName("patternBreak", .global)

.global$key <- list(
    single = 1L,
    calendar =  2L,
    origin = 3L,
    normal = 1L,
    gamma = 2L,
    standardised = 1L,
    modified = 2L,
    parametric = 3L,
    studentised = 4L,
    conditional = 1L,
    unconditional = 2L
)