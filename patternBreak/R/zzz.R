.onLoad <- function(libname, pkgname) {

    if (is.na(Sys.getenv("COLUMNS", unset = NA))) {
    
        Sys.setenv(COLUMNS = getOption("width"))
    
    }
}
