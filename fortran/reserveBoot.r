reserveBootFort <- function(triangle, ndev) {
    dyn.load("fortran/reserve_boot.so")
    triangle[is.na(triangle)] <- 0
    boot.triangle <- matrix(0, nrow = ndev, ncol = ndev)
    res <- .Fortran("reserve_boot", ndev, triangle, boot.triangle)
    return(res[[3]])
}
