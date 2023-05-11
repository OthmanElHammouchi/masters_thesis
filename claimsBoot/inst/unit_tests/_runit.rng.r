test.rng <- function() {
    checkEquals(claimsBoot:::validate_rng(1e3), "Success")
}