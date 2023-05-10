test.rng <- function() {
    checkEquals(patternBreak:::validate_rng(1e3), "Success")
}