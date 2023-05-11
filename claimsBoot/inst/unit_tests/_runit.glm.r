test.glmBoot <- function() {
    reserve <- glmBoot(test.triangle, 1e3)
    checkTrue(!any(is.na(reserve)))
}

test.glmSim <- function() {
    ndev <- dim(test.triangle)[1]
    config <- claimsBoot:::glmConfig(ndev, factors = seq(0.5, 1.5, by = 0.25))
    test.triangle <- cum2incr(test.triangle)
    results <- glmSim(test.triangle, 1e3, config, "single")
    checkTrue(!any(is.na(results$reserve)))
}