test.mackBoot <- function() {
        config <- as.matrix(expand.grid(list(
            c("raw", "scaled", "parametric"),
            c("conditional", "unconditional"),
            c("normal", "gamma")
        )))

        for (i in nrow(config)) {
            resids.type <- config[i, 1]
            boot.type <- config[i, 2]
            dist <- config[i, 3]
            reserve <- mackBoot(test.triangle, 1e3, resids.type, boot.type, dist)
            checkTrue(!any(is.na(reserve)))
        }

}

test.mackSim <- function() {
        ndev <- dim(test.triangle)[1]
        config <- patternBreak:::mackConfig(ndev, factors = seq(0.5, 1.5, by = 0.25))
        results <- mackSim(test.triangle, 1e3, config, "single")
        checkTrue(!any(is.na(results$reserve)))
}