test.mackBoot <- function() {

        reserve <- mackBoot(test.triangle, 1e3, "parametric", "conditional", "normal")

        checkTrue(!any(is.na(reserve)))
}

# test.mackSim <- function() {

#         ndev <- dim(test.triangle)[1]
#         config <- patternBreak:::mackConfig(ndev, factors = seq(0.5, 1.5, by = 0.25))

#         results <- mackSim(test.triangle, 1e3, config, 1)

# }