load(file.path(system.file(package = "patternBreak"), "data", "test_data.RData"))

test.mackBoot <- function() {
    for (triangle in datasets) {
        reserve <- mackBoot(triangle, 1e3, "parametric", "conditional", "normal")

        checkTrue(!any(is.na(reserve)))
    }
}

test.mackSim <- function() {

        triangle <- datasets[[1]]

        ndev <- dim(triangle)[1]
        config <- mackConfig(ndev, factors = seq(0.5, 1.5, by = 0.25))

        results <- mackSim(triangle, 1e3, config, 1)

}