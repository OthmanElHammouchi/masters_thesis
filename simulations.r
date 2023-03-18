library(patternBreak)
load(file.path(system.file(package = "patternBreak"), "data", "test_triangle.RData"))

triangle <- test.triangle
ndev <- nrow(triangle)
nboot <- 1e3

# Single outlier

config <- patternBreak:::mackConfig(ndev, factors = seq(0.5, 1.5, by = 0.25), type = "single")

results <- mackSim(triangle, nboot, config, "single")

# saveRDS(results, "results/data_objects/single_outlier.RDS")

# # Calendar year outlier

# config <- patternBreak:::mackConfig(ndev, factors = seq(0.5, 1.5, by = 0.25), type = "calendar")

# results <- mackSim(triangle, nboot, config, "calendar")

# # saveRDS(results, "results/data_objects/calendar_outlier.RDS")

# # Origin year outlier

# config <- patternBreak:::mackConfig(ndev, factors = seq(0.5, 1.5, by = 0.25), type = "origin")

# results <- mackSim(triangle, nboot, config, "origin")

# # saveRDS(results, "results/data_objects/origin_outlier.RDS")
