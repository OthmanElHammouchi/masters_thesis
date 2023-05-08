library(patternBreak)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
boot.types <- c("residuals", "parametric", "pairs")

# Single outlier
mean.factors <- c(0.5, 0.75, 1.25, 1.5)
sd.factors <- c(0.75, 1.25)
single.res <- patternBreak:::mackSim(triangle, "single", 5e2, mean.factors, sd.factors, boot.types)
saveRDS(single.res, file = "results/single.RDS")

# Calendar year outlier
mean.factors <- c(0.8, 0.9, 1.1, 1.2)
sd.factors <- c(0.8, 0.9, 1.1, 1.2)
calendar.res <- patternBreak:::mackSim(triangle, "calendar", 1e3, mean.factors, sd.factors, boot.types)
saveRDS(calendar.res, file = "results/calendar.RDS")

# Origin year outlier
mean.factors <- c(0.8, 0.9, 1.1, 1.2)
sd.factors <- c(0.8, 0.9, 1.1, 1.2)
origin.res <- patternBreak:::mackSim(triangle, "origin", 1e3, mean.factors, sd.factors, boot.types)
saveRDS(origin.res, file = "results/origin.RDS")
