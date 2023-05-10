library(claimsBoot)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- cum2incr(UKMotor)
boot.types <- c("residuals", "parametric")

# Single outlier
factors <- c(0.75, 1.25)
single.res <- glmSim(triangle, "single", 1e3, factors, boot.types)
saveRDS(single.res, file = file.path("results", "glm_single.RDS"))

# Calendar year outlier
factors <- c(0.75, 1.25)
calendar.res <- glmSim(triangle, "calendar", 1e3, factors, boot.types)
saveRDS(calendar.res, file = file.path("results", "glm_calendar.RDS"))

# Origin year outlier
factors <- c(0.75, 1.25)
origin.res <- glmSim(triangle, "origin", 1e3, factors, boot.types)
saveRDS(origin.res, file = file.path("results", "glm_origin.RDS"))