library(patternBreak)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
ndev <- nrow(triangle)
nboot <- 1e3

factors <- seq(0.5, 1.5, by = 0.25)
boot_types <- c("parametric", "residuals", "pairs")
proc_dists <- c("normal", "gamma")
conds <- c(TRUE, FALSE)
resids_types <- c("standardised", "modified", "studentised", "lognormal")

# # Single outlier
# single.res <- mackSim(triangle, "single", 1e3, factors, boot_types, proc_dists, conds, resids_types, show_progress = TRUE)
# save(single.res, file = "results/single.rda")

# Calendar year outlier
calendar.res <- mackSim(triangle, "calendar", 1e3, factors, boot_types, proc_dists, conds, resids_types, show_progress = TRUE)
save(calendar.res, file = "results/calendar.rda")

# # Origin year outlier
# origin.res <- mackSim(triangle, "origin", 1e3, factors, boot_types, proc_dists, conds, resids_types, show_progress = TRUE)
# save(origin.res, file = "results/origin.rda")
