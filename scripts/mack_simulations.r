library(claimsBoot)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
boot.types <- c("residuals", "parametric", "pairs")

# Single outlier
mean.factors <- c(0.5, 1, 2)
sd.factors <- c(0.5, 1, 2)
single.res <- mackSim(triangle, "single", 1e3, mean.factors, sd.factors, boot.types)
saveRDS(single.res, file = file.path("results", "mack_single.RDS"))

# Calendar year outlier
mean.factors <- c(0.5, 1, 2)
sd.factors <- c(0.5, 1, 2)
calendar.res <- mackSim(triangle, "calendar", 1e3, mean.factors, sd.factors, boot.types)
saveRDS(calendar.res, file = file.path("results", "mack_calendar.RDS"))

# Origin year outlier
mean.factors <- c(0.5, 1, 2)
sd.factors <- c(0.5, 1, 2)
origin.res <- mackSim(triangle, "origin", 1e3, mean.factors, sd.factors, boot.types)
saveRDS(origin.res, file = file.path("results", "mack_origin.RDS"))

##############################################################################

library(ggplot2)
single.res <- readRDS("results/mack_single.RDS")

contaminated <- single.res[
  (outlier.rowidx != excl.rowidx | outlier.colidx != excl.colidx) &
    boot.type == "residuals" &
    opt == "log-normal" &
    cond == FALSE &
    mean.factor == 2 &
    sd.factor == 1
]
clean <- single.res[
  (outlier.rowidx == excl.rowidx & outlier.colidx == excl.colidx) &
    boot.type == "residuals" &
    opt == "log-normal" &
    cond == FALSE &
    mean.factor == 2 &
    sd.factor == 1
]

ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(aes(reserve), clean, colour = "red") +
  facet_wrap(vars(factor(outlier.rowidx), factor(outlier.colidx)), scales = "free", labeller = label_wrap_gen(multi_line = FALSE))
