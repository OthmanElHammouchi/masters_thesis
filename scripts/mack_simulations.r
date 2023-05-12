library(claimsBoot)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
boot.types <- c("residuals", "parametric", "pairs")

# Single outlier
mean.factors <- c(100)
sd.factors <- c(1)
single.res <- mackSim(triangle, "single", 1e3, mean.factors, sd.factors, boot.types)
# saveRDS(single.res, file = file.path("results", "mack_single.RDS"))

# # Calendar year outlier
# mean.factors <- c(0.25, 0.5, 2, 4)
# sd.factors <- c(0.25, 0.5, 2, 4)
# calendar.res <- mackSim(triangle, "calendar", 1e3, mean.factors, sd.factors, boot.types)
# saveRDS(calendar.res, file = file.path("results", "mack_calendar.RDS"))

# # Origin year outlier
# mean.factors <- c(0.25, 0.5, 2, 4)
# sd.factors <- c(0.25, 0.5, 2, 4)
# origin.res <- mackSim(triangle, "origin", 1e3, mean.factors, sd.factors, boot.types)
# saveRDS(origin.res, file = file.path("results", "mack_origin.RDS"))

# library(ggplot2)
# library(formattable)

# i <- 1
# j <- 2

# contaminated <- single.res[
#   outlier.rowidx == i &
#   outlier.colidx == j &
#   (outlier.rowidx != excl.rowidx | outlier.colidx != excl.colidx) &
#     boot.type == "residuals" &
#     opt == "standardised" &
#     cond == TRUE
# ]
# clean <- single.res[
#   outlier.rowidx == i &
#   outlier.colidx == j &
#   (outlier.rowidx == excl.rowidx & outlier.colidx == excl.colidx) &
#     boot.type == "residuals" &
#     opt == "standardised" &
#     cond == TRUE]

# sds <- sapply(split(contaminated, list(contaminated$excl.rowidx, contaminated$excl.colidx), drop = TRUE), function(df) { sd(df$reserve) })

contaminated <- single.res[
  (outlier.rowidx != excl.rowidx | outlier.colidx != excl.colidx) &
    boot.type == "residuals" &
    opt == "log-normal" &
    cond == FALSE
]
clean <- single.res[
  (outlier.rowidx == excl.rowidx & outlier.colidx == excl.colidx) &
    boot.type == "residuals" &
    opt == "log-normal" &
    cond == FALSE]

contaminated <- contaminated[reserve < quantile(reserve, 0.975)]
clean <- clean[reserve < quantile(reserve, 0.975)]

ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(aes(reserve), clean, colour = "red") +
  facet_wrap(vars(outlier.rowidx, outlier.colidx), scales = "free", labeller = label_wrap_gen(multi_line = FALSE))
