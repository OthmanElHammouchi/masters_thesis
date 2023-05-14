library(claimsBoot)
suppressPackageStartupMessages(library(ChainLadder))
library(ggplot2)

calendar.res <- mackSim(UKMotor, "calendar", 1e3, c(2), c(1), "pairs")

single.res <- readRDS("results/mack_single.RDS")
calendar.res <- readRDS("results/mack_calendar.RDS")
origin.res <- readRDS("results/mack_origin.RDS")

contaminated <- single.res[
  (outlier.rowidx != excl.rowidx | outlier.colidx != excl.colidx) &
    boot.type == "pairs"
]
clean <- single.res[
  (outlier.rowidx == excl.rowidx & outlier.colidx == excl.colidx) &
    boot.type == "pairs"
]

ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(aes(reserve), clean, colour = "red") +
  facet_wrap(
    vars(factor(outlier.rowidx), factor(outlier.colidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE)
  )

#######################################################



library(formattable)
min_formatter <- formatter("span", style = x ~ style(color = ifelse(x == min(x), "red", "black")))
formattable(mean.table, lapply(seq_len(nrow(mean.table)), function(row) {
  area(row) ~ min_formatter
}))