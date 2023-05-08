library(ggplot2)
library(ggforce)
library(patternBreak)
library(data.table)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
ndev <- ncol(triangle)

# Single outlier
single.res <- as.data.table(readRDS("results/single.RDS"))
mean.factors <- unique(single.res$mean.factor)
sd.factors <- unique(single.res$sd.factor)
boot.types <- unique(single.res$boot.type)

# in mm
width <- 418
height <- 591

mm.per.pt <- 0.3528
width.mm <- mm.per.pt * width
height.mm <- mm.per.pt * height

theme_set(theme_bw())

plot.dir <- "plots"

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm")
)

# for (boot.type in c("residuals")) {
#   opts <- unique(single.res[boot.type == boot.type]$opt)

#   for (opt in opts) {
#     # Conditional
#     contaminated <- single.res[
#       boot.type == boot.type &
#         excl.rowidx != outlier.rowidx &
#         excl.colidx != outlier.colidx
#     ]

#     uncontaminated <- single.res[
#       boot.type == boot.type &
#         excl.rowidx == outlier.rowidx &
#         excl.colidx == outlier.colidx
#     ]

#     contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
#     uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

#     p <- ggplot() +
#       geom_density(aes(reserve, group = interaction(outlier.colidx, outlier.rowidx)), contaminated[cond == TRUE]) +
#       geom_density(mapping = aes(reserve), uncontaminated[cond == TRUE], colour = "red", linewidth = 1.25) +
#       facet_grid(rows = vars(factor(mean.factor)), cols = vars(factor(sd.factor)), scales = "free")

#     ggsave(
#       file.path(plot.dir, paste0("mack_sim_densities_cond_", boot.type, opt, ".eps")),
#       p,
#       units = "mm",
#       height = width.mm, # landscape
#       width = height.mm
#     )

#     # Unconditional
#     p <- ggplot() +
#       geom_density(aes(reserve, group = interaction(outlier.colidx, outlier.rowidx)), contaminated[cond == FALSE]) +
#       geom_density(mapping = aes(reserve), uncontaminated[cond == FALSE], colour = "red", linewidth = 1.25) +
#       facet_grid(rows = vars(factor(mean.factor)), cols = vars(factor(sd.factor)), scales = "free")

#     ggsave(
#       file.path(plot.dir, paste0("mack_sim_densities_uncond_", boot.type, opt, ".eps")),
#       p,
#       units = "mm",
#       height = width.mm, # landscape
#       width = height.mm
#     )

#   }
# }

# calendar.res <- as.data.table(readRDS("results/calendar.RDS"))
# origin.res <- as.data.table(readRDS("results/origin.RDS"))

###################

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm"),
  strip.text.x = element_text(size = 6)
)

width <- 410
height <- 630

mm.per.pt <- 0.3528
width.mm <- mm.per.pt * width
height.mm <- mm.per.pt * height

mf <- 0.75
sf <- 0.75

for (bt in c("residuals", "parametric")) {
  opts <- unique(single.res[boot.type == bt]$opt)

  for (o in opts) {

    contaminated <- single.res[
      boot.type == bt &
      opt == o &
        (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
        mean.factor == mf &
        sd.factor == sf
    ]

    uncontaminated <- single.res[
      boot.type == boot.type &
      opt == opt &
        excl.rowidx == outlier.rowidx &
        excl.colidx == outlier.colidx &
        mean.factor == mf &
        sd.factor == sf
    ]

    contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
    uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

    for (i_page in 1:2) {

      p <- ggplot() +
        geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated[cond == TRUE]) +
        geom_density(mapping = aes(reserve), uncontaminated[cond == TRUE], colour = "red") +
        facet_wrap_paginate(
          vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
          scales = "free",
          labeller = label_wrap_gen(multi_line = FALSE),
          ncol = 4,
          nrow = 3,
          page = i_page
        ) +
        xlab("Reserve") +
        ylab("Density")

      path <- file.path(plot.dir, paste0(paste(
        "mack_sim_densities_cond", bt, o, i_page, sep = "_"),
        ".eps")
      )

      ggsave(
        path,
        p,
        units = "mm",
        height = width.mm, # landscape
        width = height.mm
      )

      p <- ggplot() +
        geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated[cond == FALSE]) +
        geom_density(mapping = aes(reserve), uncontaminated[cond == FALSE], colour = "red") +
        facet_wrap_paginate(
          vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
          scales = "free",
          labeller = label_wrap_gen(multi_line = FALSE),
          ncol = 4,
          nrow = 3,
          page = i_page
        ) +
        xlab("Reserve") +
        ylab("Density")

      path <- file.path(plot.dir, paste0(paste(
        "mack_sim_densities_uncond", bt, o, i_page, sep = "_"),
        ".eps")
      )

      ggsave(
        path,
        p,
        units = "mm",
        height = width.mm, # landscape
        width = height.mm
      )

    }
  }
}
