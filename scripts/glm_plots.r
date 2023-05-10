library(ggplot2)
library(ggforce)
library(claimsBoot)
library(data.table)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
ndev <- ncol(triangle)

# Single outlier
single.res <- as.data.table(readRDS("results/glm_single.RDS"))
factors <- unique(single.res$factor)
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
  legend.key.size = unit(0.6, "cm"),
  strip.text.x = element_text(size = 6)
)

width <- 410
height <- 630

mm.per.pt <- 0.3528
width.mm <- mm.per.pt * width
height.mm <- mm.per.pt * height

f <- 0.75

bt <- "parametric"
opts <- unique(single.res[boot.type == bt]$opt)
for (o in opts[!is.na(opts)]) {
  contaminated <- single.res[
    boot.type == bt &
      opt == o &
      (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
      factor == f
  ]

  uncontaminated <- single.res[
    boot.type == bt &
      opt == o &
      excl.rowidx == outlier.rowidx &
      excl.colidx == outlier.colidx &
      factor == f
  ]

  contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
  uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

  for (i_page in 1:2) {
    p <- ggplot() +
      geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
      geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
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
      "glm_sim_densities", bt, o, i_page, sep = "_"),
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

bt <- "residuals"
contaminated <- single.res[
  boot.type == bt &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    factor == f
]

uncontaminated <- single.res[
  boot.type == bt &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    factor == f
]

contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

for (i_page in 1:2) {

  p <- ggplot() +
    geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
    geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
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
    "glm_sim_densities", bt, i_page, sep = "_"),
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

####################################
bt <- "parametric"
f <- 1.25
o <- "normal"
i <- 3
j <- 3

contaminated <- single.res[
  boot.type == bt &
    opt == o &
    outlier.rowidx == i &
    outlier.colidx == j &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    factor == f
]

uncontaminated <- single.res[
  boot.type == bt &
    opt == o &
    outlier.rowidx == i &
    outlier.colidx == j &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    factor == f
]

ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
  xlab("Reserve") +
  ylab("Density")

temp <- contaminated[, sd(reserve), by = .(excl.rowidx, excl.colidx)]
