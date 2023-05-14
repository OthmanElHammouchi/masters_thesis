library(ggplot2)
library(ggforce)
library(claimsBoot)
library(data.table)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
ndev <- ncol(triangle)

plot.dir <- "plots"

theme_set(theme_bw())

# in mm
width <- 410
height <- 630

mm.per.pt <- 0.3528
width.mm <- mm.per.pt * width
height.mm <- mm.per.pt * height

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm"),
  strip.text.x = element_text(size = 6)
)

## Single outlier:
###############################################################################
single.res <- as.data.table(readRDS("results/glm_single.RDS"))
factors <- unique(single.res$factor)
boot.types <- unique(single.res$boot.type)

f <- 2

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
      "glm_single_densities", bt, o, i_page, sep = "_"),
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
    "glm_single_densities", bt, i_page, sep = "_"),
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

## Calendar outlier:
###############################################################################
calendar.res <- as.data.table(readRDS("results/glm_calendar.RDS"))
factors <- unique(calendar.res$factor)
boot.types <- unique(calendar.res$boot.type)

f <- 2

bt <- "parametric"
opts <- unique(calendar.res[boot.type == bt]$opt)
for (o in opts[!is.na(opts)]) {
  contaminated <- calendar.res[
    boot.type == bt &
      opt == o &
      excl.diagidx != outlier.diagidx &
      factor == f
  ]

  uncontaminated <- calendar.res[
    boot.type == bt &
      opt == o &
      excl.diagidx == outlier.diagidx &
      factor == f
  ]

  contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
  uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

  p <- ggplot() +
    geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
    geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
    facet_wrap(
      vars(factor(outlier.diagidx)),
      scales = "free",
      labeller = label_wrap_gen(multi_line = FALSE),
      ncol = 4,
      nrow = 2
    ) +
    xlab("Reserve") +
    ylab("Density")

  path <- file.path(plot.dir, paste0(paste(
    "glm_calendar_densities", bt, o, sep = "_"),
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


bt <- "residuals"
contaminated <- calendar.res[
  boot.type == bt &
    excl.diagidx != outlier.diagidx &
    factor == f
]

uncontaminated <- calendar.res[
  boot.type == bt &
    excl.diagidx == outlier.diagidx &
    factor == f
]

contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
  facet_wrap(
    vars(factor(outlier.diagidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 4,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "glm_calendar_densities", bt, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

## Origin outlier:
###############################################################################
origin.res <- as.data.table(readRDS("results/glm_origin.RDS"))
factors <- unique(origin.res$factor)
boot.types <- unique(origin.res$boot.type)

f <- 2

bt <- "parametric"
opts <- unique(origin.res[boot.type == bt]$opt)
for (o in opts[!is.na(opts)]) {
  contaminated <- origin.res[
    boot.type == bt &
      opt == o &
      excl.rowidx != outlier.rowidx &
      factor == f
  ]

  uncontaminated <- origin.res[
    boot.type == bt &
      opt == o &
      excl.rowidx == outlier.rowidx &
      factor == f
  ]

  contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
  uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]


  p <- ggplot() +
    geom_density(aes(reserve, group = factor(excl.rowidx)), contaminated) +
    geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
    facet_wrap(
      vars(factor(outlier.rowidx)),
      scales = "free",
      labeller = label_wrap_gen(multi_line = FALSE),
      ncol = 3,
      nrow = 2
    ) +
    xlab("Reserve") +
    ylab("Density")

  path <- file.path(plot.dir, paste0(paste(
    "glm_origin_densities", bt, o, sep = "_"),
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


bt <- "residuals"
contaminated <- origin.res[
  boot.type == bt &
    excl.rowidx != outlier.rowidx &
    factor == f
]

uncontaminated <- origin.res[
  boot.type == bt &
    excl.rowidx == outlier.rowidx &
    factor == f
]

contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
  facet_wrap(
    vars(factor(outlier.rowidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 3,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "glm_origin_densities", bt, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)